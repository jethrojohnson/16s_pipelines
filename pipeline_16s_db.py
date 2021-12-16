###############################################################################
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
############################################################################### 
"""
==================
pipeline_16s_db.py
==================
A pipeline to create a 16s sequence reference database from NCBI genome
assemblies.
"""


import os,sys,re
import sqlite3
import pandas as pd
import wget
import shutil

import rpy2
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
pandas2ri.activate()
from rpy2.robjects import r as R

from ruffus import *
from cgatcore import pipeline as P
from cgatcore import iotools as IOTools
from cgatcore import experiment as E
import cgat.FastaIterator as FastaIterator

###############################################################################
# Utility Settings
###############################################################################
# load options from the config file
PARAMS = P.get_parameters(
    ["%s/pipeline.yml" % os.path.splitext(__file__)[0],
     "pipeline.yml"])


# Utility functions
def connect():
    '''Returns an sqlite3 database handle.
    '''

    dbh = sqlite3.connect(PARAMS["database"])

    return dbh


###############################################################################
# Pipeline steps
###############################################################################
## 1. Details of genomes to be downloaded 
###############################################################################
@originate('ncbi_assembly_summary.tsv.gz')
def fetchAssemblySummary(outfile):
    '''Fetch the NCBI assembly summary file'''

    tmpf = P.get_temp_filename('.')
    wget.download(PARAMS['location_assembly_summary'], tmpf)

    # UTF8 encoding in file is causing problems downstream
    statement = ("iconv -f UTF8 -t ascii//translit %(tmpf)s | gzip > %(outfile)s")
    P.run(statement, to_cluster=False)

    os.unlink(tmpf)


@transform(fetchAssemblySummary, suffix('.tsv.gz'), '.load')
def loadAssemblySummary(infile, outfile):
    '''Load current assembly summary'''
    
    # fetch assembly summary table and edit
    df = pd.read_table(infile, skiprows=1, index_col=0)
    df.index.name = 'assembly_accession'

    # load... will throw error if table already exists
    table_name = P.snip(outfile, '.load', strip_path=True)
    df.to_sql(name=table_name, con=connect())

    open(outfile, 'w').close()


@files(loadAssemblySummary,
       'ncbi_accessions.tsv.gz')
def fetchAccessions(infile, outfile):
    '''Fetch details of genomes to be fetched with rsync
    Option to specify WHERE statement in pipeline configuration file
    '''
    
    infile = P.snip(infile[0], '.load')
    outf = IOTools.open_file(outfile, 'w')
    outf_failed = IOTools.open_file(P.snip(outfile, '.tsv.gz') \
                                    + '_failed.tsv.gz', 'w')
    
    outf.write('\t'.join(["assembly_accession",
                          "bioproject",
                          "biosample",
                          "taxid",
                          "species_taxid",
                          "organism_name",
                          "ftp_path"]) + '\n')
        
    statement = ("SELECT"
                 "  assembly_accession,"
                 "  bioproject,"
                 "  biosample,"
                 "  taxid,"
                 "  species_taxid,"
                 "  organism_name,"
                 "  ftp_path"
                 " FROM {infile}".format(**locals()))

    if PARAMS['genomes_selection']: 
        statement = statement + " WHERE " + PARAMS['genomes_selection']
    
    df = pd.read_sql(sql=statement, con=connect())

    for row in df.iterrows():
        row = row[1]
        genome_id = os.path.basename(row['ftp_path'])
        genome_id = genome_id + '_genomic.fna.gz'
        ftp_path = row['ftp_path'] + '/' + genome_id

        line_out = '\t'.join(map(str, [row["assembly_accession"],
                                       row["bioproject"],
                                       row["biosample"],
                                       row["taxid"],
                                       row["species_taxid"],
                                       row["organism_name"],
                                       ftp_path])) + '\n'

        # There are entries with 'na' in ftp_path
        if not ftp_path.startswith('http'):
            outf_failed.write(line_out)
        else:
            outf.write(line_out)
    outf.close()


@transform(fetchAccessions, suffix('_accessions.tsv.gz'), '_taxonomy.tsv.gz')
def fetchTaxonomy(infile, outfile):
    '''Fetch full taxonomy for each ncbi taxid'''
     
    assert IOTools.open_file(infile).readline().split()[3] == 'taxid'

    names_dmp = os.path.join(PARAMS['location_ncbi_taxonomy'], 'names.dmp')
    nodes_dmp = os.path.join(PARAMS['location_ncbi_taxonomy'], 'nodes.dmp')
    assert os.path.exists(names_dmp), names_dmp
    assert os.path.exists(nodes_dmp), nodes_dmp

    statement = ("zcat %(infile)s | cut -f4 |"
                 " python %(mbscriptsdir)s/taxid2taxonomy.py"
                 "  --names-dmp %(names_dmp)s"
                 "  --nodes-dmp %(nodes_dmp)s"
                 "  --header"
                 "  --header-out"
                 "  --log %(outfile)s.log |"
                 " gzip > %(outfile)s")
    P.run(statement)


@transform(fetchTaxonomy, suffix('.tsv.gz'), '.load')
def loadTaxonomy(infile, outfile):

    # fetch assembly summary table and edit
    df = pd.read_table(infile)

    # load... will throw error if table already exists
    table_name = P.snip(outfile, '.load', strip_path=True)
    df.to_sql(name=table_name, con=connect(), if_exists='replace', index=False)

    open(outfile, 'w').close()


###############################################################################
# Fetch genomes and genome annotations
###############################################################################
@follows(mkdir('ncbi_assemblies.dir'))
@subdivide(fetchAccessions,
           regex('(.+).tsv.gz'),
           r'ncbi_assemblies.dir/\1_*.tsv.gz')
def splitTaxonAccessions(infile, outfiles):
    '''Split the download list into groups of 10 for download'''

    # Handle the number of accessions appropriately when chunking files
    l = len(str(sum(1 for l in IOTools.open_file(infile))))

    outf_stub = os.path.join('ncbi_assemblies.dir',
                             P.snip(os.path.basename(infile), '.tsv.gz'))
    outf = IOTools.open_file(outf_stub + '_' + '0'*l + '.tsv.gz', 'w')
    for n, line in enumerate(IOTools.open_file(infile)):
        if n == 0:
            continue
        if n % 10 == 0:
            outf.close()
            outf = IOTools.open_file(outf_stub + '_' + str(n).zfill(l) + '.tsv.gz',
                                    'w') 
        outf.write(line)
    outf.close()


@subdivide(splitTaxonAccessions,
           regex('(.+)/(.+).tsv.gz'),
           r'\1/\2.dir/*gz')
def fetchGenomeAssemblyAndAnnotations(infile, outfiles):
    '''Download all assemblies using rsync, as per ncbi recommendation'''
    
    outdir = P.snip(infile, '.tsv.gz') + '.dir'
    if os.path.exists(outdir):
        shutil.rmtree(outdir)
    os.mkdir(outdir)

    for genome in IOTools.open_file(infile):
        genome_path = genome.split().pop()
        genome_path = re.sub('ftp', 'rsync', genome_path)
        genome_path = re.sub('https', 'rsync', genome_path)
        gff_path = P.snip(genome_path, '.fna.gz') + '.gff.gz'

        statement = (" rsync --copy-links --quiet"
                     "  %(genome_path)s"
                     "  %(outdir)s &&"
                     " rsync --copy-links --quiet"
                     "  %(gff_path)s"
                     "  %(outdir)s")
        P.run(statement)


###############################################################################
# Extract and optionally extended 16S genes from microbial genomes. 
###############################################################################
@transform(fetchGenomeAssemblyAndAnnotations, suffix('.fna.gz'), '.fasta')
def indexGenomeFastas(infile, outfile):
    '''Index each genome fasta file'''

    indexed_fasta = P.snip(infile, '.fna.gz')
    statement = ("cgat index_fasta"
                 " --log=%(indexed_fasta)s.log"
                 " %(indexed_fasta)s"
                 " %(infile)s")
    P.run(statement)

@follows(

    
###############################################################################
# Generic pipeline tasks
###############################################################################
@follows(loadTaxonomy)
def full():
    pass


def main(argv=None):
    if argv is None:
        argv = sys.argv
    P.main(argv)


if __name__ == "__main__":
    sys.exit(P.main(sys.argv))    
