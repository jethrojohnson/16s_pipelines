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
===========================
pipeline_extend_16s.py
===========================
For every genome in refseq, download and extract the area surrounding
each 16S gene. 

"""

from ruffus import *

import sys
import os
import re
import time
import glob
import math
import pickle
import gzip
import sqlite3
import collections
import shutil
import pandas as pd
import numpy as np
import random
import logging as L
import wget
import MySQLdb
import ftplib

import rpy2
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
pandas2ri.activate()
from rpy2.robjects import r as R

import CGAT.Experiment as E
import CGATPipelines.Pipeline as P
import CGAT.FastaIterator as FastaIterator
import CGAT.IOTools as IOTools


###############################################################################
# Specify parameters and set utility functions
###############################################################################

# load options from the config file
PARAMS = P.getParameters(
    ["%s/pipeline.ini" % os.path.splitext(__file__)[0],
     "../pipeline.ini",
     "pipeline.ini"])

# Utility functions
def connect():
    '''utility function to connect to database.

    Use this method to connect to the pipeline database.
    Additional databases can be attached here as well.

    Returns an sqlite3 database handle.
    '''

    dbh = sqlite3.connect(PARAMS["database"])

    return dbh


def fetchBioSQLTaxonomy(ncbi_tax_id,
                        host='weinstock_dbs',
                        user='guest',
                        passwd='',
                        db='biosql'):
    '''Fetch the full parent taxonomy (to Kingdom) for
       a given NCBI tax id number'''

    phylogeny = collections.OrderedDict()

    # connect to BioSQL database
    db = MySQLdb.connect(host=host, user=user, passwd=passwd, db=db)
    cur = db.cursor()

    # break if more than 50 iterations...
    safety = 50
    first = True
    node_rank = None
    while node_rank != 'superkingdom':
        safety -= 1
        if safety < 0:
            break

        # internal BioSQL taxon_ids don't correspond to ncbi taxon IDs
        # WARNING: There is overlap between the two. 
        field = 'taxon_id'
        if first:
            tax_id = ncbi_tax_id
            field = 'ncbi_taxon_id'
            first = False
            
        # fetch tax_id, parent_tax_id node_rank, scientific_name
        statement = ("SELECT "
                     "  taxon.taxon_id,"
                     "  taxon.parent_taxon_id,"
                     "  taxon.node_rank,"
                     "  taxon_name.name"
                     " FROM taxon JOIN taxon_name"
                     "  ON taxon.taxon_id = taxon_name.taxon_id"
                     " WHERE taxon.{} = {} "
                     "AND name_class = 'scientific name'".format(field, tax_id))
        cur.execute(statement)
        taxonomy = cur.fetchall()
        if len(taxonomy) == 1:
            old_tax_id, tax_id, node_rank, scientific_name = taxonomy[0]
            phylogeny[node_rank] = scientific_name
        elif len(taxonomy) == 0:
            E.warn('No taxonomy information for taxon ID %s beyond %s' \
                   % (ncbi_tax_id, node_rank))
            phylogeny['genus'] = 'unavailable'
            phylogeny['species'] = 'unavailable'
            
            return phylogeny

        else:
            raise ValueError('multiple entries for tax id %s: %s' \
                             % (ncbi_tax_id, str(taxonomy)))
            
    return phylogeny


###############################################################################
## Section: Fetch NCBI RefSeq/Genbank Genome Assemblies For Taxa of Interest
###############################################################################
@originate('ncbi_assembly_summary.txt')
def fetchAssemblySummary(outfile):
    '''Fetch the NCBI assembly summary file... and remove UTF8 chars'''

    # wget assumes filename extension
    tmpf = P.getTempFilename('.')`
    os.unlink(tmpf)
    wget.download(PARAMS['location_assembly_summary'], out=tmpf)

    
    statement = ("iconv -f UTF8 -t ascii//translit %(tmpf)s > %(outfile)s"
    to_cluster = False
    P.run()

    os.unlink(tmpf)


@transform(fetchAssemblySummary, suffix('.txt'), '.load')
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
       'ncbi_accessions.txt')
def fetchTaxonAccessions(infile, outfile):
    '''Find the accession details of all genomes belonging to taxon
    of interest'''
    
    infile = P.snip(infile[0], '.load')
    outf = IOTools.openFile(outfile, 'w')
    
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
           
    df = pd.read_sql(sql=statement, con=connect())

    for row in df.iterrows():
        row = row[1]
        genome_id = os.path.basename(row['ftp_path'])
        genome_id = genome_id + '_genomic.fna.gz'
        ftp_path = row['ftp_path'] + '/' + genome_id
        
        outf.write('\t'.join(map(str, [row["assembly_accession"],
                                       row["bioproject"],
                                       row["biosample"],
                                       row["taxid"],
                                       row["species_taxid"],
                                       row["organism_name"],
                                       ftp_path])) + '\n')
    outf.close()


@follows(mkdir('ncbi_assemblies.dir'))
@subdivide(fetchTaxonAccessions,
           regex('(.+).txt'),
           r'ncbi_assemblies.dir/\1_*.txt')
def splitTaxonAccessions(infile, outfiles):
    '''Split the download list into groups of 10 for download'''

    outf_stub = os.path.join('ncbi_assemblies.dir',
                             P.snip(os.path.basename(infile), '.txt'))
    outf = IOTools.openFile(outf_stub + '_00000.txt', 'w')
    for n, line in enumerate(IOTools.openFile(infile)):
        if n == 0:
            continue
        if n % 10 == 0:
            outf.close()
            outf = IOTools.openFile(outf_stub + '_' + str(n).zfill(5) + '.txt',
                                    'w') 
        outf.write(line)
    outf.close()


@subdivide(splitTaxonAccessions,
           regex('(.+)/(.+).txt'),
           r'\1/\2.dir/*gz')
def fetchGenomeAssemblyAndAnnotations(infile, outfiles):
    '''Download all assemblies using rsync, as per ncbi recommendation'''
    
    outdir = P.snip(infile, '.txt') + '.dir'
    if os.path.exists(outdir):
        shutil.rmtree(outdir)
    os.mkdir(outdir)

    # out_fail = os.path.join(outdir, 'failed.txt')
    # out_fail = IOTools.openFile(out_fail, 'w')
    
    for genome in IOTools.openFile(infile):
        genome_path = genome.split().pop()
        genome_path = re.sub('ftp', 'rsync', genome_path)
        gff_path = P.snip(genome_path, '.fna.gz') + '.gff.gz'

        ## Something about ftplib hangs after the first few hundred queries...
        # # not all entries have an '*_rna_from_genomic.fna.gz' file
        # subdir = os.path.dirname(genome_path)\
        #     [len('ftp://ftp.ncbi.nlm.nih.gov/'):]
    
        # server = 'ftp.ncbi.nlm.nih.gov'
        # user = 'anonymous'
        # password = 'jethro.johnson@jax.org'
    
        # ftp = ftplib.FTP(server)
        # ftp.login(user,password)
        # l = [os.path.basename(x) for x in ftp.nlst(subdir)]
    
        # if os.path.basename(genome_path) in l \
        #    and os.path.basename(gff_path) in l:

        statement = (" rsync --copy-links --quiet"
                     "  %(genome_path)s"
                     "  %(outdir)s &&"
                     " rsync --copy-links --quiet"
                     "  %(gff_path)s"
                     "  %(outdir)s")
            # to_cluster = False
        P.run()
    #     else:
    #         out_fail.write(genome)
        
    # out_fail.close()


###############################################################################
# Fetch full taxonomy for downloaded assemblies
###############################################################################
@transform(fetchTaxonAccessions,
           regex('(.+)_accessions.txt'),
           r'\1_taxonomy.txt')
def fetchAssemblyTaxonomy(infile, outfile):
    '''Output the full taxonomy for the downloaded assemblies, WARNING:
    empty dictionaries are replaced with genus "unknown", species "unknown"'''

    output_dict = {}
    
    header = True
    for line in IOTools.openFile(infile):
        line = line.split('\t')
        if header:
            assert line[4] == 'species_taxid'
            header = False
            continue
        
        ncbi_taxon_id = line[4]

        # key = assembly_accession, value = full taxonomy
        output_dict[line[0]] = fetchBioSQLTaxonomy(ncbi_taxon_id)
        E.info('Fetched taxonomy for %s' % ncbi_taxon_id)
        
    output_df = pd.DataFrame.from_dict(output_dict, orient='index')
    output_df.to_csv(IOTools.openFile(outfile, 'w'),
                     sep='\t',
                     index_label='accession_id')


@transform(fetchAssemblyTaxonomy, suffix('.txt'), '.load')
def loadAssemblyTaxonomy(infile, outfile):
    df = pd.read_table(infile)
    table = P.snip(outfile, '.load', strip_path=True)
    df.to_sql(name=table, con=connect(), if_exists='replace')
    open(outfile, 'w').close()


###############################################################################
# Slop gtf intervals and extract extended 16S gene sequences. 
###############################################################################
@transform(fetchGenomeAssemblyAndAnnotations, suffix('.fna.gz'), '.fasta')
def indexGenomeFastas(infile, outfile):
    '''Create a simple index of the genome fasta files'''

    indexed_fasta = P.snip(infile, '.fna.gz')
    statement = ("python %(scriptsdir)s/index_fasta.py"
                 " --log=%(indexed_fasta)s.log"
                 " %(indexed_fasta)s"
                 " %(infile)s")
    P.run()
    

@follows(mkdir('ncbi_16s_genes.dir'))
@transform(indexGenomeFastas,
           regex('.+/(.+)_genomic.fasta'),
           r'ncbi_16s_genes.dir/\1_16s_genes.fasta.gz')
def fetch16Sgenes(infile, outfile):
    '''Fetch fasta sequences for 16S genes annotated in gff'''

    in_gff = P.snip(infile, '.fasta') + '.gff.gz'

    # There are some assemblies with no 16S gene annotations
    matched = False
    annotation = re.compile(PARAMS['search_regex'], re.IGNORECASE)
    for line in IOTools.openFile(in_gff):
        if annotation.search(line):
            matched = True
            break

    if matched:
        statement = ("zcat %(in_gff)s |"
                     " grep -i '%(search_regex)s' |"
                     " python %(scriptsdir)s/gff2gff.py"
                     "  --genome-file=%(infile)s"
                     "  --method=add-flank"
                     "  --flank-method=extend"
                     "  --extension-upstream=%(slop_upstream)s"
                     "  --extension-downstream=%(slop_downstream)s"
                     "  --log=%(outfile)s.log |"
                     " python %(scriptsdir)s/gff2fasta.py"
                     "  --genome-file=%(infile)s"
                     "  --log=%(outfile)s.log |"
                     " gzip > %(outfile)s")
        job_options = '-l walltime=00:15:00,mem=1GB,nodes=1:ppn=1'
        
        P.run()
    else:
        IOTools.openFile(outfile, 'w').close()


@jobs_limit(10)
@transform(fetch16Sgenes, suffix('.fasta.gz'), '_renamed.fasta.gz')
def rename16SGenes(infile, outfile):
    '''Add the file name to the fasta headers'''

    filename = P.snip(os.path.basename(infile), '_16s_genes.fasta.gz')

    outf = IOTools.openFile(outfile, 'w')
    
    n = 0
    for fasta in FastaIterator.FastaIterator(IOTools.openFile(infile)):
        n += 1
        header = filename + '.' + str(n)
        outf.write('\n'.join(['>' + header + ' ' + fasta.title,
                              fasta.sequence]) + '\n')

    outf.close()               

    
@collate(rename16SGenes,
         regex('(.+)/.+_16s_genes_renamed.fasta.gz'),
         r'\1/all_16s_genes_stats.tsv.gz')
def fetch16SGeneStats(infiles, outfile):
    '''Fetch the basic statistics for all 16S genes downloaded from NCBI
    '''

    infiles = ' '.join(infiles)
    statement = ("zcat %(infiles)s |"
                 " python %(scriptsdir)s/fasta2table.py"
                 "  --split-fasta-identifier"
                 "  --section=gaps,na,length"
                 "  --log=%(outfile)s.log |"
                 " gzip > %(outfile)s")
    P.run()
    
    
@collate(rename16SGenes,
         regex('(.+)/.+_16s_genes_renamed.fasta.gz'),
         r'\1/all_16s_genes_filtered.fasta.gz')
def collate16SGenes(infiles, outfile):
    '''Collate all 16S gene sequences, at the same time filter to remove
    those above/below specified length thresholds.
    '''

    # Too many files for zcat... 
    # infiles = ' '.join(infiles)
    rgx = os.path.dirname(outfile) + '/*_renamed.fasta.gz'
    tmpf = P.getTempFilename('.')
    statement = ("for FILE in %(rgx)s; do cat $FILE >> %(tmpf)s; done &&" 
                 " zcat %(tmpf)s |"
                 " python %(scriptsdir)s/fasta2fasta.py"
                 "  --method=filter"
                 "  --filter-method=min-length=%(filter_min_length)s"
                 "  --filter-method=max-length=%(filter_max_length)s"
                 "  --log=%(outfile)s.log |"
                 " gzip > %(outfile)s")
    to_cluster = False
    P.run()



###############################################################################
# Filter 16S genes using HMMR3
###############################################################################
@follows(mkdir('ncbi_16s_genes_filtered.dir'))
@transform(collate16SGenes,
           regex('.+/all_16s_genes_filtered.fasta.gz'),
           r'ncbi_16s_genes_filtered.dir/16s_genes_hmmr_hits.tsv.gz')
def filter16SGenesWithHmmr(infile, outfile):
    '''Use hmmr3 trained on mock community sequences to filter 
    16S genes'''

    outfile = P.snip(outfile, '.gz')
    outf_main = P.snip(outfile, '.tsv') + '_main.tsv'
    
    statement = ("nhmmscan"
                 " -o %(outf_main)s"
                 " --tblout %(outfile)s"
                 " %(filter_nhmmscan)s &&"
                 " gzip %(outfile)s"
                 " gzip %(outf_main)s")
    to_cluster = False
    P.run()
    
    

    
###############################################################################
# Align 16S genes and perform entropy calculations
###############################################################################
@transform(collate16SGenes, suffix('_filtered.fasta.gz'), '_base.fasta')
def extractBaselineSequence(infile, outfile):
    '''Pull out baseline sequence along which to calculate entropy.'''

    base_seq_id = PARAMS['entropy_base']
    check = False
    with IOTools.openFile(outfile, 'w') as outf:
        for fasta in FastaIterator.FastaIterator(IOTools.openFile(infile)):
            if fasta.title.startswith(base_seq_id):
                outf.write('>BASE\n' + fasta.sequence + '\n')
                check = True
            else:
                continue
    assert check, "Couldn't find base sequence in ref db"


@transform(collate16SGenes,
           suffix('.fasta.gz'),
           add_inputs(extractBaselineSequence),
           '_aligned.fasta.gz')
def align16SGenes(infiles, outfile):
    '''Align the extracted sequences using muscle'''

    seq_fasta, base_fasta = infiles
    tmpf = P.getTempFilename('.')
    tmpf_out = P.snip(outfile, '.gz')

    statement = ("zcat %(seq_fasta)s > %(tmpf)s &&"
                 " cat %(base_fasta)s >> %(tmpf)s &&"
                 " muscle -in %(tmpf)s -out %(tmpf_out)s"
                 " -maxhours %(muscle_run_time)s &&"
                 " gzip %(tmpf_out)s")
    cluster_options=PARAMS['muscle_run_options']
    P.run()
    
    os.unlink(tmpf)


@transform(align16SGenes, suffix('.fasta.gz'), '_entropy.tsv')
def  calculateEntropy(infile, outfile):
    '''Calculate the shannon entropy for multiple sequence alignments'''

    entropy_scr = os.path.join(os.path.dirname(__file__),
                               'alignment2entropy.py')
    statement = ("python {entropy_scr}"
                 " -a {infile}"
                 " -b BASE"
                 " -o {outfile}")
    P.run()

    
if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
