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
import numpy as np
import wget
import shutil

from ruffus import *
from cgatcore import pipeline as P
from cgatcore import iotools as IOTools
from cgatcore import experiment as E

import cgat.FastaIterator as FastaIterator
import cgat.GTF as GTF

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
#@originate('ncbi_assembly_summary.tsv.gz')
@files('assembly_summary.txt', 'ncbi_assembly_summary.tsv.gz')
def fetchAssemblySummary(infile, outfile):
    '''Fetch the NCBI assembly summary file'''

    # tmpf = P.get_temp_filename('.')
    # wget.download(PARAMS['location_assembly_summary'], tmpf)

    # UTF8 encoding in file is causing problems downstream
    statement = ("iconv -f UTF8 -t ascii//translit %(infile)s | gzip > %(outfile)s")
    P.run(statement, to_cluster=False)

    #os.unlink(tmpf)


@transform(fetchAssemblySummary, suffix('.tsv.gz'), '.load')
def loadAssemblySummary(infile, outfile):
    '''Load current assembly summary'''
    
    # fetch assembly summary table and edit
    df = pd.read_table(infile, skiprows=1, index_col=0, sep="\t")
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

    infile = P.snip(infile, '.load')
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
## 2.  Fetch genomes and genome annotations
###############################################################################
@follows(mkdir('01_ncbi_assemblies.dir'))
@subdivide(fetchAccessions,
           regex('(.+).tsv.gz'),
           r'01_ncbi_assemblies.dir/\1_*.tsv.gz')
def splitTaxonAccessions(infile, outfiles):
    '''Split the download list into groups of 10 for download'''

    # Handle the number of accessions appropriately when chunking files
    l = len(str(sum(1 for l in IOTools.open_file(infile))))

    outf_stub = os.path.join('01_ncbi_assemblies.dir',
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
#           regex('(.+)/(ncbi_accessions_14380).tsv.gz'),
           r'\1/\2.dir/*.gff.gz')
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
        rna_path = P.snip(genome_path, '_genomic.fna.gz') + '_rna_from_genomic.fna.gz' 
        
        statement = (" rsync --copy-links --quiet"
                     "  %(genome_path)s"
                     "  %(outdir)s &&"
                     " rsync --copy-links --quiet"
                     "  %(gff_path)s"
                     "  %(outdir)s &&"
                     " rsync --copy-links --quiet"
                     "  %(rna_path)s"
                     "  %(outdir)s")
        P.run(statement)


###############################################################################
## 3. Extract and optionally extended 16S genes from microbial genomes. 
###############################################################################
@active_if(PARAMS['target_from'] == 'genome')
@jobs_limit(PARAMS['local_n'])
@transform(fetchGenomeAssemblyAndAnnotations,
           regex('(.+)/(.+)_genomic.gff.gz'),
           r'\1/\2_to_keep.gtf.gz')
def selectTargetFeatures(infile, outfile):
    '''Iterate over feature file and select features to be kept
    based on criteria specified in configuration file'''

    locations = []
    sources = [x.lower() for x in PARAMS['target_source'].split(',')]
    features = [x.lower() for x in PARAMS['target_feature'].split(',')]

    with IOTools.open_file(outfile, 'w') as outf:
        for gtf in GTF.iterator(IOTools.open_file(infile)):
            if gtf.source.lower() in sources \
               and gtf.feature.lower() in features \
               and re.search(PARAMS['target_attribute'], gtf.attributes, re.IGNORECASE):
                loc = '_'.join(map(str, [gtf.start, gtf.end]))
                if loc in locations:
                    continue
                else:
                    outf.write(str(gtf) + '\n')
            else:
                continue
    

@transform(selectTargetFeatures,
           suffix('_to_keep.gtf.gz'),
           '_genomic.fasta')
def indexGenomeFastas(infile, outfile):
    '''Index each genome fasta file'''

    indexed_fasta = P.snip(infile, '_to_keep.gtf.gz') + '_genomic'
    input_fasta = indexed_fasta + '.fna.gz'
    statement = ("cgat index_fasta"
                 " --log=%(indexed_fasta)s.log"
                 " %(indexed_fasta)s"
                 " %(input_fasta)s")
    P.run(statement)


@follows(mkdir('02_ncbi_target_regions.dir'))
@transform(indexGenomeFastas,
           regex('.+/(.+)_genomic.fasta'),
           r'02_ncbi_target_regions.dir/\1_targets.fasta.gz')
def fetchTargetRegion(infile, outfile):
    '''Fetch fasta sequences for 16S genes annotated in gff3'''

    in_gff = P.snip(infile, '_genomic.fasta') + '_to_keep.gtf.gz'

    # Possibility of assemblies with no 16S gene annotations
    matched = False
    if IOTools.open_file(in_gff).readlines():
        matched = True

    if matched:
        statement = ("zcat %(in_gff)s |"
                     " cgat gff2gff"
                     "  --genome-file=%(infile)s"
                     "  --method=add-flank"
                     "  --flank-method=extend"
                     "  --extension-upstream=%(target_slop_up)s"
                     "  --extension-downstream=%(target_slop_down)s"
                     "  --log=%(outfile)s.log |"
                     " cgat gff2fasta"
                     "  --genome-file=%(infile)s"
                     "  --log=%(outfile)s.log |"
                     " gzip > %(outfile)s")        
        P.run(statement)
    else:
        IOTools.open_file(outfile, 'w').close()

###############################################################################
# Extract targets from rna feature file
@active_if(PARAMS['target_from'] == 'rna')
@follows(mkdir('02_ncbi_target_regions.dir'))
@transform(fetchGenomeAssemblyAndAnnotations,
           regex('.+/(.+)_genomic.gff.gz'),
           r'ncbi_target_regions.dir/\1_targets.fasta.gz')
def fetchTargetRegionFromRNAFeatureFile(infile, outfile):
    '''As an alternative to fetching target region from genome sequence,
    fetch it directly from the rna fasta file of genome features. 
    '''

    infile = P.snip(infile, '_genomic.gff.gz') + '_rna_from_genomic.fna.gz'

    statement = ("zcat %(infile)s |"
                 " awk '/^>/ {printf(\"\\n%%s\\n\",$0);next; }"
                 "  { printf(\"%%s\",$0);}  END {printf(\"\\n\");}' |"
                 " grep -A1 '%(target_regex)s' |"
                 " sed '/^--$/d' |" # remove grep separator for context lines
                 " gzip > %(outfile)s")
    P.run(statement)

    
###############################################################################
# Summarize and filter out partial annotations
@follows(fetchTargetRegionFromRNAFeatureFile, fetchTargetRegion)
@merge('02_ncbi_target_regions.dir/*fasta.gz',
       '02_ncbi_target_regions_passed.tsv.gz')
def summarizeAndFilterTargetRegions(infiles, outfile):
    '''Summarise the number and length of target regions, 
    additionally filter out genomes that have one or more target region
    that's smaller than expected (hardcoded for 16S). 
    As not all genomes are taken forward...@merge is used.
    '''
    passed = []
    failed = []
    
    for infile in infiles:
        if not IOTools.open_file(infile).readlines():
            failed.append([infile, 'no_genes'])
            continue

        lengths = []
        for fasta in FastaIterator.FastaIterator(IOTools.open_file(infile)):
            lengths.append(len(fasta.sequence))

        max_l = max(lengths)
        min_l = min(lengths)
        diff = max_l - min_l

        ## Some rudimentary filtering...
        # dropping genomes with one or more partial 16s, based on length
        if max_l < 1450:
            failed.append([infile, 'no_complete_annotation'])
        elif diff > 0.01*max_l:
            failed.append([infile, 'one_or_more_incomplete_annotation'])
        else:
            passed.append('\t'.join(map(str, [infile,
                                              len(lengths),
                                              np.mean(lengths),
                                              np.median(lengths),
                                              max_l,
                                              min_l,
                                              diff,
                                              ','.join([str(x) for x in lengths])])))

    out_failed = P.snip(outfile, '_passed.tsv.gz') + '_failed.tsv.gz'
    with IOTools.open_file(out_failed, 'w') as outf:
        outf.write('\n'.join(['\t'.join(x) for x in failed]))

    with IOTools.open_file(outfile, 'w') as outf:
        outf.write('\t'.join(['accession', 'n_16s_genes',
                              'mean_length', 'median_length',
                              'max_length', 'min_length', 'diff']) + '\n')
        outf.write('\n'.join(passed))    


###############################################################################
# 4. Filter and align target sequences for each genome
###############################################################################
# Filter target sequences
# Perhaps add filtering step based on hmmscan (see refseq hmm profiles for
# annotated genes: https://ftp.ncbi.nih.gov/hmm/current/) or BLAST. This seems
# a bit circular as both hmm search and blastrules are already in the RefSeq
# functional annotation pipeline from which these 16S genes originate.
# (https://ftp.ncbi.nih.gov/pub/blastrules/README)
@follows(mkdir('03_align_target_regions.dir'))
@split(summarizeAndFilterTargetRegions, '03_align_target_regions.dir/*_targets.fasta.gz')
def fetchFilteredTargetSequences(infile, outfiles):
    '''A place holder'''

    header = True
    for line in IOTools.open_file(infile):
        if header:
            header = False
            continue
        infile = line.split()[0]
        assert os.path.exists(infile), infile

        outfile = os.path.join('03_align_target_regions.dir',
                               os.path.basename(infile))

        if os.path.exists(outfile):
            os.unlink(outfile)
        shutil.copyfile(infile, outfile)
    

@transform(fetchFilteredTargetSequences,
           suffix('_targets.fasta.gz'),
           '_targets_aligned.fasta.gz')
def alignFilteredTargetSequences(infile, outfile):
    '''Align the extracted target sequences for each genome using MUSCLE'''

    statement = ("zcat %(infile)s |"
                 " muscle -quiet |"
                 " awk '/^>/ {printf(\"\\n%%s\\n\",$0);next; }"
                 "  { printf(\"%%s\",$0);}  END {printf(\"\\n\");}' |"
                 " sed '/^$/d' |"
                 " gzip > %(outfile)s")
    P.run(statement)


@transform(alignFilteredTargetSequences,
           suffix('_aligned.fasta.gz'),
           '_entropy.tsv.gz')
def calculatePerBaseEntropy(infile, outfile):
    '''Calculate the entropy for each base position in sequence alignments'''

    statement = ("python %(mbscriptsdir)s/alignment2entropy.py"
                 " --alignment-file %(infile)s"
                 " --min-occurrence 1"
                 " --out-table %(outfile)s"
                 " --log %(outfile)s.log")
    P.run(statement)


@transform(alignFilteredTargetSequences,
           suffix('_aligned.fasta.gz'),
           '_variants.tsv.gz')
def calculatePerBaseSubstitutions(infile, outfile):
    '''Calculate the entropy for each base position in sequence alignments'''

    statement = ("python %(mbscriptsdir)s/alignment2entropy.py"
                 " --alignment-file %(infile)s"
                 " --variants"
                 " --min-occurrence 1"
                 " --out-table %(outfile)s"
                 " --log %(outfile)s.log")
    P.run(statement)


@follows(calculatePerBaseEntropy, calculatePerBaseSubstitutions)
def calculateEntropy():
    pass

###############################################################################
# Summarise
###############################################################################
@merge(fetchFilteredTargetSequences, '04_feature_summary.tsv.gz')
def countFeatures(infiles, outfile):
    '''Iterate across all filtered target sequences, for each genome, 
    summarise the number of features and the number of unique features
    '''

    out = []

    for infile in infiles:
        file_name = P.snip(os.path.basename(infile), '_targets.fasta.gz')

        n = 0
        unique = set()
        for fasta in FastaIterator.FastaIterator(IOTools.open_file(infile)):
            n += 1
            unique.add(fasta.sequence)

        unique = len(unique)
        out.append([file_name, str(n), str(unique), str(unique/float(n))])

    outf = IOTools.open_file(outfile, 'w')
    outf.write('\t'.join(['genome', 'n_features', 'n_unique', 'prop_unique']) + '\n')
    outf.write('\n'.join(['\t'.join(x) for x in out]))
    outf.close()


@merge(calculatePerBaseSubstitutions, '04_feature_indel_summary.tsv.gz')
def countIndels(infiles, outfile):
    '''Iterate across summaries of subs/indels at each base position, count the
    number (and proportion of subs/indels
    '''

    out = []

    for infile in infiles:
        file_name = P.snip(os.path.basename(infile), '_targets_variants.tsv.gz')

        df = pd.read_csv(IOTools.open_file(infile), header=0, sep='\t')
        l = len(df.index)
        n = sum(df.notna()['Entropy'])
        p = n/float(l)*100
        
        out.append([file_name, str(l), str(n), str(round(p, 3))])

    outf = IOTools.open_file(outfile, 'w')
    outf.write('\t'.join(['genome', 'aligned_len', 'n_in_del_sub', 'perc_in_del_sub']) + '\n')
    outf.write('\n'.join(['\t'.join(x) for x in out]))
    outf.close()


@follows(countFeatures, countIndels)
def summarizeFeatures():
    pass
    
###############################################################################
# Generic pipeline tasks
###############################################################################
@follows(loadTaxonomy, calculateEntropy, summarizeFeatures)
def full():
    pass


def main(argv=None):
    if argv is None:
        argv = sys.argv
    P.main(argv)


if __name__ == "__main__":
    sys.exit(P.main(sys.argv))    
