## tools dependencies: fastqc, AdapterRemoval, STAR, picard, Subread, GATK.
import sys
import pathlib
import argparse
import shlex, subprocess
import glob
import logging


################## checking dependencies ##################








################################################################################################
################################################################################################
##--------------------------------------------------------------------------------------------##
##                          RNA-seq data processing functions                                 ##
##--------------------------------------------------------------------------------------------##
################################################################################################
################################################################################################
'''
Fastqc is sometimes funny in differenet platform, so we don't implement it in the pipeline.
'''


##--------------------------------------------------------------------------------------------##
## Fastqc
##--------------------------------------------------------------------------------------------##

#def fastqc(data,outdir):
##    commands = shlex.split(commandline)
#
#    process = subprocess.Popen(commands, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
#    outlog, errlog = process.communicate()
#
#    logging.info(outlog.decode("utf-8"))

##--------------------------------------------------------------------------------------------##
## Trimming sequencing adapters off raw data
##--------------------------------------------------------------------------------------------##

def trimming(data1,data2,name,outdir):
    commandline = 'AdapterRemoval --file1 ' + data1 + ' --file2 ' + data2 + ' --basename ' + TRIMDATA + '/' + name + ' --trimns --trimqualities --collapse --gzip --threads ' + num_threads
    commands = shlex.split(commandline)

    process = subprocess.Popen(commands, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    outlog, errlog = process.communicate()

    logging.info(outlog.decode("utf-8"))

##--------------------------------------------------------------------------------------------##
## Aligning trimmed data to reference genome
##--------------------------------------------------------------------------------------------##

def star(data1,data2,prefix,mapref):
    commandline = 'STAR --outFileNamePrefix ' + prefix + ' --readFilesIn ' + data1 + ' ' + data2 + ' --readFilesCommand zcat --genomeDir ' + mapref + ' --outFilterMultimapNmax 5 --outFilterMismatchNmax 8 --outSAMunmapped Within --outSAMmapqUnique 60 --outSAMstrandField intronMotif --outSAMtype BAM SortedByCoordinate --runThreadN ' + num_threads

    commands = shlex.split(commandline)

    process = subprocess.Popen(commands, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    outlog, errlog = process.communicate()

    logging.info(outlog.decode("utf-8"))

##--------------------------------------------------------------------------------------------##
## AddOrReplaceReadGroups
##--------------------------------------------------------------------------------------------##

def addreadgroups(data,outdir,name,picard):
    commandline = 'java -jar ' + picard + ' AddOrReplaceReadGroups I=' + data + ' O=' + outdir + '/' + name + '.rg.bam SO=coordinate RGID=' + name + ' RGLB=' + name + ' RGPL=Illumina RGPU=GenomeAnalyzerII RGSM=' + name

    commands = shlex.split(commandline)

    process = subprocess.Popen(commands, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    outlog, errlog = process.communicate()

    logging.info(outlog.decode("utf-8"))

##--------------------------------------------------------------------------------------------##
## MarkDuplicates
##--------------------------------------------------------------------------------------------##

def markduplicates(data,outdir,name,picard):
    commandline = 'java -jar ' + picard + ' MarkDuplicates I=' + data + ' O=' + outdir + '/' + name + '.markdup.bam M=' + outdir + '/' + name + '.markdup.metrics.txt CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT'

    commands = shlex.split(commandline)

    process = subprocess.Popen(commands, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    outlog, errlog = process.communicate()

    logging.info(outlog.decode("utf-8"))

##--------------------------------------------------------------------------------------------##
## SplitNCigarReads
##--------------------------------------------------------------------------------------------##

def splitncigar(data,ref,outdir,name,gatk):
    commandline = gatk + ' SplitNCigarReads -R ' + ref + ' -I ' + data + ' --output ' + outdir + '/' + name + '.markdup.split.bam'

    commands = shlex.split(commandline)

    process = subprocess.Popen(commands, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    outlog, errlog = process.communicate()

    logging.info(outlog.decode("utf-8"))

##--------------------------------------------------------------------------------------------##
## featureCounts
##--------------------------------------------------------------------------------------------##

def featurecount(annotation,outdir,name,data):
    commandline = 'featureCounts -t exon -g gene_id -a ' + annotation + ' -o ' + outdir + '/' + name + '.counts.txt' + data + ' -T ' + num_threads

    commands = shlex.split(commandline)

    process = subprocess.Popen(commands, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    outlog, errlog = process.communicate()

    logging.info(outlog.decode("utf-8"))

################################################################################################
################################################################################################
##--------------------------------------------------------------------------------------------##
##                               Genotyping functions                                         ##
##--------------------------------------------------------------------------------------------##
################################################################################################
################################################################################################

##--------------------------------------------------------------------------------------------##
## HaplotypeCaller
##--------------------------------------------------------------------------------------------##

def haplotypecaller(data,ref,outdir,dbsnp,name,gatk):
    commandline = gatk + ' HaplotypeCaller -R ' + ref + ' -I ' + data + ' --output ' + outdir + '/' + name + '.markdup.split.vcf.gz --dbsnp ' + dbsnp + ' --dont-use-soft-clipped-bases --standard-min-confidence-threshold-for-calling 30.0 --native-pair-hmm-threads ' + num_threads

    commands = shlex.split(commandline)

    process = subprocess.Popen(commands, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    outlog, errlog = process.communicate()

    logging.info(outlog.decode("utf-8"))

##--------------------------------------------------------------------------------------------##
## VariantFiltration
##--------------------------------------------------------------------------------------------##

def genotypefilter(data,ref,outdir,name,gatk):
    commandline = gatk + ' VariantFiltration -R ' + ref + ' -V ' + data + ' -window 35 -cluster 3 --filter-expression "FS > 30.0 || QD < 2.0"  --filter-name "RNASeqFilters_FS_QD" --filter-expression "QUAL < 30.0 || MQ < 30.0 || DP < 10.0 "  --filter-name "LowQualFilter" --output ' + outdir + '/' + name + '.markdup.split.filtered.vcf.gz'

    commands = shlex.split(commandline)

    process = subprocess.Popen(commands, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    outlog, errlog = process.communicate()

    logging.info(outlog.decode("utf-8"))

##--------------------------------------------------------------------------------------------##
## VariantEval
##--------------------------------------------------------------------------------------------##

#def varianteval(data,ref,dbsnp,outdir,name,gatk):
#    commandline = gatk + ' -T VariantEval -R ' + ref + ' -eval ' + data + ' -D ' + dbsnp + ' -noEV -EV CompOverlap -EV IndelSummary -EV TiTvVariantEvaluator -EV CountVariants -EV MultiallelicSummary -o ' + outdir + '/' + name + '.markdup.split.filtered.vcf.eval.grp'

#    commands = shlex.split(commandline)
#    p = subprocess.Popen(commands)


##--------------------------------------------------------------------------------------------##
## ASEReadCounter
##--------------------------------------------------------------------------------------------##

def allelespecificexpression(bam,vcf,ref,outdir,name,gatk):
    commandline = gatk + ' ASEReadCounter -R ' + ref + ' --output ' + outdir + '/' + name + '.ASE.csv' + ' -I ' + bam + ' --variant ' + vcf + ' -min-depth 10 --min-mapping-quality 10 --min-base-quality 2 --disable-read-filter NotDuplicateReadFilter'

    commands = shlex.split(commandline)
    process = subprocess.Popen(commands, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    outlog, errlog = process.communicate()

    logging.info(outlog.decode("utf-8"))




##--------------------------------------------------------------------------------------------##
## Argument
##--------------------------------------------------------------------------------------------##

parser = argparse.ArgumentParser(description='Pipeline for RNA-seq analysis',
                                 epilog='Ahhhhhhhhhh life')
parser.add_argument('-data', nargs='?', help='path to raw sequencing data in fastq format and gzipped. Please enter the _1.fastq.gz')
parser.add_argument('-annotation', nargs='?', help='path to annotation file in gtf format. Please download the right format from Gencode and unzip it')
parser.add_argument('-out', default='Pipeline-output', help='The output directory, by default it is ./Pipeline-output')
parser.add_argument('-ref', nargs='?', help='reference genome to use, please download from the GATK resource bundle, should including a .fasta file, a .fai file and a .dict file, please unzip them')
parser.add_argument('-mapref', nargs='?', help='alignment index, please use STAR to make the index, should be generated from the -ref file')
parser.add_argument('-dbsnp', nargs='?', help='the dbsnp file download from GATK resource bundle, please unzip')
parser.add_argument('-picard', nargs='?', help='where is the picard.jar')
parser.add_argument('-GATK', default='gatk', help='If you have executable gatk globally, levave it as default, otherwise please specify the code to execute gatk (e.g. "path/to/gatk")')
parser.add_argument('-threads', default='16', help='the number of threads to run the pipeline, default is 16')





## read arguments
args = vars(parser.parse_args())


##--------------------------------------------------------------------------------------------##
## Directories
##--------------------------------------------------------------------------------------------##

print('Start setting up the output directories.')

OUTPUT = args['out']
TRIMDATA = OUTPUT+'/Trimmed_data'
ALIGNDATA = OUTPUT+'/Mapping'
RGDATA = OUTPUT+'/AddOrReplaceReadGroups'
MDDATA = OUTPUT+'/MarkDuplicates'
SPDATA = OUTPUT+'/SplitNCigarReads'
EXPLEVEL = OUTPUT+'/featureCounts'
#FASTQC1 = OUTPUT+'/Raw_fastqc'
#FASTQC2 = TRIMDATA+'/Trim_fastqc'
GENOTYPING = OUTPUT+'/HaplotypeCaller'
VARFILTER = OUTPUT+'/VariantFiltration'
ASEREADCOUNTER = OUTPUT+'/ASEReadCounter'

Dirlist = [OUTPUT, TRIMDATA, ALIGNDATA, RGDATA, MDDATA, SPDATA, EXPLEVEL, GENOTYPING, VARFILTER, ASEREADCOUNTER]

for directory in Dirlist:
    pathlib.Path(directory).mkdir(parents=True, exist_ok=True)


##--------------------------------------------------------------------------------------------##
## Data
##--------------------------------------------------------------------------------------------##

data1 = args['data']

if data1.endswith('_1.fastq.gz'):
    name = data1.split('/')[-1].rstrip('_1.fastq.gz')
    data2 = data1.rstrip(data1.split('/')[-1]) + name + '_2.fastq.gz'

elif data1.endswith('_R1.fastq.gz'):
    name = data1.split('/')[-1].rstrip('_R1.fastq.gz')
    data2 = data1.rstrip(data1.split('/')[-1]) + name + '_R2.fastq.gz'

else:
    print('Please use the standard file name of raw data (e.g. either *_R1.fastq.gz or *_1.fastq.gz)')


print('The sample name is :', name)

annotation = args['annotation']
dbsnp = args['dbsnp']
ref = args['ref']
mapref = args['mapref']

##--------------------------------------------------------------------------------------------##
## Other settings
##--------------------------------------------------------------------------------------------##

num_threads = args['threads'] ## number of threads to run

################################################################################################
################################################################################################
##--------------------------------------------------------------------------------------------##
##                                    START RUNNING                                           ##
##--------------------------------------------------------------------------------------------##
################################################################################################
################################################################################################

logging.basicConfig(format='%(asctime)s %(message)s', filename="\"name\"-RNAseq-pipeline.log", level=logging.INFO)

##--------------------------------------------------------------------------------------------##
## Fastqc on raw data
##--------------------------------------------------------------------------------------------##
#print('Starting to QC raw data')

#fastqc(data1, FASTQC1)
#fastqc(data2, FASTQC1)

#print('Finished raw data QC, the results are in', FASTQC1)

##--------------------------------------------------------------------------------------------##
## Trimming sequencing adapters off raw data and fastqc on trimmed data
##--------------------------------------------------------------------------------------------##
print('Starting to trim raw data')

trimming(data1,data2,name,TRIMDATA)

trimdata1 = glob.glob(TRIMDATA + '/' + name + '.pair1.truncated.gz')[0]
trimdata2 = glob.glob(TRIMDATA + '/' + name + '.pair2.truncated.gz')[0]

#fastqc(trimdata1, FASTQC2)
#fastqc(trimdata2, FASTQC2)

print('Finished trimming, so far so good, the results are in', TRIMDATA)

##--------------------------------------------------------------------------------------------##
## Aligning trimmed data to reference genome
##--------------------------------------------------------------------------------------------##

print('Starting to alignment')

star(trimdata1,trimdata2,ALIGNDATA + '/' + name + '_',mapref)

aligndata = glob.glob(ALIGNDATA + '/' + name + '_Aligned.sortedByCoord.out.bam')[0]

print('Finished alignment, hopefully nothing went wrong, the results are in', ALIGNDATA)

##--------------------------------------------------------------------------------------------##
## AddOrReplaceReadGroups
##--------------------------------------------------------------------------------------------##

print('Starting to add read group to bam files')

addreadgroups(aligndata,RGDATA,name,args['picard'])

rgdata = glob.glob(RGDATA + '/' + name + '.rg.bam')[0]

print('Finished adding the read groups, it is reuired for lots of the downstream analysis, the results are in', RGDATA)

##--------------------------------------------------------------------------------------------##
## MarkDuplicates
##--------------------------------------------------------------------------------------------##

print('Starting to mark duplicates')

markduplicates(rgdata,MDDATA,name,args['picard'])

mddata = glob.glob(MDDATA + '/' + name + '.markdup.bam')[0]

print('Finished marking the duplicates, hate the duplicates, the results are in', MDDATA)

##--------------------------------------------------------------------------------------------##
## SplitNCigarReads
##--------------------------------------------------------------------------------------------##

print('Starting to splitncigar')

splitncigar(mddata,ref,SPDATA,name,args['GATK'])

spdata = glob.glob(SPDATA + '/' + name + '.markdup.split.bam')[0]

print('Finished splitncigar, it is recommended in the best practise tho, the results are in', SPDATA)

##--------------------------------------------------------------------------------------------##
## featureCounts
##--------------------------------------------------------------------------------------------##
print('Starting to calculate expression level')

featurecount(annotation,EXPLEVEL,name,spdata)

print('Finished calculate the expression level, featureCounts is better than htseq-count, the results are in', EXPLEVEL)

################################################################################################
################################################################################################
##--------------------------------------------------------------------------------------------##
##                                 START Genotyping                                           ##
##--------------------------------------------------------------------------------------------##
################################################################################################
################################################################################################

##--------------------------------------------------------------------------------------------##
## HaplotypeCaller
##--------------------------------------------------------------------------------------------##
print('Starting to call genotype')

haplotypecaller(spdata,ref,GENOTYPING,dbsnp,name,args['GATK'])

genotype = glob.glob(GENOTYPING + '/' + name + '.markdup.split.vcf.gz')[0]

print('Finished calling the genotype from RNAseq data, well done, the results are in', GENOTYPING)

##--------------------------------------------------------------------------------------------##
## VariantFiltration
##--------------------------------------------------------------------------------------------##

print('Starting to filter genotype')

genotypefilter(genotype,ref,VARFILTER,name,args['GATK'])

filteredgenotype = glob.glob(VARFILTER + '/' + name + '.markdup.split.filtered.vcf.gz')[0]

print('Finished setting filters in the vcf, just cannot trust the data without filters, the results are in', VARFILTER)

##--------------------------------------------------------------------------------------------##
## VariantEval
##--------------------------------------------------------------------------------------------##

#varianteval(filteredgenotype,ref,dbsnp,VAREVAL,name,args['GATK'])


##--------------------------------------------------------------------------------------------##
## ASEReadCounter
##--------------------------------------------------------------------------------------------##
print('Starting to do ASE analysis')

allelespecificexpression(mddata,filteredgenotype,ref,ASEREADCOUNTER,name,args['GATK'])

print('Finally! Finished the allele specific analysis, farewell bro, the results are in', ASEREADCOUNTER)
