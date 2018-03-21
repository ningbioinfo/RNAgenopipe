## Author: Ning Liu
## Version: 0.1.0 beta

import sys
import pathlib
import argparse
import shlex, subprocess
import glob
import logging


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
    commandline = 'featureCounts -t exon -g gene_id -a ' + annotation + ' -o ' + outdir + '/' + name + '.counts.txt ' + data + ' -T ' + num_threads

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
## SelectVariant
##--------------------------------------------------------------------------------------------##

def selctvaiant(data,ref,outdir,name,gatk):
    commandline = gatk + ' SelectVariants -R ' + ref + ' -V ' + data + ' --output ' + outdir + '/' + name + '.markdup.split.filtered.biallelic.vcf.gz --select-type-to-include SNP --restrict-alleles-to BIALLELIC'

    commands = shlex.split(commandline)

    process = subprocess.Popen(commands, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    outlog, errlog = process.communicate()

    logging.info(outlog.decode("utf-8"))


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
parser.add_argument('-genotypemode', default=1, help='by default, the genotype mode is on, meaning the script will require "-dbsnp" option to be specified, if you want to turn off this mode and just use this script to do simple RNA-seq analysis, speicify this option to anything but 1')
parser.add_argument('-dbsnp', nargs='?', help='the dbsnp file download from GATK resource bundle, please unzip')
parser.add_argument('-picard', nargs='?', help='where is the picard.jar')
parser.add_argument('-GATK', default='gatk', help='If you have executable gatk globally, levave it as default, otherwise please specify the code to execute gatk (e.g. "path/to/gatk")')
parser.add_argument('-threads', default='16', help='the number of threads to run the pipeline, default is 16')





## read arguments
args = vars(parser.parse_args())


##--------------------------------------------------------------------------------------------##
## Directories
##--------------------------------------------------------------------------------------------##

print('Start setting up the output directories.','\n')

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
    name = data1.split('/')[-1].rstrip('1.fastq.gz').rstrip('_')
    data2 = data1.rstrip(data1.split('/')[-1]) + name + '_2.fastq.gz'

elif data1.endswith('_R1.fastq.gz'):
    name = data1.split('/')[-1].rstrip('R1.fastq.gz').rstrip('_')
    data2 = data1.rstrip(data1.split('/')[-1]) + name + '_R2.fastq.gz'

else:
    print('Please use the standard file name of raw data (e.g. either *_R1.fastq.gz or *_1.fastq.gz)','\n')


print('The sample name is :', name, '\n')

annotation = args['annotation']
dbsnp = args['dbsnp']
ref = args['ref']
mapref = args['mapref']

##--------------------------------------------------------------------------------------------##
## Other settings
##--------------------------------------------------------------------------------------------##

num_threads = args['threads'] ## number of threads to run



################## checking dependencies ##################

print('Checking dependencies....','\n')

try:
    subprocess.Popen(['AdapterRemoval','-version'],stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
except:
    print('AdapterRemoval has not been installed.')
    sys.exit(1)
try:
    subprocess.Popen(['STAR','-h','|','grep','^versionSTAR'],stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
except:
    print('STAR has not been installed.')
    sys.exit(1)
try:
    subprocess.Popen(['java','-jar',args['picard'],'-h'],stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
except:
    print('picard has not been detected.')
    sys.exit(1)
try:
    subprocess.Popen([args['GATK'],'-h'],stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
except:
    print('GATK has not been installed.')
    sys.exit(1)

print('Finished checking dependencies.','\n')

################################################################################################
################################################################################################
##--------------------------------------------------------------------------------------------##
##                                    START RUNNING                                           ##
##--------------------------------------------------------------------------------------------##
################################################################################################
################################################################################################

logging.basicConfig(format='%(asctime)s %(message)s', filename='%s-RNAseq-pipeline.log' % name, level=logging.INFO)

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
print('Starting to trim raw data','\n')

trimming(data1,data2,name,TRIMDATA)

trimdata1 = glob.glob(TRIMDATA + '/' + name + '.pair1.truncated.gz')[0]
trimdata2 = glob.glob(TRIMDATA + '/' + name + '.pair2.truncated.gz')[0]

#fastqc(trimdata1, FASTQC2)
#fastqc(trimdata2, FASTQC2)

print('Finished trimming, so far so good, the results are in', TRIMDATA,'\n')

##--------------------------------------------------------------------------------------------##
## Aligning trimmed data to reference genome
##--------------------------------------------------------------------------------------------##

print('Starting to alignment','\n')

star(trimdata1,trimdata2,ALIGNDATA + '/' + name + '_',mapref)

aligndata = glob.glob(ALIGNDATA + '/' + name + '_Aligned.sortedByCoord.out.bam')[0]

print('Finished alignment, hopefully nothing went wrong, the results are in', ALIGNDATA,'\n')

##--------------------------------------------------------------------------------------------##
## AddOrReplaceReadGroups
##--------------------------------------------------------------------------------------------##

print('Starting to add read group to bam files','\n')

addreadgroups(aligndata,RGDATA,name,args['picard'])

rgdata = glob.glob(RGDATA + '/' + name + '.rg.bam')[0]

print('Finished adding the read groups, it is reuired for lots of the downstream analysis, the results are in', RGDATA,'\n')

##--------------------------------------------------------------------------------------------##
## MarkDuplicates
##--------------------------------------------------------------------------------------------##

print('Starting to mark duplicates')

markduplicates(rgdata,MDDATA,name,args['picard'])

mddata = glob.glob(MDDATA + '/' + name + '.markdup.bam')[0]

print('Finished marking the duplicates, hate the duplicates, the results are in', MDDATA,'\n')

##--------------------------------------------------------------------------------------------##
## SplitNCigarReads
##--------------------------------------------------------------------------------------------##

print('Starting to splitncigar','\n')

splitncigar(mddata,ref,SPDATA,name,args['GATK'])

spdata = glob.glob(SPDATA + '/' + name + '.markdup.split.bam')[0]

print('Finished splitncigar, it is recommended in the best practise tho, the results are in', SPDATA,'\n')

##--------------------------------------------------------------------------------------------##
## featureCounts
##--------------------------------------------------------------------------------------------##
print('Starting to calculate expression level','\n')

featurecount(annotation,EXPLEVEL,name,spdata)

print('Finished calculate the expression level, featureCounts is better than htseq-count, the results are in', EXPLEVEL,'\n')

################################################################################################
################################################################################################
##--------------------------------------------------------------------------------------------##
##                                 START Genotyping                                           ##
##--------------------------------------------------------------------------------------------##
################################################################################################
################################################################################################

if args['genotypemode']==1:
    print('OK now you choose the genotype mode, it will start to do the genotyping steps.','\n')


##--------------------------------------------------------------------------------------------##
## HaplotypeCaller
##--------------------------------------------------------------------------------------------##
    print('Starting to call genotype','\n')

    haplotypecaller(spdata,ref,GENOTYPING,dbsnp,name,args['GATK'])

    genotype = glob.glob(GENOTYPING + '/' + name + '.markdup.split.vcf.gz')[0]

    print('Finished calling the genotype from RNAseq data, well done, the results are in', GENOTYPING,'\n')

##--------------------------------------------------------------------------------------------##
## VariantFiltration
##--------------------------------------------------------------------------------------------##

    print('Starting to filter genotype','\n')

    genotypefilter(genotype,ref,VARFILTER,name,args['GATK'])

    filteredgenotype = glob.glob(VARFILTER + '/' + name + '.markdup.split.filtered.vcf.gz')[0]

    print('Finished setting filters in the vcf, just cannot trust the data without filters, the results are in', VARFILTER,'\n')

##--------------------------------------------------------------------------------------------##
## VariantEval
##--------------------------------------------------------------------------------------------##

#varianteval(filteredgenotype,ref,dbsnp,VAREVAL,name,args['GATK'])

##--------------------------------------------------------------------------------------------##
## SelectVariant
##--------------------------------------------------------------------------------------------##

    print('Starting to restrict to bialleic','\n')

    selctvaiant(filteredgenotype,ref,VARFILTER,name,args['GATK'])

    bialleicvcf = glob.glob(VARFILTER + '/' + name + '.markdup.split.filtered.biallelic.vcf.gz')[0]

    print('Finsished restrict variant to bialleic because it is required for ASE analysis, the results are in', VARFILTER,'\n')


##--------------------------------------------------------------------------------------------##
## ASEReadCounter
##--------------------------------------------------------------------------------------------##
    print('Starting to do ASE analysis','\n')

    allelespecificexpression(mddata,bialleicvcf,ref,ASEREADCOUNTER,name,args['GATK'])

    print('Finally! Finished the allele specific analysis, farewell bro, the results are in', ASEREADCOUNTER,'\n')

else:
    print('Since you did not use the genotype mode, your analysis is finished now!','\n')
