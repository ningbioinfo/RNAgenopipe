#!/bin/bash
#SBATCH -p batch
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --time=24:00:00
#SBATCH --mem=100GB
#SBATCH -o /fast/users/a1692215/RNA-seq_Montgomery/out.log
#SBATCH -e /fast/users/a1692215/RNA-seq_Montgomery/err.log
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=ning.liu@student.adelaide.edu.au

## script for RNAseq & Genotyping

## Modules
module load fastqc
module load AdapterRemoval
module load STAR
module load picard/2.2.4-Java-1.8.0_71
module load Subread
module load GATK
module load Python/3.6.1-foss-2016b


##--------------------------------------------------------------------------------------------##
## Make directories
##--------------------------------------------------------------------------------------------##

## required directories and files
REF=/fast/users/a1692215/RNA-eqtl/Ref/human_g1k_v37.fasta
MAPREF=/fast/users/a1692215/RNA-eqtl/Ref/STAR_index #mapping index
SCRIPTS=/fast/users/a1692215/RNA-eqtl/


## Dirctories to make
OUTPUT=/fast/users/a1692215/RNA-eqtl/Pipeline-output
TRIMDATA=${OUTPUT}/Trimmed_data
ALIGNDATA=${OUTPUT}/Mapping
RGDATA=${OUTPUT}/AddOrReplaceReadGroups
MDDATA=${OUTPUT}/MarkDuplicates
SPDATA=${OUTPUT}/SplitNCigarReads
EXPLEVEL=${OUTPUT}/featureCounts
#FASTQC1=${OUTPUT}/Raw_fastqc
FASTQC2=${TRIMDATA}/Trim_fastqc
GENOTYPING=${OUTPUT}/HaplotypeCaller
VARFILTER=${OUTPUT}/VariantFiltration
VAREVAL=${OUTPUT}/VariantEval
ASEREADCOUNTER=${OUTPUT}/ASEReadCounter

## make array of directories to make
declare -a arr=(${OUTPUT} ${TRIMDATA} ${ALIGNDATA} ${RGDATA} ${MDDATA} ${SPDATA} \
                ${EXPLEVEL} ${FASTQC2}\
                ${GENOTYPING} ${VARFILTER} ${VAREVAL} ${ASEREADCOUNTER})

## Additional python scripts
if [ -e ${SCRIPTS}/Directories_security.py ]
then
    DS=${SCRIPTS}/Directories_security.py
else
    echo "There is no directories_security scripts in the location."
fi


## making directory if not exist
for DIRECTORY in "${arr[@]}"
  do
    python3 ${DS} ${DIRECTORY}
  done



##--------------------------------------------------------------------------------------------##
## DATA
##--------------------------------------------------------------------------------------------##

DATA=$1 #Run each sample in parallel, *** use _1 as input ***
Annotation_gtf=fast/users/a1692215/RNA-eqtl/Ref/gencode.v27lift37.annotation.gtf
DBSNP=/fast/users/a1692215/RNA-eqtl/Ref/dbsnp_138.b37.vcf
NAME=$(basename ${DATA} _1.fastq.gz) # change the suffix if needed


##--------------------------------------------------------------------------------------------##
## Program
##--------------------------------------------------------------------------------------------##

PICARD=$EBROOTPICARD/picard.jar
GATK=$EBROOTGATK/GenomeAnalysisTK.jar

################################################################################################
################################################################################################
##--------------------------------------------------------------------------------------------##
##                                   RNA-seq data processing                                  ##
##--------------------------------------------------------------------------------------------##
################################################################################################
################################################################################################

##--------------------------------------------------------------------------------------------##
## Fastqc on raw data
##--------------------------------------------------------------------------------------------##


#fastqc -t 16 -f fastq -o ${FASTQC1} ${DATA}
#fastqc -t 16 -f fastq -o ${FASTQC1} ${DATA/_1/_2}


##--------------------------------------------------------------------------------------------##
## Trimming sequencing adapters off raw data
##--------------------------------------------------------------------------------------------##

AdapterRemoval --file1 ${DATA} --file2 ${DATA/_1/_2} \
    --basename ${TRIMDATA}/${NAME} \
    --trimns \
    --trimqualities \
    --collapse \
    --threads 16 \
    --gzip

## QC again


#fastqc -t 16 -f fastq -o ${FASTQC2} ${TRIMDATA}/${NAME}.pair1.truncated.gz
#fastqc -t 16 -f fastq -o ${FASTQC2} ${TRIMDATA}/${NAME}.pair2.truncated.gz



##--------------------------------------------------------------------------------------------##
## Aligning trimmed data to reference genome
##--------------------------------------------------------------------------------------------##


STAR --outFileNamePrefix ${ALIGNDATA}/${NAME}_ \
		--readFilesIn ${TRIMDATA}/${NAME}.pair1.truncated.gz ${TRIMDATA}/${NAME}.pair2.truncated.gz \
		--readFilesCommand zcat \
		--genomeDir ${MAPREF} \
		--runThreadN 16 \
		--outFilterMultimapNmax 5 \
		--outFilterMismatchNmax 8 \
		--outSAMunmapped Within \
    --outSAMmapqUnique 60 \
		--outSAMstrandField intronMotif \
		--outSAMtype BAM SortedByCoordinate



##--------------------------------------------------------------------------------------------##
## AddOrReplaceReadGroups
##--------------------------------------------------------------------------------------------##

java -jar ${PICARD} AddOrReplaceReadGroups \
   I=${ALIGNDATA}/${NAME}_Aligned.sortedByCoord.out.bam \
   O=${RGDATA}/${NAME}.rg.bam \
   SO=coordinate RGID=${NAME} RGLB=${NAME} RGPL=Illumina RGPU=GenomeAnalyzerII RGSM=${NAME}

##--------------------------------------------------------------------------------------------##
## MarkDuplicates
##--------------------------------------------------------------------------------------------##

java -jar ${PICARD} MarkDuplicates \
   I=${RGDATA}/${NAME}.rg.bam \
   O=${MDDATA}/${NAME}.markdup.bam \
   M=${MDDATA}/${NAME}.markdup.metrics.txt CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT

##--------------------------------------------------------------------------------------------##
## SplitNCigarReads
##--------------------------------------------------------------------------------------------##

gatk SplitNCigarReads \
   -R ${REF} -I ${MDDATA}/${NAME}.markdup.bam \
   --output ${SPDATA}/${NAME}.markdup.split.bam

##--------------------------------------------------------------------------------------------##
## featureCounts
##--------------------------------------------------------------------------------------------##

featureCounts -T 16 -t exon -g gene_id \
    -a ${Annotation_gtf} \
    -o ${EXPLEVEL}/${NAME}.counts.txt \
    ${SPDATA}/${NAME}.markdup.split.bam


################################################################################################
################################################################################################
##--------------------------------------------------------------------------------------------##
##                                       Genotyping                                           ##
##--------------------------------------------------------------------------------------------##
################################################################################################
################################################################################################


##--------------------------------------------------------------------------------------------##
## HaplotypeCaller
##--------------------------------------------------------------------------------------------##

gatk HaplotypeCaller \
   -R ${REF} -I ${SPDATA}/${NAME}.markdup.split.bam \
   --output ${GENOTYPING}/${NAME}.markdup.split.vcf.gz --dbsnp ${DBSNP} --dont-use-soft-clipped-bases \
   --standard-min-confidence-threshold-for-calling 30.0 -nct 16


##--------------------------------------------------------------------------------------------##
## VariantFiltration
##--------------------------------------------------------------------------------------------##

gatk VariantFiltration -R ${REF} \
   -V ${GENOTYPING}/${NAME}.markdup.split.vcf.gz \
   -window 35 -cluster 3 \
   --filter-expression "FS > 30.0 || QD < 2.0"  --filter-name "RNASeqFilters_FS_QD" \
   --filter-expression "QUAL < 30.0 || MQ < 30.0 || DP < 10.0 "  --filter-name "LowQualFilter" \
   --output ${VARFILTER}/${NAME}.markdup.split.filtered.vcf.gz


##--------------------------------------------------------------------------------------------##
## VariantEval
##--------------------------------------------------------------------------------------------##

#gatk VariantEval -R ${REF} \
#   -eval ${VARFILTER}/${NAME}.markdup.split.filtered.vcf.gz -D ${DBSNP} \
#   -noEV -EV CompOverlap -EV IndelSummary -EV TiTvVariantEvaluator -EV CountVariants \
#   -EV MultiallelicSummary --output ${VAREVAL}/${NAME}.markdup.split.filtered.vcf.eval.grp

##--------------------------------------------------------------------------------------------##
## ASEReadCounter
##--------------------------------------------------------------------------------------------##

gatk ASEReadCounter -R ${REF} \
   --output ${ASEREADCOUNTER}/${NAME}.ASE.csv \
   -I ${MDDATA}/${NAME}.markdup.bam \
   --variant ${VARFILTER}/${NAME}.markdup.split.filtered.vcf.gz \
   -min-depth 10 --min-mapping-quality 10 --min-base-quality 2 --disable-read-filter NotDuplicateReadFilter
