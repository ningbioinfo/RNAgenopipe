#!/bin/bash
#SBATCH -p batch
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --time=36:00:00
#SBATCH --mem=32GB
#SBATCH -o /fast/users/a1692215/RNA-seq_Montgomery/out.log
#SBATCH -e /fast/users/a1692215/RNA-seq_Montgomery/err.log
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=ning.liu@student.adelaide.edu.au

## Preparation of the expression data (RNAseq)

## Modules
module load fastqc
module load AdapterRemoval
module load STAR
module load picard/2.2.4-Java-1.8.0_71
module load Subread
module load R/3.4.2-foss-2016b
#module load GATK
module load Python/3.6.1-foss-2016b

## goto working directory
cd /fast/users/a1692215/RNA-seq_Montgomery/

## Directories
DATA=/fast/users/a1692215/RNA-seq_Montgomery/raw_data/
REFs=/fast/users/a1692215/RNA-seq_Montgomery/Ref/
MAPREF=${REFs}/STAR_index/
OUTPUT=/fast/users/a1692215/RNA-seq_Montgomery/OUTPUT/
ROOT=$(pwd)


TRIMDATA=${OUTPUT}/trimmed_data
ALIGNDATA=${OUTPUT}/mapping
MDDATA=${OUTPUT}/after_mapping/mark_dup
EXPLEVEL=${OUTPUT}/after_mapping/featureCounts
RNA_METRICS=${OUTPUT}/after_mapping/CollectRnaSeqMetrics
IS_METRICS=${OUTPUT}/after_mapping/CollectInsertSizeMetrics
SPLITNCIGAR=${OUTPUT}/after_mapping/SplitNCigarReads
ADDREADG=${OUTPUT}/after_mapping/AddOrReplaceReadGroups
FASTQC1=${DATA}/fastqc
FASTQC2=${TRIMDATA}/fastqc

TOMAKEDI=(${OUTPUT} ${TRIMDATA} ${ALIGNDATA} ${MDDATA} ${EXPLEVEL} ${RNA_METRICS} ${IS_METRICS} ${SPLITNCIGAR} ${ADDREADG} ${FASTQC1} ${FASTQC2})

## Additional python scripts
if [ -e ${ROOT}/Directories_security.py ]
then
    DS=${ROOT}/Directories_security.py
else
    echo "There is no directories_security scripts in the location."
fi


## making directory if not exist
for DIRECTORY in ${TOMAKEDI}
  do
    python3 ${DS} ${DIRECTORY}
  done



## Data files
#GATK_ref=${REFs}/Homo_sapiens_assembly38.fasta
Annotation_gtf=${REFs}/hg38_annotation.gtf
xenoRefFlat=${REFs}/hg38_xenoRefFlat.txt


##--------------------------------------------------------------------------------------------##
## Fastqc on raw data
##--------------------------------------------------------------------------------------------##

for FQGZ in ${DATA}/*.fastq.gz
  do
    fastqc -t 16 -f fastq -o ${FASTQC1} ${FQGZ}
  done

##--------------------------------------------------------------------------------------------##
## Trimming sequencing adapters off raw data
##--------------------------------------------------------------------------------------------##
cd ${TRIMDATA}

for FQGZ in ${DATA}/*_1.fastq.gz
	do
		Basename=$(basename ${FQGZ} _1.fastq.gz)

		AdapterRemoval --file1 ${FQGZ} --file2 ${FQGZ/_1/_2} \
    --basename ${Basename} \
    --trimns \
    --trimqualities \
    --collapse \
    --threads 16
    --gzip
	done

## QC again

for FQGZ in ${TRIMDATA}/*.fastq.gz
  do
    fastqc -t 16 -f fastq -o ${FASTQC2} ${FQGZ}
  done


##--------------------------------------------------------------------------------------------##
## Aligning trimmed data to reference genome
##--------------------------------------------------------------------------------------------##

cd ${ALIGNDATA}

for FQGZ in ${TRIMDATA}/*pair1.truncated.fastq.gz
	do
		Basename=$(basename ${FQGZ} .pair1.truncated.fastq.gz)

		STAR --outFileNamePrefix ${OUTPUT}/${Basename} \
		--readFilesIn ${FQGZ} ${FQGZ/pair1/pair2} \
		--readFilesCommand zcat \
		--genomeDir ${MAPREF} \
		--runThreadN 16 \
		--outFilterMultimapNmax 5 \
		--outFilterMismatchNmax 8 \
		--outSAMunmapped Within \
		--outSAMstrandField intronMotif \
		--outSAMtype BAM SortedByCoordinate
	done


##--------------------------------------------------------------------------------------------##
## MarkDuplicates
##--------------------------------------------------------------------------------------------##

cd ${MDDATA}

for BAM in ${ALIGNDATA}/*Aligned.sortedByCoord.out.bam
	do
		Basename=$(basename ${BAM} Aligned.sortedByCoord.out.bam)

		java -jar /apps/software/picard/2.2.4-Java-1.8.0_71/picard.jar MarkDuplicates \
		I=${BAM} \
		O=${MDDATA}/${Basename}_markdup.bam \
		M=${MDDATA}/${Basename}_md_metrics.txt \
		AS=TRUE
	done


##--------------------------------------------------------------------------------------------##
## featureCounts
##--------------------------------------------------------------------------------------------##
cd ${EXPLEVEL}

for BAM in ${MDDATA}/*_markdup.bam
	do
		Basename=$(basename ${BAM} _markdup.bam)

		featureCounts -T 16 -t exon -g gene_id \
    -a ${Annotation_gtf} \
    -o ${Basename}.counts.txt \
    ${BAM}
	done


##--------------------------------------------------------------------------------------------##
## CollectRnaSeqMetrics
##--------------------------------------------------------------------------------------------##

cd ${RNA_METRICS}

for BAM in ${MDDATA}/*_markdup.bam
	do
		Basename=$(basename ${BAM} _markdup.bam)

		java -jar /apps/software/picard/2.2.4-Java-1.8.0_71/picard.jar CollectRnaSeqMetrics \
		I=${BAM} \
		O=${RNA_METRICS}/${Basename}.RNA_Metrics \
		REF_FLAT=${xenoRefFlat} \
		STRAND=NONE
	done


##--------------------------------------------------------------------------------------------##
## CollectInsertSizeMetrics
##--------------------------------------------------------------------------------------------##

cd ${IS_METRICS}


for BAM in ${MDDATA}/*_markdup.bam
	do
		Basename=$(basename ${BAM} _markdup.bam)

		java -jar /apps/software/picard/2.2.4-Java-1.8.0_71/picard.jar CollectInsertSizeMetrics \
		I=${BAM} \
		O=${IS_METRICS}/${Basename}_isizes.txt \
		H=${IS_METRICS}/${Basename}_isizes_histogram.pdf
	done

##--------------------------------------------------------------------------------------------##
## AddOrReplaceReadGroups
##--------------------------------------------------------------------------------------------##

##--------------------------------------------------------------------------------------------##
## SplitNCigarReads
##--------------------------------------------------------------------------------------------##
