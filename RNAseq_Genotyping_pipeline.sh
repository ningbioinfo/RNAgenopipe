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
module load R/3.4.2-foss-2016b
module load GATK
module load Python/3.6.1-foss-2016b


##--------------------------------------------------------------------------------------------##
## Make directories
##--------------------------------------------------------------------------------------------##

## required directories and files
MAPREF=/fast/users/a1692215/RNA-seq_Montgomery/Ref/STAR_index #mapping index



## Dirctories to make
OUTPUT=/fast/users/a1692215/RNA-seq_Montgomery/Output
TRIMDATA=${OUTPUT}/Trimmed_data
ALIGNDATA=${OUTPUT}/Mapping
RGDATA=${OUTPUT}/AddOrReplaceReadGroups
MDDATA=${OUTPUT}/MarkDuplicates
SPDATA=${OUTPUT}/SplitNCigarReads
EXPLEVEL=${OUTPUT}/after_mapping/featureCounts
RNA_METRICS=${OUTPUT}/CollectRnaSeqMetrics
IS_METRICS=${OUTPUT}/CollectInsertSizeMetrics
FASTQC1=${OUTPUT}/Raw_fastqc
FASTQC2=${TRIMDATA}/Trim_fastqc
GENOTYPING=${OUTPUT}/HaplotypeCaller
VARFILTER=${OUTPUT}/VariantFiltration
VAREVAL=${OUTPUT}/VariantEval
ASEREADCOUNTER=${OUTPUT}/ASEReadCounter

## make array of directories to make
declare -a arr=(${OUTPUT} ${TRIMDATA} ${ALIGNDATA} ${RGDATA} ${MDDATA} ${SPDATA} \
                ${EXPLEVEL} ${RNA_METRICS} ${IS_METRICS} ${FASTQC1} ${FASTQC2}\
                ${GENOTYPING} ${VARFILTER} ${VAREVAL} ${ASEREADCOUNTER})

## Additional python scripts
if [ -e ${ROOT}/Directories_security.py ]
then
    DS=${ROOT}/Directories_security.py
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
REF=/fast/users/a1692215/RNA-seq_Montgomery/Ref/Homo_sapiens_assembly38.fasta
Annotation_gtf=fast/users/a1692215/RNA-seq_Montgomery/Ref/hg38_annotation.gtf
xenoRefFlat=fast/users/a1692215/RNA-seq_Montgomery/Ref/hg38_xenoRefFlat.txt
DBSNP=/fast/users/a1692215/RNA-seq_Montgomery/Ref/Homo_sapiens_assembly38.dbsnp138.vcf
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


fastqc -t 16 -f fastq -o ${FASTQC1} ${DATA}
fastqc -t 16 -f fastq -o ${FASTQC1} ${DATA/_1/_2}


##--------------------------------------------------------------------------------------------##
## Trimming sequencing adapters off raw data
##--------------------------------------------------------------------------------------------##

AdapterRemoval --file1 ${DATA} --file2 ${DATA/_1/_2} \
    --basename ${TRIMDATA}/${NAME} \
    --trimns \
    --trimqualities \
    --collapse \
    --threads 16 \
    #--qualitymax 93 \ #sometimes is needed, if found 0b sample in ${TRIMDATA}
    --gzip

## QC again


fastqc -t 16 -f fastq -o ${FASTQC2} ${DATA}
fastqc -t 16 -f fastq -o ${FASTQC2} ${DATA/_1/_2}



##--------------------------------------------------------------------------------------------##
## Aligning trimmed data to reference genome
##--------------------------------------------------------------------------------------------##


STAR --outFileNamePrefix ${ALIGNDATA}/${NAME} \
		--readFilesIn ${TRIMDATA}/{NAME}.pair1.truncated.gz ${TRIMDATA}/{NAME}.pair2.truncated.gz \
		--readFilesCommand zcat \
		--genomeDir ${MAPREF} \
		--runThreadN 16 \
		--outFilterMultimapNmax 5 \
		--outFilterMismatchNmax 8 \
		--outSAMunmapped Within \
		--outSAMstrandField intronMotif \
		--outSAMtype BAM SortedByCoordinate



##--------------------------------------------------------------------------------------------##
## AddOrReplaceReadGroups
##--------------------------------------------------------------------------------------------##

java -jar ${PICARD} AddOrReplaceReadGroups \
   I=${ALIGNDATA}/${NAME}.Aligned.sortedByCoord.out.bam \
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

java -jar ${GATK} -T SplitNCigarReads \
   -R ${REF} -I ${MDDATA}/${NAME}.markdup.bam \
   -o ${SPDATA}/${NAME}.markdup.split.bam -rf ReassignOneMappingQuality \
   -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS


##--------------------------------------------------------------------------------------------##
## featureCounts
##--------------------------------------------------------------------------------------------##

featureCounts -T 16 -t exon -g gene_id \
    -a ${Annotation_gtf} \
    -o ${EXPLEVEL}/${NAME}.counts.txt \
    ${SPDATA}/${NAME}.markdup.split.bam


##--------------------------------------------------------------------------------------------##
## CollectRnaSeqMetrics
##--------------------------------------------------------------------------------------------##

java -jar ${PICARD} CollectRnaSeqMetrics \
		I=${SPDATA}/${NAME}.markdup.split.bam \
		O=${RNA_METRICS}/${NAME}.RNA_Metrics \
		REF_FLAT=${xenoRefFlat} \
		STRAND=NONE


##--------------------------------------------------------------------------------------------##
## CollectInsertSizeMetrics
##--------------------------------------------------------------------------------------------##

java -jar ${PICARD} CollectInsertSizeMetrics \
		I=${SPDATA}/${NAME}.markdup.split.bam \
		O=${IS_METRICS}/${NAME}_isizes.txt \
		H=${IS_METRICS}/${NAME}_isizes_histogram.pdf
	done

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

java -jar ${GATK} -T HaplotypeCaller \
   -R ${REF} -I ${SPDATA}/${NAME}.markdup.split.bam \
   -o ${GENOTYPING}/${NAME}.markdup.split.vcf.gz --dbsnp ${DBSNP} -dontUseSoftClippedBases \
   -stand_call_conf 30 -nct 16


##--------------------------------------------------------------------------------------------##
## VariantFiltration
##--------------------------------------------------------------------------------------------##

java -jar ${GATK} -T VariantFiltration -R ${REF} \
   -V ${GENOTYPING}/${NAME}.markdup.split.vcf.gz \
   -window 35 -cluster 3 \
   --filterExpression "FS > 30.0 || QD < 2.0"  -filterName "RNASeqFilters_FS_QD" \
   --filterExpression "QUAL < 30.0 || MQ < 30.0 || DP < 10.0 "  -filterName "LowQualFilter" \
   -o ${VARFILTER}/${NAME}.markdup.split.filtered.vcf.gz


##--------------------------------------------------------------------------------------------##
## VariantEval
##--------------------------------------------------------------------------------------------##

java -jar ${GATK} -T VariantEval -R ${REF} \
   -eval ${VARFILTER}/${NAME}.markdup.split.filtered.vcf.gz -D ${DBSNP} \
   -noEV -EV CompOverlap -EV IndelSummary -EV TiTvVariantEvaluator -EV CountVariants \
   -EV MultiallelicSummary -o ${VAREVAL}/${NAME}.markdup.split.filtered.vcf.eval.grp

##--------------------------------------------------------------------------------------------##
## ASEReadCounter
##--------------------------------------------------------------------------------------------##

java -jar ${GATK} -T ASEReadCounter -R ${REF} \
   -o ${ASEREADCOUNTER}/${NAME}.ASE.csv \
   -I ${MDDATA}/${NAME}.markdup.bam \
   -sites ${VARFILTER}/${NAME}.markdup.split.filtered.vcf.gz -U ALLOW_N_CIGAR_READS \
   -minDepth 10 --minMappingQuality 10 --minBaseQuality 2 -drf DuplicateRead
