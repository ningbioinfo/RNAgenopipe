# This is the script to do RNA-seq analysis and Genotyping

## This is a preparation pipeline for RNA-seq expression data and genotyping data before eQTL mapping

## The pipeline is specifically designed to run in [Phoenix](https://www.adelaide.edu.au/phoenix/)

### To execute the pipeline

1. Login to Phoenix

2. Go to your data directory

3. do

  ```
  for file in ./*_1.fastq.gz; do sbatch /path/to/the/scripts/RNAseq_Genotyping_pipeline.sh; done
  ```

**Change the suffix of the raw input data if needed**
