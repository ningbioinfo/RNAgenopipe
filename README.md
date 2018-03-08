# This is a pipeline to do RNA-seq analysis and impute genotype from RNA-seq data

## The script is required python (version > 3) (tested with 3.6.1)

## Tools dependencies

1. AdapterRemoval

  - <https://github.com/MikkelSchubert/adapterremoval>

  To install AdapterRemoval: `conda install -c maxibor adapterremoval2` (if you have conda), otherwise go the the github page and follow the instruction.

2. STAR

  - <https://github.com/alexdobin/STAR>

  To install STAR: **Linux**

  ```
   git clone https://github.com/alexdobin/STAR.git
   cd STAR/source

   # Build STAR
   make STAR
  ```

  **Mac**

  ```
   git clone https://github.com/alexdobin/STAR.git
   cd STAR/source

   # Build STAR
   make STARforMacStatic
  ```

1. picard

  - <https://broadinstitute.github.io/picard/>

  To install picard:

  - Go to <https://github.com/broadinstitute/picard/releases/tag/2.17.11>
  - Download the .jar file
  - make sure you have Java, to test: `java -version`, otherwise go <http://www.oracle.com/technetwork/java/javase/downloads/index.html> to download Java
  - test the picard tool by `java -jar path/to/picard.jar -h`, you should see the manual of picard

1. GATK (>4.0)

  - <https://software.broadinstitute.org/gatk/download/>

  To install GATK:

  ```
  wget https://github.com/broadinstitute/gatk/releases/download/4.0.2.1/gatk-4.0.2.1.zip
  unzip gatk-4.0.2.1.zip
  ```

Tools are preferred to be the lastest version.

## If you are working in [Phoenix](https://www.adelaide.edu.au/phoenix/) the HPC

Just

```
module load AdapterRemoval
module load STAR
module load picard
module load GATK
module load Python/3.6.1-foss-2016b
```

# How to use the script

1. This script including two parts:

  1. RNA-seq analysis, including:

    - Trimming
    - Alignment
    - MarkDuplicates
    - Calculate expression level

  2. Genotype stage, including:

    - Impute genotype from RNA-seq
    - Filtering
    - Allele specific expression analysis **To disable the second stage, specific the argument '-genotypemode'**

2. To see the manual and run the script properly, see `python RNA-seq-pipeline.py -h`
