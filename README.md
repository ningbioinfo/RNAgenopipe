# This is a pipeline to do RNA-seq analysis and impute genotype from RNA-seq data

## Tools dependencies

 1. AdapterRemoval
  - https://github.com/MikkelSchubert/adapterremoval
  
  To install AdapterRemoval: ```conda install -c maxibor adapterremoval2``` (if you have conda), otherwise go the the github page and follow the instruction.
  
 2. STAR
  - https://github.com/alexdobin/STAR
  
  To install STAR:
    **Linux**
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
    

 3. picard
  - https://broadinstitute.github.io/picard/
  
  To install picard:
   - Go to https://github.com/broadinstitute/picard/releases/tag/2.17.11
   - Download the .jar file
   - make sure you have Java, to test: `java -version`, otherwise go http://www.oracle.com/technetwork/java/javase/downloads/index.html to download Java
   - test the picard tool by `java -jar path/to/picard.jar -h`, you should see the manual of picard
   
   
 4. GATK (>4.0)
  - https://software.broadinstitute.org/gatk/download/
  
  To install GATK:
   ```
   wget https://github.com/broadinstitute/gatk/releases/download/4.0.2.1/gatk-4.0.2.1.zip
   unzip gatk-4.0.2.1.zip
   ```
 
  
 Tools are preferred to be the lastest version.

