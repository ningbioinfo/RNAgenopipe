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
  
   
 4. GATK (>4.0)
  - https://software.broadinstitute.org/gatk/download/
  
 Tools are preferred to be the lastest version.

