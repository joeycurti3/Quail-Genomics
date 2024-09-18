'Scripts' directory is broked down into four subdirectories: 

1) step01_sequencing: scripts for checking md5 sums of sequence data, downloading and preparing the reference genome, etc.
2) step02_preprocessing: all steps in the process up to joint genotype calling, including aligning FASTQ to refernece, SAM to BAM conversion, GATK MarkDuplicates, and GATK HaplotypeCaller
3) step03_genotyping: joint genotyping and filtering in GATK
4) step04_analyses: all analyses presented in the manuscript 
