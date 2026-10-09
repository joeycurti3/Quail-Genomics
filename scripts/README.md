'Scripts' directory is broken down into four subdirectories: 

1) step01_sequencing: scripts for checking md5 sums of sequence data, downloading and preparing the reference genome, etc.
2) step02_preprocessing: all steps in the process up to joint genotype calling, including aligning FASTQ to reference, SAM to BAM conversion, GATK MarkDuplicates, and GATK HaplotypeCaller
3) step03_genotyping: joint genotyping and filtering in GATK
4) step04_analyses: all analyses presented in the manuscript

Each of these four folders have folder-level README files that indicate the purpose of each script, the order they should be executed, and helpful notes on their usage. Typically, larger processes are numbered and should be executed in numerical order (e.g., complete "step03_genotyping" before moving to "step04_analyses"). Within-step subprocesses are typically designated by letters and should be executed in alphabetical order (e.g., run "step03_a_CAQU_GenotypeGVCF_ExcludeList_20240512.sh" then "step03_b_CAQU_TrimAlternates_VariantAnnotator_ExcludeList_20240531.sh"). 

Occasionally, wrapper scripts are needed to submit jobs in the bash environment, and they are typically denoted by the "run_" script header.


