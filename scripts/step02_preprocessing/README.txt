The following files/scripts are in this directory and should be executed in this order: 

1. step02_a_CAQU_fastqc_20240404.sh - Run Fastqc on all fastq files to assess quality 
2. step02_b_CAQU_FastqToSam_MarkAdapters_20240423.sh - Mark illumina adapters and create an unmapped BAM file needed in downstream steps
3. step02_c_CAQU_Align_Ref2fq_20240418.sh - Align sequence data to reference genome using BWA-MEM
4. step02_d_CAQU_mergealign_qualimap_20240418.sh - Merge unmapped and mapped BAM files, assess quality of alignment using Qualimap
5. step02_e_CAQU_markduplicates_20240422.sh - Mark duplicate reads
6. step02_f_CAQU_Haplotypecaller_ExcludeList_20240422.sh - Call haplotypes, needed for joint genotyping 
