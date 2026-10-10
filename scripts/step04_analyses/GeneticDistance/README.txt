The following files/scripts are in this directory and should be executed in this order: 

1. run_step04_a_CAQU_SNPRelate_LDandMAFpruning_20240821.sh - Wrapper script needed to submit Step 2
2. step04_a_CAQU_SNPRelate_LDandMAFpruning_20240821.R - Perform linkage disequilibrium (LD) and minor allele frequency (MAF) filtering on VCF, output GDS file 
3. run_step04_b_CAQU_SeqArray2VCF_20240821.sh - Wrapper script needed to submit Step 4
4. step04_b_CAQU_SeqArray2VCF_20240821.R - Take the output of Step 2 and convert to a VCF
5. run_step04_c_CAQU_IndGeneticDist_20240820.sh - Wrapper script needed to submit Step 6
6. step04_c_CAQU_IndGeneticDist_20240820.R - Calculate individual-level genetic distance metrics
