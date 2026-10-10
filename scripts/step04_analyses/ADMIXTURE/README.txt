The following files/scripts are in this directory and should be executed in this order: 

1. run_step04_a_CAQU_SNPRelate_LDandMAFpruning_20240624.sh - Wrapper script required to submit Step 2
2. step04_a_CAQU_SNPRelate_LDandMAFpruning_20240624.R - Perform linkage and MAF pruning on GFS file, output new GDS file (input to Step 4)
3. run_step04_b_CAQU_SeqArray2VCF_20240624.sh - Wrapper script needed to run Step 4
4. step04_b_CAQU_SeqArray2VCF_20240624.R - Take filtered GDS from Step 2 and convert to a VCF
5. step04_c_CAQU_VCF2PLINK_20240624.sh - Take VCF from Step 4 and convert to a PLINK file
6. step04_d_CAQU_ADMIXTURE_SAMO_20240624.sh - Run ADMIXTURE analysis on a subset of samples from the Santa Monica Mountains
7. step04_d_CAQU_ADMIXTURE_unrelated_20240619.sh - Run ADMIXTURE analysis on a subset of all unrelated samples
8. step04_e_CAQU_ADMIXTURE_pophelper_20240823.R - Visualize the outputs of Steps 6-7
