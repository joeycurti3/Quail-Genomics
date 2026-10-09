The following files/scripts are in this directory and should be executed in this order: 

1. step03_a_CAQU_GenotypeGVCF_ExcludeList_20240512.sh - Perform joint genotyping on samples
2. step03_b_CAQU_TrimAlternates_VariantAnnotator_ExcludeList_20240531.sh - Trim alternative alleles and add standard annotation to VCF
3. step03_c_CAQU_VariantFiltration_20240610.sh - Filter joint VCF at the site-level, executes Step 4 "step03_c_CAQU_CustomPyFilter_20240612.py"
4. step03_c_CAQU_CustomPyFilter_20240612.py - Submitted as part of Step 3, additional individual-level VCF filters 
5. step03_d_CAQU_FilterTally_20240612.sh - Generate counts of sites that were caught in various VCF filters
6. step03_d_CAQU_FilterTallyPlot_20240613.R - Visualize the output of Step 6 "step03_d_CAQU_FilterTally_20240612.sh"
7. step03_e_CAQU_CatVariants_mergeScaffolds_20240613.sh - Concatenate scaffold-level VCFs, all-sites VCF
8. step03_f_CAQU_SelectVariantsSNPS_GetPassSites_20240613.sh - Select biallelic SNPs for final VCF
9. Annotations and Coverage/ - Directory containing scripts needed to query VCF annotations to get filtering cutoffs and thresholds
10. Intervals/ - Directory of genomic intervals that are used for "scatter-gather" parallelization in GATK
11. Masking/ - Directory of files needed to create a masking file of repetitive regions across the genome needed in joint genotyping 
