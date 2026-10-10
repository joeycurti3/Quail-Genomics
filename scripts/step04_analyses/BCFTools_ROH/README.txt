The following files/scripts are in this directory and should be executed in this order: 

1. step04_a_CAQU_BCFTools_ROH_20240822.sh - Calculate length of runs of homozygosity (ROH) for an allsites VCF
2. step04_b_CAQU_BCFTools_subsetROH_20240626.sh - Subset the output of Step 1 to the individual-level
3. step04_c_CAQU_BCFTools_ROH_plot_20250401.R - Bin the output of Step 2 and visualize. Output a CSV of these binned, individual-level ROH values for input into Step 4
4. step04_d_CAQU_NHST_20240814.R - Analyze difference in ROH age estimates between road sampling locations using bootstrapping
5. samples.txt - Text file of sample names needed to subset ROH ouput in Step 2. Full column-level metadata below

samples.txt Metadata
Column Name / Data Type / Description / Comments 
Column 1 / Character / Sample names of 46 unrelated quail from the Santa Monica Mountains and Simi Hills / NA
