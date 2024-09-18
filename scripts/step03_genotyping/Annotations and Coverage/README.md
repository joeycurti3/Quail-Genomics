This directory contains scripts needed to perform filtering steps. Specifically, it contains the following scripts:

1) step03_a_CAQU_getCombDP_20240602.sh - Get all SNPs in a VCF that has already had VariantAnnotator run on it and export to a table
2) step03_a_CAQU_getCombDP_20240602.R - take the table from #1 and query the depth field to determine the 99th percentile of depth across all samples
3) step03_b_CAQU_getINDcoverage_20240602.py - get individual depth annotations per sample
4) step03_c_CAQU_maxINDcoverage_99perct_20240602.sh - take the output from #3 and get 99th percentile of depth per sample
5) step03_d_CAQU_SummarizeSNPTable_20240610.R - take output from #1 and determine hard filtering thresholds for other GATK recommended filters 
