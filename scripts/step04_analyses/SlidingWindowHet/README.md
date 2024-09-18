This directory contains several scripts, including: 

1) run_step04_a_CAQU_SlidingWindowHet_20240623 - wrapper script for step04_a_CAQU_SlidingWindowHet_20240623.py
2) step04_a_CAQU_SlidingWindowHet_20240623.py - script to calculate sliding window heterozygosity per sample
3) step04_b_SlidingWindowHet_plot_20240625.R - script to visualize output of #2
4) step04_c_CAQU_hetsignif_20240902 - script to calculate stats related to site-level mean het

Some of these scripts rely on input, including chromosome lengths. I generated that using a modified one-liner from Chris Kyriazis:

## one liner for creating chrom_lengths.txt from vcf file:
zcat my.vcf.gz | grep "contig=<ID" | cut -d'=' -f3,4 | sed 's/,length=/\t/g' | sed 's/>//g' > chrom_lengths.txt

## filter the resulting chrom_lengths.txt file by the scaffholds of interest in chrom_lengths_filter.txt file that is one column of scaffold names that exactly match those in first file
awk 'NR==FNR{a[$1][$0];next} $0 in a {for (i in a[$0]) print i}' chrom_lengths.txt chrom_lengths_filter.txt > subset_chrom_lengths.txt




