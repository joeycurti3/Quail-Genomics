### <ins>The following files/scripts are in this directory and should be executed in this order:</ins> 

1. run_step04_a_CAQU_SlidingWindowHet_20240623 - wrapper script for step04_a_CAQU_SlidingWindowHet_20240623.py
2. step04_a_CAQU_SlidingWindowHet_20240623.py - script to calculate sliding window heterozygosity per sample
3. step04_b_SlidingWindowHet_plot_20240625.R - script to visualize output of #2
4. step04_c_CAQU_hetsignif_20240902 - script to calculate stats related to site-level mean het
5. mean_het.csv - Output of mean per basepair heterozygosity by sample
6. subset_chrom_lengths_2024.txt - Text file of chromosomes of interest and their associated length, based on the NCBI Assembly Report for the reference genome 
(https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/023/055/505/GCA_023055505.1_bCalCai1.0.p/GCA_023055505.1_bCalCai1.0.p_assembly_report.txt)

### <ins>Some of these scripts rely on input files, including chromosome lengths. I generated that input and modified it using the following code from Chris Kyriazis:</ins> 

#### *one liner for creating chrom_lengths.txt from vcf file:*
zcat my.vcf.gz | grep "contig=<ID" | cut -d'=' -f3,4 | sed 's/,length=/\t/g' | sed 's/>//g' > chrom_lengths.txt

#### *filter the resulting chrom_lengths.txt file by the scaffholds of interest in chrom_lengths_filter.txt file that is one column of scaffold names that exactly match those in first file:*
awk 'NR==FNR{a[$1][$0];next} $0 in a {for (i in a[$0]) print i}' chrom_lengths.txt chrom_lengths_filter.txt > subset_chrom_lengths.txt

### <ins>Metadata</ins> 
#### *mean_het.csv*
| Column Name | Data Type | Description | Comments |
| :---: | :---: | :---: | :---: | 
| individuals | character | Name of the sample | NA
| mean_het | numeric | Value of per base pair heterozygosity | Units: heterozygote/basepair of sequence

#### *subset_chrom_lengths_2024.txt*
| Column Name | Data Type | Description | Comments |
| :---: | :---: | :---: | :---: | 
| Column 1 | character | GenBanc accession number for scaffold | NA
| Column 2 | numeric | Sequence length of scaffold | NA
