## Sort MASK .bed file and remove header for use with GATK ## 

# I first sorted the BED file by the first column (scaffold name) and then by the second column (start position):

sort -k 1,1 -k2,2n GCA_023055505.1_bCalCai1.0.p_genomic_mask.bed > GCA_023055505.1_bCalCai1.0.p_genomic_mask_NoHead_sorted.bed

# Then I remove the first two rows that had weirdly named scaffolds:

tail -n +2 GCA_023055505.1_bCalCai1.0.p_genomic_mask_NoHead_sorted.bed > GCA_023055505.1_bCalCai1.0.p_genomic_mask_FINAL.bed

# The resulting .bed file was accepted by GATK 
