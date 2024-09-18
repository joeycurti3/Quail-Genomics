This directory contains a .gff file outputted from RepeatMasker of all repeat regions in the CAQU reference genome (GCA_023055505.1_bCalCai1.0.p). To use this in GATK, we need to run the following scripts:

1) step03_a_CAQU_Gff2Bed_20240918.py - this converts the .gff file to a .bed file
2) step03_b_CAQU_RepeatMaskerBed_20240918.sh - this sorts the .bed file and removes the header 
