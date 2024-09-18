#! /usr/bin/env python3

import pandas as pd

#script to convert gff to bedfile
#import tab delimited file into pandas dataframe
Gff = pd.read_csv('GCA_023055505.1_bCalCai1.0.p_genomic.fna.out.gff', sep='\t', header=None)
print(Gff.head())
#Keep columns for bed file
Bed = Gff[[0,3,4]]

#add headers (comment out if you do not want a bed file with headers)
Bed.rename(columns={0: 'Chromosome', 3: 'Start', 4:'End'}, inplace=True)
print(Bed.head())

Bed['Start'] = Bed['Start'].apply(lambda x: x-1)
print(Bed.head())

#Output bed file
Bed.to_csv('GCA_023055505.1_bCalCai1.0.p_genomic_mask.bed', sep='\t', index=False)



'''
##gff-version 2
##date 2022-12-05
##sequence-region GCA_023055505.1_bCalCai1.0.p_genomic.fna
'''