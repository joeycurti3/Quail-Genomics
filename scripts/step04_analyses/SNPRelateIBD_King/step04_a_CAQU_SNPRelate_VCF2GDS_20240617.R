###########################################################################
###### R Script to Convert VCF Files to GDS for Input into SNPRelate#######
###########################################################################

# Author: Joey Curti (jcurti3@g.ucla.edu) 
# Revised from Script by: Chris Kyriazis
# Date: MON JUN 17 2024

## Clean Workspace 

rm(list = ls())

## Dependencies 

library(gdsfmt)
library(SNPRelate)

## Define Variables

todaysdate=format(Sys.Date(),format="%Y%m%d")
INDIR=<insert path to VCF containing only biallelic PASS SNPs> # Where your snp vcf file is 
INFILE_NAME="GCA_023055505.1_bCalCai1.0.p_PassSNPs" # exclude the .vcf.gz suffix 
OUTDIR=<insert path to output directory>
vcf.fn = paste(INDIR,INFILE_NAME,".vcf.gz",sep="")

## Main 

# VCF2GDS

snpgdsVCF2GDS(vcf.fn, paste(OUTDIR,INFILE_NAME,".gds",sep=""), method="biallelic.only")

# summary 

# Open an output file

sink(file=paste(OUTDIR,INFILE_NAME,"_NoLDNoMAF.gds.summary.txt",sep=""))
snpgdsSummary(paste(OUTDIR,INFILE_NAME,".gds",sep=""))

# close summary file

sink()

## Cleanup

closeAllConnections()
