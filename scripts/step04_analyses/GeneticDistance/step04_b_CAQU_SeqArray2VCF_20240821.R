####################################################################
###### Converting SNPRelate GDS file to VCF File in SNPRelate ######
####################################################################

# Author: Joseph Curti (jcurti3@g.ucla.edu)
# Date: WED AUG 21 2024

## References: 
# http://bioconductor.org/packages/release/bioc/vignettes/SNPRelate/inst/doc/SNPRelate.html
# https://uw-gac.github.io/SISG_2021/ancestry-and-relatedness-inference.html

## Dependencies

library(gdsfmt)
library(SNPRelate)
library(SeqArray)

## Read in Data, Define Variables

infile_pruned <- "GCA_023055505.1_bCalCai1.0.p_PassSNPs_maf01_ld2_SAMO.gds"
outfile1_pruned <- "GCA_023055505.1_bCalCai1.0.p_PassSNPs_maf01_ld2_SAMO_SeqArrayIntFile.gds"
outfile2_pruned <- "GCA_023055505.1_bCalCai1.0.p_PassSNPs_maf01_ld2_SAMO_GDS2VCF.vcf"

infile_notpruned <- "GCA_023055505.1_bCalCai1.0.p_PassSNPs_SAMO_SubsetNoPrune.gds"
outfile1_notpruned <- "GCA_023055505.1_bCalCai1.0.p_PassSNPs_SAMO_SubsetNoPrune_SeqArrayIntFile.gds"
outfile2_notpruned <- "GCA_023055505.1_bCalCai1.0.p_PassSNPs_SAMO_SubsetNoPrune_GDS2VCF.vcf"

majorref <- FALSE

## Main 

# note this: https://github.com/zhengxwen/SeqArray/issues/22
# default to: major.ref	= TRUE
# if TRUE, use the major allele as a reference allele; otherwise, use A allele in SNP GDS file as a reference allele
# SA_mrF: SeqArray format; majorref = FALSE

# First on all SAMO samples except Hollywood Res with pruning

SeqArray::seqSNP2GDS(infile_pruned, outfile1_pruned, verbose=TRUE, optimize=FALSE, major.ref = majorref)

# open the file

SeqArray::seqSummary(outfile1_pruned)
genofile_pruned <- SeqArray::seqOpen(outfile1_pruned)

# convert to vcf (don't gzip)

SeqArray::seqGDS2VCF(genofile_pruned, vcf.fn = outfile2_pruned)

# close the GDS file

seqClose(genofile_pruned)

# Now just on SAMO samples

SeqArray::seqSNP2GDS(infile_notpruned, outfile1_notpruned, verbose=TRUE, optimize=FALSE, major.ref = majorref)

# open the file

SeqArray::seqSummary(outfile1_notpruned)
genofile_notpruned <- SeqArray::seqOpen(outfile1_notpruned)

# convert to vcf (don't gzip)

SeqArray::seqGDS2VCF(genofile_notpruned, vcf.fn = outfile2_notpruned)

# close the GDS file

seqClose(genofile_notpruned)

## Clean up 
write("Finished converting GDS to VCF on:", stdout())
date()
closeAllConnections()


