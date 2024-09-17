#####################################################################
###### Estimating Kinship Using SNPRelate - King Robust Method ######
#####################################################################

# Author: Joseph Curti (jcurti3@g.ucla.edu)
# Date: MON JUN 17 2024

## References: 
# https://www.rdocumentation.org/packages/SNPRelate/versions/1.6.4/topics/snpgdsIBDMoM
# https://www.bioconductor.org/packages/devel/bioc/vignettes/SNPRelate/inst/doc/SNPRelate.html#installation-of-the-package-snprelate

## Dependencies

library(gdsfmt)
library(SNPRelate)
library(SeqArray)

## Read in Data, Define Variables

gdsfile <- "GCA_023055505.1_bCalCai1.0.p_PassSNPs.gds"

# Note:  KING Robust website says: "Please do not prune or filter any "good" SNPs that pass QC prior to any KING inference, unless the number of variants is too many to fit
# the computer memory, e.g., > 100,000,000 as in a WGS study, in which case rare variants can be filtered out. LD pruning is not recommended in KING."
# https://www.kingrelatedness.com/manual.shtml

## Main 

# open GDS file - no LD Pruning
genofile <- SNPRelate::snpgdsOpen(gdsfile)

# Run IBD Analysis using King Robust

ibd <- snpgdsIBDKING(genofile, verbose = TRUE,autosome.only = FALSE,remove.monosnp = T) # could add in maf here but I don't think it is needed
ibd.coeff <- snpgdsIBDSelection(ibd)
print(ibd.coeff)
write.csv(ibd.coeff,file="step04_b_CAQU_SNPRelateIBD_KING_NoLDNoMAF_table_20240617.csv",sep=",",col.names=T)

# close the GDS file
closefn.gds(genofile)

## Clean up 
write("Finished running SNPRelate IBD Analysis using King Robust Method", stdout())
date()
closeAllConnections()


