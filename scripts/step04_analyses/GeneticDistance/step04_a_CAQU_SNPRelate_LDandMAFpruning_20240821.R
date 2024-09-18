####################################################################################################
###### LD Pruning and subsetting on all SAMO CAQU samples before calculating genetic distance ######
###################################################################################################

# Author: Joseph Curti (jcurti3@g.ucla.edu)
# Date: MON JUN 24 2024

## References: 
# http://bioconductor.org/packages/release/bioc/vignettes/SNPRelate/inst/doc/SNPRelate.html
# https://uw-gac.github.io/SISG_2021/ancestry-and-relatedness-inference.html
# https://github.com/snigenda/FinWhale_PopGenomics_2022/blob/eb6fc440af3357771216d8bb6f45490e6e7d9499/population_structure/pca_relatedness_fst_admixture/step1_gdsPASSLDPruning_all50_20210316.R

## Dependencies

library(gdsfmt)
library(SNPRelate)

## Read in Data

gdsfmt::showfile.gds(closeall=TRUE) # make sure file is not already open
gdsfile <- <insert path to .gds file 'GCA_023055505.1_bCalCai1.0.p_PassSNPs.gds'>
genofile <- snpgdsOpen(gdsfile, readonly = TRUE)
all_samples <- c("LACM107541","LACM112287","MVZCCGP-CaOr104_I-H10","MVZCCGP-CaOr30_I-F03","MVZCCGP-CaOr33_I-A04","MVZCCGP-CaOr35_I-C04","T05B038","T11B098","T11B111","T1B050","T1B052","T1B054","T1B075","T1B081","T1B082","T2B001","T2B029","T2B060","T2B084","T2B085","T2B086","T3B024","T3B032","T3B066","T3B091","T3B092","T3B109","T4B004","T4B009","T4B057","T4B073","T4B074","T4B096","T4B110","T5B093","T6B055","T6B067","T6B107","T7B013","T7B036","T8B041","T8B070","T8B094","T9B090","WFVZ52698","WFVZ53206")

## LD Pruning 
set.seed(1000)

# Subset and prune samples

snpset_maf01_ld2_SAMO <- snpgdsLDpruning(genofile, ld.threshold = 0.2, maf = .01, autosome.only = FALSE, method = 'composite', sample.id=all_samples)

# Get the SNP ID names
snpset_id_maf01_ld2_SAMO <- unlist(unname(snpset_maf01_ld2_SAMO))

# Output new GDS files
snpgdsCreateGenoSet(gdsfile,paste0(<insert path to .gds file without .gds extention>,"_maf01_ld2_SAMO.gds"),snp.id=snpset_id_maf01_ld2_SAMO, sample.id=all_samples)

# Confirm that subset worked
gdsfmt::showfile.gds(closeall=TRUE)
gdsfile_pruned_SAMO <- <insert path to new .gds file 'GCA_023055505.1_bCalCai1.0.p_PassSNPs_maf01_ld2_SAMO.gds'>
genofile_pruned_SAMO <- snpgdsOpen(gdsfile_pruned_SAMO, readonly = TRUE)
sample_id_pruned_SAMO <- read.gdsn(index.gdsn(genofile_pruned_SAMO, "sample.id"))
print("Sample order after subset and pruning, unrelated SAMO samples")
print(sample_id_pruned_SAMO)

## another GDS but without pruning

gdsfmt::showfile.gds(closeall=TRUE) # make sure file is not already open
gdsfile <- <insert path to .gds file 'GCA_023055505.1_bCalCai1.0.p_PassSNPs.gds'>
genofile <- snpgdsOpen(gdsfile, readonly = TRUE)

# Output new GDS files
snpgdsCreateGenoSet(gdsfile,paste0(<insert path to .gds file without .gds extention>,"_SAMO_SubsetNoPrune.gds"),sample.id=all_samples)

# Confirm that subset worked
gdsfile_SAMO_NoPrune <- <insert path to new .gds file 'GCA_023055505.1_bCalCai1.0.p_PassSNPs_maf01_ld2_SAMO_SubsetNoPrune.gds'>
genofile_SAMO_NoPrune <- snpgdsOpen(gdsfile_SAMO_NoPrune, readonly = TRUE)
sample_id_SAMO_NoPrune <- read.gdsn(index.gdsn(genofile_SAMO_NoPrune, "sample.id"))
print("Sample order after subset and pruning, unrelated SAMO samples")
print(sample_id_SAMO_NoPrune)
gdsfmt::showfile.gds(closeall=TRUE)

