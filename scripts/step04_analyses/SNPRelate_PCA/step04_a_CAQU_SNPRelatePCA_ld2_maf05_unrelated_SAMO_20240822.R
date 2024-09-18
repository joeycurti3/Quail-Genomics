#####################################################################
###### Plotting PCA with Preferred Settings: LD = .2, maf = .05 #####
#####################################################################

# Author: Joseph Curti (jcurti3@g.ucla.edu)
# Date: THU AUG 22 2024
# Version: 3 - removing hollywood res samples due to high missingness

## References: 
# http://bioconductor.org/packages/release/bioc/vignettes/SNPRelate/inst/doc/SNPRelate.html
# https://uw-gac.github.io/SISG_2021/ancestry-and-relatedness-inference.html
# https://github.com/snigenda/FinWhale_PopGenomics_2022/blob/eb6fc440af3357771216d8bb6f45490e6e7d9499/population_structure/pca_relatedness_fst_admixture/step1_gdsPASSLDPruning_all50_20210316.R

## Set up workspace

rm(list = ls())
setwd(<insert working directory path>)

## Dependencies

library(gdsfmt)
library(SNPRelate)
library(ggpubr)
library(ggplot2)
library(RColorBrewer)
library(ggrepel)

## Read in Data

gdsfmt::showfile.gds(closeall=TRUE) # make sure file is not already open
gdsfile <- <insert path to GDS file of all passing SNPs>
genofile <- snpgdsOpen(gdsfile, readonly = TRUE)
pop_code2 <- read.csv(<insert path to pop code file called 'SNPRelate_PopID_unrelated_SAMO_v2_20240822.csv'>, stringsAsFactors = F)
side_code <- read.csv(<insert path to side code file called 'SNPRelate_SideID_unrelated_SAMO_20240822.csv'>, stringsAsFactors = F)
samples <- c("LACM107541","LACM112287","MVZCCGP-CaOr104_I-H10","MVZCCGP-CaOr30_I-F03","MVZCCGP-CaOr33_I-A04","MVZCCGP-CaOr35_I-C04","T05B038","T11B098","T11B111","T1B050","T1B052","T1B054","T1B075","T1B081","T1B082","T2B001","T2B029","T2B060","T2B084","T2B085","T2B086","T3B024","T3B032","T3B066","T3B091","T3B092","T3B109","T4B004","T4B009","T4B057","T4B073","T4B074","T4B096","T4B110","T5B093","T6B055","T6B067","T6B107","T7B013","T7B036","T8B041","T8B070","T8B094","T9B090","WFVZ52698","WFVZ53206")

## Define Functions 
# From Meixi Lin

format_pca <- function(pca) {
  pc.percent <- pca$varprop*100
  tab <- data.frame(sample.id = pca$sample.id,
		    pop = pop_code2[,2],
		    side = side_code[,2],
                    PC1 = pca$eigenvect[,1]*(-1), # the first eigenvector, rotated
                    PC2 = pca$eigenvect[,2]*(-1), # the second eigenvector, rotated
                    PC3 = pca$eigenvect[,3]*(-1), # the third eigenvector, rotated
                    PC4 = pca$eigenvect[,4]*(-1), # the fourth eigenvector, rotated
                    stringsAsFactors = FALSE)
  lbls <- paste0("PC", 1:4, "(", format(pc.percent[1:4], digits=2), "%", ")")
  lims <- apply(pca$eigenvect[,1:4], 2, range)
  return(list(tab, lbls, lims))
}

## LD and MAF Pruning 
set.seed(1000)

# Perform SNP LD Pruning and set Minor Allele Frequency (MAF) thresholds

snpset_maf05_ld2 <- snpgdsLDpruning(genofile, ld.threshold = 0.2, maf = .05, autosome.only = FALSE, method = 'composite')

# Get the SNP ID names for the PCA

snpset_id_maf05_ld2 <- unlist(unname(snpset_maf05_ld2))

## PCA Analysis

# Define additional variables

my_colors <- c("#88CCEE","#CC6677","#DDCC77","#117733","#332288","#AA4499","#44AA99","#999933","#882255","#661100","#6699CC")
my_shapes <- c(1,19)

# Confirm sample order before subset
sample_id <- read.gdsn(index.gdsn(genofile, "sample.id"))
print("Sample order before subset")
print(sample_id)

# Create PCA Object using LD.threshold = .2, maf = .05 and subset of samples from unrelated SAMO CAQU

pcaout_LDpruned_maf05_ld2 <- snpgdsPCA(genofile, sample.id = samples, snp.id = snpset_id_maf05_ld2, num.thread=8, autosome.only=FALSE, verbose = FALSE)
pcaformatted_LDpruned_maf05_ld2 <- format_pca(pcaout_LDpruned_maf05_ld2)
tab_LDpruned_maf05_ld2 = pcaformatted_LDpruned_maf05_ld2[[1]]
lbls = pcaformatted_LDpruned_maf05_ld2[[2]]
lims = pcaformatted_LDpruned_maf05_ld2[[3]]

# Get sample order after subsetting
print("Sample order after subset")
print(tab_LDpruned_maf05_ld2$sample.id)

print("Sample order I specified for plotting (confirm this matches with above")
print(samples)

## Plot it! 

# Plot PC1xPC2 pcaout_LDpruned_maf05_ld2

pca1v2_maf05_ld2 <- ggplot(data=tab_LDpruned_maf05_ld2, aes(x = PC1, y = PC2)) +
    geom_point(aes(colour=factor(pop),shape=factor(side)),size=5, stroke = 2) +
    labs(x = lbls[1], y = lbls[2]) +
    scale_colour_manual(values=my_colors, name="Location") +
    scale_shape_manual(values = my_shapes, name="Side",labels=c("North/West","South/East/NA")) + 
    theme(panel.background=element_blank(),panel.border=element_rect(color="black",fill=NA), axis.text=element_text(size=15),
        axis.title=element_text(size=15)) 

#ggsave(filename="SNPRelate_PCA_maf05_ld2_unrelated_SAMO_PC1v2_nolabs_20240822.png",plot=pca1v2_maf05_ld2, device="png", height = 6, width = 8)

# Plot PC1xPC2 pcaout_LDpruned_maf05_ld2 with sample labels

pca1v2_maf05_ld2_labs <- ggplot(data=tab_LDpruned_maf05_ld2, aes(x = PC1, y = PC2)) +
    geom_point(aes(colour=factor(pop),shape=factor(side)),size=5, stroke = 2) +
    labs(x = lbls[1], y = lbls[2]) +
    geom_label_repel(aes(label = sample.id),
                  box.padding   = 0.35, 
                  point.padding = 0.5,
                  segment.color = 'grey50',
		  min.segment.length = 0,
		  max.overlaps = Inf,
		  size=1.5) +
    scale_colour_manual(values=my_colors, name="Location") +
    scale_shape_manual(values = my_shapes,name="Side",labels=c("North/West","South/East/NA")) +
    theme(panel.background=element_blank(),panel.border=element_rect(color="black",fill=NA)) +
    coord_cartesian(xlim = c(-.42,.35), ylim = c(-.5,.3))

#ggsave(filename="SNPRelate_PCA_maf05_ld2_unrelated_SAMO_PC1v2_labs_20240811.png", plot=pca1v2_maf05_ld2_labs, device="png", height = 6, width = 8)

## Plot other PCs for this dataset without labels

# Legend

#leg_pca_LDpruned_maf05_ld2 = as_ggplot(get_legend(pca1v2_maf05_ld2))
#pca1v2_maf05_ld2_noleg  <- pca1v2_maf05_ld2 + theme(legend.position = "none")

# PC1 v PC3

pca1v3_maf05_ld2 <- ggplot(data=tab_LDpruned_maf05_ld2, aes(x = PC1, y = PC3)) +
    geom_point(aes(colour=factor(pop),shape=factor(side)),size=5, stroke = 2) +
    labs(x = lbls[1], y = lbls[3]) +
    scale_colour_manual(values=my_colors, name="Location") +
    scale_shape_manual(values = my_shapes, name="Side",labels=c("North/West","South/East/NA")) +
    theme(panel.background=element_blank(),
    panel.border=element_rect(color="black",fill=NA),
    axis.text=element_text(size=15),
    axis.title=element_text(size=15))

# PC2 v PC3

pca2v3_maf05_ld2 <- ggplot(data=tab_LDpruned_maf05_ld2, aes(x = PC2, y = PC3)) +
    geom_point(aes(colour=factor(pop),shape=factor(side)),size=5, stroke = 2) +
    labs(x = lbls[2], y = lbls[3]) +
    scale_colour_manual(values=my_colors, name="Location") +
    scale_shape_manual(values = my_shapes, name="Side",labels=c("North/West","South/East/NA")) +
    theme(panel.background=element_blank(),
    panel.border=element_rect(color="black",fill=NA),
    axis.text=element_text(size=15),
    axis.title=element_text(size=15))

pp <- ggarrange(pca1v2_maf05_ld2,pca1v3_maf05_ld2,pca2v3_maf05_ld2,
	  common.legend = TRUE, 
	  legend="right",
          labels = c("A", "B", "C"),
          ncol = 2, nrow = 2)

ggsave(filename="SNPRelate_PCA_maf05_ld2_unrelated_SAMO_PC1to3_nolabs_20240822.png", plot = pp, height = 10, width = 13,device="png",dpi=600,bg="white")

## Plot other PCs for this dataset with labels

# Legend

#leg_pca_LDpruned_maf05_ld2_labs = as_ggplot(get_legend(pca1v2_maf05_ld2_labs))
#pca1v2_maf05_ld2_labs_noleg  <- pca1v2_maf05_ld2_labs + theme(legend.position = "none")

# PC1 v PC3

pca1v3_maf05_ld2_labs <- ggplot(data=tab_LDpruned_maf05_ld2, aes(x = PC1, y = PC3)) +
    geom_point(aes(colour=factor(pop),shape=factor(side)),size=5, stroke = 2) +
    labs(x = lbls[1], y = lbls[3]) +
    geom_label_repel(aes(label = sample.id),
                  box.padding   = 0.35,
                  point.padding = 0.5,
                  segment.color = 'grey50',
                  min.segment.length = 0,
                  max.overlaps = Inf,
                  size=1.5) +
    scale_colour_manual(values=my_colors, name="Location") +
    scale_shape_manual(values = my_shapes,name="Side",labels=c("North/West","South/East/NA")) +
    theme(panel.background=element_blank(),panel.border=element_rect(color="black",fill=NA)) +
        coord_cartesian(xlim = c(-.42,.35), ylim = c(-.42,.4))


# PC2 v PC3

pca2v3_maf05_ld2_labs <- ggplot(data=tab_LDpruned_maf05_ld2, aes(x = PC2, y = PC3)) +
    geom_point(aes(colour=factor(pop),shape=factor(side)),size=5, stroke = 2) +
    labs(x = lbls[2], y = lbls[3]) +
    geom_label_repel(aes(label = sample.id),
                  box.padding   = 0.35,
                  point.padding = 0.5,
                  segment.color = 'grey50',
                  min.segment.length = 0,
                  max.overlaps = Inf,
                  size=1.5) +
    scale_colour_manual(values=my_colors, name="Location") +
    scale_shape_manual(values = my_shapes,name="Side",labels=c("North/West","South/East/NA")) +
    theme(panel.background=element_blank(),panel.border=element_rect(color="black",fill=NA)) +
    coord_cartesian(xlim = c(-.45,.3), ylim = c(-.42,.4))

pp2 <- ggarrange(pca1v2_maf05_ld2_labs,pca1v3_maf05_ld2_labs,pca2v3_maf05_ld2_labs,
          common.legend = TRUE, 
          legend="right",
          labels = c("A", "B", "C"),
          ncol = 2, nrow = 2)

#ggsave(filename="SNPRelate_PCA_maf05_ld2_unrelated_SAMO_PC1to3_labs_20240811.png", plot = pp2, height = 6, width = 8,device="png",dpi=300,bg="white")
