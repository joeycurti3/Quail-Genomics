###################################################################
###### Calculating individual-based genetic distance matrices #####
###################################################################

# Author: Zac MacDonald (zmacdonald@ioes.ucla.edu)
# Adapted by: Joseph Curti (jcurti3@g.ucla.edu)
# Date: TUE AUG 20 2024
# Version: V1

## References: 

## Set up workspace

rm(list = ls())
setwd(<insert working directory path>)

library(vcfR)
library(adegenet)
library(StAMPP)

## Import vcf and convert to genlight object

vcf_caqu <-read.vcfR("GCA_023055505.1_bCalCai1.0.p_PassSNPs_maf01_ld2_SAMO_GDS2VCF.vcf")
genlight_caqu <- vcfR2genlight(vcf_caqu)
ind_names_genlight_caqu <- as.data.frame(indNames(genlight_caqu))
colnames(ind_names_genlight_caqu) = "individual"
#print("check that there are no populations assigned")
#genlight_disstria@pop # no populations assigned

### Euclidean genetic distance (dist {adegenet})
gen_dist_matrix_eucl <- as.matrix(dist(genlight_caqu, method = "euclidean", diag = FALSE, upper = FALSE, p = 2))
hist(gen_dist_matrix_eucl, breaks=1000) # looks good
print("Dimensions of euclidean distance matrix")
dim(gen_dist_matrix_eucl)
print("Range of euclidean distance matrix")
range(gen_dist_matrix_eucl)
write.table(gen_dist_matrix_eucl, file="outputs/euclidean_genetic_distance_CAQU_46inds.txt", 
            sep = " ", dec = ".", row.names = TRUE, col.names = TRUE, quote = FALSE)

### Bray–Curtis using ECODIST
#library("ecodist")
#bcdist(genlight_caqu)

#install.packages("gstudio")
#library("gstudio")
# import as genepop file
#genepop_file <- "vcf_str_files/Malacosoma_maxmiss_0.95_maf_0.05_disstria_east_genepop.txt"
#gstudio_file <-  read_population(genepop_file, type="genepop")

### proportion of shared alleles (should be the same as Bray-Curtis according to Shirk et al. 2017)
pop_names <- c("PacificPallisades","E_KananRd","E_N23","CorralCyn","NE_TopangaCyn","MonteNido","S_US101","W_MalibuCyn","S_US101","N_US101","N_US101","N_US101","E_MalibuCyn","W_I405","W_I405","S_StuntRd","N_US101","SE_TopangaCyn","W_I405","W_I405","W_I405","N_US101","N_US101","N_US101","E_I405","E_I405","S_US101","S_StuntRd","N_US101","SE_TopangaCyn","E_MalibuCyn","E_MalibuCyn","E_I405","S_US101","E_I405","SE_TopangaCyn","N_US101","W_MalibuCyn","N_US101","S_US101","S_US101","N_US101","E_I405","W_I405","NewburyPark","W_KananRd")
genind_caqu <- vcfR2genind(vcf_caqu)
genind_caqu@pop <- as.factor(pop_names) # needs to be null
#rm(gen_dist_matrix_Dps)
gen_dist_matrix_Dps <- as.matrix(propShared(genind_caqu))
length(unique(gen_dist_matrix_Dps))
write.table(gen_dist_matrix_Dps, file="outputs/Dps_genetic_distance_CAQU_46inds.txt", 
            sep = " ", dec = ".", row.names = TRUE, col.names = TRUE, quote = FALSE)
dim(as.matrix(genind_caqu@tab))


### Nei's 1972 distance among indivudals using StAMPP (not sure if this GD even makes sense... is this the same as Bray-Curtis?)
# need to use genlight object for this --> genlight_disstria
# assign locations (called populations here), as was done for genind above
# even though we don't use pop, need to stipulate (and put "pop = FALSE") for stamppNeisD() to run

genlight_caqu@pop <- as.factor(pop_names)
gen_dist_matrix_nei <- as.matrix(stamppNeisD(genlight_caqu, pop = FALSE))
colnames(gen_dist_matrix_nei) <- rownames(gen_dist_matrix_nei)
hist(gen_dist_matrix_nei, breaks=1000) # looks good
print("Dimensions of Nei's distance matrix")
dim(gen_dist_matrix_nei)
print("Range of Nei's distance matrix")
range(gen_dist_matrix_nei)
write.table(gen_dist_matrix_nei, file="outputs/nei_genetic_distance_CAQU_46inds.txt", 
            sep = " ", dec = ".", row.names = TRUE, col.names = TRUE, quote = FALSE)

