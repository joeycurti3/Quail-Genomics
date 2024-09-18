####################################################################
###### R script to visualize pairwise Fst values from Plink2 #######
####################################################################

# Author: Joey Curti 
# Date: WED JUL 31 2024
# Description: Make heatmaps of pairwise Fst values for quail samples

## Clean Workspace and set working directory

setwd(<insert working directory path>)
rm(list = ls())

## Dependencies
library(ggplot2)
library(corrplot)
library(graph4lg)
library(circlize)
library(vegan)

## Read in data 

# F_st data 

hudson_fst <- read.csv(<insert path to file 'GCA_023055505.1_bCalCai1.0.p_PassSNPs_Plink_Hudson.fst.summary'>, sep = "\t", header = T, stringsAsFactors = F)
hudson_fst_filtered <- hudson_fst[hudson_fst$X.POP1 %in% c("NewburyPark","E_N23","W_KananRd","E_KananRd","N1_US101","N2_US101","N3_US101","N5_US101",
                                                        "S3_US101","S4_US101","W_MalibuCyn","E_MalibuCyn","E_TopangaCyn","MonteNido","S_StuntRd","CorralCyn",
                                                        "PacificPallisades","W_I405","E_I405","HollywoodReservoir") & hudson_fst$POP2 %in% 
                                 c("NewburyPark","E_N23","W_KananRd","E_KananRd","N1_US101","N2_US101","N3_US101","N5_US101",
                                   "S3_US101","S4_US101","W_MalibuCyn","E_MalibuCyn","E_TopangaCyn","MonteNido","S_StuntRd","CorralCyn",
                                   "PacificPallisades","W_I405","E_I405","HollywoodReservoir"),]

# geospatial data

geodata <- read.csv(<insert path to file '20240731_CAQU_Locations.csv'>, sep=",",header = T,stringsAsFactors = F)

## Clean up data 

# Calculate linearized F_st

hudson_fst_filtered$HUDSON_FST_LINEAR <- (hudson_fst_filtered$HUDSON_FST/(1-hudson_fst_filtered$HUDSON_FST))

# Edgelist to matrix
hudson_fst_matrix <- df_to_pw_mat(hudson_fst_filtered, from="X.POP1", to="POP2", value="HUDSON_FST_LINEAR")
names_order <- c("NewburyPark","E_N23","W_KananRd","E_KananRd","N1_US101","N2_US101","N3_US101","N5_US101",
                 "S3_US101","S4_US101","W_MalibuCyn","E_MalibuCyn","E_TopangaCyn","MonteNido","S_StuntRd","CorralCyn",
                 "PacificPallisades","W_I405","E_I405","HollywoodReservoir")
hudson_fst_matrix_ordered <- reorder_mat(hudson_fst_matrix,names_order)

geodata_mat <- mat_geo_dist(data=geodata,ID="X",x="Longitude",y="Latitude",crds_type = "polar", gc_formula = "vicenty")
geodata_mat_ordered <- reorder_mat(geodata_mat,names_order)

## Set self comparison in F_st data == 1

hudson_fst_matrix[hudson_fst_matrix == 0.00000000] <- 1

# Sanity check that the rows and columns match up for both symmetrical matrices

rownames(hudson_fst_matrix_ordered) == rownames(geodata_mat_ordered) 
colnames(hudson_fst_matrix_ordered) == colnames(geodata_mat_ordered)

## Plot it

# Set color bins

min(hudson_fst_filtered$HUDSON_FST_LINEAR) #-0.008104697
max(hudson_fst_filtered$HUDSON_FST_LINEAR) #0.1597027

offset = 1e-10
col_fun_hudson = colorRamp2(c(-0.008104697, -0.008104697+offset, .01, .0499, .05, .05+offset, .1, .1+offset, 1), 
                            c("white","white","white","white","lightyellow","lightyellow","coral","coral","grey30"))

# Heatmap of Fst values

png("CAQU_Plink_Fst_plot_20240731.png",res=300,units = "in", height=10, width = 10)
hudson_fst_matrix_plot <- Heatmap(hudson_fst_matrix, col = col_fun_hudson,
                                  heatmap_legend_param = list(title=expression("Pairwise F"[ST])),
                                  rect_gp = gpar(col = "white", lwd = 2),
                                  column_order=c("NewburyPark","E_N23","W_KananRd","E_KananRd","N1_US101","N2_US101","N3_US101","N5_US101",
                                                 "S3_US101","S4_US101","W_MalibuCyn","E_MalibuCyn","E_TopangaCyn","MonteNido","S_StuntRd","CorralCyn",
                                                 "PacificPallisades","W_I405","E_I405","HollywoodReservoir"),
                                  row_order=c("NewburyPark","E_N23","W_KananRd","E_KananRd","N1_US101","N2_US101","N3_US101","N5_US101",
                                                 "S3_US101","S4_US101","W_MalibuCyn","E_MalibuCyn","E_TopangaCyn","MonteNido","S_StuntRd","CorralCyn",
                                                 "PacificPallisades","W_I405","E_I405","HollywoodReservoir"),
                                  show_column_dend = FALSE, show_row_dend = FALSE,
                                  cell_fun = function(j, i, x, y, width, height, fill) {
                                    if(hudson_fst_matrix[i, j] > -.1)
                                      grid.text(sprintf("%.3f", hudson_fst_matrix[i, j]), x, y, gp = gpar(fontsize = 8))
                                  })
draw(hudson_fst_matrix_plot)
dev.off()

# Scatterplot of euclidean distance v genetic distance

genoXgeo_plot_hudson <- scatter_dist(mat_gd = as.matrix(hudson_fst_matrix_ordered), mat_ld = geodata_mat_ordered, method="lm", pts_col = "black") + 
  geom_point(size=2) + 
  labs(y=expression(frac(F["ST"],(1-F["ST"]))),x="Geographic Distance (m)") +
  ggpubr::theme_pubr()

ggsave("CAQU_Plink_Fst_IBDPlot_20240731.png",plot=genoXgeo_plot_hudson,device = "png",dpi = 300)

# test to see if significant relationship between euclidean distance and genetic distance 

vegan::mantel(hudson_fst_matrix_ordered,geodata_mat_ordered, na.rm = T)

# Mantel statistic based on Pearson's product-moment correlation 
# 
# Call:
# vegan::mantel(xdis = hudson_fst_matrix_ordered, ydis = geodata_mat_ordered,      na.rm = T) 
# 
# Mantel statistic r: 0.4868 
#       Significance: 0.017 
# 
# Upper quantiles of permutations (null model):
#   90%   95% 97.5%   99% 
# 0.240 0.357 0.443 0.514 
# Permutation: free
# Number of permutations: 999

