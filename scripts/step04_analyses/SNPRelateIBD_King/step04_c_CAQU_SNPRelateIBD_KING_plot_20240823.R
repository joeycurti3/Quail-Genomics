##########################################################################
###### R Script to Read SNPRelate Pairwise IBD Table and Visualize #######
##########################################################################

# Author: Joey Curti 
# Date: FRI AUG 23 2024
# Description: Make IBD plots of California Quail WGS Data to select unrelated individuals for analysis

## Clean Workspace 

rm(list = ls())

## Dependencies
library(ggplot2)
library(corrplot)
library(RColorBrewer)
library(graph4lg)
install.packages("colorRamp2")
library(colorRamp2)

## Read in data 

data_king_unfiltered <- read.csv(<insert path to output table from SNPRelate King IBD analysis>, header = T, stringsAsFactors = F, sep = ",")
data_king_unfiltered_tidy <- data_king_unfiltered[,c(2,3,5)]
data_king_unfiltered_tidy_NoHollywoodRes <- data_king_unfiltered_tidy[data_king_unfiltered_tidy$ID1 != "LACM107363",]
data_king_unfiltered_tidy_NoHollywoodRes <- data_king_unfiltered_tidy_NoHollywoodRes[data_king_unfiltered_tidy_NoHollywoodRes$ID2 != "LACM107363",]

# Check that worked
unique(data_king_unfiltered_tidy$ID1)
unique(data_king_unfiltered_tidy$ID2)
unique(data_king_unfiltered_tidy_NoHollywoodRes$ID1)
unique(data_king_unfiltered_tidy_NoHollywoodRes$ID2)

## edgelist to matrix
IBDmatrix_King_unfiltered <- df_to_pw_mat(data_king_unfiltered_tidy_NoHollywoodRes, from="ID1", to="ID2", value="kinship")

# make self-comparison 1 instead of 0

IBDmatrix_King_unfiltered[IBDmatrix_King_unfiltered<0] <- -.1 # anything that is negative, set to same value
IBDmatrix_King_unfiltered[IBDmatrix_King_unfiltered==0.0000] <- 1 # set self comparison = 1


## Plot it

# Note: since these are kinship estimates, we need to change the color binning. See here for list of kinship coefficients: 
# https://jvanderw.une.edu.au/Mod8Lectures_Relatedness.pdf
offset = 1e-10
col_fun_King = colorRamp2(c(-.1,.0625, .0625+offset, .125, .125+offset, .25, .25+offset, 1), 
                          c("grey","grey","yellow","orange","red","darkred","darkred","black"))

pdf(file="SNPRelateIBD_King_NoLDNoMAF_plot_20240823.pdf", width = 15, height = 15) # can't save this with ggsave, since not ggplot object 
pp1 <- Heatmap(IBDmatrix_King_unfiltered, col = col_fun_King, 
                                   heatmap_legend_param = list(title="Kinship \nCoefficient"),
                                   rect_gp = gpar(col = "white", lwd = 2),
                                   cell_fun = function(j, i, x, y, width, height, fill) {
                                     if(IBDmatrix_King_unfiltered[i, j] > .009)
                                       grid.text(sprintf("%.2f", IBDmatrix_King_unfiltered[i, j]), x, y, gp = gpar(fontsize = 4))
                                   })
draw(pp1)
dev.off() # saves the plot to wd 

## Cleanup

closeAllConnections()
