#############################################################
###### R Script to plot sliding window heterozygosity #######
#############################################################

#Author: Joey Curti(jcurti3@g.ucla.edu)
#Adapted From: Chris Kyriazis and Jacqueline Robinson
#Date: MON JUN 24 2024
#Description: Script for plotting heterozygosity in sliding windows across the genome
#Uses text files containing the following columns (for example):
#scaff	window_start		sites_total	sites_unmasked	sites_passing	sites_variant	calls_sample1	calls_sample2	calls_sample3	hets_sample1	hets_sample2		hets_sample3

## References

#https://github.com/ckyriazis/moose_WGS_project/blob/master/plotting/swHet_moose_20200903.R

## Clean Workspace
rm(list=ls())
setwd(<insert path to output directory from step04_a_CAQU_SlidingWindowHet_20240623.py>)

# Dependencies
library(ggplot2)
library(gtools)
library(ggpubr)

# Define Variables
fileprefix="GCA_023055505.1_bCalCai1.0.p_het_1000000win_1000000step"
winfiles=list.files(pattern="step.txt")
winfiles=gtools::mixedsort(winfiles)
window_size=1000000
subsetscaff=c(1,2,3,4,5,6,7,8,10,11,12,13,14,15,16,17,18,19,20,21,22,24,26,27,28,29,32,33,34,35,36,37,39,40,42,43,44,45,46,47,48,50,55,56,61,63,64,65,68,70,73,74,77,78,84,86,89,102,108,109,113,114,115,117)
n_scaffolds=length(winfiles)

# Get colors for scaffolds
mycols=rep(c("#2171b5","#6baed6"),n_scaffolds/2)

# Read in Data
hetdata=read.table(winfiles[1], header=T, sep='\t')
for(i in 2:n_scaffolds){
  temp=read.table(winfiles[i],header=T,sep='\t')
  hetdata=rbind(hetdata,temp)
}

# Debug: Inspect hetdata
print(head(hetdata))
print(colnames(hetdata))

# Filter out windows with fewer than x% of sites called
threshold=.8
hetdata=subset(hetdata,hetdata$sites_total>window_size*threshold)
hist(hetdata$sites_total)
rownames(hetdata) <- 1:nrow(hetdata)

# Get scaffold position info needed for plotting
pos = as.numeric(rownames(unique(data.frame(hetdata$scaff)[1])))
pos = append(pos, length(hetdata$scaff))

# Debug: Inspect pos
print(pos)

# Get the midpoints to center scaffolds labels on x-axis
numpos = NULL
for(i in 1:(length(pos) - 1)) {
  numpos[i] = (pos[i] + pos[i + 1]) / 2
}

# Debug: Inspect numpos
print(numpos)

# Which columns contain the numbers of calls and the numbers of hets?
colnames(hetdata)
callcols = c(4:65)
hetcols = c(66:127)

# Get sample names for plot
sampnames = gsub("calls_", "", names(hetdata)[callcols])

# Which individuals to plot?
results <- character(length(sampnames)) 
for (i in seq_along(sampnames)) {
  var_name <- paste0("samp", i)
  assign(var_name, sampnames[i])
}

# Make a function to plot

# Initialize an empty data frame to store the results

mean_het_df <- data.frame(sampname = character(), mean_het = numeric(), stringsAsFactors = FALSE)

# Initialize an empty list to store the plots

het_list <- list()

# Define the hetplot function

hetplot <- function(sampname, ymax, title, xlab, ylab, size){
  # Create a blank plot
  plot(0, 0, type="n", xlim=c(0, pos[length(pos)]), ylim=c(0, ymax), axes=F, xlab="", ylab=ylab, main=title, cex.lab=size, cex.main=size)
  
  aa <- which(sampnames == sampname)
  for (i in 1:n_scaffolds){
    temp <- hetdata[which(hetdata$scaff == unique(hetdata$scaff)[i]),]
    
    # Check if temp is not empty
    if(nrow(temp) == 0) {
      next
    }
    
    x <- as.numeric(rownames(temp))
    y <- temp[, hetcols[aa]] / temp[, callcols[aa]]
    
    # Debug: Check lengths of x and y
    print(paste("Scaffold:", i))
    print(paste("Length of x:", length(x)))
    print(paste("Length of y:", length(y)))
    
    if(length(x) != length(y)) {
      stop("Lengths of x and y differ")
    }
    
    for(j in 1:length(x)){
      lines(x = c(x[j], x[j]), c(0, y[j]), col = mycols[i], lwd = 1.1)
    }
    
    lines(as.numeric(rownames(temp)), y, col = mycols[i], lwd = 1.5)
  }
  
  # Draw line of mean genome-wide heterozygosity
  mean_het <- mean(hetdata[, hetcols[aa]] / hetdata[, callcols[aa]], na.rm = TRUE)
  abline(h = mean_het, col = "red", lty = 2)
  mean_het_df <<- cat(sampname, ",", mean_het, "\n") ## Need to debug this. Some samples are repeated, some are missing
  
  # Add y-axis
  axis(2)
  
  # Add x-axis and labels
  title(xlab = xlab, line = 2, cex.lab = size)
  axis(side = 1, at = pos, labels = FALSE)
  axis(side = 1, at = numpos, tick = FALSE, labels = subsetscaff, las = 3, cex.axis = .8, line = -.2)
  
  # Capture the plot and return it
  plot_obj <- recordPlot()
  return(plot_obj)
}

## Plot

ymax=0.01
size=1.5

# East I-405

pdf(paste(fileprefix,"_E405_20240625.pdf",sep=""),height=8,width=10)

# 3 columns, 2 rows (6 samples)
par(mfrow=c(2,3))

samples_E405 <- c("T3B091","T3B092","T8B094","T5B095","T5B093","T4B096")
  
for (i in samples_E405) {
  hetplot(i,ymax,i,"Scaffold","Heterozygosity",size)
}

#Close figure file
dev.off()

# West I-405

pdf(paste(fileprefix,"_W405_20240625.pdf",sep=""),height=8,width=8)

# 3 columns, 2 rows (6 samples)
par(mfrow=c(4,3))

samples_W405 <- c("T1B081","T1B083","T2B088","T2B086","T9B090","T1B082","T2B085","T2B084","T2B089","T2B087")

for (i in samples_W405) {
  hetplot(i,ymax,i,"Scaffold","Heterozygosity",size)
}

#Close figure file
dev.off()

# East Malibu Cyn.

pdf(paste(fileprefix,"EMC_202406256.pdf",sep=""),height=5,width=8)

#  columns,  rows (3 samples)
par(mfrow=c(1,3))

samples_EMC <- c("T4B073","T4B074","T1B075")

for (i in samples_EMC) {
  hetplot(i,ymax,i,"Scaffold","Heterozygosity",size)
}

#Close figure file
dev.off()

# West Malibu Cyn.

pdf(paste(fileprefix,"WMC_202406256.pdf",sep=""),height=5,width=8)

# 3 columns, 1 rows (3 samples)
par(mfrow=c(1,3))

samples_WMC <- c("T11B098","T3B103","T6B107")

for (i in samples_WMC) {
  hetplot(i,ymax,i,"Scaffold","Heterozygosity",size)
}

#Close figure file
dev.off()

# North US101

pdf(paste(fileprefix,"N101_202406256.pdf",sep=""),height=8,width=10)

# 4 columns, 3 rows (12 samples)
par(mfrow=c(3,4))

samples_N101 <- c("T4B006","T4B009","T7B013","T3B024","T2B029","T3B032","T1B050","T1B052","T1B054",
                  "T3B066","T6B067","T8B070")

for (i in samples_N101) {
  hetplot(i,ymax,i,"Scaffold","Heterozygosity",size)
}

#Close figure file
dev.off()

# South US101

pdf(paste(fileprefix,"S101_202406256.pdf",sep=""),height=8,width=8)

# 3 columns, 3 rows (7 samples)
par(mfrow=c(3,3))

samples_S101 <- c("T4B035","T7B036","T8B041","T05B038","T3B109","T4B110","T11B111")

for (i in samples_S101) {
  hetplot(i,ymax,i,"Scaffold","Heterozygosity",size)
}

#Close figure file
dev.off()

# East and West of N9 (Kanan Rd)

pdf(paste(fileprefix,"EN9WN9_202406256.pdf",sep=""),height=4,width=8)

# 3 columns, 3 rows (2 samples)
par(mfrow=c(1,2))

samples_EN9WN9 <- c("LACM112287","WFVZ53206")

for (i in samples_EN9WN9) {
  hetplot(i,ymax,i,"Scaffold","Heterozygosity",size)
}

#Close figure file
dev.off()

# Northern CA/Southern CA

pdf(paste(fileprefix,"NCASCA_202406256.pdf",sep=""),height=8,width=8)

# 3 columns, 2 rows (6 samples)
par(mfrow=c(2,3))

samples_NCASCA <- c("MVZCCGP.CaOr96_I.A10","MVZCCGP.CaOr45_I.C05","MVZCCGP.CaOr37_II.A01","MVZCCGP.CaOr81_II.B02","MVZCCGP.CaOr78_II.H01","MVZCCGP.CaOr71_II.F01")

for (i in samples_NCASCA) {
  hetplot(i,ymax,i,"Scaffold","Heterozygosity",size)
}

#Close figure file
dev.off()

# Create plots of per population het

SampleID <- c(
  "T11B098", "T1B075", "T1B081", "T1B083", "T2B088", "T3B066", "T3B091", "T3B092", 
  "T3B103", "T4B035", "T4B073", "T4B074", "T6B067", "T6B107", "T7B036", "T8B041", 
  "T8B070", "T8B094", "T1B050", "T2B086", "T5B095", "T4B006", "T9B090", "T4B110", 
  "T4B005", "T1B082", "T1B052", "T4B009", "T1B054", "T2B085", "T05B038", "T2B084", 
  "T2B029", "T4B057", "T11B111", "T5B093", "T2B001", "T2B089", "T4B096", "T4B004", 
  "T6B055", "T2B087", "T2B060", "T3B024", "T3B109", "T3B032", "T7B013", 
  "LACM107363", "LACM107541", "LACM112287", "WFVZ52698", "WFVZ53206", 
  "MVZCCGP-CaOr35_I-C04", "MVZCCGP-CaOr33_I-A04", "MVZCCGP-CaOr30_I-F03", 
  "MVZCCGP-CaOr104_I-H10", "MVZCCGP-CaOr81_II-B02", "MVZCCGP-CaOr78_II-H01", 
  "MVZCCGP-CaOr71_II-F01", "MVZCCGP-CaOr96_I-A10", "MVZCCGP-CaOr45_I-C05", 
  "MVZCCGP-CaOr37_II-A01"
)

locality <- c(
  "W_MalibuCyn", "E_MalibuCyn", "W_I405", "W_I405", "W_I405", "N_US101", "E_I405", "E_I405", 
  "W_MalibuCyn", "S_US101", "E_MalibuCyn", "E_MalibuCyn", "N_US101", "W_MalibuCyn", "S_US101", 
  "S_US101", "N_US101", "E_I405", "N_US101", "W_I405", "E_I405", "N_US101", "W_I405", "S_US101", 
  "S_StuntRd", "W_I405", "N_US101", "N_US101", "N_US101", "W_I405", "S_US101", "W_I405", "N_US101", 
  "E_TopangaCyn", "S_US101", "E_I405", "S_StuntRd", "W_I405", "E_I405", "S_StuntRd", "E_TopangaCyn", "W_I405", 
  "E_TopangaCyn", "N_US101", "S_US101", "N_US101", "N_US101", "HollywoodReservoir", "PacificPallisades", 
  "E_KananRd", "NewburyPark", "W_KananRd", "MonteNido", "E_TopangaCyn", "CorralCyn", "E_N23", 
  "SouthernClade", "SouthernClade", "SouthernClade", "NorthernClade", "NorthernClade", "NorthernClade"
)

mean_het <- c(
  0.0051871100, 0.0052869000, 0.0050062900, 0.0052669900, 0.0048924500, 0.0049995800, 
  0.0050836600, 0.0047434000, 0.0051745600, 0.0051066000, 0.0051980500, 0.0049590200, 
  0.0052288300, 0.0049459700, 0.0051248500, 0.0051025900, 0.0050236700, 0.0050776500, 
  0.0050555100, 0.0051286700, 0.0049667400, 0.0050299500, 0.0051687400, 0.0051332100, 
  0.0049900240, 0.0051073100, 0.0050212300, 0.0051544000, 0.0050938700, 0.0046404000, 
  0.0049016300, 0.0051624300, 0.0048744400, 0.0050636850, 0.0051642800, 0.0050044200, 
  0.0051801530, 0.0051090000, 0.0051056500, 0.0050535160, 0.0050495370, 0.0050414900, 
  0.0050428310, 0.0050816700, 0.0050556100, 0.0047922800, 0.0051659500, 0.0038109760, 
  0.0052684250, 0.0052507000, 0.0052875480, 0.0052266700, 0.0046048490, 0.0048077170, 
  0.0048123130, 0.0047826340, 0.0048287900, 0.0048667100, 0.0050249600, 0.0042482900, 
  0.0047067000, 0.0043493300
)

# Create the data frame
df <- data.frame(SampleID = SampleID, locality = locality, mean_het = mean_het)


# plot dataframe by population

df$locality <- factor(df$locality, levels=c("NorthernClade","NewburyPark",
                                                          "E_N23","W_KananRd","E_KananRd",
                                                          "N_US101","S_US101","W_MalibuCyn",
                                                          "E_MalibuCyn","E_TopangaCyn","MonteNido",
                                                          "S_StuntRd", "CorralCyn","PacificPallisades",
                                                          "W_I405","E_I405","HollywoodReservoir",
                                                          "SouthernClade"), ordered = T)

df_subset <- df[df$locality %in% c("NorthernClade",
                          "N_US101","S_US101","W_MalibuCyn",
                          "E_MalibuCyn","E_TopangaCyn",
                          "S_StuntRd","W_I405","E_I405","SouthernClade"),]

het_pop_plot <- ggplot(df_subset, aes(x=locality, y=mean_het)) + 
  geom_boxplot() + 
  geom_hline(yintercept=mean(df_subset$mean_het),linetype="dashed", color = "red", linewidth=1) +
  theme_pubclean() + 
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  labs(y = expression(paste("Heterozygosity (", bp^{-1},")")),x="Population")

setwd("~/Downloads/")

ggsave("CAQU_HetByPop_20240626.png",dpi=300,device = "png")                                         

