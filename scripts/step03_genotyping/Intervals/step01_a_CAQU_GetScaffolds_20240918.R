### Script for getting long scaffolds from hog deer reference genome ###

# Author: Chris Kyriazis
# Adapted by: Joseph Curti (jcurti3@g.ucla.edu)
# Description: get individual bed files for each scaffold >1Mb
# Date: WED SEP 18 2024

## Setup Workspace

setwd("~/Downloads/")

table <- read.csv("~/Downloads/GCA_023055505.1_bCalCai1.0.p_assembly_report.txt", skip=32, header=T,sep = "\t")
table$Start.pos <- rep(0,nrow(table)) # add this for making the bed files later
head(table)
hist(table$Sequence.Length, breaks = 100000, xlim = c(0,max(table$Sequence.Length)))

max(table$Sequence.Length)
num_scaffolds <- dim(table)[1]

## get number of scaffolds greater than some length
length <- 1000000
total_len = sum(table$Sequence.Length)
n_scaf_long <- sum(table$Sequence.Length>length)
n_scaf_long

## calculate percent of reference covered by scaffolds greater than length 
sum(head(sort(table$Sequence.Length, decreasing = T), n = n_scaf_long))/total_len

## length of reference not included in large scaffs
total_len-sum(head(sort(table$Sequence.Length, decreasing = T), n = n_scaf_long))

# create new data frame of the largest n_scaf scaffolds
table_sorted <- table[order(-table$Sequence.Length),]
df_large_scafs <- head(table_sorted, n=n_scaf_long)
dim(df_large_scafs)

write.csv(x = df_large_scafs, file = "hogdeer_scaffs_1Mb.csv")

plot(table_sorted$Sequence.Length, log="y")
plot(df_large_scafs$Sequence.Length)

## write bed files for long scaffolds

cols <- c("Assigned.Molecule.Location.Type","Start.pos","Sequence.Length")

for(file in 1:(n_scaf_long)){
  #cat(start_row, end_row, "\n")
  file_contents <- df_large_scafs[file,]
  file_contents <- file_contents[cols]
  #print(head(file_contents))
  
  ### Write file
  write.table(x=file_contents, file = paste("scaffold_",file,".bed", sep = ""), row.names = F, col.names = F, quote = F)
  
}





