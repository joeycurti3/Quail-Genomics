##############################################################
## Plotting Output of BCFTools ROH Analysis on CAQU Genomes ##
##############################################################

# Author: Joseph Curti (jcurti3@g.ucla.edu)
# Adapted From: Chris Kyriazis
# Description: Bin ROH outputs and visualize
# Version: V3 - Removing Hollywood Res sample
# Date: WED SEP 18 2024

## References
#https://github.com/ckyriazis/moose_WGS_project/blob/master/plotting/bcftools_roh_plot_moose.R

## Clean Workspace, Set WD 

rm(list = ls())
setwd(<insert working directory with ROH output files>)

## Dependencies 

library(ggplot2)
library(tidyr)
library(dplyr)
library(ggpubr)

## Define Functions

classify_roh <- function(roh_dataframe, min_roh_length){
  short_roh <- subset(roh_dataframe,length>min_roh_length & length<1000000) 
  med_roh <- subset(roh_dataframe, length>1000000 & length<10000000)
  long_roh <-  subset(roh_dataframe, length>10000000 & length<100000000)
  
  #sum each class and divide by 1000000 to convert to Mb
  sum_short_Mb <- sum(short_roh$length)/1000000
  sum_med_Mb <- sum(med_roh$length)/1000000
  sum_long_Mb <- sum(long_roh$length)/1000000
  
  return(c(sum_short_Mb, sum_med_Mb, sum_long_Mb))

}

## Define function to read in output files for each individual and filter out ROHs less than min_roh_length
read_filter_roh <- function(data, min_roh_length){
  output <- read.table(paste(data,".gz",sep=""), col.names=c("row_type","sample","chrom","start","end","length","num_markers","qual"), fill=T)
  output1 <- subset(output, row_type == "RG")
  output_class_sums <- classify_roh(output1, min_roh_length=min_roh_length)
  return(output_class_sums)
}

## Main

# Apply function to individual ROH files 

# set filename excluding individual and ".out.gz"
file <- "GCA_023055505.1_bCalCai1.0.p_BCFTools_ROH.out_"
individuals <- c(
  "T11B098", "T1B075", "T1B081", "T1B083", "T2B088", "T3B066", "T3B091", "T3B092", 
  "T3B103", "T4B035", "T4B073", "T4B074", "T6B067", "T6B107", "T7B036", "T8B041", 
  "T8B070", "T8B094", "T1B050", "T2B086", "T5B095", "T4B006", "T9B090", "T4B110", 
  "T4B005", "T1B082", "T1B052", "T4B009", "T1B054", "T2B085", "T05B038", "T2B084", 
  "T2B029", "T4B057", "T11B111", "T5B093", "T2B001", "T2B089", "T4B096", "T4B004", 
  "T6B055", "T2B087", "T2B060", "T3B024", "T3B109", "T3B032", "T7B013", "LACM107541", "LACM112287", 
  "WFVZ52698", "WFVZ53206", 
  "MVZCCGP-CaOr35_I-C04", "MVZCCGP-CaOr33_I-A04", "MVZCCGP-CaOr30_I-F03", 
  "MVZCCGP-CaOr104_I-H10", "MVZCCGP-CaOr81_II-B02", "MVZCCGP-CaOr78_II-H01", 
  "MVZCCGP-CaOr71_II-F01", "MVZCCGP-CaOr96_I-A10", "MVZCCGP-CaOr45_I-C05", 
  "MVZCCGP-CaOr37_II-A01"
)

genome_length <- 940.401282 # length of 108 scaffolds used in analysis that are
# bigger than 1Mb and do not map to sex chromosomes


# min size allowable for ROHs - typically 100kb or 300kb
min_roh_length=100000

## Read in data, filter, and classify ROHs ##

## initialize data frame
roh_size_df <- data.frame(matrix(nrow=3, ncol=length(individuals)))
colnames(roh_size_df) <- individuals
froh_total <- c()
froh_medium <- c()
froh_large <- c()

## read in data for each individual
## note that this can be VERY slow - takes several min per individual
for(i in 1:length(individuals)){
  roh_size_df[,i] <- read_filter_roh(paste(file,individuals[i], sep=""), min_roh_length=min_roh_length)
  froh_total = c(froh_total,sum(roh_size_df[1:3,i])/genome_length) # sum ROHs >min_roh_length and divide by genome length to estimate Froh
  froh_medium = c(froh_medium,sum(roh_size_df[2:3,i])/genome_length)
  froh_large = c(froh_large,sum(roh_size_df[3,i])/genome_length)
}

froh_df <- data.frame(as.numeric(froh_total),as.numeric(froh_medium),as.numeric(froh_large),individuals)

## Plotting

## F_ROH distribution

max(froh_df$as.numeric.froh_total.) # 0.08997032
min(froh_df$as.numeric.froh_total.) # 0.008722219
mean(froh_df$as.numeric.froh_total.) # 0.02659859

png(filename = "CAQU_ROHdist_20240627.png",width=8,height=6,units = "in",res=400)
hist(froh_df$as.numeric.froh_total.,breaks=(length(froh_df$individuals)/3),xlim = c(0,.1),main = "",xlab = "F_ROH")
abline(v = mean(froh_df$as.numeric.froh_total.), col = 'red', lwd = 2, lty = 'dashed')
dev.off()

# FROH data is skewed

## binned ROH data

# Need to get the data in a usable format...

roh_size_df_tidy <- roh_size_df %>%  
  pivot_longer(cols = everything(), names_to = "Sample_name", values_to = "value") %>%
  group_by(Sample_name) %>%
  mutate(id = row_number()) %>%
  pivot_wider(names_from = id, values_from = value)

colnames(roh_size_df_tidy) <- c("Sample_name", "short_roh", "med_roh", "long_roh")

# oops, barplot needs long data... 

roh_size_df_tidy_long <- roh_size_df_tidy %>%
  pivot_longer(cols = c(short_roh, med_roh, long_roh), 
               names_to = "ROH_type", 
               values_to = "value")

#Reorder so that you move from Western SAMO to Eastern SAMO

roh_size_df_tidy_long$Sample_name <- factor(roh_size_df_tidy_long$Sample_name, 
                                         levels=c("MVZCCGP-CaOr37_II-A01","MVZCCGP-CaOr45_I-C05","MVZCCGP-CaOr96_I-A10",
                                                  "WFVZ52698","MVZCCGP-CaOr104_I-H10","WFVZ53206","LACM112287",
                                                  "T4B006","T4B009","T7B013","T3B024","T2B029","T3B032","T1B050",
                                                  "T1B052","T1B054","T3B066","T6B067","T8B070","T4B035","T7B036",
                                                  "T05B038","T8B041","T3B109","T4B110","T11B111","T11B098","T3B103",
                                                  "T6B107","T4B073","T4B074","T1B075","T6B055","T4B057","T2B060",
                                                  "MVZCCGP-CaOr33_I-A04","MVZCCGP-CaOr35_I-C04","T2B001","T4B004",
                                                  "T4B005","MVZCCGP-CaOr30_I-F03","LACM107541","T1B081","T1B082",
                                                  "T1B083","T2B084","T2B085","T2B086","T2B087","T2B088","T2B089",
                                                  "T9B090","T3B091","T3B092","T5B093","T8B094","T5B095","T4B096",
                                                  "MVZCCGP-CaOr71_II-F01","MVZCCGP-CaOr78_II-H01","MVZCCGP-CaOr81_II-B02"), 
                                         ordered = T)

# Reorder roh labels so the legend has them in the right order
roh_size_df_tidy_long$ROH_type <- factor(roh_size_df_tidy_long$ROH_type, 
                                         levels=c("long_roh","med_roh","short_roh"), 
                                         ordered = T)

# Add in population data 

individuals <- c(
  "T11B098", "T1B075", "T1B081", "T1B083", "T2B088", "T3B066", "T3B091", "T3B092", 
  "T3B103", "T4B035", "T4B073", "T4B074", "T6B067", "T6B107", "T7B036", "T8B041", 
  "T8B070", "T8B094", "T1B050", "T2B086", "T5B095", "T4B006", "T9B090", "T4B110", 
  "T4B005", "T1B082", "T1B052", "T4B009", "T1B054", "T2B085", "T05B038", "T2B084", 
  "T2B029", "T4B057", "T11B111", "T5B093", "T2B001", "T2B089", "T4B096", "T4B004", 
  "T6B055", "T2B087", "T2B060", "T3B024", "T3B109", "T3B032", "T7B013", 
  "LACM107541", "LACM112287", "WFVZ52698", "WFVZ53206", 
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
  "E_TopangaCyn", "N_US101", "S_US101", "N_US101", "N_US101", "PacificPallisades", 
  "E_KananRd", "NewburyPark", "W_KananRd", "MonteNido", "E_TopangaCyn", "CorralCyn", "E_N23", 
  "SouthernClade", "SouthernClade", "SouthernClade", "NorthernClade", "NorthernClade", "NorthernClade"
)

pop_data <- data.frame(individuals,locality)

# Merge ROH size df with the sample pop map

roh_size_df_pop <- merge(roh_size_df_tidy_long,pop_data,by.x='Sample_name',by.y='individuals')

# Reorder so that you move from Western SAMO to Eastern SAMO

roh_size_df_pop$locality <- factor(roh_size_df_pop$locality, levels=c("NorthernClade","NewburyPark",
                                                                    "E_N23","W_KananRd","E_KananRd",
                                                                    "N_US101","S_US101","W_MalibuCyn",
                                                                    "E_MalibuCyn","E_TopangaCyn","MonteNido",
                                                                    "S_StuntRd", "CorralCyn","PacificPallisades",
                                                                    "W_I405","E_I405",
                                                                    "SouthernClade"), ordered = T)

# plot it

setwd("~/Downloads/")
library(ggh4x)

roh_size_plot <- ggplot(roh_size_df_pop, aes(x = interaction(Sample_name,locality, sep = "!"), y = value, fill=ROH_type)) + 
  geom_bar(stat = "identity") +
  scale_fill_manual(breaks = c("short_roh","med_roh","long_roh"),
                    values = c("seashell","seashell3", "seashell4"),
                    labels = c(paste(min_roh_length/1000000,"-1 Mb",sep=""), "1-10 Mb", "10-100 Mb"),
                    name = "ROH Size") +
  theme_pubclean() +
  theme(axis.text.x = element_text(angle = 90,hjust=1.02,vjust=.5,size=9),
        ggh4x.axis.nesttext.x = element_text(angle = 90,hjust=1)) + 
  scale_x_discrete(guide = guide_axis_nested(delim = "!",extend = .8), 
                   name = "Locality") +
  labs(x="Sample", y="Summed ROH length (Mb)")

ggsave(plot = roh_size_plot,"roh_size_plot_20240903.png",device="png",dpi=400)

## Make bar plots of FROH and ROH by length for each populations


# Merge F_ROH df with the sample pop map

froh_merged_df <- merge(froh_df,pop_data,by.x='individuals',by.y='individuals')

# Reorder so that you move from Western SAMO to Eastern SAMO

froh_merged_df$locality <- factor(froh_merged_df$locality, levels=c("NorthernClade","NewburyPark",
                                                          "E_N23","W_KananRd","E_KananRd",
                                                          "N_US101","S_US101","W_MalibuCyn",
                                                          "E_MalibuCyn","E_TopangaCyn","MonteNido",
                                                          "S_StuntRd", "CorralCyn","PacificPallisades",
                                                          "W_I405","E_I405",
                                                          "SouthernClade"), ordered = T)

# Plot it

Froh_pop_plot <- ggplot(froh_merged_df, aes(x=locality, y=as.numeric.froh_total.)) + 
  geom_boxplot() + 
  geom_hline(yintercept=mean(froh_merged_df$as.numeric.froh_total.),linetype="dashed", color = "red", linewidth=1) +
  theme_pubclean() + 
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  labs(y = bquote(F[ROH]),x="Population")

ggsave(plot=Froh_pop_plot,"CAQU_FROHbyPop_20240627.PNG",device="png",dpi=300)

## New boxplot between NORCAL, SAMO, and SOCAL, with significance testing

froh_df_state <- read.csv("~/Downloads/20240902_CAQU_seq_metadata.csv",sep=",",header = T,stringsAsFactors = F)
froh_df_state <- froh_df_state[-1,] # Remove Hollywood Res sample

# Test assumptions of ANOVA
library(afex)
library(performance)
aov_roh_state <- aov_ez("individuals", "froh", froh_df_state, between = "state_site")

# Anova Table (Type 3 tests)
# 
# Response: froh
# Effect    df  MSE    F  ges p.value
# 1 state_site 2, 58 0.00 1.45 .048    .244
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘+’ 0.1 ‘ ’ 1

check_homogeneity(aov_roh_state)
#OK: There is not clear evidence for different variances across groups (Levene's Test, p = 0.926)
is_norm <- check_normality(aov_roh_state)
plot(is_norm)
#Warning: Non-normality of residuals detected (p < .001)

kw_test <- kruskal.test(froh ~ state_site, data = froh_df_state)

# Kruskal-Wallis rank sum test
# 
# data:  froh by state_site
# Kruskal-Wallis chi-squared = 3.801, df = 2, p-value = 0.1495

# Plot it

Froh_state_barplot <- ggplot(froh_df_state, aes(x=state_site,y=froh,fill=state_site)) +
  geom_boxplot() + 
  labs(y=expression(F[ROH]),x="\nSampling Location") +
  scale_fill_manual(values=c("seashell","seashell2","seashell3")) +
  scale_x_discrete(labels=c("NORCAL","SAMO","SOCAL")) +
  theme_pubclean() +
  theme(legend.position = "none")

ggsave("FROHbyStateSite_CAQU_20240903.png",plot = Froh_state_barplot, device="png",dpi=400,heigh=5, width = 6)

## New boxplot of just the 405 and 101 samples against everything else 

froh_df_final <- read.csv("~/Downloads/20240902_CAQU_seq_metadata.csv",sep=",",header = T,stringsAsFactors = F)
froh_df_final <- froh_df_final[-1,] # Remove Hollywood Res sample
froh_df_final$groups <- NA

froh_df_final %>% #Add in grouping variable for ggplot
  mutate(groups = case_when((locality == "S_US101") ~ "a",
                             (locality == "N_US101") ~ "a",
                             (locality == "W_I405") ~ "a",
                             (locality == "E_I405") ~ "a",
                             .default = "b")) %>%
  group_by(groups) %>% summarise(n = n()) # Test to see if the counts make sense

# groups     n
# <chr>  <int>
# 1 a         35
# 2 b         26
# They do! 

froh_df_final_plot <- froh_df_final %>%
  mutate(groups = case_when((locality == "S_US101") ~ "a",
                            (locality == "N_US101") ~ "a",
                            (locality == "W_I405") ~ "a",
                            (locality == "E_I405") ~ "a",
                            .default = "b"))


# Plot it

Froh_2pop_barplot <- ggplot(froh_df_final_plot, aes(x=groups,y=froh,fill=groups)) +
  geom_boxplot() + 
  labs(y=expression(F[ROH]),x="") +
  scale_fill_manual(values=c("seashell3","seashell")) +
  scale_x_discrete(labels=c("US101 & I405","All Other Roads")) +
  theme_pubr() +
  geom_signif(comparisons = list(c("a", "b")), 
              map_signif_level=TRUE,
              annotations = "n.s.",
              y_position = 0.1) +
  theme(legend.position = "none")

ggsave("20240718_FROH_405and101comparisonplot.png",plot=Froh_2pop_barplot,device = "png",dpi=300)

## Now make boxplot with the three binned sizes b/w 405/101 and everything else

froh_merged_df_groups <- merge(roh_size_df_tidy_long,pop_data,by.x='Sample_name',by.y='individuals')

froh_merged_df_groups$groups <- NA

froh_merged_df_groups <- froh_merged_df_groups %>% #Add in grouping variable for ggplot
  mutate(groups = case_when((locality == "S_US101") ~ "a",
                            (locality == "N_US101") ~ "a",
                            (locality == "W_I405") ~ "a",
                            (locality == "E_I405") ~ "a",
                            .default = "b")) 

# Plot it

Froh_2pop_barplot_bysize <- ggplot(froh_merged_df_groups, aes(x=groups,y=value,fill=ROH_type)) +
  geom_boxplot(position = "dodge") + 
  geom_text(label="n=0",x=1.75, y=2.7) +
  labs(y=expression(F[ROH]),x="") +
  scale_fill_manual(values=c("seashell3","seashell2","seashell")) +
  scale_x_discrete(labels=c("US101 & I405","All Other Roads")) +
  theme_pubr() +
  geom_signif(map_signif_level=TRUE,
              annotations = c("n.s.","n.s."),
              y_position = c(80, 90),
              xmin=c(1,1.3),
              xmax=c(2,2.3)) +
  theme(legend.position = "none") 

ggsave(plot=Froh_2pop_barplot_bysize, "20240718_FROH_405and101comparisonplot_bysize.png",device="png",dpi=300)

## The last plot I want to make uses bins based on generation time

# Set genome-wide recombination rate - This is for the chicken following Robinson et al. 2021, 
# "Genome-wide diversity in the California condor tracks its prehistoric abundance and decline" and based on a 
# formula in Stoffel et al. 2021, "Genetic architecture and lifetime dynamics of inbreeding depression in a wild mammal"

# Set WD 

setwd("~/Downloads/ROH_2024/")

# Define Variables 

recomb = 2.7

# Determine bin sizes in Mb

g1 = (100/(2*1))/recomb
#18.51852

g2 = (100/(2*2))/recomb
#9.259259

g4 = (100/(2*4))/recomb
#4.62963

g8 = (100/(2*8))/recomb
#2.314815                 

g16 = (100/(2*16))/recomb
#1.157407

g32 = (100/(2*32))/recomb
#0.5787037

# Re-define Functions

read_filter_roh_g <- function(data, min_roh_length){
  output_g <- read.table(paste(data,".gz",sep=""), col.names=c("row_type","sample","chrom","start","end","length","num_markers","qual"), fill=T)
  output1_g <- subset(output_g, row_type == "RG")
  output_class_sums_g <- classify_roh_g(output1_g, min_roh_length=min_roh_length)
  return(output_class_sums_g)
}

classify_roh_g <- function(roh_dataframe, min_roh_length){
  roh_g2 <- subset(roh_dataframe, length>=9259259 & length<18518520)
  roh_g2_g4 <- subset(roh_dataframe, length>=4629630 & length<9259259)
  roh_g4_g8 <-  subset(roh_dataframe, length>=2314815 & length<4629630)
  roh_g8_g16 <-  subset(roh_dataframe, length>=1157407 & length<2314815)
  roh_g16_g32 <-  subset(roh_dataframe, length>=578704 & length<1157407)
  roh_g32_min <- subset(roh_dataframe, length>min_roh_length & length<578704)
  
  #sum each class and divide by 1000000 to get Mb, then divide by genome size to get % genome 
  Froh_g2 <- (sum(roh_g2$length)/1000000)/genome_length
  Froh_g2_g4 <- (sum(roh_g2_g4$length)/1000000)/genome_length
  Froh_g4_g8 <- (sum(roh_g4_g8$length)/1000000)/genome_length
  Froh_g8_g16 <- (sum(roh_g8_g16$length)/1000000)/genome_length
  Froh_g16_g32 <- (sum(roh_g16_g32$length)/1000000)/genome_length
  Froh_g32_min <- (sum(roh_g32_min$length)/1000000)/genome_length
  
  return(c(Froh_g2, Froh_g2_g4, Froh_g4_g8, Froh_g8_g16, Froh_g16_g32, Froh_g32_min))
  
}

# initialize a new data frame
roh_size_df_g <- data.frame(matrix(nrow=6, ncol=length(individuals)))
colnames(roh_size_df_g) <- individuals

for(i in 1:length(individuals)){
  roh_size_df_g[,i] <- read_filter_roh_g(paste(file,individuals[i], sep=""), min_roh_length=min_roh_length)
}


# Tidy up the data frame 
roh_size_df_g_tidy <- roh_size_df_g %>%  
  pivot_longer(cols = everything(), names_to = "Sample_name", values_to = "value") %>%
  group_by(Sample_name) %>%
  mutate(id = row_number()) %>%
  pivot_wider(names_from = id, values_from = value)

colnames(roh_size_df_g_tidy) <- c("Sample_name", "roh_g2", "roh_g2_g4", "roh_g4_g8","roh_g8_g16","roh_g16_g32","roh_g32_min")

# oops, barplot needs long data... 

roh_size_df_g_tidy_long <- roh_size_df_g_tidy %>%
  pivot_longer(cols = c("roh_g2", "roh_g2_g4", "roh_g4_g8","roh_g8_g16","roh_g16_g32","roh_g32_min"), 
               names_to = "ROH_time", 
               values_to = "value")

# reorder the levels for the ROH_time variable

roh_size_df_g_tidy_long$ROH_time <- factor(roh_size_df_g_tidy_long$ROH_time, 
                                         levels=c("roh_g2", "roh_g2_g4", "roh_g4_g8","roh_g8_g16","roh_g16_g32","roh_g32_min"), 
                                         ordered = T)

# Get counts per ROH_time

ROH_size_summary <- roh_size_df_g_tidy_long %>%
  group_by(ROH_time, value) %>%
  summarise(n=n())


# Plot it

roh_age_plot <- ggplot(roh_size_df_g_tidy_long, aes(x=ROH_time, y=value*100, fill = ROH_time)) +
  geom_violin(alpha=.7) +
  geom_point(position = "jitter", color=ifelse(roh_size_df_g_tidy_long$ROH_time == "roh_g1", "white","grey60"),alpha=.7) +
  theme_pubr() +
  theme(axis.text.x = element_text(angle=45, vjust = .55),
        legend.position = "none") +
  stat_summary(fun = "mean",
               geom = "crossbar", 
               colour = "grey30",
               width=c(0,.34,.27,.43,.65,.78)) +
  scale_x_discrete(labels=c(">=18.52\n(2g)","9.26-4.63\n(2-4g)","4.63-2.31\n(4-8g)","2.31-1.16\n(8-16g)", "1.16-0.58\n(16-32g)","<=0.58\n(>=32g)")) +
  labs(x="ROH length in Mb (~ generations to MRCA)", y="% genome")

ggsave(plot=roh_age_plot,"20240721_CAQU_ROHagePlot.png",device="png",dpi=300)

## now out of curiosity, facet by 405/101 and everything else

roh_size_df_g_tidy_long_merged <- merge(roh_size_df_g_tidy_long,pop_data,by.x='Sample_name',by.y='individuals') # get population information

roh_size_df_g_tidy_long_merged$groups <- NA # initialize empty column

roh_size_df_g_tidy_long_merged_groups <- roh_size_df_g_tidy_long_merged %>% #Add in grouping variable for ggplot
  mutate(groups = case_when((locality == "S_US101") ~ "I405 & US101",
                            (locality == "N_US101") ~ "I405 & US101",
                            (locality == "W_I405") ~ "I405 & US101",
                            (locality == "E_I405") ~ "I405 & US101",
                            .default = "All Other Roads")) 

roh_size_df_g_tidy_long_merged_groups$groups <- factor(roh_size_df_g_tidy_long_merged_groups$groups,
                                                       levels=c("I405 & US101","All Other Roads"), ordered = T)


# plot it

roh_age_plot_groups <- ggplot(roh_size_df_g_tidy_long_merged_groups, aes(x=ROH_time, y=value*100,fill=groups)) +
  geom_violin(alpha=.7, position="dodge",show.legend = F) +
  geom_point(position = position_jitterdodge(dodge.width =.9), color="grey60",alpha=.7) +
  stat_summary(fun = "mean",
               geom = "crossbar", 
               colour = c("grey30","white","grey30","grey30","grey30","grey30","grey30","grey30","grey30","grey30","grey30","grey30"),
               width=c(.1,.1,.3,.3,.18,.18,.24,.24,.45,.45,.5,.5),
               position = position_dodge(width = .9)) +
  geom_signif(map_signif_level=TRUE,
              annotations = c("*"),
              y_position = c(3.2),
              xmin = c(5.75),
              xmax = c(6.2),
              fontface="bold",
              tip_length = .025,
              textsize = 6,
              vjust = .5) +
  scale_x_discrete(labels=c("≥18.52\n(≤2g)","9.26-4.63\n(2-4g]","4.63-2.31\n(4-8g]","2.31-1.16\n(8-16g]", "1.16-0.58\n(16-32g]","<0.58\n(>32g)")) +
  scale_fill_manual(values=c("seashell4","seashell2"),name="") +
  labs(x="\nROH length in Mb (~ generations to MRCA)", y="% genome\n") +
  theme_pubclean() +
  theme(axis.text.x = element_text(angle=45, vjust = .55),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.position = "right",
        legend.key = element_rect(fill = NA, color = NA, linewidth = .1)) 


ggsave(plot=roh_age_plot_groups,file="CAQU_roh_age_plot_groups_20240903.png",device="png",dpi=300, height=5,width = 8)


## Now output CSV files for samples for I405 and US101 samples seperately. Will use these in next script
# to test for differences between these groups.

#2-4g
ttest_g2_g4_405and101 <- roh_size_df_g_tidy_long_merged_groups %>%
  filter(groups == "I405 & US101" & ROH_time == "roh_g2_g4") %>%
  pull(value)

write.csv(ttest_g2_g4_405and101,file="g2_g4_405and101.csv",row.names = F)

ttest_g2_g4_405and101_df <- data.frame(ttest_g2_g4_405and101)
ttest_g2_g4_405and101_df$group <- "I405 & US101"
colnames(ttest_g2_g4_405and101_df) <- c("roh","location")

ttest_g2_g4_otherroads <- roh_size_df_g_tidy_long_merged_groups %>%
  filter(groups == "All Other Roads" & ROH_time == "roh_g2_g4") %>%
  pull(value)

write.csv(ttest_g2_g4_otherroads,file="g2_g4_otherroads.csv",row.names = F)

#4-8g
ttest_g4_g8_405and101 <- roh_size_df_g_tidy_long_merged_groups %>%
  filter(groups == "I405 & US101" & ROH_time == "roh_g4_g8") %>%
  pull(value)

write.csv(ttest_g4_g8_405and101,file="g4_g8_405and101.csv",row.names = F)


ttest_g4_g8_otherroads <- roh_size_df_g_tidy_long_merged_groups %>%
  filter(groups == "All Other Roads" & ROH_time == "roh_g4_g8") %>%
  pull(value)

write.csv(ttest_g4_g8_otherroads,file="g4_g8_otherroads.csv",row.names = F)

#8-16g

ttest_g8_g16_405and101 <- roh_size_df_g_tidy_long_merged_groups %>%
  filter(groups == "I405 & US101" & ROH_time == "roh_g8_g16") %>%
  pull(value)

write.csv(ttest_g8_g16_405and101,file="g8g16_405and101.csv",row.names = F)

ttest_g8_g16_otherroads <- roh_size_df_g_tidy_long_merged_groups %>%
  filter(groups == "All Other Roads" & ROH_time == "roh_g8_g16") %>%
  pull(value)

write.csv(ttest_g8_g16_otherroads,file="g8g16_otherroads.csv",row.names = F)

#16-32g

ttest_g16_g32_405and101 <- roh_size_df_g_tidy_long_merged_groups %>%
  filter(groups == "I405 & US101" & ROH_time == "roh_g16_g32") %>%
  pull(value)

write.csv(ttest_g16_g32_405and101,file="g16g32_405and101.csv",row.names = F)

ttest_g16_g32_otherroads <- roh_size_df_g_tidy_long_merged_groups %>%
  filter(groups == "All Other Roads" & ROH_time == "roh_g16_g32") %>%
  pull(value)

write.csv(ttest_g16_g32_otherroads,file="g16g32_otherroads.csv",row.names = F)

#>32g

ttest_g32_min_405and101 <- roh_size_df_g_tidy_long_merged_groups %>%
  filter(groups == "I405 & US101" & ROH_time == "roh_g32_min") %>%
  pull(value)

write.csv(ttest_g32_min_405and101,file="g32_405and101.csv",row.names = F)

ttest_g32_min_otherroads <- roh_size_df_g_tidy_long_merged_groups %>%
  filter(groups == "All Other Roads" & ROH_time == "roh_g32_min") %>%
  pull(value)

write.csv(ttest_g32_min_otherroads,file="g32_otherroads.csv",row.names = F)

  
