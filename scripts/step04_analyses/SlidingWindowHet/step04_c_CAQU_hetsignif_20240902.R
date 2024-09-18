#######################################################################################
###### R script to test for differences in mean het between sampled populations #######
#######################################################################################

# Author: Joey Curti 
# Date: MON SEP 02 2024
# Description: test for significance in mean_het between populations
# References:
# https://cran.r-project.org/web/packages/afex/vignettes/assumptions_of_ANOVAs.html
# http://www.sthda.com/english/wiki/kruskal-wallis-test-in-r

## Clean Workspace and set working directory

setwd("~/Downloads/")
rm(list = ls())

## Dependencies

library(afex)
library(performance)
library(qqplotr)
library(rstatix)


## Read in data 

data <- read.csv(<insert path to mean_het.csv>,stringsAsFactors = F, header = T, sep = ",")
data_tidy <- data[-1,c(-28:-39)]

## Main 

# Quick plot to visualize differences 

pp1 <- ggplot(data=data_tidy, aes(x=state_site, y=mean_het, fill=state_site)) + 
  geom_boxplot() +
  theme_pubclean() +
  scale_fill_manual(values=c("seashell","seashell2","seashell3")) +
  labs(x="Sampling Location", y = expression(paste("Heterozygosity (", bp^{-1},")")))

# Specify model
o1 <- aov_ez("individuals", "mean_het", data_tidy, between = "state_site")

# Response: mean_het
# Effect    df  MSE         F  ges p.value
# 1 state_site 2, 58 0.00 22.23 *** .434   <.001
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘+’ 0.1 ‘ ’ 1

# Check assumptions
check_homogeneity(o1)
#OK: There is not clear evidence for different variances across groups (Levene's Test, p = 0.439).
is_norm <- check_normality(o1)
plot(is_norm)
#Warning: Non-normality of residuals detected (p = 0.036).

# Given this violation of assumptions, I am going to go with a Kruskal-Wallace instead

kw_test <- kruskal.test(mean_het ~ state_site, data = data_tidy)

# Kruskal-Wallis rank sum test
# 
# data:  mean_het by state_site
# Kruskal-Wallis chi-squared = 10.758, df = 2, p-value = 0.004612

pwc <-  dunn_test(data=data_tidy, mean_het ~ state_site, p.adjust.method = "bonferroni") 

# # A tibble: 3 × 9
# .y.      group1 group2    n1    n2 statistic       p  p.adj p.adj.signif
# * <chr>    <chr>  <chr>  <int> <int>     <dbl>   <dbl>  <dbl> <chr>       
#   1 mean_het NORCAL SAMO       3    55     2.92  0.00354 0.0106 *           
#   2 mean_het NORCAL SOCAL      3     3     0.920 0.358   1      ns          
# 3 mean_het SAMO   SOCAL     55     3    -1.65  0.0990  0.297  ns 

## Finalize figure with significance values

pp2 <- ggplot(data=data_tidy, aes(x=state_site, y=mean_het, fill=state_site)) + 
  geom_boxplot(show.legend = F) +
  theme_pubclean() +
  scale_fill_manual(values=c("seashell","seashell2","seashell3"),name="") +
  labs(x="\nSampling Location", y = expression(paste("Heterozygosity (", bp^{-1},")"))) +
  geom_signif(annotations = c("*"),
              y_position = .00551,
              xmin = 1,
              xmax =2, 
              textsize = 6,
              fontface="bold",
              vjust = .1)

ggsave("20240902_CAQU_hetsignif.png",plot=pp2,dpi=400,width = 6,height = 5)

