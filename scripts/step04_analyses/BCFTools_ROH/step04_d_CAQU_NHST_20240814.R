####################################################################################
## Ranked NHST with and without replacement to calculate significance of ROH data ##
####################################################################################

# Author: Joseph Curti (jcurti3@g.ucla.edu)
# Adapted From: Morgan Tingley (mtingley@g.ucla.edu)
# Description: Bootstrapping p-values using ranked data
# Version: V1
# Date: WED AUG 24 2024
# Note: All files that are loaded into this script are generated in the script step04_c_CAQU_BCFTools_ROH_plot_20240918

## References

## Clean Workspace, Set WD 

rm(list = ls())
setwd(<insert working directory path>)

## Dependencies
library(patchwork)

### Run for g2 - g4

## Load in data

g2g4_101and405 <- unlist(read.csv(<insert path to file 'g2_g4_405and101.csv'>,header = T, row.names = NULL, stringsAsFactors = F),use.names=FALSE)
g2g4_other <- unlist(read.csv(<insert path to file 'g2_g4_otherroads.csv'>,header = T, row.names = NULL, stringsAsFactors = F),use.names=FALSE)

## Convert to ranks

g2g4_101and405.rank <- rank(c(g2g4_101and405, g2g4_other))[1:length(g2g4_101and405)]
g2g4_other.rank <- rank(c(g2g4_101and405, g2g4_other))[(length(g2g4_101and405) + 1):(length(c(g2g4_101and405,g2g4_other)))]
g2g4.obs <- mean(g2g4_other.rank) - mean(g2g4_101and405.rank)

## Big box
g2g4_big.box <- rank(c(g2g4_101and405.rank, g2g4_other.rank))

## Randomly sample *with replacement* both groups from the box box
g2g4.replace <- NA
for(i in 1:10000) {
  sim.g2g4_101and405 <- sample(x = g2g4_big.box, size = length(g2g4_101and405), replace = TRUE)
  sim.g2g4_other <- sample(x = g2g4_big.box, size = length(g2g4_other), replace = TRUE) 
  g2g4.replace[i] <- mean(sim.g2g4_other) - mean(sim.g2g4_101and405)
}
g2g4_p.val <- round(sum(abs(g2g4.replace) >= abs(g2g4.obs)) / 10000, 3)
print(g2g4_p.val)

## Plot
g2g4_breakpt <- round(max(abs(range(g2g4.replace))), 0) + 5
hist_g2g4 <- hist(g2g4.replace, border = F, main = "2-4g", col = "#00AA20AA",
     breaks = seq(-1*g2g4_breakpt, g2g4_breakpt, by = 1), xlab = "difference in mean rank", xlim = c(-1*(g2g4_breakpt-5), g2g4_breakpt - 5))
abline(v = quantile(g2g4.replace, probs = c(0.025, 0.975)), col = "dark green", lwd = 2, lty = 2)
abline(v = g2g4.obs, lwd = 3, lty = 1, col = "dark blue")
text((g2g4_breakpt - 7), min(hist_g2g4$counts) + 100, labels = paste0("p = ", g2g4_p.val))

### Repeat for other binned categories 

## g4 - g8

## Load in data

g4g8_101and405 <- unlist(read.csv(<insert path to file 'g4_g8_405and101.csv'>,header = T, row.names = NULL, stringsAsFactors = F),use.names=FALSE)
g4g8_other <- unlist(read.csv(<insert path to file 'g4_g8_otherroads.csv'>,header = T, row.names = NULL, stringsAsFactors = F),use.names=FALSE)

## Convert to ranks

g4g8_101and405.rank <- rank(c(g4g8_101and405, g4g8_other))[1:length(g4g8_101and405)]
g4g8_other.rank <- rank(c(g4g8_101and405, g4g8_other))[(length(g4g8_101and405) + 1):(length(c(g4g8_101and405,g4g8_other)))]
g4g8.obs <- mean(g4g8_other.rank) - mean(g4g8_101and405.rank)

## Big box
g4g8_big.box <- rank(c(g4g8_101and405.rank, g4g8_other.rank))

## Randomly sample *with replacement* both groups from the box box
g4g8.replace <- NA
for(i in 1:10000) {
  sim.g4g8_101and405 <- sample(x = g4g8_big.box, size = length(g4g8_101and405), replace = TRUE)
  sim.g4g8_other <- sample(x = g4g8_big.box, size = length(g4g8_other), replace = TRUE) 
  g4g8.replace[i] <- mean(sim.g4g8_other) - mean(sim.g4g8_101and405)
}
g4g8_p.val <- round(sum(abs(g4g8.replace) >= abs(g4g8.obs)) / 10000, 3)
print(g4g8_p.val)

## Plot
g4g8_breakpt <- round(max(abs(range(g4g8.replace))), 0) + 5
hist_g4g8 <- hist(g4g8.replace, border = F, main = "4-8g", col = "#00AA20AA",
                  breaks = seq(-1*g4g8_breakpt, g4g8_breakpt, by = 1), xlab = "difference in mean rank", xlim = c(-1*(g4g8_breakpt-5), g4g8_breakpt - 5))
abline(v = quantile(g4g8.replace, probs = c(0.025, 0.975)), col = "dark green", lwd = 2, lty = 2)
abline(v = g4g8.obs, lwd = 3, lty = 1, col = "dark blue")
text((g4g8_breakpt - 7), min(hist_g4g8$counts) + 100, labels = paste0("p = ", g4g8_p.val))

## 8 - 16 g

## Load in data

g8g16_101and405 <- unlist(read.csv(<insert path to file 'g8g16_405and101.csv'>,header = T, row.names = NULL, stringsAsFactors = F),use.names=FALSE)
g8g16_other <- unlist(read.csv(<insert path to file 'g8g16_otherroads.csv'>,header = T, row.names = NULL, stringsAsFactors = F),use.names=FALSE)

## Convert to ranks

g8g16_101and405.rank <- rank(c(g8g16_101and405, g8g16_other))[1:length(g8g16_101and405)]
g8g16_other.rank <- rank(c(g8g16_101and405, g8g16_other))[(length(g8g16_101and405) + 1):(length(c(g8g16_101and405,g8g16_other)))]
g8g16.obs <- mean(g8g16_other.rank) - mean(g8g16_101and405.rank)

## Big box
g8g16_big.box <- rank(c(g8g16_101and405.rank, g8g16_other.rank))

## Randomly sample *with replacement* both groups from the box box
g8g16.replace <- NA
for(i in 1:10000) {
  sim.g8g16_101and405 <- sample(x = g8g16_big.box, size = length(g8g16_101and405), replace = TRUE)
  sim.g8g16_other <- sample(x = g8g16_big.box, size = length(g8g16_other), replace = TRUE) 
  g8g16.replace[i] <- mean(sim.g8g16_other) - mean(sim.g8g16_101and405)
}
g8g16_p.val <- round(sum(abs(g8g16.replace) >= abs(g8g16.obs)) / 10000, 3)
print(g8g16_p.val)

## Plot
g8g16_breakpt <- round(max(abs(range(g8g16.replace))), 0) + 5
hist_g8g16 <- hist(g8g16.replace, border = F, main = "8-16g", col = "#00AA20AA",
                 breaks = seq(-1*g8g16_breakpt, g8g16_breakpt, by = 1), xlab = "difference in mean rank", xlim = c(-1*(g8g16_breakpt-5), g8g16_breakpt - 5))
abline(v = quantile(g8g16.replace, probs = c(0.025, 0.975)), col = "dark green", lwd = 2, lty = 2)
abline(v = g8g16.obs, lwd = 3, lty = 1, col = "dark blue")
text((g8g16_breakpt - 7), min(hist_g8g16$counts) + 100, labels = paste0("p = ", g8g16_p.val))

## 16 - 32 g

## Load in data

g16g32_101and405 <- unlist(read.csv(<insert path to file 'g16g32_405and101.csv'>,header = T, row.names = NULL, stringsAsFactors = F),use.names=FALSE)
g16g32_other <- unlist(read.csv(<insert path to file 'g16g32_otherroads.csv'>,header = T, row.names = NULL, stringsAsFactors = F),use.names=FALSE)

## Convert to ranks

g16g32_101and405.rank <- rank(c(g16g32_101and405, g16g32_other))[1:length(g16g32_101and405)]
g16g32_other.rank <- rank(c(g16g32_101and405, g16g32_other))[(length(g16g32_101and405) + 1):(length(c(g16g32_101and405,g16g32_other)))]
g16g32.obs <- mean(g16g32_other.rank) - mean(g16g32_101and405.rank)

## Big box
g16g32_big.box <- rank(c(g16g32_101and405.rank, g16g32_other.rank))

## Randomly sample *with replacement* both groups from the box box
g16g32.replace <- NA
for(i in 1:10000) {
  sim.g16g32_101and405 <- sample(x = g16g32_big.box, size = length(g16g32_101and405), replace = TRUE)
  sim.g16g32_other <- sample(x = g16g32_big.box, size = length(g16g32_other), replace = TRUE) 
  g16g32.replace[i] <- mean(sim.g16g32_other) - mean(sim.g16g32_101and405)
}
g16g32_p.val <- round(sum(abs(g16g32.replace) >= abs(g16g32.obs)) / 10000, 3)
print(g16g32_p.val)

## Plot

g16g32_breakpt <- round(max(abs(range(g16g32.replace))), 0) + 5
hist_g16g32 <- hist(g16g32.replace, border = F, main = "16-32g", col = "#00AA20AA",breaks = seq(-1*g16g32_breakpt, g16g32_breakpt, by = 1), xlab = "difference in mean rank", xlim = c(-1*(g16g32_breakpt-5), g16g32_breakpt - 5))
hist_g16g32_2 <- abline(v = quantile(g16g32.replace, probs = c(0.025, 0.975)), col = "dark green", lwd = 2, lty = 2)
hist_g16g32_2 <- abline(v = g16g32.obs, lwd = 3, lty = 1, col = "dark blue") 
hist_g16g32_2 <- text((g16g32_breakpt - 7), min(hist_g16g32$counts) + 100, labels = paste0("p = ", g16g32_p.val))

## g32 and greater

## Load in data

g32_101and405 <- unlist(read.csv(<insert path to file 'g32_405and101.csv'>,header = T, row.names = NULL, stringsAsFactors = F),use.names=FALSE)
g32_other <- unlist(read.csv(<insert path to file 'g32_otherroads.csv'>,header = T, row.names = NULL, stringsAsFactors = F),use.names=FALSE)

## Convert to ranks

g32_101and405.rank <- rank(c(g32_101and405, g32_other))[1:length(g32_101and405)]
g32_other.rank <- rank(c(g32_101and405, g32_other))[(length(g32_101and405) + 1):(length(c(g32_101and405,g32_other)))]
g32.obs <- mean(g32_other.rank) - mean(g32_101and405.rank)

## Big box
g32_big.box <- rank(c(g32_101and405.rank, g32_other.rank))

## Randomly sample *with replacement* both groups from the box box
g32.replace <- NA
for(i in 1:10000) {
  sim.g32_101and405 <- sample(x = g32_big.box, size = length(g32_101and405), replace = TRUE)
  sim.g32_other <- sample(x = g32_big.box, size = length(g32_other), replace = TRUE) 
  g32.replace[i] <- mean(sim.g32_other) - mean(sim.g32_101and405)
}
g32_p.val <- round(sum(abs(g32.replace) >= abs(g32.obs)) / 10000, 3)
print(g32_p.val)

## Plot

g32_breakpt <- round(max(abs(range(g32.replace))), 0) + 5
hist_g32 <- hist(g32.replace, border = F, main = "≥32g", col = "#00AA20AA",breaks = seq(-1*g32_breakpt, g32_breakpt, by = 1), xlab = "difference in mean rank", xlim = c(-1*(g32_breakpt-5), g32_breakpt - 5))
hist_g32_2 <- abline(v = quantile(g32.replace, probs = c(0.025, 0.975)), col = "dark green", lwd = 2, lty = 2)
hist_g32_2 <- abline(v = g32.obs, lwd = 3, lty = 1, col = "dark blue") 
hist_g32_2 <- text((g32_breakpt - 7), min(hist_g32$counts) + 100, labels = paste0("p = ", g32_p.val))

# save multiplot
par(mar=c(1,1,1,1))
tiff(filename="CAQU_nhst_20240838.tiff",res = 400,width = 10, height = 10, units = 'in')
par(mfrow=c(3,2))
# Redo all plots, then:
dev.off()

## Plot
g2g4_breakpt <- round(max(abs(range(g2g4.replace))), 0) + 5
hist_g2g4 <- hist(g2g4.replace, border = F, main = "2-4g", col = "#00AA20AA",
                  breaks = seq(-1*g2g4_breakpt, g2g4_breakpt, by = 1), xlab = "difference in mean rank", xlim = c(-1*(g2g4_breakpt-5), g2g4_breakpt - 5))
abline(v = quantile(g2g4.replace, probs = c(0.025, 0.975)), col = "dark green", lwd = 2, lty = 2)
abline(v = g2g4.obs, lwd = 3, lty = 1, col = "dark blue")
text((g2g4_breakpt - 7), min(hist_g2g4$counts) + 100, labels = paste0("p = ", g2g4_p.val))

g4g8_breakpt <- round(max(abs(range(g4g8.replace))), 0) + 5
hist_g4g8 <- hist(g4g8.replace, border = F, main = "4-8g", col = "#00AA20AA",
                  breaks = seq(-1*g4g8_breakpt, g4g8_breakpt, by = 1), xlab = "difference in mean rank", xlim = c(-1*(g4g8_breakpt-5), g4g8_breakpt - 5))
abline(v = quantile(g4g8.replace, probs = c(0.025, 0.975)), col = "dark green", lwd = 2, lty = 2)
abline(v = g4g8.obs, lwd = 3, lty = 1, col = "dark blue")
text((g4g8_breakpt - 7), min(hist_g4g8$counts) + 100, labels = paste0("p = ", g4g8_p.val))

g8g16_breakpt <- round(max(abs(range(g8g16.replace))), 0) + 5
hist_g8g16 <- hist(g8g16.replace, border = F, main = "8-16g", col = "#00AA20AA",
                   breaks = seq(-1*g8g16_breakpt, g8g16_breakpt, by = 1), xlab = "difference in mean rank", xlim = c(-1*(g8g16_breakpt-5), g8g16_breakpt - 5))
abline(v = quantile(g8g16.replace, probs = c(0.025, 0.975)), col = "dark green", lwd = 2, lty = 2)
abline(v = g8g16.obs, lwd = 3, lty = 1, col = "dark blue")
text((g8g16_breakpt - 7), min(hist_g8g16$counts) + 100, labels = paste0("p = ", g8g16_p.val))

g16g32_breakpt <- round(max(abs(range(g16g32.replace))), 0) + 5
hist_g16g32 <- hist(g16g32.replace, border = F, main = "16-32g", col = "#00AA20AA",breaks = seq(-1*g16g32_breakpt, g16g32_breakpt, by = 1), xlab = "difference in mean rank", xlim = c(-1*(g16g32_breakpt-5), g16g32_breakpt - 5))
hist_g16g32_2 <- abline(v = quantile(g16g32.replace, probs = c(0.025, 0.975)), col = "dark green", lwd = 2, lty = 2)
hist_g16g32_2 <- abline(v = g16g32.obs, lwd = 3, lty = 1, col = "dark blue") 
hist_g16g32_2 <- text((g16g32_breakpt - 7), min(hist_g16g32$counts) + 100, labels = paste0("p = ", g16g32_p.val))

g32_breakpt <- round(max(abs(range(g32.replace))), 0) + 5
hist_g32 <- hist(g32.replace, border = F, main = "≥32g", col = "#00AA20AA",breaks = seq(-1*g32_breakpt, g32_breakpt, by = 1), xlab = "difference in mean rank", xlim = c(-1*(g32_breakpt-5), g32_breakpt - 5))
hist_g32_2 <- abline(v = quantile(g32.replace, probs = c(0.025, 0.975)), col = "dark green", lwd = 2, lty = 2)
hist_g32_2 <- abline(v = g32.obs, lwd = 3, lty = 1, col = "dark blue") 
hist_g32_2 <- text((g32_breakpt - 7), min(hist_g32$counts) + 100, labels = paste0("p = ", g32_p.val))

