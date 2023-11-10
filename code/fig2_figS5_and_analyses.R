# Description ----
# Fig. 2 and Fig. S5: Mortality rate vs. size

# Citation:
# Borrego, I., Perez, T.M., Bentley, L.P., Bison, N.N., Byrnes, L, Candido, H.G., Chmurzynski, A. Durán, S.M., 
# Fox, T.S., Gaitan, M.A., Garen, J.C., Orwig, D.A., Pau, S. Scott, J.L., Simovic, M. Swenson, N.G., Wieczynski, D.J.,
# Buzzard, V. Enquist, B.J., & Michaletz, S.T. (in review) From tropics to treeline: extending and assessing 
# metabolic scaling theory for global variation in plant mortality rates.

# July 2023, updated November 2023
# Sean Michaletz, sean.michaletz@gmail.com
# Based on earlier code from Isaac Borrego, Sean Michaletz


# Initialize ----
#--Libraries
library(ggpmisc)
library(ggplot2)
library(ggpubr)
library(grid)
library(lmodel2)

# Load data for mortality rate, size, temperature, and traits (file from Isaac, JUL 2023)
mort_data <- read.csv("https://raw.githubusercontent.com/MichaletzLab/mortality-scaling/main/data/mortality_data.csv")
# Load site climate data
climate <- read.csv("https://raw.githubusercontent.com/MichaletzLab/mortality-scaling/main/data/climate.csv",header=T)

# Rename DBH size class from 'size' to 'dbh_cm'
names(mort_data)[names(mort_data) == 'size'] <- 'dbh_cm'
# Rename leaf P and N to avoid confusion with C
names(mort_data)[names(mort_data) == 'lpc'] <- 'P_percent'
names(mort_data)[names(mort_data) == 'lnc'] <- 'N_percent'

# Calculate leaf P:N
mort_data$PN_ratio <- mort_data$P_percent / mort_data$N_percent

# Add growing season length, MAT, and GST data to mort_data
mort_data$L_gs_mo <- climate$L_gs_mo[match(mort_data$site, climate$Site)]
mort_data$MAT_C <- climate$MAT_C[match(mort_data$site, climate$Site)]
mort_data$GST_C <- climate$GST_C[match(mort_data$site, climate$Site)]

# Calculate monthly mortality rate (correct annual rates for growing season length)
mort_data$mort_rate_mo <- mort_data$avg_mort / mort_data$L_gs_mo


# Individual sites ----
# Plotting RMA and 95% CI from lmodel2 in ggplot2
# Updated using stat_ma_line() from package ggpmisc, Sean Michaletz, 6 Dec 2022

# create subsets for each site
AND <- subset(mort_data, site == "AND" & !(is.na(avg_mort)))
BCI <- subset(mort_data, site == "BCI" & !(is.na(avg_mort)))
CAP <- subset(mort_data, site == "CAP" & !(is.na(avg_mort)))
CWT <- subset(mort_data, site == "CWT" & !(is.na(avg_mort)))
HFR <- subset(mort_data, site == "HFR" & !(is.na(avg_mort)))
LUQ <- subset(mort_data, site == "LUQ" & !(is.na(avg_mort)))
NWT <- subset(mort_data, site == "NWT" & !(is.na(avg_mort)))
PSR <- subset(mort_data, site == "PSR" & !(is.na(avg_mort)))

## AND ----
sma_size_AND <- lmodel2(log10(mort_rate_mo) ~ log10(dbh_cm), "interval", "interval", data=AND, nperm=99)
sma_size_AND
# report p-value
sma_size_AND$P.param

# Plot the data, SMA, and 95% CI using plot.lmodel2(). (Note that this function is called using just plot() below). 
# The code for this function is at https://github.com/vegandevs/lmodel2/blob/master/R/lmodel2.R#L229.
# NOTE: The plotted 95% CI do not simply use the 2.5% and 97.5% slopes and intercepts given in *$confidence.intervals.  
# Rather, the reported 2.5% and 97.5% intercepts are corrected using the mean x and y values from
# the original fitted data, along with the reported 2.5% and 97.5% slopes. For example, the plotted lower intercept
# is calculated as intercept_lower_corrected = mean(y) + slope_lower * mean(x). For details, see 
# https://github.com/vegandevs/lmodel2/blob/master/R/lmodel2.R#L287 and https://github.com/vegandevs/lmodel2/blob/master/R/lmodel2.R#L289.

# Plot the data, SMA, and 95% CI with ggplot2 using stat_ma_line() from ggpmisc package.
# See https://rdrr.io/cran/ggpmisc/man/stat_ma_line.html
plot.AND <- ggplot(AND) + 
  stat_ma_line(aes(x = dbh_cm, y = mort_rate_mo), method="SMA", color = "black", linewidth = 0) + 
  geom_errorbar(aes(x = dbh_cm, ymin = lower, ymax = upper), width = 0.04) + 
  geom_point(aes(x = dbh_cm, y = mort_rate_mo), fill = "#92C5DE", pch=21, size=3.5) +
  stat_ma_line(aes(x = dbh_cm, y = mort_rate_mo), method="SMA", color = "black", linewidth = 1, se=F) + 
  xlab("Diameter at breast height (cm)") +
  ylab(expression(paste('Mortality rate (ind ', ~ind^-1, ~mo^-1, ')'))) +
  coord_cartesian(xlim = c(0.3, 25), ylim = c(0.00001, 0.11)) +
  scale_x_continuous(trans="log10", breaks = c(1, 10, 30)) +
  scale_y_continuous(trans = "log10", breaks = c(0.0001, 0.001, 0.010, 0.100)) +
  theme_bw()
plot.AND

## BCI ----
sma_size_BCI <- lmodel2(log10(mort_rate_mo) ~ log10(dbh_cm), "interval", "interval", data=BCI, nperm=99)
sma_size_BCI
# report p-value
sma_size_BCI$P.param

plot.BCI <- ggplot(BCI) + 
  stat_ma_line(aes(x = dbh_cm, y = mort_rate_mo), method="SMA", color = "black", linewidth = 0) + 
  geom_errorbar(aes(x = dbh_cm, ymin = lower, ymax = upper), width = 0.04) + 
  geom_point(aes(x = dbh_cm, y = mort_rate_mo), fill = "#B2182B", pch=21, size=3.5) +
  stat_ma_line(aes(x = dbh_cm, y = mort_rate_mo), method="SMA", color = "black", linewidth = 1, se=F) + 
  xlab("Diameter at breast height (cm)") +
  ylab(expression(paste('Mortality rate (ind ', ~ind^-1, ~mo^-1, ')'))) +
  coord_cartesian(xlim = c(0.3, 25), ylim = c(0.00001, 0.11)) + 
  scale_x_continuous(trans="log10", breaks = c(1, 10, 30)) +
  scale_y_continuous(trans = "log10", breaks = c(0.0001, 0.001, 0.010, 0.100)) +
  theme_bw()
plot.BCI

## CAP ----
sma_size_CAP <- lmodel2(log10(mort_rate_mo) ~ log10(dbh_cm), "interval", "interval", data=CAP, nperm=99) 
sma_size_CAP
# report p-value
sma_size_CAP$P.param

plot.CAP <- ggplot(CAP) + 
  stat_ma_line(aes(x = dbh_cm, y = mort_rate_mo), method="SMA", color = "black", linewidth = 0) + 
  geom_errorbar(aes(x = dbh_cm, ymin = lower, ymax = upper), width = 0.04) + 
  geom_point(aes(x = dbh_cm, y = mort_rate_mo), fill = "#4393C3", pch=21, size=3.5) +
  stat_ma_line(aes(x = dbh_cm, y = mort_rate_mo), method="SMA", color = "black", linewidth = 1, se=F) + 
  xlab("Diameter at breast height (cm)") +
  ylab(expression(paste('Mortality rate (ind ', ~ind^-1, ~mo^-1, ')'))) +
  coord_cartesian(xlim = c(0.3, 25), ylim = c(0.00001, 0.11)) + 
  scale_x_continuous(trans="log10", breaks = c(1, 10, 30)) +
  scale_y_continuous(trans = "log10", breaks = c(0.0001, 0.001, 0.010, 0.100)) +
  theme_bw()
plot.CAP

## CWT ----
sma_size_CWT <- lmodel2(log10(mort_rate_mo) ~ log10(dbh_cm), "interval", "interval", data=CWT, nperm=99) 
sma_size_CWT
# report p-value
sma_size_CWT$P.param

plot.CWT <- ggplot(CWT) + 
  stat_ma_line(aes(x = dbh_cm, y = mort_rate_mo), method="SMA", color = "black", linewidth = 0) + 
  geom_errorbar(aes(x = dbh_cm, ymin = lower, ymax = upper), width = 0.04) + 
  geom_point(aes(x = dbh_cm, y = mort_rate_mo), fill = "#FDDBC7", pch=21, size=3.5) +
  stat_ma_line(aes(x = dbh_cm, y = mort_rate_mo), method="SMA", color = "black", linewidth = 1, se=F) + 
  xlab("Diameter at breast height (cm)") +
  ylab(expression(paste('Mortality rate (ind ', ~ind^-1, ~mo^-1, ')'))) +
  coord_cartesian(xlim = c(0.3, 25), ylim = c(0.00001, 0.11)) + 
  scale_x_continuous(trans="log10", breaks = c(1, 10, 30)) +
  scale_y_continuous(trans = "log10", breaks = c(0.0001, 0.001, 0.010, 0.100)) +
  theme_bw()
plot.CWT

## HFR ----
sma_size_HFR <- lmodel2(log10(mort_rate_mo) ~ log10(dbh_cm), "interval", "interval", data=HFR, nperm=99) 
sma_size_HFR
# report p-value
sma_size_HFR$P.param

plot.HFR <- ggplot(HFR) + 
  stat_ma_line(aes(x = dbh_cm, y = mort_rate_mo), method="SMA", color = "black", linewidth = 0) + 
  geom_errorbar(aes(x = dbh_cm, ymin = lower, ymax = upper), width = 0.04) + 
  geom_point(aes(x = dbh_cm, y = mort_rate_mo), fill = "#D1E5F0", pch=21, size=3.5) +
  stat_ma_line(aes(x = dbh_cm, y = mort_rate_mo), method="SMA", color = "black", linewidth = 1, se=F) + 
  xlab("Diameter at breast height (cm)") +
  ylab(expression(paste('Mortality rate (ind ', ~ind^-1, ~mo^-1, ')'))) +
  coord_cartesian(xlim = c(0.3, 25), ylim = c(0.00001, 0.11)) + 
  scale_x_continuous(trans="log10", breaks = c(1, 10, 30)) +
  scale_y_continuous(trans = "log10", breaks = c(0.0001, 0.001, 0.010, 0.100)) +
  theme_bw()
plot.HFR

## LUQ ----
sma_size_LUQ <- lmodel2(log10(mort_rate_mo) ~ log10(dbh_cm), "interval", "interval", data=LUQ, nperm=99) 
sma_size_LUQ
# report p-value
sma_size_LUQ$P.param

plot.LUQ <- ggplot(LUQ) + 
  stat_ma_line(aes(x = dbh_cm, y = mort_rate_mo), method="SMA", color = "black", linewidth = 0) + 
  geom_errorbar(aes(x = dbh_cm, ymin = lower, ymax = upper), width = 0.04) + 
  geom_point(aes(x = dbh_cm, y = mort_rate_mo), fill = "#F4A582", pch=21, size=3.5) +
  stat_ma_line(aes(x = dbh_cm, y = mort_rate_mo), method="SMA", color = "black", linewidth = 1, se=F) + 
  xlab("Diameter at breast height (cm)") +
  ylab(expression(paste('Mortality rate (ind ', ~ind^-1, ~mo^-1, ')'))) +
  coord_cartesian(xlim = c(0.3, 25), ylim = c(0.00001, 0.11)) + 
  scale_x_continuous(trans="log10", breaks = c(1, 10, 30)) +
  scale_y_continuous(trans = "log10", breaks = c(0.0001, 0.001, 0.010, 0.100)) +
  theme_bw()
plot.LUQ

## NWT ----
sma_size_NWT <- lmodel2(log10(mort_rate_mo) ~ log10(dbh_cm), "interval", "interval", data=NWT, nperm=99)
sma_size_NWT
# report p-value
sma_size_NWT$P.param

plot.NWT <- ggplot(NWT) + 
  stat_ma_line(aes(x = dbh_cm, y = mort_rate_mo), method="SMA", color = "black", linewidth = 0) + 
  geom_errorbar(aes(x = dbh_cm, ymin = lower, ymax = upper), width = 0.04) + 
  geom_point(aes(x = dbh_cm, y = mort_rate_mo), fill = "#2166AC", pch=21, size=3.5) +
  stat_ma_line(aes(x = dbh_cm, y = mort_rate_mo), method="SMA", color = "black", linewidth = 1, se=F) + 
  xlab("Diameter at breast height (cm)") +
  ylab(expression(paste('Mortality rate (ind ', ~ind^-1, ~mo^-1, ')'))) +
  coord_cartesian(xlim = c(0.3, 25), ylim = c(0.00001, 0.11)) +
  scale_x_continuous(trans = "log10", breaks = c(1, 10, 30)) +
  scale_y_continuous(trans = "log10", breaks = c(0.0001, 0.001, 0.010, 0.100)) +
  theme_bw()
plot.NWT

## PSR ----
sma_size_PSR <- lmodel2(log10(mort_rate_mo) ~ log10(dbh_cm), "interval", "interval", data=PSR, nperm=99) 
sma_size_PSR
# report p-value
sma_size_PSR$P.param

plot.PSR <- ggplot(PSR) + 
  stat_ma_line(aes(x = dbh_cm, y = mort_rate_mo), method="SMA", color = "black", linewidth = 0) + 
  geom_errorbar(aes(x = dbh_cm, ymin = lower, ymax = upper), width = 0.04) + 
  geom_point(aes(x = dbh_cm, y = mort_rate_mo), fill = "#D6604D", pch=21, size=3.5) +
  stat_ma_line(aes(x = dbh_cm, y = mort_rate_mo), method="SMA", color = "black", linewidth = 1, se=F) + 
  xlab("Diameter at breast height (cm)") +
  ylab(expression(paste('Mortality rate (ind ', ~ind^-1, ~mo^-1, ')'))) +
  coord_cartesian(xlim = c(0.3, 25), ylim = c(0.00001, 0.11)) + 
  scale_x_continuous(trans = "log10", breaks = c(1, 10, 30)) +
  scale_y_continuous(trans = "log10", breaks = c(0.0001, 0.001, 0.010, 0.100)) +
  theme_bw()
plot.PSR


# Fig. 2a: Plot sites together ----
Fig_2a <- ggarrange(plot.BCI + rremove("xlab") + rremove("ylab"), 
                    plot.PSR + rremove("xlab") + rremove("ylab"), 
                    plot.LUQ + rremove("xlab") + rremove("ylab"), 
                    plot.CWT + rremove("xlab") + rremove("ylab"), 
                    plot.HFR + rremove("xlab") + rremove("ylab"), 
                    plot.AND + rremove("xlab") + rremove("ylab"), 
                    plot.CAP + rremove("xlab") + rremove("ylab"), 
                    plot.NWT + rremove("xlab") + rremove("ylab"))

## arrange plots together
Fig_2a <- annotate_figure(Fig_2a, left = textGrob(expression(paste('Mortality rate (ind ', ~ind^-1, ~mo^-1, ')')), rot = 90, gp = gpar(cex = 1.3)),
                          bottom = textGrob("Diameter at breast height (cm)", gp = gpar(cex = 1.3, face = "bold")))
## add axis titles to arranged plots
Fig_2a


# Fig. 2b: Fitted exponents by site ----

# Calculate mean (95% CI) rate-dbh scaling exponent
calc_mean_ci <- function(numbers) {
  n <- length(numbers)
  mean_val <- mean(numbers)
  std_err <- sd(numbers) / sqrt(n)
  margin <- qt(0.975, df = n - 1) * std_err
  lower_interval <- mean_val - margin
  upper_interval <- mean_val + margin
  result <- list(mean = mean_val, lower_ci = lower_interval, upper_ci = upper_interval)
  return(result)
}
calc_mean_ci(c(sma_size_AND$regression.results$Slope[3],
               sma_size_BCI$regression.results$Slope[3], 
               sma_size_CAP$regression.results$Slope[3], 
               sma_size_CWT$regression.results$Slope[3], 
               sma_size_HFR$regression.results$Slope[3], 
               sma_size_LUQ$regression.results$Slope[3], 
               sma_size_NWT$regression.results$Slope[3], 
               sma_size_PSR$regression.results$Slope[3]))

# Make dataframe with estimated slopes and confidence intervals
exp_sma <- data.frame(site=c("AND", "BCI", "CAP", "CWT",  "HFR", "LUQ", "NWT", "PSR"), slope=NA, CI_low=NA, CI_high=NA)
# Add GST data to map dataframe
exp_sma$GST_C <- climate$GST_C[match(exp_sma$site, climate$Site)]
# Order points by GST to assign colors for plotting
exp_sma$site=factor(exp_sma$site, levels=unique(exp_sma$site[order(-exp_sma$GST_C)]), ordered=TRUE)

# Fill in fitted slopes and 95% CI by iterating over each row
for (i in 1:nrow(exp_sma)) {
  # Get the current site
  current_site <- exp_sma$site[i]
  
  # Check if the current site matches any of the dataframes
  if (current_site == "AND") {
    exp_sma$slope[i] <- sma_size_AND$regression.results[3, 3]
    exp_sma$CI_low[i] <- sma_size_AND$confidence.intervals[3, 4]
    exp_sma$CI_high[i] <- sma_size_AND$confidence.intervals[3, 5]
  } else if (current_site == "BCI") {
    exp_sma$slope[i] <- sma_size_BCI$regression.results[3, 3]
    exp_sma$CI_low[i] <- sma_size_BCI$confidence.intervals[3, 4]
    exp_sma$CI_high[i] <- sma_size_BCI$confidence.intervals[3, 5]
  } else if (current_site == "CAP") {
    exp_sma$slope[i] <- sma_size_CAP$regression.results[3, 3]
    exp_sma$CI_low[i] <- sma_size_CAP$confidence.intervals[3, 4]
    exp_sma$CI_high[i] <- sma_size_CAP$confidence.intervals[3, 5]
  } else if (current_site == "CWT") {
    exp_sma$slope[i] <- sma_size_CWT$regression.results[3, 3]
    exp_sma$CI_low[i] <- sma_size_CWT$confidence.intervals[3, 4]
    exp_sma$CI_high[i] <- sma_size_CWT$confidence.intervals[3, 5]
  } else if (current_site == "HFR") {
    exp_sma$slope[i] <- sma_size_HFR$regression.results[3, 3]
    exp_sma$CI_low[i] <- sma_size_HFR$confidence.intervals[3, 4]
    exp_sma$CI_high[i] <- sma_size_HFR$confidence.intervals[3, 5]
  } else if (current_site == "LUQ") {
    exp_sma$slope[i] <- sma_size_LUQ$regression.results[3, 3]
    exp_sma$CI_low[i] <- sma_size_LUQ$confidence.intervals[3, 4]
    exp_sma$CI_high[i] <- sma_size_LUQ$confidence.intervals[3, 5]
  } else if (current_site == "NWT") {
    exp_sma$slope[i] <- sma_size_NWT$regression.results[3, 3]
    exp_sma$CI_low[i] <- sma_size_NWT$confidence.intervals[3, 4]
    exp_sma$CI_high[i] <- sma_size_NWT$confidence.intervals[3, 5]
  } else if (current_site == "PSR") {
    exp_sma$slope[i] <- sma_size_PSR$regression.results[3, 3]
    exp_sma$CI_low[i] <- sma_size_PSR$confidence.intervals[3, 4]
    exp_sma$CI_high[i] <- sma_size_PSR$confidence.intervals[3, 5]
  }
}

Fig_2b <- ggplot(exp_sma, aes(x=site, y=slope, fill = site)) +
  geom_point(size = 3.5) + 
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high), width = 0.3) +
  coord_cartesian(ylim = c(-2.1, 0.2)) +
  labs(y = "Diameter-scaling exponent (dimensionless)", x = "Site") + 
  scale_y_continuous(breaks = seq(-2, 0, by = 1)) +
  geom_abline(intercept = -0.66666667, slope = 0, color = "black", size = 0.75, lty = "dashed") + theme_bw(base_size = 15) +
  scale_fill_manual(values = c("BCI" = "#B2182B","PSR" ="#D6604D","LUQ" = "#F4A582","CWT" = "#FDDBC7","HFR" = "#D1E5F0","AND" = "#92C5DE","CAP" = "#4393C3", "NWT" = "#2166AC")) +
  geom_point(pch=21, size=3.5) +
  scale_color_discrete(limits=c("BCI","PSR","LUQ","CWT","HFR","MTB","AND","CAP","NWT")) +
  theme(axis.text.y=element_text(size=15), axis.title=element_text(size=15),legend.text=element_text(size=15),legend.title=element_blank())
Fig_2b


# Fig. 2: Mortality rate vs. size ----
Fig_2 <- ggarrange(Fig_2a, Fig_2b, labels = 'auto')
Fig_2


# Fig. S5: SMA fit to all sites combined ----
sma_size_all <- lmodel2(log10(mort_rate_mo) ~ log10(dbh_cm), "interval", "interval", data=mort_data, nperm=99) # Note this uses base 10 log rather than base e
sma_size_all
# report p-value
sma_size_all$P.param

Fig_S5 <- ggplot(mort_data) + 
  stat_ma_line(aes(x = dbh_cm, y = mort_rate_mo), method="SMA", color = "black", linewidth = 0) + 
  geom_errorbar(aes(x = dbh_cm, ymin = lower, ymax = upper), width = 0.04) + 
  geom_point(aes(x = dbh_cm, y = mort_rate_mo), data=AND, fill = "#92C5DE", pch=21, size=3.5) +
  geom_point(aes(x = dbh_cm, y = mort_rate_mo), data=BCI, fill = "#B2182B", pch=21, size=3.5) +
  geom_point(aes(x = dbh_cm, y = mort_rate_mo), data=CAP, fill = "#4393C3", pch=21, size=3.5) +
  geom_point(aes(x = dbh_cm, y = mort_rate_mo), data=CWT, fill = "#FDDBC7", pch=21, size=3.5) +
  geom_point(aes(x = dbh_cm, y = mort_rate_mo), data=HFR, fill = "#D1E5F0", pch=21, size=3.5) +
  geom_point(aes(x = dbh_cm, y = mort_rate_mo), data=LUQ, fill = "#F4A582", pch=21, size=3.5) +
  geom_point(aes(x = dbh_cm, y = mort_rate_mo), data=NWT, fill = "#2166AC", pch=21, size=3.5) +
  geom_point(aes(x = dbh_cm, y = mort_rate_mo), data=PSR, fill = "#D6604D", pch=21, size=3.5) +
  stat_ma_line(aes(x = dbh_cm, y = mort_rate_mo), method="SMA", color = "black", linewidth = 1, se=F) + 
  xlab("Diameter at breast height (cm)") +
  ylab(expression(paste('Mortality rate (ind ', ~ind^-1, ~mo^-1, ')'))) +
  coord_cartesian(xlim = c(10^-0.4115726, 10^1.5803646), ylim = c(0.00001, 0.11)) +
  scale_x_continuous(trans = "log10", breaks = c(1, 10, 30)) +
  scale_y_continuous(trans = "log10", breaks = c(0.0001, 0.001, 0.010, 0.100)) +
  theme_bw()
Fig_S5
