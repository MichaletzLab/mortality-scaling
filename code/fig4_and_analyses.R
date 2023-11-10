# Description ----
# Figure 4: Mortality rate vs. traits

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
library(ggplot2)
library(ggpmisc)
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

# Add -1/kT to mort_data for modified Arrhenius analyses (Michaletz & Garen, in review)
mort_data$neg_1_kT <- -1/(0.00008617*(mort_data$GST_C + 273.15))

# Calculate monthly mortality rate (correct annual rates for growing season length)
mort_data$mort_rate_mo <- mort_data$avg_mort / mort_data$L_gs_mo

# Order sites by GST so legend appears in order of descending GST (as per Isaac's design)
mort_data$site <-factor(mort_data$site, levels=c("BCI","PSR","LUQ","CWT","HFR","AND","CAP","NWT"))


# Fig. 4: Mortality rate vs. traits ----
#--Mortality rate vs. P:N
# Fit the RMA using lmodel2
rma_PN <- lmodel2(log10(mort_rate_mo) ~ log10(PN_ratio), "interval", "interval", data = mort_data, nperm=99)
rma_PN
# report p-value
rma_PN$P.param

# Plot the data, RMA, and 95% confidence band in ggplot2 using predict.lmodel2() from package ggpmisc.
Fig_4a <- ggplot(mort_data) + 
  stat_ma_line(aes(x = PN_ratio, y = mort_rate_mo), method="RMA", range.y = "interval", range.x = "interval", color = "black", linewidth = 0) +
  geom_point(aes(x = PN_ratio, y = mort_rate_mo, fill = site), pch=21, colour="black", size=3.5) + 
  stat_ma_line(aes(x = PN_ratio, y = mort_rate_mo), method="RMA", range.y = "interval", range.x = "interval", se=F, color = "black", linewidth = 1) +
  xlab(expression(paste("Leaf P:N (g P g", ~N^-1, ")"))) + 
  ylab(expression(paste('Mortality rate (ind ', ~ind^-1, ~mo^-1, ')'))) +
  guides(fill=guide_legend(title="")) +
  scale_y_continuous(breaks = c(0.0005, 0.002, 0.005, 0.015), trans = "log10") +
  scale_x_continuous(breaks = c(0.05, 0.1, 0.2),trans = "log10") +
  scale_fill_manual(values = c("BCI" = "#B2182B","PSR" ="#D6604D","LUQ" = "#F4A582","CWT" = "#FDDBC7","HFR" = "#D1E5F0","AND" = "#92C5DE","CAP" = "#4393C3","NWT" = "#2166AC")) + 
  theme_bw(base_size = 15)
Fig_4a

#--Mortality rate vs. LMA
# Fit the RMA using lmodel2
rma_lma <- lmodel2(log10(mort_rate_mo) ~ log10(lma), "interval", "interval", data = mort_data, nperm=99)
rma_lma
# report p-value
rma_lma$P.param

# Plot the data, RMA, and 95% confidence band in ggplot2 using predict.lmodel2() from package ggpmisc.
Fig_4b <- ggplot(mort_data) +
  stat_ma_line(aes(x = lma, y = mort_rate_mo), method="RMA", range.y = "interval", range.x = "interval", color = "black", linewidth = 0) +
  geom_point(aes(x = lma, y = mort_rate_mo, fill = site), pch=21, colour="black", size=3.5) +
  stat_ma_line(aes(x = lma, y = mort_rate_mo), method="RMA", range.y = "interval", range.x = "interval", se=F, color = "black", linewidth = 1) +
  xlab(expression(paste("Leaf mass per area (g", ~cm^-2, ")"))) + 
  ylab(expression(paste('Mortality rate (ind ', ~ind^-1, ~mo^-1, ')'))) +
  guides(fill=guide_legend(title="")) +
  scale_y_continuous(breaks = c(0.0005, 0.002, 0.005, 0.015), trans = "log") + 
  scale_x_continuous(breaks = c(0.005, 0.01, 0.02), trans = "log") +
  scale_fill_manual(values = c("BCI" = "#B2182B","PSR" ="#D6604D","LUQ" = "#F4A582","CWT" = "#FDDBC7","HFR" = "#D1E5F0","AND" = "#92C5DE","CAP" = "#4393C3","NWT" = "#2166AC")) + 
  theme_bw(base_size = 15)
Fig_4b

#--Mortality rate vs. wood density
# Fit the RMA using lmodel2
rma_wd <- lmodel2(log10(mort_rate_mo) ~ log10(wood_density), "interval", "interval", data = mort_data, nperm=99)
rma_wd
# report p-value
rma_wd$P.param

# Plot the data, RMA, and 95% confidence band in ggplot2 using predict.lmodel2() from package ggpmisc.
Fig_4c <- ggplot(mort_data) + 
  stat_ma_line(aes(x = wood_density, y = mort_rate_mo), method="RMA", range.y = "interval", range.x = "interval", color = "black", linewidth = 0) + 
  geom_point(aes(x = wood_density, y = mort_rate_mo, fill = site), pch=21, colour="black", size=3.5) +
  stat_ma_line(aes(x = wood_density, y = mort_rate_mo), method="RMA", range.y = "interval", range.x = "interval", se=F, color = "black", linewidth = 1) + 
  xlab(expression(paste("Wood density (g", ~cm^-3, ")"))) + 
  ylab(expression(paste('Mortality rate (ind ', ~ind^-1, ~mo^-1, ')'))) +
  guides(fill=guide_legend(title="")) +
  scale_x_continuous(trans = "log10") +
  scale_y_continuous(breaks = c(0.0005, 0.002, 0.005, 0.015), trans = "log10") + 
  scale_fill_manual(values = c("BCI" = "#B2182B","PSR" ="#D6604D","LUQ" = "#F4A582","CWT" = "#FDDBC7","HFR" = "#D1E5F0","AND" = "#92C5DE","CAP" = "#4393C3","NWT" = "#2166AC")) + 
  theme_bw(base_size = 15)
Fig_4c

Fig_4 <- ggarrange(Fig_4a + rremove("ylab"), Fig_4b + rremove("ylab"), Fig_4c + rremove("ylab"), common.legend = TRUE, legend = 'right', nrow = 1, labels = 'auto')
Fig_4 <- annotate_figure(Fig_4, left = textGrob(expression(paste('Mortality rate (ind ', ~ind^-1, ~mo^-1, ')')), rot = 90, gp = gpar(cex = 1.3)))
Fig_4
