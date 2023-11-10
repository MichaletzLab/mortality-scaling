# Description ----
# Fig. 3 and Fig. S6: Mortality rate vs. temperature

# Citation:
# Borrego, I., Perez, T.M., Bentley, L.P., Bison, N.N., Byrnes, L, Candido, H.G., Chmurzynski, A. Durán, S.M., 
# Fox, T.S., Gaitan, M.A., Garen, J.C., Orwig, D.A., Pau, S. Scott, J.L., Simovic, M. Swenson, N.G., Wieczynski, D.J.,
# Buzzard, V. Enquist, B.J., & Michaletz, S.T. (in review) From tropics to treeline: extending and assessing 
# metabolic scaling theory for global variation in plant mortality rates.

# Plotting RMA and 95% CI from lmodel2 in ggplot2
# Sean Michaletz, 2 Dec 2022
# Updated using stat_ma_line() from package ggpmisc, Sean Michaletz, 6 Dec 2022
# Updated modified Arrhenius space following Michaletz & Garen (in review), Sean Michaletz, 11 July 2023
# Updated November 2023 for GitHub upload

# Initialize ----
#--Libraries
library(lmodel2)
library(ggplot2)
library(ggpmisc)
library(scales)

# Load data for mortality rate, size, temperature, and traits (file from Isaac, JUL 2023)
mort_data <- read.csv("https://raw.githubusercontent.com/MichaletzLab/mortality-scaling/main/data/mortality_data.csv")
# Load data for site-level mortality rates (file from Isaac, JUL 2023)
mort_data_site <- read.csv("https://raw.githubusercontent.com/MichaletzLab/mortality-scaling/main/data/macrosystems_MAT_Mort_update9.9.21.csv")
# Load site climate data
climate <- read.csv("https://raw.githubusercontent.com/MichaletzLab/mortality-scaling/main/data/climate.csv",header=T)

# Add growing season length, MAT, and GST data to mort_data
mort_data$L_gs_mo <- climate$L_gs_mo[match(mort_data$site, climate$Site)]
mort_data$MAT_C <- climate$MAT_C[match(mort_data$site, climate$Site)]
mort_data$GST_C <- climate$GST_C[match(mort_data$site, climate$Site)]
mort_data_site$L_gs_mo <- climate$L_gs_mo[match(mort_data_site$site, climate$Site)]
mort_data_site$MAT_C <- climate$MAT_C[match(mort_data_site$site, climate$Site)]
mort_data_site$GST_C <- climate$GST_C[match(mort_data_site$site, climate$Site)]

# Add -1/kT to mort_data for modified Arrhenius analyses (Michaletz & Garen, in review)
mort_data$neg_1_kT <- -1/(0.00008617*(mort_data$GST_C + 273.15))
mort_data_site$neg_1_kT <- -1/(0.00008617*(mort_data_site$GST_C + 273.15))

# Calculate monthly mortality rate (correct annual rates for growing season length)
mort_data$mort_rate_mo <- mort_data$avg_mort / mort_data$L_gs_mo
mort_data_site$mort_rate_mo <- mort_data_site$avg_mort / mort_data_site$L_gs_mo
mort_data_site$SE_mort_rate_mo <- mort_data_site$se_mort / mort_data_site$L_gs_mo

# Order sites by GST so legend appears in order of descending GST (as per Isaac's design)
mort_data$site <-factor(mort_data$site, levels=c("BCI","PSR","LUQ","CWT","HFR","AND","CAP","NWT"))
mort_data_site$site <-factor(mort_data_site$site, levels=c("BCI","PSR","LUQ","CWT","HFR","AND","CAP","NWT"))


# Mortality rate vs. temperature ----
# First, fit the RMA using lmodel2
rma_mort.temp <- lmodel2(log(mort_rate_mo) ~ neg_1_kT, "interval", "interval", data = mort_data, nperm=99)
rma_mort.temp
# report p-value
rma_mort.temp$P.param


## Mortality vs. T by size class (Fig. 3) ----
# Plot the data, RMA, and 95% confidence band in ggplot2 using predict.lmodel2() from package ggpmisc.
Fig_3 <- ggplot(mort_data) +
  stat_ma_line(aes(x = neg_1_kT, y = mort_rate_mo), method="RMA", range.y = "interval", range.x = "interval", color = "black", linewidth = 0) +
  geom_point(aes(x = neg_1_kT, y = mort_rate_mo, fill = site), pch=21, size=3.5) +
  stat_ma_line(aes(x = neg_1_kT, y = mort_rate_mo), method="RMA", range.y = "interval", range.x = "interval", se=F, color = "black", linewidth = 1) +
  xlab(expression(paste('Reciprocal thermal energy, ', '-1/',italic('k')[italic("B")], italic('T'), ' (',  eV^{-1}, ')'))) +
  ylab(expression(paste('Mortality rate (ind ', ~ind^-1, ~mo^-1, ')'))) +
  scale_x_continuous(sec.axis = sec_axis(trans = ~ (-1/(.*0.00008617))-273.15 , name = expression(paste('Temperature'~(degree*C))))) +
  scale_y_continuous(trans="log", breaks = trans_breaks("log", function(x) exp(x), n=3), labels = trans_format("log", math_format(e^.x))) +
  scale_fill_manual(values = c("BCI" = "#B2182B","PSR" ="#D6604D","LUQ" = "#F4A582","CWT" = "#FDDBC7","HFR" = "#D1E5F0","AND" = "#92C5DE","CAP" = "#4393C3","NWT" = "#2166AC")) + 
  guides(fill=guide_legend(title="")) +
  theme_bw(base_size = 15)
Fig_3


## Mortality vs. T by site (Fig. S6) ----
# First, fit the RMA using lmodel2 # E = 0.64 eV (0.21 - 1.41)
rma_mort.temp.site <- lmodel2(log(mort_rate_mo) ~ neg_1_kT, "interval", "interval", data = mort_data_site, nperm=99)
rma_mort.temp.site
# report p-value
rma_mort.temp.site$P.param

Fig_S6 <- ggplot(mort_data_site) +
  stat_ma_line(aes(x = neg_1_kT, y = mort_rate_mo), method="RMA", range.y = "interval", range.x = "interval", color = "black", linewidth = 0) +
  geom_errorbar(aes(x = neg_1_kT, ymin = mort_rate_mo - SE_mort_rate_mo, ymax = mort_rate_mo + SE_mort_rate_mo)) +
  geom_point(aes(x = neg_1_kT, y = mort_rate_mo, fill = site), pch=21, size=3.5) +
  stat_ma_line(aes(x = neg_1_kT, y = mort_rate_mo), method="RMA", range.y = "interval", range.x = "interval", se=F, color = "black", linewidth = 1) +
  xlab(expression(paste('Reciprocal thermal energy, ', '-1/',italic('k')[italic("B")], italic('T'), ' (',  eV^{-1}, ')'))) +
  ylab(expression(paste('Mortality rate (ind ', ~ind^-1, ~mo^-1, ')'))) +
  scale_x_continuous(sec.axis = sec_axis(trans = ~ (-1/(.*0.00008617))-273.15 , name = expression(paste('Temperature'~(degree*C))))) +
  scale_y_continuous(trans="log", breaks = trans_breaks("log", function(x) exp(x), n=3), labels = trans_format("log", math_format(e^.x))) +
  scale_fill_manual(values = c("BCI" = "#B2182B","PSR" ="#D6604D","LUQ" = "#F4A582","CWT" = "#FDDBC7","HFR" = "#D1E5F0","AND" = "#92C5DE","CAP" = "#4393C3","NWT" = "#2166AC")) + 
  guides(fill=guide_legend(title="")) +
  theme_bw(base_size = 15)
Fig_S6