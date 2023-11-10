# Description ----
# Fig. S7: Traits vs. MAT

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
library(ggpubr)
library(grid)

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

# Order sites by GST so legend appears in order of descending GST (as per Isaac's design)
mort_data$site <-factor(mort_data$site, levels=c("BCI","PSR","LUQ","CWT","HFR","AND","CAP","NWT"))

# Add MAT to mort_data
mort_data$MAT_C <- climate$MAT_C[match(mort_data$site, climate$Site)]


## Fig. S7: Traits vs. MAT ----
# Fig. S7a: P:N vs. temperature
model_pn.temp <- lm(PN_ratio ~ MAT_C, data = mort_data)
summary(model_pn.temp)

Fig_S7a <- ggplot(mort_data, aes(x = MAT_C, y = PN_ratio)) + 
  geom_smooth(method = "lm", se=TRUE, color="black", linewidth=0) +
  geom_point(aes(fill = site), pch=21,size = 3.5) + 
  geom_smooth(method = "lm", se=F, color="black") +
  xlab(expression(paste('Mean annual temperature'~(degree*C)))) +
  ylab(expression(paste("Leaf P:N (g P g", ~N^-1, ")"))) +
  scale_fill_manual(values = c("BCI" = "#B2182B","PSR" ="#D6604D","LUQ" = "#F4A582","CWT" = "#FDDBC7","HFR" = "#D1E5F0","AND" = "#92C5DE","CAP" = "#4393C3","NWT" = "#2166AC")) + 
  guides(fill=guide_legend(title="")) +
  theme_bw(base_size = 15)
Fig_S7a

# Fig. S7b: LMA vs. temperature
model_lma.temp <- lm(lma ~ MAT_C, data = mort_data)
summary(model_lma.temp)

Fig_S7b <- ggplot(mort_data, aes(x = MAT_C, y = lma)) + 
  geom_smooth(method = "lm", se=TRUE, color="black", linewidth=0) +
  geom_point(aes(fill = site), pch=21,size = 3.5) + 
  geom_smooth(method = "lm", se=F, color="black") + 
  xlab(expression(paste('Mean annual temperature'~(degree*C)))) +
  ylab(expression(paste("Leaf mass per area (g", ~cm^-2, ")"))) +
  scale_fill_manual(values = c("BCI" = "#B2182B","PSR" ="#D6604D","LUQ" = "#F4A582","CWT" = "#FDDBC7","HFR" = "#D1E5F0","AND" = "#92C5DE","CAP" = "#4393C3","NWT" = "#2166AC")) + 
  guides(fill=guide_legend(title="")) +
  theme_bw(base_size = 15)
Fig_S7b

# Fig. S7c: Wood density vs. temperature
model_wd.temp <- lm(wood_density ~ MAT_C, data = mort_data)
summary(model_wd.temp)

Fig_S7c <- ggplot(mort_data, aes(x = MAT_C, y = wood_density)) + 
  geom_smooth(method = "lm", se=TRUE, color="black", linewidth=0) +
  geom_point(aes(fill = site), pch=21,size = 3.5) + 
  geom_smooth(method = "lm", se=F, color="black") +
  xlab(expression(paste('Mean annual temperature'~(degree*C)))) +
  ylab(expression(paste("Wood density (g", ~cm^-3, ")"))) +
  scale_fill_manual(values = c("BCI" = "#B2182B","PSR" ="#D6604D","LUQ" = "#F4A582","CWT" = "#FDDBC7","HFR" = "#D1E5F0","AND" = "#92C5DE","CAP" = "#4393C3","NWT" = "#2166AC")) +
  guides(fill=guide_legend(title="")) +
  theme_bw(base_size = 15)
Fig_S7c

# Fig. S7d: Composite trait vs. temperature
model_traitscomp.temp <- lm(PN_ratio^0.76*lma^-1*wood_density^-0.75 ~ MAT_C, data = mort_data)
summary(model_traitscomp.temp)

Fig_S7d <- 
  ggplot(mort_data, aes(x = MAT_C, y = PN_ratio^0.76*lma^-1*wood_density^-0.75)) + 
  #geom_smooth(method = "lm", se=TRUE, color="black", linewidth=0) +
  geom_point(aes(fill = site), pch=21,size = 3.5) + 
  #geom_smooth(method = "lm", se=F, color="black") + 
  xlab(expression(paste('Mean annual temperature'~(degree*C)))) +
  ylab(expression(paste("Composite trait  (g ", ~P^0.76, "g ", ~N^-0.76, ~g^-1.75, ~cm^4.25, ")"))) +
  scale_fill_manual(values = c("BCI" = "#B2182B","PSR" ="#D6604D","LUQ" = "#F4A582","CWT" = "#FDDBC7","HFR" = "#D1E5F0","AND" = "#92C5DE","CAP" = "#4393C3","NWT" = "#2166AC")) + 
  guides(fill=guide_legend(title="")) +
  theme_bw(base_size = 15)
Fig_S7d

# Fig. S7: Traits vs. MAT
Fig_S7 <- ggarrange(Fig_S7a + rremove("xlab"), Fig_S7b + rremove("xlab"), Fig_S7c + rremove("xlab"), Fig_S7d + rremove("xlab"), common.legend = TRUE, legend = "right", labels = "auto")
Fig_S7 <- annotate_figure(Fig_S7, bottom = textGrob('Mean annual temperature'~(degree*C), gp = gpar(cex = 1.3)))
Fig_S7