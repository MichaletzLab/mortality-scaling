# Description ----
# Figure S2: Number of surviving stems versus time

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

# Load data
Nt <- read.csv("https://raw.githubusercontent.com/MichaletzLab/mortality-scaling/main/data/macrosystems_decay_curve.csv",header=T)

# Order sites by GST so legend appears in order of descending GST (as per Isaac's design)
Nt$site <-factor(Nt$site, levels=c("BCI","PSR","LUQ","CWT","HFR","AND","CAP","NWT"))

## Fig. S2: Plot ln(Nt) vs. t ----
Fig_S2 <- ggplot(Nt, aes(x = time_cumul, y = alive, fill=site)) + 
  geom_point(pch=21, colour="black", size=3.5) +
  geom_smooth(aes(color=site), method='lm', se=FALSE) +
  xlab("Time (mo)") + 
  ylab("Number of surviving stems (dimensionless)") + 
  scale_y_continuous(trans = "log", breaks = c(700, 1000, 3000, 5700)) +
  scale_fill_manual(values = c("BCI" = "#B2182B","PSR" ="#D6604D","LUQ" = "#F4A582","CWT" = "#FDDBC7","HFR" = "#D1E5F0","AND" = "#92C5DE","CAP" = "#4393C3","NWT" = "#2166AC")) + 
  scale_color_manual(values = c("BCI" = "#B2182B","PSR" ="#D6604D","LUQ" = "#F4A582","CWT" = "#FDDBC7","HFR" = "#D1E5F0","AND" = "#92C5DE","CAP" = "#4393C3","NWT" = "#2166AC")) + 
  theme_bw(base_size = 15)
Fig_S2