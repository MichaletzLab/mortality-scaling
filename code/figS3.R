# Description ----
# Figure S3: Diameter-scaling exponents vs. bin size

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
bt <- read.csv("https://raw.githubusercontent.com/MichaletzLab/mortality-scaling/main/data/macrosystems_bin_test.csv",header=T)

# Order sites by GST so legend appears in order of descending GST (as per Isaac's design)
bt$site <-factor(bt$site, levels=c("BCI","PSR","LUQ","CWT","HFR","AND","CAP","NWT"))

## Fig. S3: Plot diameter-scaling exponent vs. bin size ----
Fig_S3 <- ggplot(bt, aes(x = bin_size, y = sma_slope, fill=site)) + 
  geom_line(aes(color=site)) +
  geom_point(pch=21, colour="black", size=3.5) +
  geom_errorbar(aes(ymin = sma_ci_low, ymax = sma_ci_up, color=site), width = 0.1) +
  geom_abline(intercept = -0.666666, slope = 0, lty = "dashed") +
  xlab("Bin size (cm)") + 
  ylab("Estimated diameter-scaling exponent (dimensionless)") + 
  scale_x_continuous(limits = c(0.5, 10.5), breaks = c(1,2,3,4,5,6,7,8,9,10)) +
  scale_fill_manual(values = c("BCI" = "#B2182B","PSR" ="#D6604D","LUQ" = "#F4A582","CWT" = "#FDDBC7","HFR" = "#D1E5F0","AND" = "#92C5DE","CAP" = "#4393C3","NWT" = "#2166AC")) + 
  scale_color_manual(values = c("BCI" = "#B2182B","PSR" ="#D6604D","LUQ" = "#F4A582","CWT" = "#FDDBC7","HFR" = "#D1E5F0","AND" = "#92C5DE","CAP" = "#4393C3","NWT" = "#2166AC")) + 
  guides(fill=guide_legend(title=""), color="none") +
  theme_bw(base_size = 15)
Fig_S3
