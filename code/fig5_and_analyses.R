# Description ----
# Fig. 5: Mortality rate vs. size, temperature, and climate

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
library(car)
library(ggplot2)
library(ggpubr)
library(grid)
library(rsq)

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


# Mortality rate vs. size, temperature, and climate ----

# Multiple regression with free trait-scaling parameters
m_reg_fit <- lm(log(mort_rate_mo) ~ log(PN_ratio) + log(lma) + log(wood_density) + log(dbh_cm) + neg_1_kT, mort_data)
vif(m_reg_fit) # All variance-inflation factors < 5, so OK
summary(m_reg_fit)
rsq.partial(m_reg_fit)
confint(m_reg_fit)

# Partial regression plots
#--Use avPlots to calculate adjusted values and store as lists
panel_5a.data <- avPlots(m_reg_fit, ~log(dbh_cm), col="gray70", col.lines="black", main="", 
                         xlab=expression(atop(paste("Diameter at breast height (cm)"))),
                         ylab=expression(paste('Adjusted mortality rate (ind ', ~ind^-1, ~mo^-1, ')')))
panel_5b.data <- avPlots(m_reg_fit, ~neg_1_kT, col="gray70", col.lines="black", main="", 
                         xlab=expression(paste('Reciprocal thermal energy, ', '-1/',italic('k')[italic("B")], italic('T'), ' (',  eV^{-1}, ')')),
                         ylab=expression(paste('Adjusted mortality rate (ind ', ~ind^-1, ~mo^-1, ')')))
panel_5c.data <- avPlots(m_reg_fit, ~log(wood_density), col="gray70", col.lines="black", main="", 
                         xlab=expression(paste("Wood density (g", ~cm^-3, ")")),
                         ylab=expression(paste('Adjusted mortality rate (ind ', ~ind^-1, ~mo^-1, ')')))
panel_5d.data <- avPlots(m_reg_fit, ~log(lma), col="gray70", col.lines="black", main="", 
                         xlab=expression(paste("Leaf mass per area (g", ~cm^-2, ")")),
                         ylab=expression(paste('Adjusted mortality rate (ind ', ~ind^-1, ~mo^-1, ')')))
panel_5e.data <- avPlots(m_reg_fit, ~log(PN_ratio), col="gray70", col.lines="black", main="", 
                         xlab=expression(paste("Leaf P:N (g P g", ~N^-1, ")")),
                         ylab=expression(paste('Adjusted mortality rate (ind ', ~ind^-1, ~mo^-1, ')')))

#--Store adjusted values in a dataframe
Fig_5.data <- data.frame(panel_5a.data, panel_5b.data, panel_5c.data, panel_5d.data, panel_5e.data)
colnames(Fig_5.data) <- c("dbh", "dbh_mort", "neg_1_kT", "neg_1_kT_mort", "wood_density", "wood_density_mort", "lma", "lma_mort", "PN_ratio", "PN_ratio_mort")

#--Convert ln'ed data back to absolute numbers for plotting
Fig_5.data$dbh <- exp(Fig_5.data$dbh)
Fig_5.data$dbh_mort <- exp(Fig_5.data$dbh_mort)
Fig_5.data$neg_1_kT <- exp(Fig_5.data$neg_1_kT)
Fig_5.data$neg_1_kT_mort <- exp(Fig_5.data$neg_1_kT_mort)
Fig_5.data$wood_density <- exp(Fig_5.data$wood_density)
Fig_5.data$wood_density_mort <- exp(Fig_5.data$wood_density_mort)
Fig_5.data$lma <- exp(Fig_5.data$lma)
Fig_5.data$lma_mort <- exp(Fig_5.data$lma_mort)
Fig_5.data$PN_ratio <- exp(Fig_5.data$PN_ratio)
Fig_5.data$PN_ratio_mort <- exp(Fig_5.data$PN_ratio_mort)

# Add a column to identify what site points belong to. Do this by 
# getting the site names from mort_data, ignoring rows for which there
# are no mortality data in that size class.
Fig_5.data$site <- mort_data$site[!is.na(mort_data$mort_rate_mo)]

#--Plot panels
base_breaks <- function(n = 10){
  function(x) {
    axisTicks(log(range(x, na.rm = TRUE)), log = TRUE, n = n)
  }
}

panel_5a <- ggplot(Fig_5.data, aes(x=dbh, y=dbh_mort)) + 
  geom_smooth(method = "lm", se=TRUE, color="black", linewidth=0) +
  geom_point(aes(fill = site), pch=21, colour="black", size=3.5) +
  geom_smooth(method = "lm", se=F, color="black") + 
  xlab(expression(atop("Adjusted diameter", paste("at breast height (cm)")))) + 
  ylab(expression(paste('Adjusted mortality rate (ind ', ~ind^-1, ~mo^-1, ')'))) +
  scale_y_continuous(trans = "log", breaks = c(0.2, 0.5, 1.0, 2.0, 5.0)) + 
  theme_bw(base_size=12) + 
  theme(plot.title=element_text(hjust=0.94, vjust=-1.8)) +
  scale_fill_manual(values = c("BCI" = "#B2182B","PSR" ="#D6604D","LUQ" = "#F4A582","CWT" = "#FDDBC7","HFR" = "#D1E5F0","AND" = "#92C5DE","CAP" = "#4393C3","NWT" = "#2166AC")) +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=15),legend.text=element_text(size=15),legend.title=element_blank()) 
panel_5a

panel_5b <- ggplot(Fig_5.data, aes(x=neg_1_kT, y=neg_1_kT_mort)) + 
  geom_smooth(method = "lm", se=TRUE, color="black", linewidth=0) +
  geom_point(aes(fill = site), pch=21, colour="black", size=3.5) +
  geom_smooth(method = "lm", se=F, color="black") + 
  xlab(expression(atop("Adjusted reciprocal thermal", paste('energy, -1/',italic('k')[italic("B")], italic('T'), ' (',  eV^{-1}, ')')))) +
  ylab(expression(paste('Adjusted mortality rate (ind ', ~ind^-1, ~mo^-1, ')'))) +
  scale_y_continuous(trans = "log", breaks = c(0.2, 0.5, 1.0, 2.0, 5.0)) + 
  theme_bw(base_size=12) + 
  theme(plot.title=element_text(hjust=0.94, vjust=-1.8)) +
  scale_fill_manual(values = c("BCI" = "#B2182B","PSR" ="#D6604D","LUQ" = "#F4A582","CWT" = "#FDDBC7","HFR" = "#D1E5F0","AND" = "#92C5DE","CAP" = "#4393C3","NWT" = "#2166AC")) +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=15),legend.text=element_text(size=15),legend.title=element_blank()) 
panel_5b

panel_5c <- ggplot(Fig_5.data, aes(x=wood_density, y=wood_density_mort)) + 
  geom_smooth(method = "lm", se=TRUE, color="black", linewidth=0) +
  geom_point(aes(fill = site), pch=21, colour="black", size=3.5) +
  geom_smooth(method = "lm", se=F, color="black") + 
  xlab(expression(atop("Adjusted wood", paste("density (g", ~cm^-3, ")")))) +
  ylab(expression(paste('Adjusted mortality rate (ind ', ~ind^-1, ~mo^-1, ')'))) +
  scale_y_continuous(trans = "log", breaks = c(0.2, 0.5, 1.0, 2.0, 5.0)) + 
  theme_bw(base_size=12) + 
  theme(plot.title=element_text(hjust=0.94, vjust=-1.8)) +
  scale_fill_manual(values = c("BCI" = "#B2182B","PSR" ="#D6604D","LUQ" = "#F4A582","CWT" = "#FDDBC7","HFR" = "#D1E5F0","AND" = "#92C5DE","CAP" = "#4393C3","NWT" = "#2166AC")) +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=15),legend.text=element_text(size=15),legend.title=element_blank()) 
panel_5c

panel_5d <- ggplot(Fig_5.data, aes(x=lma, y=lma_mort)) + 
  geom_smooth(method = "lm", se=TRUE, color="black", linewidth=0) +
  geom_point(aes(fill = site), pch=21, colour="black", size=3.5) +
  geom_smooth(method = "lm", se=F, color="black") + 
  xlab(expression(atop("Adjusted leaf mass", paste(" per area (g", ~cm^-2, ")")))) +
  ylab(expression(paste('Adjusted mortality rate (ind ', ~ind^-1, ~mo^-1, ')'))) +
  scale_y_continuous(trans = "log", breaks = c(0.2, 0.5, 1.0, 2.0, 5.0)) + 
  theme_bw(base_size=12) + 
  theme(plot.title=element_text(hjust=0.94, vjust=-1.8)) +
  scale_fill_manual(values = c("BCI" = "#B2182B","PSR" ="#D6604D","LUQ" = "#F4A582","CWT" = "#FDDBC7","HFR" = "#D1E5F0","AND" = "#92C5DE","CAP" = "#4393C3","NWT" = "#2166AC")) +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=15),legend.text=element_text(size=15),legend.title=element_blank()) 
panel_5d

panel_5e <- ggplot(Fig_5.data, aes(x=PN_ratio, y=PN_ratio_mort)) + 
  geom_smooth(method = "lm", se=TRUE, color="black", linewidth=0) +
  geom_point(aes(fill = site), pch=21, colour="black", size=3.5) +
  geom_smooth(method = "lm", se=F, color="black") + 
  xlab(expression(atop("Adjusted leaf P:N", paste("(g P g", ~N^-1, ")")))) +
  ylab(expression(paste('Adjusted mortality rate (ind ', ~ind^-1, ~mo^-1, ')'))) +
  scale_y_continuous(trans = "log", breaks = c(0.2, 0.5, 1.0, 2.0, 5.0)) + 
  theme_bw(base_size=12) + 
  theme(plot.title=element_text(hjust=0.94, vjust=-1.8)) +
  scale_fill_manual(values = c("BCI" = "#B2182B","PSR" ="#D6604D","LUQ" = "#F4A582","CWT" = "#FDDBC7","HFR" = "#D1E5F0","AND" = "#92C5DE","CAP" = "#4393C3","NWT" = "#2166AC")) +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=15),legend.text=element_text(size=15),legend.title=element_blank()) 
panel_5e

# Fig. 5: Partial regression plots ----
Fig_5 <- ggarrange(panel_5a + rremove("ylab"), panel_5b + rremove("ylab"), panel_5c + rremove("ylab"), panel_5d + rremove("ylab"), panel_5e + rremove("ylab"), common.legend = TRUE, legend = "right", labels = 'auto')
Fig_5 <- annotate_figure(Fig_5, left = textGrob(expression(paste('Adjusted mortality rate (ind ', ~ind^-1, ~mo^-1, ')')), rot = 90, gp = gpar(cex = 1.3)))
Fig_5