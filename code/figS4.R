# Description ----
# Figure S4: Traits by size class

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
library(dplyr)
library(grid)
library(ggpubr)

# Load data for mortality rate, size, temperature, and traits (file from Isaac, JUL 2023)
mort_data <- read.csv("https://raw.githubusercontent.com/MichaletzLab/mortality-scaling/main/data/mortality_data.csv")

# Rename DBH size class from 'size' to 'dbh_cm'
names(mort_data)[names(mort_data) == 'size'] <- 'dbh_cm'
# Rename leaf P and N to avoid confusion with C
names(mort_data)[names(mort_data) == 'lpc'] <- 'P_percent'
names(mort_data)[names(mort_data) == 'lnc'] <- 'N_percent'

# Calculate leaf P:N
mort_data$PN_ratio <- mort_data$P_percent / mort_data$N_percent

# Order sites by GST so legend appears in order of descending GST (as per Isaac's design)
mort_data$site <-factor(mort_data$site, levels=c("BCI","PSR","LUQ","CWT","HFR","AND","CAP","NWT"))

## Fig. S4: Traits by size class ----
# This figure illustrates variation in community-weighted traits by size class

# P:N vs. DBH (Fig. S4a)
# Test for significant variation in P:N with DBH
# Group the data by the "site" variable and fit the linear model for each group
grouped_models <- mort_data %>%
  group_by(site) %>%
  do(model_pn.dbh_cm = lm(PN_ratio ~ dbh_cm, data = .))
# Access the summary of each model
summary_list <- lapply(grouped_models$model_pn.dbh_cm, summary)
# Print the summaries
for (i in seq_along(summary_list)) {
  cat("Summary for Site", i, "\n")
  print(summary_list[[i]])
  cat("\n")
}

Fig_S4a <- ggplot(mort_data, aes(x = dbh_cm, y = PN_ratio)) + 
  geom_point(aes(fill = site), pch=21,size = 3.5) + 
  geom_smooth(aes(color=site), method = "lm", se=F) +
  xlab(expression(paste("Diameter at breast height (cm)"))) +
  ylab(expression(paste("Leaf P:N (g P g", ~N^-1, ")"))) +
  scale_fill_manual(values = c("BCI" = "#B2182B","PSR" ="#D6604D","LUQ" = "#F4A582","CWT" = "#FDDBC7","HFR" = "#D1E5F0","AND" = "#92C5DE","CAP" = "#4393C3","NWT" = "#2166AC")) + 
  scale_color_manual(values = c("BCI" = "#B2182B","PSR" ="#D6604D","LUQ" = "#F4A582","CWT" = "#FDDBC7","HFR" = "#D1E5F0","AND" = "#92C5DE","CAP" = "#4393C3","NWT" = "#2166AC")) + 
  guides(fill=guide_legend(title=""), color="none") +
  theme_bw(base_size = 15)
Fig_S4a

# LMA vs. dbh_cm (Fig. S4b)
# Test for significant variation in P:N with DBH
# Group the data by the "site" variable and fit the linear model for each group
grouped_models <- mort_data %>%
  group_by(site) %>%
  do(model_lma.dbh = lm(lma ~ dbh_cm, data = .))
# Access the summary of each model
summary_list <- lapply(grouped_models$model_lma.dbh, summary)
# Print the summaries
for (i in seq_along(summary_list)) {
  cat("Summary for Site", i, "\n")
  print(summary_list[[i]])
  cat("\n")
}

Fig_S4b <- ggplot(mort_data, aes(x = dbh_cm, y = lma)) + 
  geom_point(aes(fill = site), pch=21,size = 3.5) + 
  geom_smooth(aes(color=site), method = "lm", se=F) +
  xlab(expression(paste("Diameter at breast height (cm)"))) +
  ylab(expression(paste("Leaf mass per area (g", ~cm^-2, ")"))) +
  scale_fill_manual(values = c("BCI" = "#B2182B","PSR" ="#D6604D","LUQ" = "#F4A582","CWT" = "#FDDBC7","HFR" = "#D1E5F0","AND" = "#92C5DE","CAP" = "#4393C3","NWT" = "#2166AC")) + 
  scale_color_manual(values = c("BCI" = "#B2182B","PSR" ="#D6604D","LUQ" = "#F4A582","CWT" = "#FDDBC7","HFR" = "#D1E5F0","AND" = "#92C5DE","CAP" = "#4393C3","NWT" = "#2166AC")) + 
  guides(fill=guide_legend(title=""), color="none") +
  theme_bw(base_size = 15)
Fig_S4b

# Wood density vs. dbh_cm (Fig. S4c)
# Test for significant variation in wood density with DBH
# Group the data by the "site" variable and fit the linear model for each group
grouped_models <- mort_data %>%
  group_by(site) %>%
  do(model_density.dbh = lm(wood_density ~ dbh_cm, data = .))
# Access the summary of each model
summary_list <- lapply(grouped_models$model_density.dbh, summary)
# Print the summaries
for (i in seq_along(summary_list)) {
  cat("Summary for Site", i, "\n")
  print(summary_list[[i]])
  cat("\n")
}

Fig_S4c <- ggplot(mort_data, aes(x = dbh_cm, y = wood_density)) + 
  geom_point(aes(fill = site), pch=21,size = 3.5) + 
  geom_smooth(aes(color=site), method = "lm", se=F) +
  xlab(expression(paste("Diameter at breast height (cm)"))) +
  ylab(expression(paste("Wood density (g", ~cm^-3, ")"))) +
  scale_fill_manual(values = c("BCI" = "#B2182B","PSR" ="#D6604D","LUQ" = "#F4A582","CWT" = "#FDDBC7","HFR" = "#D1E5F0","AND" = "#92C5DE","CAP" = "#4393C3","NWT" = "#2166AC")) + 
  scale_color_manual(values = c("BCI" = "#B2182B","PSR" ="#D6604D","LUQ" = "#F4A582","CWT" = "#FDDBC7","HFR" = "#D1E5F0","AND" = "#92C5DE","CAP" = "#4393C3","NWT" = "#2166AC")) + 
  guides(fill=guide_legend(title=""), color="none") +
  theme_bw(base_size = 15)
Fig_S4c

# Composite trait vs. dbh_cm (Fig. S4d)
# Test for significant variation in wood density with DBH
# Group the data by the "site" variable and fit the linear model for each group
grouped_models <- mort_data %>%
  group_by(site) %>%
  do(model_composite.dbh = lm(PN_ratio^0.76*lma^-1*wood_density^-0.75 ~ dbh_cm, data = .))
# Access the summary of each model
summary_list <- lapply(grouped_models$model_composite.dbh, summary)
# Print the summaries
for (i in seq_along(summary_list)) {
  cat("Summary for Site", i, "\n")
  print(summary_list[[i]])
  cat("\n")
}

Fig_S4d <- 
  ggplot(mort_data, aes(x = dbh_cm, y = PN_ratio^0.76*lma^-1*wood_density^-0.75)) + 
  geom_point(aes(fill = site), pch=21,size = 3.5) + 
  geom_smooth(aes(color=site), method = "lm", se=F) +
  xlab(expression(paste("Diameter at breast height (cm)"))) +
  ylab(expression(paste("Composite trait  (g ", ~P^0.76, "g ", ~N^-0.76, ~g^-1.75, ~cm^4.25, ")"))) +
  scale_fill_manual(values = c("BCI" = "#B2182B","PSR" ="#D6604D","LUQ" = "#F4A582","CWT" = "#FDDBC7","HFR" = "#D1E5F0","AND" = "#92C5DE","CAP" = "#4393C3","NWT" = "#2166AC")) + 
  scale_color_manual(values = c("BCI" = "#B2182B","PSR" ="#D6604D","LUQ" = "#F4A582","CWT" = "#FDDBC7","HFR" = "#D1E5F0","AND" = "#92C5DE","CAP" = "#4393C3","NWT" = "#2166AC")) + 
  guides(fill=guide_legend(title=""), color="none") +
  theme_bw(base_size = 15)
Fig_S4d

# Traits vs. DBH (Fig. S4)
Fig_S4 <- ggarrange(Fig_S4a + rremove("xlab"), Fig_S4b + rremove("xlab"), Fig_S4c + rremove("xlab"), Fig_S4d + rremove("xlab"), common.legend = TRUE, legend = "right", labels = "auto")
Fig_S4 <- annotate_figure(Fig_S4, bottom = textGrob("Diameter at breast height (cm)", gp = gpar(cex = 1.3)))
Fig_S4
