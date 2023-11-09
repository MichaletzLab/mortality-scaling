# Description ----
# Figure 1: Map and climate space diagram for Forest MacroSystems sites

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
library(cowplot)

# Load site coordinate data
map <- read.csv("https://raw.githubusercontent.com/MichaletzLab/mortality-scaling/main/data/FMS_map.csv",header=T)
# Load site climate data
climate <- read.csv("https://raw.githubusercontent.com/MichaletzLab/mortality-scaling/main/data/climate.csv",header=T)
# Load Whittaker biome outlines (from Michaletz et al. 2014 in Nature)
whit <- read.csv("https://raw.githubusercontent.com/MichaletzLab/mortality-scaling/main/data/whittaker.csv",header=T)

# Add site codes to map dataframe
map$site <- c("CAP", "AND", "HFR", "NWT", "CWT", "MTB", "LUQ", "PSR", "BCI")
# Remove MTB since not used in paper
map <- map[-c(6),]
# Add GST data to map dataframe
map$GST_C <- climate$GST_C[match(map$site, climate$Site)]

## Fig. 1a: Plot map ----
world_map <- map_data("world")
# Order points by GST to assign colors for plotting
map$site=factor(map$site, levels=unique(map$site[order(-map$GST_C)]), ordered=TRUE)

# Plot map
Fig_1a <- ggplot(world_map, aes(x = long, y = lat, group = group)) +
  geom_polygon(fill="gray", colour = "black")+
  geom_point(data=map, aes(x=long, y=lat, fill = site), inherit.aes = FALSE, colour="black", shape=21, size=5, alpha=0.75)+
  scale_fill_manual(values = c("BCI" = "#B2182B","PSR" ="#D6604D","LUQ" = "#F4A582","CWT" = "#FDDBC7","HFR" = "#D1E5F0","AND" = "#92C5DE","CAP" = "#4393C3","NWT" = "#2166AC"))+
  ylab(expression("Latitude"~(degree)))+
  xlab(expression("Longitude"~(degree)))+
  guides(fill=guide_legend(title=""))+
  theme_bw(base_size=13) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  coord_cartesian(xlim = c(-135,-50),ylim = c(5,55))
Fig_1a

## Fig. 1b: Plot climate space diagram ----
# Define coordinates for plotting biome labels (numbers 1-9)
biomes <- data.frame(x=c(25.8,18.7,28.1,18.5,2.2,28.4,18,-8.7,28.4), y=c(4150,2800,1425,2100,1500,675,500,250,200)) 

Fig_1b <- ggplot(climate, aes(x=MAT_C, y=MAP_mm)) + 
  geom_path(data=whit, aes(x=MAT_whittaker, y=MAP_whittaker)) +
  geom_point(data=climate, aes(x=MAT_C, y=MAP_mm, fill = Site ), inherit.aes = FALSE, colour="black", shape=21, size=5, alpha=0.75)+
  scale_fill_manual(values = c("BCI" = "#B2182B","PSR" ="#D6604D","LUQ" = "#F4A582","CWT" = "#FDDBC7","HFR" =  "#D1E5F0","AND" = "#92C5DE","CAP" = "#4393C3","NWT" = "#2166AC"))+
  xlab(expression('Mean annual temperature'~(degree*C))) + 
  ylab("Mean annual precipitation (mm)") + 
  guides(fill=guide_legend(title=""))+
  geom_text(data = biomes, aes(label = c(1:9), x=x, y=y), size=4) + 
  theme_bw(base_size=13) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none")
Fig_1b

## Fig. 1: Plot map & climate space together ----
legend <- get_legend(Fig_1a)
Fig_1a <- Fig_1a + theme(legend.position='none')

dev.new()
#tiff(filename="FMS_maps.tiff", width=8, height = 8, units="in", res=400)
pdf("Figure_1.pdf", width=10, height = 3.5)
plot_grid(Fig_1a, Fig_1b, legend, ncol=3, rel_widths = c(3,3,2))
dev.off()
