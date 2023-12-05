# Script to plot sample map

library(ggplot2)
library(dplyr)
library(sf)
library(tidyterra)
library(terra)
library(raster)
library(patchwork)
source("scripts/theme_emily.R")

#~~ Read in map data

world.map <- vect("data/map/TM_WORLD_BORDERS_SIMPL-0.3.shp")

# Crop to avoid weird lines
#world.map <- terra::crop(world.map, extent(-179.9, 179.9, -90, 84))

# Make sf object
world.map <- st_as_sf(world.map)

# Remove Antarctica
world.map <- world.map[world.map$NAME != "Antarctica",]

# Transform to robinson projection

#world.robinson <- st_transform(world.map,  crs ="+proj=robin +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs") 
#world.robinson <- terra::crop(world.robinson, extent(-179.9, 179.9, -90, 84))

#~~ Read in distribution data from the IUCN

# BIROSTRIS

birostris <- vect("data/map/redlist_species_data_632abf27-0988-40f2-aa38-d8b58d8e3e98/data_0.shp")

bir_extant <- birostris[birostris$LEGEND=='Extant (resident)',]
bir_pos_extant <- birostris[birostris$LEGEND=='Possibly Extant (resident)',]

bir_extant <- st_as_sf(bir_extant)
bir_pos_extant <- st_as_sf(bir_pos_extant)

#bir_extant_robinson <- st_transform(bir_extant,  crs ="+proj=robin +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs") 
#bir_pos_extant_robinson <- st_transform(bir_pos_extant,  crs ="+proj=robin +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs") 

# ALFREDI

alfredi <- vect("data/map/redlist_species_data_948f544a-06d7-4b7d-9f4e-f5a0e8196d15/data_0.shp")

alf_extant <- alfredi[alfredi$LEGEND=='Extant (resident)',]
alf_pos_extant <- alfredi[alfredi$LEGEND=='Possibly Extant (resident)',]

alf_extant <- st_as_sf(alf_extant)
alf_pos_extant <- st_as_sf(alf_pos_extant)

#alf_extant_robinson <- st_transform(alf_extant,  crs ="+proj=robin +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs") 
#alf_pos_extant_robinson <- st_transform(alf_pos_extant,  crs ="+proj=robin +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs") 

#~~ Read in sampling location files

sample_sites <- read.csv("data/map/sampling_sites_lat_long.csv") %>%
  filter(Species != "Mobula Sp. 1") # n = 176

colnames(sample_sites) <- c("sample", "dup", "species",
                            "country", "site", "lat", "long")

# Summarise number of samples per site

birostris_samples <- sample_sites %>%
  filter(species == "Mobula birostris") %>% # n = 85
  dplyr::group_by(species, country, lat, long) %>%
  dplyr::summarise(n = n()) %>%
  mutate(lon_samples = long,
         lat_samples = lat)

alfredi_samples <- sample_sites %>%
  filter(species == "Mobula alfredi") %>% # n = 91
  dplyr::group_by(species, country, lat, long) %>%
  dplyr::summarise(n = n()) %>%
  mutate(lon_samples = long,
         lat_samples = lat)


#~~ Plot map


pal <- c("#EAD3BF",
         "#4F1C2E",
         "#7E8AA2",
         "#C16683",
         "#B53325",
         "#F2B705")

# BIROSTRIS

#D27A15
#1D8981

m_bir <- ggplot() + 
  geom_spatvector(data = world.map, colour = "grey99", fill = "grey80", size = 0.1) +
  geom_spatvector(data = bir_pos_extant, col = "grey99", fill = "#655A6C", alpha = 0.1) +
  geom_spatvector(data = bir_extant, col = NA, fill = "#655A6C", alpha = 0.5) +
  geom_point(aes(x = lon_samples, y = lat_samples, color = country, size = n), 
             data = birostris_samples, alpha = 0.9) + # 3, 2
  scale_colour_manual(values = pal) + 
 # theme_emily() +
  guides(color = "none") +
  theme(panel.border = element_blank(),
        panel.grid = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white"),
        #axis.line = element_line(colour = "white"),
        axis.line = element_blank(),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.key = element_rect(fill = "white"),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        plot.title = element_text(hjust = 0.04),
        plot.subtitle = element_text(hjust = 0.08)) +
  scale_size(range = c(3,9)) + 
  #ggtitle("C")
  ggtitle(label = "C",
          subtitle = expression(paste("Oceanic manta ray (", italic("Mobula birostris"), ")")))
  #ggtitle(expression(paste("B - Oceanic manta ray, ", italic("(Mobula birostris)"))))

ggsave("figs/map_birostris.png", m_bir, width = 11, height = 6)
ggsave("figs/map_birostris_ppt.png", m_bir, width = 9, height = 4)

saveRDS(m_bir, "figs/map_birostris.RDS")


m_bir_ppt <- ggplot() + 
  geom_spatvector(data = world.map, colour = "grey99", fill = "grey80", size = 0.1) +
  geom_spatvector(data = bir_pos_extant, col = "grey99", fill = "#655A6C", alpha = 0.1) +
  geom_spatvector(data = bir_extant, col = NA, fill = "#655A6C", alpha = 0.5) +
  scale_colour_manual(values = pal) + 
  # theme_emily() +
  guides(color = "none") +
  theme(panel.border = element_blank(),
        panel.grid = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white"),
        #axis.line = element_line(colour = "white"),
        axis.line = element_blank(),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.key = element_rect(fill = "white"),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        plot.title = element_text(hjust = 0.04),
        plot.subtitle = element_text(hjust = 0.08)) +
  scale_size(range = c(3,9))


ggsave("figs/map_birostris_dist.png", m_bir_ppt, width = 9, height = 4)





# ALFREDI

pal <- c("#9986A5",
         "#9A8822",
         "#35274A",
         "#3B9AB2",
         "#0B775E",
         "#E2D200")

m_alf <- ggplot() + 
  geom_spatvector(data = world.map, colour = "grey99", fill = "grey80", size = 0.1) +
  geom_spatvector(data = alf_pos_extant, col = "grey99", fill = "#1D8981", alpha = 0.1) +
  geom_spatvector(data = alf_extant, col = NA, fill = "#1D8981", alpha = 0.4) +
  geom_point(aes(x = lon_samples, y = lat_samples, color = country, size = n), 
             data = alfredi_samples, alpha = 0.9) + # 3, 2
  scale_colour_manual(values = pal) + 
  # theme_emily() +
  guides(color = "none") +
  theme(panel.border = element_blank(),
        panel.grid = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white"),
        #axis.line = element_line(colour = "white"),
        axis.line = element_blank(),
        legend.position = "none",
        legend.title = element_blank(),
        legend.key = element_rect(fill = "white"),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        plot.title = element_text(hjust = 0.05),
        plot.subtitle = element_text(hjust = 0.08)) +
  scale_size(range = c(3,9)) + 
  #ggtitle("A")
  ggtitle(label = "A",
          subtitle = expression(paste("Reef manta ray (", italic("Mobula alfredi"), ")")))

ggsave("figs/map_alfredi.png", m_alf, width = 11, height = 6)
ggsave("figs/map_alfredi_ppt.png", m_alf, width = 9, height = 4)

saveRDS(m_alf, "figs/map_alfredi.RDS")

m_alf_dist <- ggplot() + 
  geom_spatvector(data = world.map, colour = "grey99", fill = "grey80", size = 0.1) +
  geom_spatvector(data = alf_pos_extant, col = "grey99", fill = "#1D8981", alpha = 0.1) +
  geom_spatvector(data = alf_extant, col = NA, fill = "#1D8981", alpha = 0.4) +
  scale_colour_manual(values = pal) + 
  # theme_emily() +
  guides(color = "none") +
  theme(panel.border = element_blank(),
        panel.grid = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white"),
        #axis.line = element_line(colour = "white"),
        axis.line = element_blank(),
        legend.position = "none",
        legend.title = element_blank(),
        legend.key = element_rect(fill = "white"),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        plot.title = element_text(hjust = 0.05),
        plot.subtitle = element_text(hjust = 0.08)) +
  scale_size(range = c(3,9))

ggsave("figs/map_alfredi_dist.png", m_alf_dist, width = 9, height = 4)
