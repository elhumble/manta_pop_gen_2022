# Oversea distance ~ Fst

library(marmap)
library(data.table)
library(dplyr)
library(adegenet)
library(ggplot2)
library(wesanderson)
library(patchwork)
library(tidyr)
library(hierfstat)
library(scales)
library(geosphere)
source("scripts/theme_emily.R")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#           Alfredi            #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~ Get distance matrix

# Sampling locations

latlong_alf <- fread("data/meta/alfredi_lat_long.txt", header = T)

latlong_alf <- latlong_alf %>%
  mutate(SITE = gsub("_", " ", SITE))

coords_alf <- latlong_alf[,c(3,2)]
rownames(coords_alf) <- latlong_alf$SITE


# import NOAA bathymetry data
# needs to be high res or else paths go over land

bat <- getNOAA.bathy(lon1 = -180,
                     lon2 = 180,
                     lat1 = -60,
                     lat2 = 30, res = 4, 
                     keep = TRUE, path = "data/marmap/")

# Get depth of points (must all be negative)

get.depth(bat, x = coords_alf$LONG, y = coords_alf$LAT, locator = F)

# Create transition object with a low negative depth constraint to avoid inaccuracies (-10)
# Takes a long time, saved and reloading

#tr <- trans.mat(bat, min.depth = -10, max.depth = NULL)
#save(tr, file = "data/marmap/transition_alfredi.Rdata")
load("data/marmap/transition_alfredi.Rdata")
tr_alf <- tr
rm(tr)

# Compute least cost paths below -10m and plot on map
# Takes a long time, saved and reloading

#lc_path <- lc.dist(tr, coords, res="path")
#save(lc_path, file = "data/marmap/lc_path_alfredi.Rdata")
load("data/marmap/lc_path_alfredi.Rdata")
lc_path_alf <- lc_path
rm(lc_path)

# Distance matrix

#lc_matrix <- lc.dist(tr, coords, res="dist")
#save(lc_matrix, file = "data/marmap/lc_matrix_alfredi.Rdata")
load("data/marmap/lc_matrix_alfredi.Rdata")
lc_matrix_alf <- lc_matrix
rm(lc_matrix)

# check distances make sense and prepare matrix

lc_matrix_alf
geodist_alf <- as.matrix(lc_matrix_alf)
rownames(geodist_alf) <- rownames(coords_alf)
colnames(geodist_alf) <- rownames(geodist_alf)

# Remove AP due to small sample size

geodist_alf <- geodist_alf[c(1:4,6),c(1:4,6)]

#~~ Prepare genetic info

meta <- read.csv("data/meta/sample_info_mastersheet_reformatted.csv") %>%
  dplyr::select(c(Sample_code, Species, Location))

gl_alf <- read.PLINK("data/plink/across_sp/alfredi_geno4_dp6_ind_geno_depth_no_dups_maf_unrelated_ld.raw")

inds_alf <- data.frame(gl_alf@ind.names) %>%
  separate(gl_alf.ind.names, c("Filename", "Dup"), sep = "_") %>%
  left_join(meta, by = c("Filename" = "Sample_code")) %>%
  unite(ind_names, c("Filename", "Dup")) %>%
  mutate(ind_names = gsub("_NA", "", ind_names))

# Assign populations
pop(gl_alf) <- inds_alf$Location

# Fsts

pw_fst_alf <- dartR::gl.fst.pop(gl_alf, nboots=1000, percent=95, nclusters=1)

# Prepare matrix

pw_fst_alf$Fsts[is.na(pw_fst_alf$Fsts)] <- t(pw_fst_alf$Fsts)[is.na(pw_fst_alf$Fsts)]
alf_melted_cor <- pw_fst_alf$Fsts

order_alf <- c("Seychelles", "Chagos", "Maldives", "Hawaii", "Fiji")
fstWC_alf <- alf_melted_cor[order_alf, order_alf]

#~~ Double check column names are the same

colnames(geodist_alf)
colnames(fstWC_alf)
colnames(geodist_alf) == colnames(fstWC_alf)

#~~ 3D space

# geodist_alf <- log(geodist_alf)
# fstWC_alf <- fstWC_alf / (1 - fstWC_alf)

#~~ Run Mantel test

#~~ Convert to dist object
geodist_alf <- as.dist(geodist_alf)  ; geodist_alf
fstWC_alf <- as.dist(fstWC_alf) ; fstWC_alf
plot(fstWC_alf ~ geodist_alf)


#~~ Mantel test using all sites
man_alf <- mantel.rtest(geodist_alf, fstWC_alf, nrepet = 1000)
man_alf

#~~ Create labels for plots

fstlab_alf <- expression(italic("F")[ST])
r1_alf <- round(man_alf$obs,2)
p1_alf <- format(round(man_alf$pvalue, 3), nsmall=3)

#~~ Create dataframe of distance matrices

df_dist_alf <- data.frame(geodistance=as.vector(geodist_alf),
                gendistance=as.vector(fstWC_alf))

#~~ Plot

ibd_alf <- ggplot(df_dist_alf) + 
  geom_point(aes(x = geodistance, y = gendistance),
             col = "#1D8981", fill = "#1D8981", # "grey69"
             shape = 21, alpha = 0.5, stroke = 0.1, size = 3) + 
  geom_smooth(aes(x = geodistance, y = gendistance), method = "lm",
              se = TRUE, level= 0.95, alpha = 0.05, size = 1, col = "#1D8981", fill = "#1D8981") +
  xlab("Geographic distance (km)") +
  ylab(fstlab_alf) +
  scale_x_continuous(labels = comma) +
  ggtitle("B") +
  theme_emily()

ibd_alf
man_alf
r1_alf
p1_alf

saveRDS(ibd_alf, "figs/ibd_alf.RDS")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#           Birostris          #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~ Get distance matrix

# Sampling locations

latlong_bir <- fread("data/meta/birostris_lat_long.txt", header = T)

latlong_bir <- latlong_bir %>%
  mutate(SITE = gsub("_", " ", SITE))

coords_bir <- latlong_bir[,c(3,2)]
rownames(coords_bir) <- latlong_bir$SITE

# import NOAA bathymetry data
# needs to be high res or else paths go over land

# bat <- getNOAA.bathy(lon1 = -180,
#                      lon2 = 180,
#                      lat1 = -60,
#                      lat2 = 30, res = 4,
#                      keep = TRUE, path = "data/marmap/")

# Get depth of points (must all be negative)

get.depth(bat, x = coords_bir$LONG, y = coords_bir$LAT, locator = F)

# Create transition object with a low negative depth constraint to avoid inaccuracies (-10)

#tr <- trans.mat(bat, min.depth = -10, max.depth = NULL)
#save(tr, file = "data/marmap/transition_birostris.Rdata")
load("data/marmap/transition_birostris.Rdata")
tr_bir <- tr
rm(tr)

# Compute least cost paths below -10m and plot on map

#lc_path <- lc.dist(tr, coords, res="path")
#save(lc_path, file = "data/marmap/lc_path_birostris.Rdata")
load("data/marmap/lc_path_birostris.Rdata")
lc_path_bir <- lc_path
rm(lc_path)

# Distance matrix

#lc_matrix <- lc.dist(tr, coords, res="dist")
#save(lc_matrix, file = "data/marmap/lc_matrix_birostris.Rdata")
load("data/marmap/lc_matrix_birostris.Rdata")
lc_matrix_bir <- lc_matrix
rm(lc_matrix)

# check distances make sense and prepare matrix

lc_matrix_bir
geodist_bir <- as.matrix(lc_matrix_bir)
rownames(geodist_bir) <- rownames(coords_bir)
colnames(geodist_bir) <- rownames(geodist_bir)

# Remove SA due to small sample size

geodist_bir <- geodist_bir[c(1:5),c(1:5)]

#~~ Prepare genetic info

meta <- read.csv("data/meta/sample_info_mastersheet_reformatted.csv") %>%
  dplyr::select(c(Sample_code, Species, Location))

gl_bir <- read.PLINK("data/plink/across_sp/birostris_geno4_dp6_ind_geno_depth_no_dups_maf_unrelated_ld.raw")
gl_bir@n.loc

inds_bir <- data.frame(gl_bir@ind.names) %>%
  separate(gl_bir.ind.names, c("Filename", "Dup"), sep = "_") %>%
  left_join(meta, by = c("Filename" = "Sample_code")) %>%
  unite(ind_names, c("Filename", "Dup")) %>%
  mutate(ind_names = gsub("_NA", "", ind_names))

# Assign populations
pop(gl_bir) <- inds_bir$Location

# Fsts

pw_fst_bir <- dartR::gl.fst.pop(gl_bir, nboots = 1000, percent = 95, nclusters = 1)

# Prepare matrix

pw_fst_bir$Fsts[is.na(pw_fst_bir$Fsts)] <- t(pw_fst_bir$Fsts)[is.na(pw_fst_bir$Fsts)]
bir_melted_cor <- pw_fst_bir$Fsts

order_bir <- c("Mexico Pacific", "Peru", "Mexico Caribbean", "Sri Lanka" ,"the Philippines", "South Africa")
fstWC_bir <- bir_melted_cor[order_bir, order_bir]

# Remove South Africa due to small sample sizes

fstWC_bir <- fstWC_bir[c(1:5), c(1:5)]

# Replace negative values for Mantel test

fstWC_bir[fstWC_bir < 0] = 0.0000001
fstWC_bir

#~~ Double check column names are the same

colnames(geodist_bir)
colnames(fstWC_bir)
colnames(geodist_bir) == colnames(fstWC_bir)

#~~ 3D space

#geodist_bir <- log(geodist_bir)
#fstWC_bir <- fstWC_bir / (1 - fstWC_bir)

#~~ Run Mantel test

#~~ Convert to dist object

geodist_bir <- as.dist(geodist_bir)  ; geodist_bir
fstWC_bir <- as.dist(fstWC_bir) ; fstWC_bir
plot(fstWC_bir ~ geodist_bir)


#~~ Mantel test using all sites
man_bir <- mantel.rtest(geodist_bir, fstWC_bir, nrepet = 1000)
man_bir

#~~ Create labels for plots

fstlab_bir <- expression(italic("F")[ST])
r1_bir <- round(man_bir$obs,2)
p1_bir <- format(round(man_bir$pvalue, 3), nsmall=3)

#~~ Create dataframe of distance matrices

df_dist_bir <- data.frame(geodistance=as.vector(geodist_bir),
                      gendistance=as.vector(fstWC_bir))

#~~ Plot

ibd_bir <- ggplot(df_dist_bir) + 
  geom_point(aes(x = geodistance, y = gendistance),
             col = "#D27A15", fill = "#D27A15",
             shape = 21, stroke = 0.1, size = 3) + 
  geom_smooth(aes(x = geodistance, y = gendistance), method = "lm",
              se = TRUE, level= 0.95, alpha = 0.05, size = 1, 
              col = "#D27A15", fill = "#D27A15") +
  xlab("Geographic distance (km)") +
  ylab(fstlab_bir) +
  #ggtitle(expression(italic("F"))) +
  scale_x_continuous(labels = comma) +
  ggtitle("D") +
  theme_emily()

ibd_bir

man_bir
r1_bir
p1_bir

saveRDS(ibd_bir, "figs/ibd_bir.RDS")


#~~ Plotting for manuscript

ibd_alf + ibd_bir

ggsave("figs/IBD.png", ibd_alf + ibd_bir, width = 7, height = 3)

man_alf
r1_alf
p1_alf

man_bir
r1_bir
p1_bir

