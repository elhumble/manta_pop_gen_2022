#~~ Population structure analysis PCA and DAPC

library(adegenet)
library(ggplot2)
library(dplyr)
library(data.table)
library(tidyr)
library(wesanderson)
library(patchwork)
library(SNPRelate)
library(SeqArray)
library(hierfstat)
library(dartR)
source("scripts/theme_emily.R")


#~~~~~~~~~~~~~~~~~~~~~~~~~~#
#      Across species      #
#~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Using polymorphic SNPs for each species from full dataset

#~~ ADEGENET

# scp ehumble@eddie.ecdf.ed.ac.uk:/exports/cmvm/eddie/eb/groups/ogden_grp/emily/manta_pop_gen_2020/data/out/stacks_PE/*ld.raw data/plink/across_sp/

#~~ Read in metadata

meta <- read.csv("data/meta/sample_info_mastersheet_reformatted.csv") %>%
  dplyr::select(c(Sample_code, Species, Location))

#~~~~~~~~~~~~~~~~~~~~~~~~~~#
#       Alfredi PCA        #
#~~~~~~~~~~~~~~~~~~~~~~~~~~#

gl <- read.PLINK("data/plink/across_sp/alfredi_geno4_dp6_ind_geno_depth_no_dups_maf_unrelated_ld.raw")
nLoc(gl)

pca1 <- glPca(gl, nf = nLoc(gl))

pc1 <- pca1$scores[,1]
pc2 <- pca1$scores[,2]
pc3 <- pca1$scores[,3]
pc4 <- pca1$scores[,4]

ind_names <- gl@ind.names

#~~ Read in metadata

ggplot_pca <- as.data.frame(cbind(pc1, pc2, pc3, pc4)) %>%
  mutate(ind_names = ind_names) %>%
  separate(ind_names, c("Filename", "Dup"), sep = "_") %>%
  left_join(meta, by = c("Filename" = "Sample_code")) %>%
  dplyr::select(c(pc1, pc2, pc3, pc4, Filename, Dup, Species, Location)) %>%
  unite(ind_names, c("Filename", "Dup")) %>%
  mutate(ind_names = gsub("_NA", "", ind_names)) %>%
  unite(pop, c("Species", "Location"), remove = F) %>%
  mutate(Location = case_when(Location == "Seychelles" ~"SEY",
                              Location == "Maldives" ~ "MAL",
                              Location == "Chagos" ~ "CHAG",
                              Location == "Australia Pacific" ~ "AP",
                              Location == "Hawaii" ~ "HAW",
                              Location == "Fiji" ~ "FIJI"))

# eig

eig <- data.frame(pca1$eig)
eig$percentage = (eig[, 1]/sum(eig$pca1.eig))*100
sum(eig$percentage)
sum(eig$percentage[1:2])

eig$percentage <- round(eig$percentage, digits = 1)
eig$percentage[1]
eig$percentage[2]
eig$percentage[3]

col_palette <- c("#9986A5",
                 "#9A8822",
                 "#35274A",
                 "#3B9AB2",
                 "#0B775E",
                 "#E2D200")

pc1_pc2 <- ggplot(ggplot_pca, aes(pc1, pc2)) + 
  geom_point(aes(col = factor(Location)), size = 4, alpha = 0.7) +
  theme_emily() +
  scale_color_manual(values = col_palette, name = "Location") +
  theme(axis.text.x = element_text(face = "plain"),
        axis.text.y = element_text(face = "plain"),
        axis.title.x = element_text(face = "plain"),
        axis.title.y = element_text(face = "plain"),
        legend.position = "none",
        plot.title = element_text(face = "italic")) +
  ylab(paste0("PC2 (",eig$percentage[2],"%)")) +
  xlab(paste0("PC1 (",eig$percentage[1],"%)")) +
 # stat_ellipse(aes(pc1, pc2, colour = Location),
  #             type="norm",
   #            lty=2, size=1) +
  ggtitle(expression(italic("Mobula alfredi")))

pc1_pc2


pc1_pc3 <- ggplot(ggplot_pca, aes(pc1, pc3)) + 
  geom_point(aes(col = factor(Location)), size = 4, alpha = 0.7) +
  theme_emily() +
  scale_color_manual(values = col_palette, name = "Location") +
  theme(axis.text.x = element_text(face = "plain"),
        axis.text.y = element_text(face = "plain"),
        axis.title.x = element_text(face = "plain"),
        axis.title.y = element_text(face = "plain")) +
  ylab(paste0("PC3 (",eig$percentage[3],"%)")) +
  xlab(paste0("PC1 (",eig$percentage[1],"%)"))
 # stat_ellipse(aes(pc1, pc3, colour = Location),
  #             type="norm",
   #            lty=2, size=1)


pc1_pc2 + pc1_pc3

ggsave("figs/PC1_PC2_alfredi.png", pc1_pc2 + pc1_pc3, height = 4, width = 9)


ggplot(ggplot_pca, aes(pc2, pc3)) + 
  geom_point(aes(col = factor(Location)), size = 4, alpha = 0.7) +
  theme_emily() +
  scale_color_manual(values = col_palette, name = "Location") +
  theme(axis.text.x = element_text(face = "plain"),
        axis.text.y = element_text(face = "plain"),
        axis.title.x = element_text(face = "plain"),
        axis.title.y = element_text(face = "plain")) +
  ylab(paste0("PC3 (",eig$percentage[3],"%)")) +
  xlab(paste0("PC1 (",eig$percentage[1],"%)"))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#               DAPC               #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Following method from Thia, code modified from Tom Jenkins
# Based on the sampling locations, K was set to 5 (peru and pacific = 1). 
# Therefore, 4 principal components were retained in the DAPC.

gl@pop <- as.factor(ggplot_pca$Location)

# Run a DAPC using site IDs as priors

dapc2 = dapc(gl,
             gl$pop,
             n.pca = 5,
             n.da = 3)

scatter(dapc2, 1, 1, col = col_palette, 
        bg = "white", scree.da = FALSE, 
        legend = TRUE, solid = .6)

scatter(dapc2, scree.pca = TRUE,
        clab = 0, posi.pca = "bottomleft")

# Percent of genetic variance explained by each axis

percent <- dapc2$eig/sum(dapc2$eig)*100
barplot(percent, ylab = "Genetic variance explained by eigenvectors (%)", ylim = c(0,60),
        names.arg = round(percent, 1))

# Create a data.frame containing individual coordinates

ind_coords <- as.data.frame(dapc2$ind.coord)

# Rename columns of dataframe

colnames(ind_coords) <- c("Axis1", "Axis2", "Axis3")

# Add individual names

ind_coords$Ind <- indNames(gl)

# Add locations

ind_coords$Site <- gl$pop

# Calculate centroid (average) position for each population

centroid <- aggregate(cbind(Axis1, Axis2, Axis3) ~ 
                        Site, data = ind_coords, FUN = mean)

centroid$Site <- factor(centroid$Site)

# Add centroid coordinates to ind_coords dataframe

ind_coords <- left_join(ind_coords, centroid, by = "Site", suffix = c("",".cen"))

# Custom x and y labels

xlab <- paste("PC1 (", format(round(percent[1], 1), nsmall=1)," %)", sep="")
ylab <- paste("PC2 (", format(round(percent[2], 1), nsmall=1)," %)", sep="")

# Plotting

alf_dapc_12 <- ggplot(data = ind_coords, aes(x = Axis1, y = Axis2)) +
  geom_segment(aes(xend = Axis1.cen, yend = Axis2.cen, colour = Site), 
               show.legend = FALSE) +
  geom_point(aes(fill = Site, col = Site), size = 3, alpha = 0.7) +
  #geom_label(data = centroid, aes(label = Site, fill = Site, 
  #                                alpha = 0.6), size = 3, show.legend = FALSE) +
  scale_fill_manual(values = col_palette, name = "Location") +
  scale_colour_manual(values = col_palette, name = "Location") +
  labs(x = xlab, y = ylab) +
  #ggtitle("Mobula alfredi") +
  # stat_ellipse(aes(Axis1, Axis2, colour = Site),
  #              type = "norm",
  #              lty = 2, size = 1) +
  theme_emily() +
  ggtitle("B")

alf_dapc_12

ylab_13 <- paste("PC3 (", format(round(percent[3], 1), nsmall=1)," %)", sep="")

alf_dapc_13 <- ggplot(data = ind_coords, aes(x = Axis1, y = Axis3)) +
  geom_segment(aes(xend = Axis1.cen, yend = Axis3.cen, colour = Site), show.legend = FALSE) +
  geom_point(aes(fill = Site, col = Site), size = 4, alpha = 0.7) +
  #geom_label(data = centroid, aes(label = Site, fill = Site, 
  #                                alpha = 0.6), size = 3, show.legend = FALSE) +
  scale_fill_manual(values = col_palette, name = "Location") +
  scale_colour_manual(values = col_palette, name = "Location") +
  labs(x = xlab, y = ylab_13) +
  ggtitle(expression(italic("Mobula alfredi"))) +
# stat_ellipse(aes(Axis1, Axis2, colour = Site),
  #              type = "norm",
  #              lty = 2, size = 1) +
  theme_emily() +
  ggtitle("A")

alf_dapc_13

# Export plot

ggsave("figs/DAPC1_2_alfredi.png", alf_dapc_12 + alf_dapc_13, height = 5, width = 12)

saveRDS(alf_dapc_12, "figs/DAPC_alf.RDS")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#       BIROSTRIS PCA        #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

meta <- read.csv("data/meta/sample_info_mastersheet_reformatted.csv") %>%
  dplyr::select(c(Sample_code, Species, Location))

gl <- read.PLINK("data/plink/across_sp/birostris_geno4_dp6_ind_geno_depth_no_dups_maf_unrelated_ld.raw")
nLoc(gl)
pca1 <- glPca(gl, nf = nLoc(gl))

pc1 <- pca1$scores[,1]
pc2 <- pca1$scores[,2]
pc3 <- pca1$scores[,3]
pc4 <- pca1$scores[,4]

ind_names <- gl@ind.names

#~~ Read in metadata

ggplot_pca <- as.data.frame(cbind(pc1, pc2, pc3, pc4)) %>%
  mutate(ind_names = ind_names) %>%
  separate(ind_names, c("Filename", "Dup"), sep = "_") %>%
  left_join(meta, by = c("Filename" = "Sample_code")) %>%
  dplyr::select(c(pc1, pc2, pc3, pc4, Filename, Dup, Species, Location)) %>%
  unite(ind_names, c("Filename", "Dup")) %>%
  mutate(ind_names = gsub("_NA", "", ind_names)) %>%
  tidyr::unite(pop, c("Species", "Location"), remove = F) %>%
  mutate(Location = case_when(Location == "the Philippines" ~ "The Philippines",
                              TRUE ~ Location)) %>%
  mutate(Location = case_when(Location == "Sri Lanka" ~ "SL",
                              Location == "The Philippines" ~ "PHI",
                              Location == "Peru" ~ "PERU",
                              Location == "Mexico Caribbean" ~ "MC",
                              Location == "Mexico Pacific" ~ "MP",
                              Location == "South Africa" ~ "SA"))
  
# eig

eig <- data.frame(pca1$eig)
eig$percentage = (eig[, 1]/sum(eig$pca1.eig))*100
sum(eig$percentage)
sum(eig$percentage[1:2])

eig$percentage <- round(eig$percentage, digits = 1)
eig$percentage[1]
eig$percentage[2]
eig$percentage[3]

col_palette <- c("#EAD3BF",
                 "#4F1C2E",
                 "#7E8AA2",
                 "#F2B705",
                 "#C16683",
                 "#B53325")

pc1_pc2 <- ggplot(ggplot_pca, aes(pc1, pc2)) + 
  geom_point(aes(col = factor(Location)), size = 4, alpha = 0.8) +
  theme_emily() +
  scale_color_manual(values = col_palette, name = "Location") +
  theme(axis.text.x = element_text(face = "plain"),
        axis.text.y = element_text(face = "plain"),
        axis.title.x = element_text(face = "plain"),
        axis.title.y = element_text(face = "plain"),
        legend.position = "none",
        plot.title = element_text(face = "italic")) +
  ylab(paste0("PC2 (",eig$percentage[2],"%)")) +
  xlab(paste0("PC1 (",eig$percentage[1],"%)")) +
 # stat_ellipse(aes(pc1, pc2, colour = Location),
 #            type="norm",
 #             lty=2, size=1) +
  ggtitle("Mobula birostris")

pc1_pc2

pc1_pc3 <- ggplot(ggplot_pca, aes(pc1, pc3)) + 
  geom_point(aes(col = factor(Location)), size = 4, alpha = 0.8) +
  theme_emily() +
  scale_color_manual(values = col_palette, name = "Location") +
  theme(axis.text.x = element_text(face = "plain"),
        axis.text.y = element_text(face = "plain"),
        axis.title.x = element_text(face = "plain"),
        axis.title.y = element_text(face = "plain")) +  
  #stat_ellipse(aes(pc1, pc2, colour = Location),
  #             type="norm",
  #             lty=2, size=1) +
  ylab(paste0("PC3 (",eig$percentage[3],"%)")) +
  xlab(paste0("PC1 (",eig$percentage[1],"%)"))

pc1_pc2

pc1_pc2 + pc1_pc3

ggsave("figs/PC1_PC2_birostris.png", pc1_pc2 + pc1_pc3, height = 4, width = 10)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#               DAPC               #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Following method from Thia using code modified from Tom Jenkins
# Based on the sampling locations, K was set to 5 (peru and pacific = 1). 
# Therefore, 4 principal components were retained in the DAPC.

gl@pop <- as.factor(ggplot_pca$Location)

# Run a DAPC using site IDs as priors

dapc2 = dapc(gl,
             gl$pop,
             n.pca = 4,
             n.da = 3)

scatter(dapc2, 1, 1, col = col_palette, 
        bg = "white", scree.da = FALSE, 
        legend = TRUE, solid = .6)

scatter(dapc2, scree.pca = TRUE,
        clab = 0, posi.pca = "bottomleft")

# Percent of genetic variance explained by each axis

percent <- dapc2$eig/sum(dapc2$eig)*100
barplot(percent, ylab = "Genetic variance explained by eigenvectors (%)", ylim = c(0,60),
        names.arg = round(percent, 1))

# Create a data.frame containing individual coordinates

ind_coords <- as.data.frame(dapc2$ind.coord)

# Rename columns of dataframe

colnames(ind_coords) <- c("Axis1", "Axis2", "Axis3")

# Add individual names

ind_coords$Ind <- indNames(gl)

# Add locations

ind_coords$Site <- gl$pop

# Calculate centroid (average) position for each population

centroid <- aggregate(cbind(Axis1, Axis2, Axis3) ~ 
                        Site, data = ind_coords, FUN = mean)

centroid$Site

ind_coords$Site <- factor(ind_coords$Site, levels = c("MC",
                                                    "MP",
                                                    "PERU",
                                                    "SA",
                                                    "SL",
                                                    "PHI"))

# Add centroid coordinates to ind_coords dataframe

ind_coords <- left_join(ind_coords, centroid, by = "Site", suffix = c("",".cen"))

# Custom x and y labels

xlab <- paste("PC1 (", format(round(percent[1], 1), nsmall=1)," %)", sep="")
ylab <- paste("PC2 (", format(round(percent[2], 1), nsmall=1)," %)", sep="")

# Plotting

col_palette <- c("#EAD3BF",
                 "#4F1C2E",
                 "#7E8AA2",
                 "#C16683",
                 "#B53325",
                 "#F2B705")

dapc_12 <- ggplot(data = ind_coords, aes(x = Axis1, y = Axis2)) +
  geom_segment(aes(xend = Axis1.cen, yend = Axis2.cen, colour = Site), 
               show.legend = FALSE) +
  geom_point(aes(fill = Site, col = Site), size = 3, alpha = 0.7) +
  #geom_label(data = centroid, aes(label = Site, fill = Site, 
  #                                alpha = 0.6), size = 3, show.legend = FALSE) +
  scale_fill_manual(values = col_palette, name = "Location") +
  scale_colour_manual(values = col_palette, name = "Location") +
  labs(x = xlab, y = ylab) +
  theme_emily() +
  ggtitle("D")

dapc_12

ylab_13 <- paste("PC3 (", format(round(percent[3], 1), nsmall=1)," %)", sep="")

dapc_13 <- ggplot(data = ind_coords, aes(x = Axis1, y = Axis3)) +
  geom_segment(aes(xend = Axis1.cen, yend = Axis3.cen, colour = Site), show.legend = FALSE) +
  geom_point(aes(fill = Site, col = Site), size = 4, alpha = 0.7) +
  #geom_label(data = centroid, aes(label = Site, fill = Site, 
  #                                alpha = 0.6), size = 3, show.legend = FALSE) +
  scale_fill_manual(values = col_palette, name = "Location") +
  scale_colour_manual(values = col_palette, name = "Location") +
  labs(x = xlab, y = ylab_13) +
 # ggtitle("Mobula birostris") +
  # stat_ellipse(aes(Axis1, Axis2, colour = Site),
  #              type = "norm",
  #              lty = 2, size = 1) +
  theme_emily() +
  ggtitle("B")

dapc_12 + dapc_13

# Export plot

ggsave("figs/DAPC1_2_birostris.png", dapc_12 + dapc_13, height = 5, width = 12)

saveRDS(dapc_12, "figs/DAPC_bir.RDS")

alf_dapc_13 + dapc_13


# PC3 plotting

ggsave("figs/DAPC1_3.png", alf_dapc_13 + dapc_13, height = 5, width = 12)

