# Heterozygosity

library(inbreedR)
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
source("scripts/theme_emily.R")
library(wesanderson)
library(gghalves)
library(forcats)
library(ggsignif)
library(RColorBrewer)
library(stringr)
library(patchwork)
library(lme4)
library(broom.mixed)
library(performance)
library(sjPlot)

#~~ Function to get MLH values from raw plink file

get_MLH_from_plinkraw <- function(file) {
  
  x <- fread(file, colClasses = "character")
  ids <- x$IID
  x <- dplyr::select(x, contains("HET"))
  #row.names(x) <- ids
  NAs <- apply(x, 1, function(x) sum(is.na(x)))
  
  MLH <- as.data.frame(MLH(x))
  MLH$ANIMAL <- ids
  MLH$NAS <- NAs
  colnames(MLH) <- c("MLH", "ANIMAL", "NAs")
  MLH <- MLH
  
}

# Run MLH on PLINK raw

# scp ehumble@eddie.ecdf.ed.ac.uk:/exports/cmvm/eddie/eb/groups/ogden_grp/emily/manta_pop_gen_2020/data/out/stacks_PE/across_sp_geno4_dp6_ind_geno_depth_no_dups_mac_geno.raw data/plink/across_sp
  
raw_file <- "data/plink/across_sp/across_sp_geno4_dp6_ind_geno_depth_no_dups_mac_geno.raw"

MLH <- get_MLH_from_plinkraw(raw_file)

#~~ Combine with metadata

meta <- read.csv("data/meta/sample_info_mastersheet_reformatted.csv") %>%
  dplyr::select(Sample_code, Species, Location) %>%
  mutate(Location = gsub("the", "The", Location))

MLH <- MLH %>%
  separate(ANIMAL, c("ID", "dup"), sep = "_") %>%
  left_join(meta, by = c("ID" = "Sample_code")) %>%
  mutate(Location = as.factor(Location)) %>%
  mutate(Location = case_when(Location == "Sri Lanka" ~ "SL",
                              Location == "The Philippines" ~ "PHI",
                              Location == "Peru" ~ "PERU",
                              Location == "Mexico Caribbean" ~ "MC",
                              Location == "Mexico Pacific" ~ "MP",
                              Location == "South Africa" ~ "SA",
                              Location == "Seychelles" ~"SEY",
                              Location == "Maldives" ~ "MAL",
                              Location == "Chagos" ~ "CHAG",
                              Location == "Australia Pacific" ~ "AP",
                              Location == "Hawaii" ~ "HAW",
                              Location == "Fiji" ~ "FIJI")) %>%
  mutate(species_abb = case_when(Species == "Mobula alfredi" ~ "M. alfredi",
                                 Species == "Mobula birostris" ~ "M. birostris"))

MLH %>% filter(Species == "Mobula birostris")

pal <- c("#53A49D", "#D27A15")

# Species level

sp_MLH <- ggplot(filter(MLH, species_abb == "M. alfredi" | species_abb == "M. birostris"), 
       aes(as.factor(species_abb), MLH, fill = species_abb)) +
  geom_half_point(side = "l", shape = 21, alpha = 0.5, stroke = 0.1, size = 3,
                  transformation_params = list(height = 0, width = 1.3, seed = 1)) +
  geom_half_boxplot(side = "r", outlier.color = NA,
                    width = 0.6, lwd = 0.3, color = "black",
                    alpha = 0.8) +
  theme_emily() +
  scale_fill_manual(values = pal) + 
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x=element_text(face=c("italic"))) +
  labs(x = "Species", y = "Multi-locus heterozygosity") +
  ggtitle("C")

sp_MLH

ggsave("figs/MLH_species.png", sp_MLH, width = 5, height = 4)

# NAs

MLH <- MLH %>%
  mutate(prop_miss = NAs/15312)

ggplot(filter(MLH, species_abb == "M. alfredi" | species_abb == "M. birostris"), 
                 aes(as.factor(species_abb), prop_miss, fill = species_abb)) +
  geom_half_point(side = "l", shape = 21, alpha = 0.5, stroke = 0.1, size = 3,
                  transformation_params = list(height = 0, width = 1.3, seed = 1)) +
  geom_half_boxplot(side = "r", outlier.color = NA,
                    width = 0.6, lwd = 0.3, color = "black",
                    alpha = 0.8) +
  theme_emily() +
  scale_fill_manual(values = pal) + 
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x=element_text(face=c("italic")))

mod <- lm(prop_miss ~ species_abb, data = filter(MLH, species_abb == "M. alfredi" | species_abb == "M. birostris"))
tidy(mod, conf.int = TRUE)
summary.lm(mod)



mean(filter(MLH, species_abb == "M. alfredi")$prop_miss)
mean(filter(MLH, species_abb == "M. birostris")$prop_miss)

#~~ Population level

pal <- c("#9986A5",
         "#9A8822",
         "#35274A",
         "#3B9AB2",
         "#0B775E",
         "#E2D200")

alf_MLH <- ggplot(filter(MLH, species_abb == "M. alfredi"), 
                  aes(fct_reorder(Location, MLH, .desc = T),
                      MLH, fill = Location)) +
  geom_half_point(side = "l", shape = 21, alpha = 0.5, stroke = 0.1, size = 3, 
                  transformation_params = list(height = 0, width = 1.3, seed = 1)) +
  geom_half_boxplot(side = "r", outlier.color = NA, 
                    width = 0.6, lwd = 0.3, color = "black",
                    alpha = 0.8) + 
  theme_emily() + 
  scale_fill_manual(values = pal) + 
  theme(legend.position = "none") +
  labs(x = "Location", y = "") +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) +
 # ggtitle(expression(italic("M. alfredi"))) +
  ggtitle("D")

alf_MLH

ggsave("figs/MLH_alfredi.png", alf_MLH, width = 5, height = 4)

bir_pal <- c("#EAD3BF",
             "#4F1C2E",
             "#7E8AA2",
             "#F2B705",
             "#C16683",
             "#B53325")

bir_MLH <- ggplot(filter(MLH, Species == "Mobula birostris"), 
       aes(fct_reorder(Location, MLH, .desc = T), 
           MLH, fill = Location)) +
  geom_half_point(side = "l", shape = 21, alpha = 0.5, stroke = 0.1, size = 3,
                  transformation_params = list(height = 0, width = 1.3, seed = 1)) +
  geom_half_boxplot(side = "r", outlier.color = NA,
                    width = 0.6, lwd = 0.3, color = "black",
                    alpha = 0.8) + 
  theme_emily() + 
  scale_fill_manual(values = bir_pal) + 
  theme(legend.position = "none") +
  labs(x = "Location", y = "") +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) +
 # ggtitle(expression(italic("M. birostris")))
  ggtitle("E")

bir_MLH

ggsave("figs/MLH_birostris.png", bir_MLH, width = 5, height = 4)

sp_MLH + bir_MLH + alf_MLH

ggsave("figs/MLH_all.png", sp_MLH + alf_MLH + bir_MLH +
         plot_layout(nrow = 1, widths = c(0.7,1.3,1.3)),
                     width = 10, height = 3.5)

ggsave("figs/MLH_all.pdf", sp_MLH + alf_MLH + bir_MLH +
         plot_layout(nrow = 1, widths = c(0.7,1.3,1.3)),
       width = 10, height = 3.5, dev = cairo_pdf)

#~~ Some stats

# Models for MLH

fit_het <- MLH %>% filter(Species == "Mobula birostris" | 
                                  Species == "Mobula alfredi") %>%
  lm(MLH ~ Species, data = .)

tidy(fit_het, conf.int = TRUE)
check_model(fit_het)
plot_model(fit_het)
summary.lm(fit_het)

# Mean het in birostris

MLH %>%
  filter(Species == "Mobula birostris") %>%
  summarise(mean = mean(MLH),
            min = min(MLH),
            max = max(MLH))

# Mean het in alfredi

MLH %>%
  filter(Species == "Mobula alfredi") %>%
  summarise(mean = mean(MLH),
            min = min(MLH),
            max = max(MLH))

# Mean het in alfredi by pop

MLH %>%
  filter(Species == "Mobula alfredi") %>%
  group_by(Location) %>%
  summarise(mean = mean(MLH),
            min = min(MLH),
            max = max(MLH))

# Min and max of mean hets in birostris

MLH %>%
  filter(Species == "Mobula birostris") %>%
  group_by(Location) %>%
  summarise(mean = mean(MLH),
            min = min(MLH),
            max = max(MLH)) %>%
  ungroup() %>%
  summarise(min_mean = min(mean),
            max_mean = max(mean))

# Mean var within pops

MLH %>%
  filter(Species == "Mobula alfredi") %>%
  group_by(Location) %>%
  summarise(var = var(MLH)) %>%
  ungroup() %>%
  summarise(mean_var = mean(var))
  
MLH %>%
  filter(Species == "Mobula birostris") %>%
  group_by(Location) %>%
  summarise(var = var(MLH)) %>%
  ungroup() %>%
  summarise(mean_var = mean(var))
