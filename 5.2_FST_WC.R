library(adegenet)
library(ggplot2)
library(dplyr)
library(data.table)
library(tidyr)
library(wesanderson)
library(patchwork)
library(ape)
library(hierfstat)
library(dartR)
library(forcats)
source("scripts/theme_emily.R")

meta <- read.csv("data/meta/sample_info_mastersheet_reformatted.csv") %>%
  dplyr::select(c(Sample_code, Species, Location))

#~~~~~~~~~~~~~~~~~~~~~~~~~~#
#      Across species      #
#~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Using polymorphic SNPs for each species from full dataset

# scp ehumble@eddie.ecdf.ed.ac.uk:/exports/cmvm/eddie/eb/groups/ogden_grp/emily/manta_pop_gen_2020/data/out/stacks_PE/*raw data/plink/across_sp/

#~~ ALFREDI

gl <- read.PLINK("data/plink/across_sp/alfredi_geno4_dp6_ind_geno_depth_no_dups_maf_unrelated_ld.raw")
nLoc(gl)

inds <- data.frame(gl@ind.names) %>%
  separate(gl.ind.names, c("Filename", "Dup"), sep = "_") %>%
  left_join(meta, by = c("Filename" = "Sample_code")) %>%
  unite(ind_names, c("Filename", "Dup")) %>%
  mutate(ind_names = gsub("_NA", "", ind_names)) %>%
  mutate(Location = case_when(Location == "Fiji" ~ "South Pacific",
                            Location == "Australia Pacific" ~ "South Pacific",
                            TRUE ~ Location)) %>%
  mutate(Location = case_when(Location == "Seychelles" ~"SEY",
                              Location == "Maldives" ~ "MAL",
                              Location == "Chagos" ~ "CHAG",
                              Location == "South Pacific" ~ "SP",
                              Location == "Hawaii" ~ "HAW"))

inds <- data.frame(gl@ind.names) %>%
  separate(gl.ind.names, c("Filename", "Dup"), sep = "_") %>%
  left_join(meta, by = c("Filename" = "Sample_code")) %>%
  unite(ind_names, c("Filename", "Dup")) %>%
  mutate(ind_names = gsub("_NA", "", ind_names)) %>%
  mutate(Location = case_when(Location == "Seychelles" ~"SEY",
                              Location == "Maldives" ~ "MAL",
                              Location == "Chagos" ~ "CHAG",
                              Location == "Australia Pacific" ~ "AP",
                              Location == "Hawaii" ~ "HAW",
                              Location == "Fiji" ~ "FIJ"))

pop(gl) <- inds$Location

pw_fst_alf <- dartR::gl.fst.pop(gl, nboots=1000, percent=95, nclusters=1)
pw_fst_boot_alf <- pw_fst_alf$Bootstraps[,c(1,2,1003:1006)]

# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

pw_fst_alf$Fsts[is.na(pw_fst_alf$Fsts)] <- t(pw_fst_alf$Fsts)[is.na(pw_fst_alf$Fsts)]

#alf_melted_cor <- pw_fst_alf$Fsts[c(5,2,4,1,3),c(5,2,4,1,3)]
#alf_melted_cor <- pw_fst_alf$Fsts[c(6,5,2,4,1,3),c(6,5,2,4,1,3)]
alf_melted_cor <- pw_fst_alf$Fsts[c(5,2,4,1,3),c(5,2,4,1,3)]

alf_upper_tri <- get_upper_tri(alf_melted_cor)

# Melt the correlation matrix
alf_melted_cormat <- melt(alf_upper_tri, na.rm = TRUE)

col_palette <- c("#9986A5",
                 "#9A8822",
                 "#35274A",
                 "#3B9AB2",
                 "#0B775E",
                 "#E2D200",
                 "#EAD3BF",
                 "#FD6467")

alf_fst_plot <- ggplot(alf_melted_cormat, aes(Var1, Var2, fill = value)) +
  geom_tile(color = "white") +
  #scale_fill_gradient2(name="Fst") +
  scale_fill_gradient(low = "#F8FCFC",
                      high = "#1D8981",
                      name = expression(italic("F")[ST])) + #limit = c(-0.0025,0.160)
  theme_emily() +
  theme(axis.line.y = element_blank(),
        axis.line.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  coord_fixed() +
  ggtitle("A")
 # ggtitle(expression(paste("Reef manta ray ", italic("(Mobula alfredi)"))))

alf_fst_plot


saveRDS(alf_fst_plot, "figs/fst_alf.RDS")


#~ Plotting CIs

alf_fst_CI <- pw_fst_boot_alf %>%
  unite(Pairwise, c("Population1", "Population2")) %>%
  mutate(Pairwise = fct_reorder(Pairwise, desc(Fst))) %>%
  mutate(sig = as.factor(case_when(`Lower bound CI limit` < 0 & `Upper bound CI limit` > 0 ~ 1,
                                   TRUE ~ 0))) %>%
  filter(!grepl("AP", Pairwise)) %>%
  ggplot(aes(Fst, Pairwise, col = sig, fill = sig)) +
  geom_pointrange(aes(xmax = `Upper bound CI limit`, xmin = `Lower bound CI limit`), 
                  size = 0.3, 
                  position=position_dodge(0.2)) +
  geom_vline(xintercept = 0) +
  scale_color_manual(values = c("black", "grey40"), name = "Significant") +
  theme_emily() +
  theme(legend.position = "none",
        plot.margin = unit(c(1,0.5,0.5,0.5), "cm"),
        axis.line.y = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(size = 12)) +
  ggtitle("A") + xlab(expression(italic("F")[ST]))

alf_fst_CI


#~~ Some stats

alf_melted_cormat %>%
  mutate(ocean = case_when(Var1 == "SP" & Var2 == "HAW" ~ "within",
                           Var1 == "MAL" & Var2 == "SEY" ~ "within",
                           Var1 == "CHAG" & Var2 == "MAL" ~ "within",
                           Var1 == "CHAG" & Var2 == "SEY" ~ "within",
                           TRUE ~ "between")) %>%
  group_by(ocean) %>%
  summarise(mean = mean(value))

alf_melted_cormat %>%
  mutate(ocean = case_when(Var1 == "FIJ" & Var2 == "HAW" ~ "within",
                           Var1 == "MAL" & Var2 == "SEY" ~ "within",
                           Var1 == "CHAG" & Var2 == "MAL" ~ "within",
                           Var1 == "CHAG" & Var2 == "SEY" ~ "within",
                           TRUE ~ "between")) %>%
  group_by(ocean) %>%
  summarise(mean = mean(value))

alf_melted_cormat %>%
  summarise(mean = mean(value),
            min = min(value),
            max = max(value))

#~~ AMOVA (Code from Kara)

library(poppr)

#~ import to genclone format
genI <- gl2gi(gl)
genC <- as.genclone(genI)

#~ extract the population column to act as strata
pops <- genI$pop
strata(genI) <- data.frame(pops)
p.amova <- poppr.amova(genI, ~pops)

#~ store the amova
amova.df <- data.frame(p.amova$results)
amova.df$Source_of_variation <- row.names(amova.df)
p.amova$componentsofcovariance$Source_of_variation <- row.names(amova.df)

#~ Calcuating p-values
amova.pvalues <- ade4::randtest(p.amova, nrepet = 999)
amova.pvalues

#~~ Combine everything
amova.df <- full_join(amova.df,p.amova$componentsofcovariance)
amova.df <- 
  amova.df %>% 
  mutate(Source = factor(c("Among populations", 
                           "Between samples within populations",
                           "Within samples",
                           "Total"))) %>% 
  
  select(Source, everything()) %>% 
  select(-Source_of_variation)

#~ Phi stat
phi <- 
  p.amova$statphi %>% 
  mutate(Source = factor(c("Within samples",
                           "Between samples within populations",
                           "Among populations"))) 
amova.df <-
  full_join(amova.df, phi) %>% 
  select(-Mean.Sq, -Sigma)

#~ pvalues
amova_df_paper <- 
  data.frame(Source = factor(c("Within samples",
                               "Between samples within populations",
                               "Among populations")),
             Obs = amova.pvalues$obs,
             Std.Obs = amova.pvalues$expvar$Std.Obs,
             Alter = amova.pvalues$alter,
             Pvalue = amova.pvalues$adj.pvalue) %>% 
  full_join(amova.df) %>% 
  select(Source, Df, Sum.Sq, `%`, Phi, Pvalue)

amova_df_paper %>% 
  write.table("figs/amova_alfredi_out.txt",
              row.names=F, quote=F, sep="\t")

#~~~~~~~~~~~~~~~~~~~~~~~~~~#
#        BIROSTRIS         #
#~~~~~~~~~~~~~~~~~~~~~~~~~~#

gl_bir <- read.PLINK("data/plink/across_sp/birostris_geno4_dp6_ind_geno_depth_no_dups_maf_unrelated_ld.raw")
nLoc(gl_bir)

inds_bir <- data.frame(gl_bir@ind.names) %>%
  separate(gl_bir.ind.names, c("Filename", "Dup"), sep = "_") %>%
  left_join(meta, by = c("Filename" = "Sample_code")) %>%
  unite(ind_names, c("Filename", "Dup")) %>%
  mutate(ind_names = gsub("_NA", "", ind_names)) %>%
  mutate(Location = case_when(Location == "Sri Lanka" ~ "SL",
                              Location == "the Philippines" ~ "PHI",
                              Location == "Peru" ~ "PERU",
                              Location == "Mexico Caribbean" ~ "MC",
                              Location == "Mexico Pacific" ~ "MP",
                              Location == "South Africa" ~ "SA"))

pop(gl_bir)
pop(gl_bir) <- inds_bir$Location
pop(gl_bir)

pw_fst_bir <- dartR::gl.fst.pop(gl_bir, nboots=1000, percent=95, nclusters=1)
pw_fst_boot_bir <- pw_fst_bir$Bootstraps[,c(1,2,1003:1006)]

pw_fst_bir$Fsts[is.na(pw_fst_bir$Fsts)] <- t(pw_fst_bir$Fsts)[is.na(pw_fst_bir$Fsts)]

bir_melted_cor <- pw_fst_bir$Fsts[c(1,3,4,2,5),c(1,3,4,2,5)]
#bir_melted_cor <- pw_fst_bir$Fsts[c(1,3,6,4,2,5),c(1,3,6,4,2,5)]

get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

bir_upper_tri <- get_upper_tri(bir_melted_cor)

# Melt the correlation matrix
bir_melted_cor <- melt(bir_upper_tri, na.rm = TRUE)

#~~ Plot

bir_fst_plot <- ggplot(bir_melted_cor, aes(Var1, Var2, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "#FFF2E9",
                      high = "#D27A15",
                      name = expression(italic("F")[ST])) + #limit = c(-0.0025,0.160)
  theme_emily() +
  theme(axis.line.y = element_blank(),
        axis.line.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  coord_fixed() +
  ggtitle("C")
#  ggtitle(expression(paste("Oceanic manta ray ", italic("(Mobula birostris)"))))

bir_fst_plot

saveRDS(bir_fst_plot, "figs/fst_bir.RDS")


#~ Plotting CIs

bir_fst_CI <- pw_fst_boot_bir %>%
  filter(Population2 != "SA") %>%
  unite(Pairwise, c("Population1", "Population2")) %>%
  mutate(Pairwise = fct_reorder(Pairwise, desc(Fst))) %>%
  #mutate(sig = as.factor(case_when(`Lower bound CI limit` < 0 & `Upper bound CI limit` > 0 ~ 1,
  #                       TRUE ~ 0))) %>%
  mutate(sig = as.factor(case_when(`p-value` < 0.05 ~ 1,
                                   TRUE ~ 0))) %>%
  ggplot(aes(Fst, Pairwise, col = sig, fill = sig)) +
  geom_pointrange(aes(xmax = `Upper bound CI limit`, xmin = `Lower bound CI limit`), 
                  size = 0.3, 
                  position=position_dodge(0.2)) +
  geom_vline(xintercept = 0) +
  scale_color_manual(values = c("grey", "black"), name = "Significant") +
  theme_emily() +
  theme(legend.position = "none",
        plot.margin = unit(c(1,0.5,0.5,0.5), "cm"),
        axis.line.y = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(size = 12)) +
  ggtitle("B") + xlab(expression(italic("F")[ST]))

bir_fst_CI

alf_fst_plot + bir_fst_plot

ggsave("figs/FST_mat_WC.png", alf_fst_plot / bir_fst_plot, width = 10, height = 6)

alf_fst_CI + bir_fst_CI

ggsave("figs/FST_CI_WC.png", alf_fst_CI + bir_fst_CI, width = 7, height = 5)


#~~ Some stats

bir_melted_cor %>%
  summarise(mean = mean(value),
            min = min(value),
            max = max(value))


#~~ AMOVA

#~ import to genclone format

genI <- gl2gi(gl_bir)
genC <- as.genclone(genI)

#~ extract the population column to act as strata

pops <- genI$pop
strata(genI) <- data.frame(pops)
p.amova <- poppr.amova(genI, ~pops)

#~ store the amova

amova.df <- data.frame(p.amova$results)
amova.df$Source_of_variation <- row.names(amova.df)
p.amova$componentsofcovariance$Source_of_variation <- row.names(amova.df)

#~ Calcuating p-values

amova.pvalues <- ade4::randtest(p.amova, nrepet = 999)
amova.pvalues

#~~ Combine everything

amova.df <- full_join(amova.df,p.amova$componentsofcovariance)
amova.df <- 
  amova.df %>% 
  mutate(Source = factor(c("Among populations", 
                           "Between samples within populations",
                           "Within samples",
                           "Total"))) %>% 
  
  select(Source, everything()) %>% 
  select(-Source_of_variation)

#~ Phi stat

phi <- 
  p.amova$statphi %>% 
  mutate(Source = factor(c("Within samples",
                           "Between samples within populations",
                           "Among populations"))) 
amova.df <-
  full_join(amova.df, phi) %>% 
  select(-Mean.Sq, -Sigma)

#~ pvalues

amova_df_paper <- 
  data.frame(Source = factor(c("Within samples",
                               "Between samples within populations",
                               "Among populations")),
             Obs = amova.pvalues$obs,
             Std.Obs = amova.pvalues$expvar$Std.Obs,
             Alter = amova.pvalues$alter,
             Pvalue = amova.pvalues$adj.pvalue) %>% 
  full_join(amova.df) %>% 
  select(Source, Df, Sum.Sq, `%`, Phi, Pvalue)

amova_df_paper %>% 
  write.table("figs/amova_birostris_out.txt",
              row.names=F, quote=F, sep="\t")

