# Plotting admixture results

library(purrr)
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
library(wesanderson)
library(patchwork)
library(readxl)
source("scripts/theme_emily.R")

#scp ehumble@eddie.ecdf.ed.ac.uk:/exports/cmvm/eddie/eb/groups/ogden_grp/emily/manta_pop_gen_2020/data/out/stacks_PE/admixture/*admixture.fam data/admixture/
#scp ehumble@eddie.ecdf.ed.ac.uk:/exports/cmvm/eddie/eb/groups/ogden_grp/emily/manta_pop_gen_2020/data/out/stacks_PE/admixture/*admixture.*.Q data/admixture/
#scp ehumble@eddie.ecdf.ed.ac.uk:/exports/cmvm/eddie/eb/groups/ogden_grp/emily/manta_pop_gen_2020/data/out/stacks_PE/admixture/*error data/admixture/
  
#~~~~~~~~~~~~~~~~~~~~~#
#      Birostris      #
#~~~~~~~~~~~~~~~~~~~~~#
  
data_path <- "data/admixture/"
files <- dir(data_path, pattern = "*Q$")

data <- tibble(filename = files[9:16]) %>%
  dplyr::mutate(file_contents = purrr::map(filename,
                                    ~ fread(file.path(data_path, .))))

# add individual names

inds <- fread("data/admixture/birostris_admixture.fam") %>%
  dplyr::select(V1) %>%
  rename(ID = V1)

# Get population info

meta <- read.csv("data/meta/sample_info_mastersheet_reformatted.csv") %>%
  dplyr::select(Sample_code, Species, Location, Site)

data <- unnest(data, cols = c(file_contents)) %>%
  mutate(ID = rep(inds$ID, 8)) %>% # number of Ks
  pivot_longer(cols = c(V1, V2, V3, V4, V5, V6, V7, V8)) %>%
  drop_na() %>%
  separate(ID, c("Filename", "Dup"), sep = "_") %>%
  left_join(meta, by = c("Filename" = "Sample_code")) %>%
  unite(ID, c("Filename", "Dup")) %>%
  mutate(ID = gsub("_NA", "", ID))  %>%
  mutate(Location = case_when(Location == "Sri Lanka" ~ "SL",
                              Location == "the Philippines" ~ "PHI",
                              Location == "Peru" ~ "PERU",
                              Location == "Mexico Caribbean" ~ "MC",
                              Location == "Mexico Pacific" ~ "MP",
                              Location == "South Africa" ~ "SA")) %>%
  mutate(Location = factor(Location, levels = c("PERU", "MP", "MC", "SA", "SL", "PHI")))

K1 <- filter(data, grepl("1.Q$", filename)) %>%
  mutate(ID = as.factor(ID),
         name = as.factor(name))

K2 <- filter(data, grepl("2.Q$", filename)) %>%
  mutate(ID = as.factor(ID),
         name = as.factor(name))

K3 <- filter(data, grepl("3.Q$", filename)) %>%
  mutate(ID = as.factor(ID),
         name = as.factor(name))

K4 <- filter(data, grepl("4.Q$", filename)) %>%
  mutate(ID = as.factor(ID),
         name = as.factor(name))

K5 <- filter(data, grepl("5.Q$", filename)) %>%
  mutate(ID = as.factor(ID),
         name = as.factor(name))

K6 <- filter(data, grepl("6.Q$", filename)) %>%
  mutate(ID = as.factor(ID),
         name = as.factor(name))

K7 <- filter(data, grepl("7.Q$", filename)) %>%
  mutate(ID = as.factor(ID),
         name = as.factor(name))

K8 <- filter(data, grepl("8.Q$", filename)) %>%
  mutate(ID = as.factor(ID),
         name = as.factor(name))
# Plot

# https://luisdva.github.io/rstats/model-cluster-plots/

col_palette <- c("#7E8AA2",
                 "#F2B705",
                 "#C16683",
                 "#EAD3BF",
                 "#4F1C2E",
                 "#B53325",
                 "#E2D200",
                 "#3B9AB2")

k1_plot <- ggplot(K1, aes(factor(ID), value, fill = factor(name))) +
  geom_col(color = "gray", size = 0.1) +
  facet_grid(~Location, switch = "x", scales = "free", space = "free") +
  theme_minimal() + labs(x = "Individuals", title = "K=1", y = "Ancestry") +
  scale_fill_manual(values = col_palette) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = expansion(add = 1)) +
  theme(panel.spacing.x = unit(0.1, "lines"),
        axis.text.x =  element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        strip.text.x = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none")

k2_plot <- ggplot(K2, aes(factor(ID), value, fill = factor(name))) +
  geom_col(color = "gray", size = 0.1) +
  facet_grid(~Location, switch = "x", scales = "free", space = "free") +
  theme_minimal() + labs(x = "Individuals", title = "K=2", y = "Ancestry") +
  scale_fill_manual(values = col_palette) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = expansion(add = 1)) +
  theme(panel.spacing.x = unit(0.1, "lines"),
        axis.text.x =  element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        strip.text.x = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none")


k3_plot <- ggplot(K3, aes(factor(ID), value, fill = factor(name))) +
  geom_col(color = "gray", size = 0.1) +
  facet_grid(~Location, switch = "x", scales = "free", space = "free") +
  theme_minimal() + labs(x = "Individuals", title = "K=3", y = "Ancestry") +
  scale_fill_manual(values = col_palette) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = expansion(add = 1)) +
  theme(panel.spacing.x = unit(0.1, "lines"),
        axis.text.x =  element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        strip.text.x = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none")

k4_plot <- ggplot(K4, aes(factor(ID), value, fill = factor(name))) +
  geom_col(color = "gray", size = 0.1) +
  facet_grid(~Location, switch = "x", scales = "free", space = "free") +
  theme_minimal() + labs(x = "Individuals", title = "K=4", y = "Ancestry") +
  scale_fill_manual(values = col_palette) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = expansion(add = 1)) +
  theme(panel.spacing.x = unit(0.1, "lines"),
        axis.text.x =  element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        strip.text.x = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none")

k5_plot <- ggplot(K5, aes(factor(ID), value, fill = factor(name))) +
  geom_col(color = "gray", size = 0.1) +
  facet_grid(~Location, switch = "x", scales = "free", space = "free") +
  theme_minimal() + labs(x = "Individuals", title = "K=5", y = "Ancestry") +
  scale_fill_manual(values = col_palette) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = expansion(add = 1)) +
  theme(panel.spacing.x = unit(0.1, "lines"),
        axis.text.x =  element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        strip.text.x = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none")

k6_plot <- ggplot(K6, aes(factor(ID), value, fill = factor(name))) +
  geom_col(color = "gray", size = 0.1) +
  facet_grid(~Location, switch = "x", scales = "free", space = "free") +
  theme_minimal() + labs(x = "Individuals", title = "K=6", y = "Ancestry") +
  scale_fill_manual(values = col_palette) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = expansion(add = 1)) +
  theme(panel.spacing.x = unit(0.1, "lines"),
        axis.text.x =  element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        strip.text.x = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none")


k7_plot <- ggplot(K7, aes(factor(ID), value, fill = factor(name))) +
  geom_col(color = "gray", size = 0.1) +
  facet_grid(~fct_inorder(Location), switch = "x", scales = "free", space = "free") +
  theme_minimal() + labs(x = "Individuals", title = "K=7", y = "Ancestry") +
  scale_fill_manual(values = col_palette) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = expansion(add = 1)) +
  theme(panel.spacing.x = unit(0.1, "lines"),
        axis.text.x =  element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        strip.text.x = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none")

k8_plot <- ggplot(K8, aes(factor(ID), value, fill = factor(name))) +
  geom_col(color = "gray", size = 0.1) +
  facet_grid(~Location, switch = "x", scales = "free", space = "free") +
  theme_minimal() + labs(x = "Individuals", title = "K=8", y = "Ancestry") +
  scale_fill_manual(values = col_palette) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = expansion(add = 1)) +
  theme(panel.spacing.x = unit(0.1, "lines"),
        #axis.text.x =  element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none")

k1_plot + k2_plot + k3_plot + k4_plot + k5_plot + k6_plot + k7_plot + k8_plot + plot_layout(ncol = 1)


ggsave("figs/admixture_birostris.png",
       k1_plot + k2_plot + k3_plot + k4_plot + k5_plot + k6_plot + k7_plot + k8_plot + plot_layout(ncol = 1),
       width = 5, height = 9)

# Cross validation

cv_err <- fread("data/admixture/birostris_cv_error")

cv_error <- ggplot(cv_err, aes(V1, V2)) +
  geom_point() +
  geom_line() +
  ylab("Cross-validation error") +
  xlab("K") +
  theme_emily() + ggtitle("B")

cv_error

ggsave("figs/admixture_birostris_cv.png",
       cv_error,
       width = 5, height = 4)

#~~~~~~~~~~~~~~~~~~~#
#      Alfredi      #
#~~~~~~~~~~~~~~~~~~~#

data <- tibble(filename = files[1:8]) %>%
  dplyr::mutate(file_contents = purrr::map(filename,
                                           ~ fread(file.path(data_path, .))))

# add individual names

inds <- fread("data/admixture/alfredi_admixture.fam") %>%
  dplyr::select(V1) %>%
  rename(ID = V1)

# Get population info

data <- unnest(data, cols = c(file_contents)) %>%
  mutate(ID = rep(inds$ID, 8)) %>% # number of Ks
  pivot_longer(cols = c(V1, V2, V3, V4, V5, V6, V7, V8)) %>%
  drop_na() %>%
  separate(ID, c("Filename", "Dup"), sep = "_") %>%
  left_join(meta, by = c("Filename" = "Sample_code")) %>%
  unite(ID, c("Filename", "Dup")) %>%
  mutate(ID = gsub("_NA", "", ID)) %>%
  mutate(Location = case_when(Location == "Seychelles" ~"SEY",
                              Location == "Maldives" ~ "MAL",
                              Location == "Chagos" ~ "CHA",
                              Location == "Australia Pacific" ~ "AP",
                              Location == "Hawaii" ~ "HAW",
                              Location == "Fiji" ~ "FIJI")) %>%
  mutate(Location = factor(Location, levels = c("MAL", "SEY", "CHA", "AP", "FIJI", "HAW")))

K1 <- filter(data, grepl("1.Q$", filename)) %>%
  mutate(ID = as.factor(ID),
         name = as.factor(name))

K2 <- filter(data, grepl("2.Q$", filename)) %>%
  mutate(ID = as.factor(ID),
         name = as.factor(name))

K3 <- filter(data, grepl("3.Q$", filename)) %>%
  mutate(ID = as.factor(ID),
         name = as.factor(name))

K4 <- filter(data, grepl("4.Q$", filename)) %>%
  mutate(ID = as.factor(ID),
         name = as.factor(name))

K5 <- filter(data, grepl("5.Q$", filename)) %>%
  mutate(ID = as.factor(ID),
         name = as.factor(name))

K6 <- filter(data, grepl("6.Q$", filename)) %>%
  mutate(ID = as.factor(ID),
         name = as.factor(name))

K7 <- filter(data, grepl("7.Q$", filename)) %>%
  mutate(ID = as.factor(ID),
         name = as.factor(name))

K8 <- filter(data, grepl("8.Q$", filename)) %>%
  mutate(ID = as.factor(ID),
         name = as.factor(name))
# Plot

# https://luisdva.github.io/rstats/model-cluster-plots/

col_palette <- c("#9986A5",
                 "#9A8822",
                 "#35274A",
                 "#3B9AB2",
                 "#0B775E",
                 "#E2D200",
                 "#EAD3BF",
                 "#FD6467")

k1_plot <- ggplot(K1, aes(factor(ID), value, fill = factor(name))) +
  geom_col(color = "gray", size = 0.1) +
  facet_grid(~Location, switch = "x", scales = "free", space = "free") +
  theme_minimal() + labs(x = "Individuals", title = "K=1", y = "Ancestry") +
  scale_fill_manual(values = col_palette) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = expansion(add = 1)) +
  theme(panel.spacing.x = unit(0.1, "lines"),
        axis.text.x =  element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        strip.text.x = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none")

k2_plot <- ggplot(K2, aes(factor(ID), value, fill = factor(name))) +
  geom_col(color = "gray", size = 0.1) +
  facet_grid(~Location, switch = "x", scales = "free", space = "free") +
  theme_minimal() + labs(x = "Individuals", title = "K=2", y = "Ancestry") +
  scale_fill_manual(values = col_palette) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = expansion(add = 1)) +
  theme(panel.spacing.x = unit(0.1, "lines"),
        axis.text.x =  element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        strip.text.x = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none")

k3_plot <- ggplot(K3, aes(factor(ID), value, fill = factor(name))) +
  geom_col(color = "gray", size = 0.1) +
  facet_grid(~Location, switch = "x", scales = "free", space = "free") +
  theme_minimal() + labs(x = "Individuals", title = "K=3", y = "Ancestry") +
  scale_fill_manual(values = col_palette) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = expansion(add = 1)) +
  theme(panel.spacing.x = unit(0.1, "lines"),
        axis.text.x =  element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        strip.text.x = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none")

k4_plot <- ggplot(K4, aes(factor(ID), value, fill = factor(name))) +
  geom_col(color = "gray", size = 0.1) +
  facet_grid(~Location, switch = "x", scales = "free", space = "free") +
  theme_minimal() + labs(x = "Individuals", title = "K=4", y = "Ancestry") +
  scale_fill_manual(values = col_palette) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = expansion(add = 1)) +
  theme(panel.spacing.x = unit(0.1, "lines"),
        axis.text.x =  element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        strip.text.x = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none")

k5_plot <- ggplot(K5, aes(factor(ID), value, fill = factor(name))) +
  geom_col(color = "gray", size = 0.1) +
  facet_grid(~Location, switch = "x", scales = "free", space = "free") +
  theme_minimal() + labs(x = "Individuals", title = "K=5", y = "Ancestry") +
  scale_fill_manual(values = col_palette) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = expansion(add = 1)) +
  theme(panel.spacing.x = unit(0.1, "lines"),
        axis.text.x =  element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        strip.text.x = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none")

k6_plot <- ggplot(K6, aes(factor(ID), value, fill = factor(name))) +
  geom_col(color = "gray", size = 0.1) +
  facet_grid(~Location, switch = "x", scales = "free", space = "free") +
  theme_minimal() + labs(x = "Individuals", title = "K=6", y = "Ancestry") +
  scale_fill_manual(values = col_palette) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = expansion(add = 1)) +
  theme(panel.spacing.x = unit(0.1, "lines"),
        axis.text.x =  element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        strip.text.x = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none")


k7_plot <- ggplot(K7, aes(factor(ID), value, fill = factor(name))) +
  geom_col(color = "gray", size = 0.1) +
  facet_grid(~Location, switch = "x", scales = "free", space = "free") +
  theme_minimal() + labs(x = "Individuals", title = "K=7", y = "Ancestry") +
  scale_fill_manual(values = col_palette) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = expansion(add = 1)) +
  theme(panel.spacing.x = unit(0.1, "lines"),
        axis.text.x =  element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        strip.text.x = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none")

k8_plot <- ggplot(K8, aes(factor(ID), value, fill = factor(name))) +
  geom_col(color = "gray", size = 0.1) +
  facet_grid(~Location, switch = "x", scales = "free", space = "free") +
  theme_minimal() + labs(x = "Individuals", title = "K=8", y = "Ancestry") +
  scale_fill_manual(values = col_palette) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = expansion(add = 1)) +
  theme(panel.spacing.x = unit(0.1, "lines"),
        #axis.text.x =  element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none")


k1_plot + k2_plot + k3_plot + k4_plot + k5_plot + k6_plot + k7_plot + k8_plot + plot_layout(ncol = 1)


ggsave("figs/admixture_alfredi.png",
       k1_plot + k2_plot + k3_plot + k4_plot + k5_plot + k6_plot + k7_plot + k8_plot + plot_layout(ncol = 1),
       width = 5, height = 9)

# Cross validation

cv_err <- fread("data/admixture/alfredi_cv_error")

cv_error_alf <- ggplot(cv_err, aes(V1, V2)) +
  geom_point() +
  geom_line() +
  ylab("Cross-validation error") +
  xlab("K") +
  theme_emily() + ggtitle("A")

cv_error_alf

ggsave("figs/admixture_alfredi_cv.png",
       cv_error_alf,
       width = 5, height = 4)

