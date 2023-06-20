library(data.table)
library(dplyr)
library(purrr)
library(tidyr)
library(ggplot2)
library(RColorBrewer)
library(readxl)
source("scripts/theme_emily.R")
library(wesanderson)
library(stringr)
library(patchwork)


# Script to detect related individuals

#~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#       Across species      #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~ Read in genome and NGSrelate output

# scp ehumble@eddie.ecdf.ed.ac.uk:/exports/cmvm/eddie/eb/groups/ogden_grp/emily/manta_pop_gen_2020/data/out/stacks_PE/relate/*res data/relate/
# scp ehumble@eddie.ecdf.ed.ac.uk:/exports/cmvm/eddie/eb/groups/ogden_grp/emily/manta_pop_gen_2020/data/out/stacks_PE/relate/*genome data/relate/
  
# ALFREDI
  
gen <- fread("data/relate/alfredi_relatedness.genome", header = T)

summary(gen$PI_HAT)
hist(gen$PI_HAT)

#~~ Combine with ngsrelate output

ngsrel <- fread("data/relate/alfredi_relatedness.res")

# individual order is the same

gen$R1 <- ngsrel$R1
gen$R0 <- ngsrel$R0
gen$KING <- ngsrel$KING


relate <- ggplot(gen, aes(R1, KING)) +
  geom_point(size = 2, alpha = 0.3, shape = 21,
             aes(fill = ifelse(KING > 1/2^(5/2) & R1 > 0.9, "red", "black"),
                 col = ifelse(KING > 1/2^(5/2) & R1 > 0.9, "red", "black"))) +
  scale_colour_manual(values=c("black", "red")) + 
  scale_fill_manual(values=c("black", "red")) + 
  theme_emily() +
  scale_x_log10(limit = c(-0.1,10)) +
  theme(legend.title = element_blank(),
        legend.position = "none") +
  guides(colour = guide_legend(override.aes = list(alpha=1),
                               keyheight = 0.1,
                               default.unit = "inch")) +
  geom_hline(yintercept = 1/2^(5/2), linetype = "dashed", alpha = 0.3) +
  geom_vline(xintercept = 0.9, linetype = "dashed", alpha = 0.3) +
  ggtitle("A")

relate

#~~~ Individuals to remove

filter(gen, R1 > 0.9 & KING > 1/2^(5/2))

# Average relatedness in remaining sample
filter(gen, R1 < 0.9 & KING < 1/2^(5/2)) %>%
  summarise(mean(PI_HAT))

# 1352 (Hawaii) / 1355 (Hawaii)
# 1352 (Hawaii) / 1356 (Hawaii)
# 1355 (Hawaii) / 1356 (Hawaii)
# 1293 (Maldives) / 1292 (Maldives)

# Remove the following, manually created file on EDDIE:

# 1355
# 1356
# 1293

ggsave("figs/R1_KING_across_sp_alfredi.png", relate, height = 4, width = 5)


# BIROSTRIS

gen <- fread("data/relate/birostris_relatedness.genome", header = T)

summary(gen$PI_HAT)
hist(gen$PI_HAT)

#~~ Combine with ngsrelate output

ngsrel <- fread("data/relate/birostris_relatedness.res")

# individual order is the same

gen$R1 <- ngsrel$R1
gen$R0 <- ngsrel$R0
gen$KING <- ngsrel$KING


relate <- ggplot(gen, aes(R1, KING)) +
  geom_point(size = 2, alpha = 0.3, shape = 21,
             aes(fill = ifelse(KING > 1/2^(5/2) & R1 > 0.9, "red", "black"),
                 col = ifelse(KING > 1/2^(5/2) & R1 > 0.9, "red", "black"))) +
  scale_colour_manual(values=c("black", "red")) + 
  scale_fill_manual(values=c("black", "red")) + 
  theme_emily() +
  scale_x_log10() +
  theme(legend.title = element_blank(),
        legend.position = "none") +
  guides(colour = guide_legend(override.aes = list(alpha=1),
                               keyheight = 0.1,
                               default.unit = "inch")) +
  geom_hline(yintercept = 1/2^(5/2), linetype = "dashed", alpha = 0.3) +
  geom_vline(xintercept = 0.9, linetype = "dashed", alpha = 0.3) +
  ggtitle("B")

relate

#~~~ Individuals to remove

filter(gen, R1 > 0.9 & KING > 1/2^(5/2))
filter(gen,KING > 1/2^(5/2)) # removing all outliers

# Average relatedness in remaining sample
filter(gen, KING < 1/2^(5/2)) %>%
  summarise(mean(PI_HAT))

# 1061 (Peru) / 1060 (Peru)
# 0728 (SL) / 0729 (SL)
# 1121 (Philippines) / 1122_lib3 (Philippines)
# 1077 (SL) / 1078 (SL)
# 1139 (Philippines) / 1137 (Philippines)

# Remove the following, manually created file on EDDIE:

# 0728
# 1077
# 1061
# 1121
# 1137

ggsave("figs/R1_KING_across_sp_birostris.png", relate, height = 4, width = 5)


