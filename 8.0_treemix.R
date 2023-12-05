library(RColorBrewer)
library(R.utils)
library(OptM)
library(ggplot2)
library(dplyr)
library(purrr)
library(tidyr)
library(patchwork)
source("scripts/plotting_funcs.R")
source("scripts/treemix_bootstrap.R")
source("scripts/newick_split.R")
source("scripts/theme_emily.R")


#~~ Get optimal number of migration edges from initial treemix runs:

# With package optM

# Used to be that you stopped adding migration edges when 99.8% variation in data explained
# optM automates the process using an ad hoc statistic based on the second order rate of 
# change in the loglik

# ALFREDI

# scp ehumble@eddie.ecdf.ed.ac.uk:/exports/cmvm/eddie/eb/groups/ogden_grp/emily/manta_pop_gen_2020/data/out/stacks_PE/treemix/initial_runs/alfredi* data/treemix/alfredi/initial_runs/

alfredi <- "data/treemix/alfredi/initial_runs"
alfredi.optM <- optM(alfredi)

plot_optM(alfredi.optM, method = "Evanno")

alfredi_test.linear <- optM(alfredi, method = "linear", tsv="linear.txt")
plot_optM(alfredi_test.linear, method = "linear")

plot_tree(cex=0.9,"data/treemix/alfredi/initial_runs/alfredi_1_m2")

# Variance explained only. Determined using treemix_fraction.pl script

# scp ehumble@eddie.ecdf.ed.ac.uk:/exports/cmvm/eddie/eb/groups/ogden_grp/emily/manta_pop_gen_2020/data/out/stacks_PE/treemix/initial_runs/var_alfredi* data/treemix/alfredi/initial_runs/

#~~ Alfredi

data_path <- "data/treemix/alfredi/initial_runs/"
files <- dir(data_path, pattern = "var*")

options(scipen=999)

var_alf <- tibble(filename = files) %>%
  mutate(file_contents = map(filename, ~ scan(file.path(data_path, .), what = "",
                                              sep = "\n",
                                              quiet = TRUE))) %>%
  unnest(cols = c(file_contents)) %>%
  separate(file_contents, c("file_contents", "variance"), sep = "variation. ") %>%
  mutate(variance = as.numeric(variance)) %>%
  separate(filename, c("summary", "species", "run", "migrations")) %>%
  mutate(migrations = gsub("m", "", migrations)) %>%
  mutate(run = as.numeric(run),
         migrations = as.numeric(migrations)) %>%
  dplyr::select(run, migrations, variance) %>%
  mutate(variance = variance * 100) %>%
  group_by(migrations) %>%
  summarise(mean = mean(variance),
            min = min(variance),
            max = max(variance))

head(var_alf)

plot_var_alf <- ggplot(var_alf, aes(as.factor(migrations), mean)) +
  geom_point() +
  geom_line() +
  geom_linerange(aes(ymin = min, ymax = max)) +
  theme_emily() +
  geom_hline(yintercept = 99.8, linetype = "dashed") +
  xlab("Migration edges") + ylab("Variance explained") +
  ggtitle("A")

plot_var_alf

# Residuals

stem_alf <- "data/treemix/alfredi/initial_runs/alfredi_1_m0"
pop_order_alf <- "data/meta/treemix_pop_order_alf.txt"

# modified the below function to include different colour scheme and label abbreviations

png("figs/residuals_alfredi.png", units = "in", res = 300, width = 5, height = 5)

plot_resid(stem_alf, pop_order_alf, min = -0.009, max = 0.009, 
           cex = 1, usemax = T, wcols = "r", sp ="alfredi")
title("A", cex.main = 1.5, font.main = 1, adj=0, line = 0)

dev.off()

# The best number of migration events was estimated to M=1, using an ad hoc analysis of the
# second order rate of change in the log-likelihood (Evanno) method (fig. S4). 99.8% of the
# variance was explained at M=1, and variance only increases marginally at M>1.

# Birostris

# scp ehumble@eddie.ecdf.ed.ac.uk:/exports/cmvm/eddie/eb/groups/ogden_grp/emily/manta_pop_gen_2020/data/out/stacks_PE/treemix/initial_runs/birostris* data/treemix/birostris/initial_runs/
  
birostris <- "data/treemix/birostris/initial_runs"
birostris.optM <- optM(birostris)

plot_optM(birostris.optM, method = "Evanno")

birostris_test.linear <- optM(birostris, method = "linear")
plot_optM(birostris_test.linear, method = "linear")

plot_tree(cex = 0.9, "data/treemix/birostris/initial_runs/birostris_1_m3")

# 3 migration edges

# Variance explained only. Determined using treemix_fraction.pl script

#scp ehumble@eddie.ecdf.ed.ac.uk:/exports/cmvm/eddie/eb/groups/ogden_grp/emily/manta_pop_gen_2020/data/out/stacks_PE/treemix/initial_runs/var_birostris* data/treemix/birostris/initial_runs/

#~~ Birostris

data_path <- "data/treemix/birostris/initial_runs/"
files <- dir(data_path, pattern = "var*")

options(scipen=999)

var_bir <- tibble(filename = files) %>%
  mutate(file_contents = map(filename, ~ scan(file.path(data_path, .), what = "",
                                              sep = "\n",
                                              quiet = TRUE))) %>%
  unnest(cols = c(file_contents)) %>%
  separate(file_contents, c("file_contents", "variance"), sep = "variation. ") %>%
  mutate(variance = as.numeric(variance)) %>%
  separate(filename, c("summary", "species", "run", "migrations")) %>%
  mutate(migrations = gsub("m", "", migrations)) %>%
  mutate(run = as.numeric(run),
         migrations = as.numeric(migrations)) %>%
  dplyr::select(run, migrations, variance) %>%
  mutate(variance = variance * 100) %>%
  group_by(migrations) %>%
  summarise(mean = mean(variance),
            min = min(variance),
            max = max(variance))

head(var_bir)

plot_var_bir <- ggplot(var_bir, aes(as.factor(migrations), mean)) +
  geom_point() +
  geom_line() +
  geom_linerange(aes(ymin = min, ymax = max)) +
  theme_emily() +
  geom_hline(yintercept = 99.8, linetype = "dashed") +
  xlab("Migration edges") + ylab("Variance explained") +
  ggtitle("B") +
  ylim(99.2, 100)

plot_var_bir

plot_var_alf + plot_var_bir

ggsave("figs/var_treemix.png", plot_var_alf + plot_var_bir, 
       width = 8, height = 4)

# Residuals

stem_bir <- "data/treemix/birostris/initial_runs/birostris_1_m0"
pop_order_bir <- "data/meta/treemix_pop_order_bir.txt"

# modified the below function to include different colour scheme and label abbreviations

png("figs/residuals_birostris.png", units = "in", res = 300, width = 5, height = 5)

plot_resid(stem_bir, pop_order_bir, min = -0.009, max = 0.009, 
           cex = 1, usemax = T, wcols = "r", sp ="birostris")
title("B", cex.main = 1.5, font.main = 1, adj=0, line = 0)

dev.off()

# Residuals above zero represent populations that are more
# closely related to each other in the data than in the best-fit tree, and are candidates for
# admixture.



#~~ Plot final run

# Cant install genabel dependency on new versions of R. Simply loading bite functions instead.

# scp ehumble@eddie.ecdf.ed.ac.uk:/exports/cmvm/eddie/eb/groups/ogden_grp/emily/manta_pop_gen_2020/data/out/stacks_PE/treemix/final_runs/alfredi_m* data/treemix/alfredi/final_runs/
# scp ehumble@eddie.ecdf.ed.ac.uk:/exports/cmvm/eddie/eb/groups/ogden_grp/emily/manta_pop_gen_2020/data/out/stacks_PE/treemix/final_runs/alfredi_m*_outtree.newick data/treemix/alfredi/final_runs/

title_alf <- "A"
sub_alf <- expression(paste("Reef manta ray (", italic("Mobula alfredi"), ")"))
  
pdf("figs/treemix_alfredi_m2.pdf", width = 6, height = 5)

treemix.bootstrap(in.file = "data/treemix/alfredi/final_runs/alfredi_m2_treemix", 
                  out.file = "data/treemix/alfredi/final_runs/alfredi_m2",
                  plus =  0.05, ybar = 0.3, lwd = 1.5,
                  phylip.file = "data/treemix/alfredi/final_runs/alfredi_m2_outtree.newick", 
                  nboot = 100, mbar = T, plotboot = T,
                  cex=0.9, xmin = -0.005, disp = 0.001, 
                  boot.legend.location = "topleft", xbar = 0)
mtext(side = 3, line = 3, col = "#333333", at = -0.02, adj = 0, padj = 1.2, cex = 1.2, title_alf)
mtext(side = 3, line = 2, col = "#333333", at = -0.02, adj = 0, padj = 1.3, cex = 0.9, sub_alf)

dev.off()

plot_resid("data/treemix/alfredi/final_runs/alfredi_m1_treemix", "data/meta/treemix_pop_order_alf.txt")


#scp ehumble@eddie.ecdf.ed.ac.uk:/exports/cmvm/eddie/eb/groups/ogden_grp/emily/manta_pop_gen_2020/data/out/stacks_PE/treemix/final_runs/birostris_m* data/treemix/birostris/final_runs/
#scp ehumble@eddie.ecdf.ed.ac.uk:/exports/cmvm/eddie/eb/groups/ogden_grp/emily/manta_pop_gen_2020/data/out/stacks_PE/treemix/final_runs/birostris_m*_outtree.newick data/treemix/birostris/final_runs/

#scp ehumble@eddie.ecdf.ed.ac.uk:/exports/cmvm/eddie/eb/groups/ogden_grp/emily/manta_pop_gen_2020/data/out/stacks_PE/treemix/final_runs/birostris_m0* data/treemix/birostris/final_runs/
#scp ehumble@eddie.ecdf.ed.ac.uk:/exports/cmvm/eddie/eb/groups/ogden_grp/emily/manta_pop_gen_2020/data/out/stacks_PE/treemix/final_runs/birostris_m0*_outtree.newick data/treemix/birostris/final_runs/
  
title_bir <- "B"
sub_bir <- expression(paste("Oceanic manta ray (", italic("Mobula birostris"), ")"))

pdf("figs/treemix_birostris_m0.pdf", width = 6, height = 4.5)

treemix.bootstrap(in.file = "data/treemix/birostris/final_runs/birostris_m0_treemix", 
                  out.file = "data/treemix/birostris/final_runs/birostris_m0",
                  plus =  0.05, lwd = 1.5,
                  phylip.file = "data/treemix/birostris/final_runs/birostris_m0_outtree.newick", 
                  nboot = 100, mbar = F, plotboot = T,
                  cex=0.9, xmin = -0.005, disp = 0.001, 
                  boot.legend.location = "topleft", xbar = 0)
mtext(side = 3, line = 3, col = "#333333", at = -0.005, adj = 0, padj = 1.2, cex = 1.2, title_bir)
mtext(side = 3, line = 2, col = "#333333", at = -0.005, adj = 0, padj = 1.3, cex = 0.9, sub_bir)

dev.off() 

plot_resid("data/treemix/birostris/final_runs/birostris_m0_treemix", "data/meta/treemix_pop_order_bir.txt")

