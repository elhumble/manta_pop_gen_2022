# Determine optimum set of denovomap params for across species analysis

library(tidyverse)
library(data.table)
library(stringr)
library(scales)

# stacks_file <- "data/denovomap_optimisation/m3_M1_n1/populations.log.distribs"

# Function to read in population log files 

get_table <- function(stacks_file) {
  file_temp <- fread(stacks_file, fill = T) %>% as_tibble()
  index1 <- which(file_temp$V1 == "BEGIN" & file_temp$V2 == "snps_per_loc_postfilters")  
  index2 <- which(file_temp$V1 == "END" & file_temp$V2 == "snps_per_loc_postfilters")  
  df <- file_temp[(index1+3):(index2-1), ] %>% 
    select(V1) %>%
    separate(V1, c("n_snps", "n_loci"))
}

path <- "data/"

files <- list.files(path = "data/", 
           pattern="populations.log.distribs",
           recursive = T)

data <- tibble(filename = files) %>% # create a data frame holding the file names
  mutate(file_contents = map(filename, ~ get_table(file.path(path, .)))) %>% # a new data column
  unnest(cols = c(file_contents)) %>% # unnest into a df
  mutate(n_snps = as.numeric(n_snps),
         n_loci = as.numeric(n_loci)) %>%
  separate(filename, c("job", "run", "params", "populations"), sep = "/") %>%
  separate(params, c("m", "M", "n"), remove = F) %>%
  mutate(populations = case_when(grepl("r80", populations) ~ "r80",
                                 TRUE ~ "full")) %>%
  mutate(m = as.numeric(gsub("m", "", m)),
         M = as.numeric(gsub("M", "", M)),
         n = as.numeric(gsub("n", "", n)))

# summarise data

df <- data %>%
  group_by(run,params, populations, m, M, n) %>%
  summarise(assembled_loci = sum(n_loci),
            polymorphic_loci = sum(n_loci[n_snps != 0]),
            n_snps = sum(n_loci * n_snps)) %>%
  pivot_longer(cols = assembled_loci:n_snps,
               names_to = "statistic",
               values_to = "total")
#%>%
 # pivot_longer(M:n)


# Both r80 and full

ggplot(df, aes(params, total, col = populations)) +
  geom_point() +
  facet_wrap(run ~ statistic)

# r80 only
ggplot(filter(df, populations == "r80"), aes(params, total, col = populations)) +
  geom_point() +
  facet_wrap(run ~ statistic, scales = "free") +
  theme(axis.text.x = element_text(angle=90))

df %>%
  filter(populations == "r80") %>%
  group_by(run, statistic) %>%
  slice(which.max(total))



pal <- c("#35274A",
         "#3B9AB2")

param_plot <- df %>%
  filter(run == "across_sp_opt") %>%
  filter(statistic == "polymorphic_loci") %>%
  filter(populations == "r80") %>%
  mutate(test = case_when(n == M ~ "n = M",
                       n > M ~ "n = M + 1")) %>%
  ggplot(aes(as.factor(M), total, col = test)) +
  theme_emily() +
  geom_point() +
  scale_color_manual(values = pal) +
  theme(legend.position = c(0.8, 0.1),
        legend.title = element_blank()) +
  xlab("M") + ylab(paste0("Number of polymorphic loci", "\n", "present in 80% of individuals")) +
  scale_y_continuous(labels = comma)

param_plot

ggsave("figs/param_plot.png", param_plot, height = 4, width = 5)

#~~ Across species

# m3 M3 n4

#~~ Alfredi

# m3 M3 n4

#~~ Birostris

# m3 M3 n4

