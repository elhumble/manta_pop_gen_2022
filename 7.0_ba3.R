# Script to visualise in BA3-SNPs output

library(stringr)
library(dplyr)
library(tidyr)
library(data.table)
library(wesanderson)
source("scripts/theme_emily.R")
library(circlize)
library(ggplot2)
library(scales)
library(purrr)

# scp ehumble@eddie.ecdf.ed.ac.uk:/exports/cmvm/eddie/eb/groups/ogden_grp/emily/manta_pop_gen_2020/data/out/stacks_PE/bayesass/alf_seed_1*/alfredi_bayesass3_* data/ba3/alf_mig_rates/
# scp ehumble@eddie.ecdf.ed.ac.uk:/exports/cmvm/eddie/eb/groups/ogden_grp/emily/manta_pop_gen_2020/data/out/stacks_PE/bayesass/*bayesass3.trace.txt data/ba3/

#~~~~~~~~~~~~~~~~~~~~~~~#
#       ALFREDI         #
#~~~~~~~~~~~~~~~~~~~~~~~#

#~ Assess mixing

alf_trace <- fread("data/ba3/alfredi_bayesass3.trace.txt")

trace_fig <- ggplot(filter(alf_trace, LogProb > -66000), aes(State, LogProb)) +
  geom_line(linewidth = 0.4) +
  theme_emily() + 
  scale_x_continuous(labels = comma) +
  scale_y_continuous(labels = comma)

ggsave("figs/trace_ba3.png", trace_fig, height = 4, width = 5)


# Function to extract migration rates from BA3 output file

get_ba3_out <- function(file) {
  
  X <- scan(
    file = file,
    what = "",
    sep = "\n",
    quiet = TRUE
  )
  
  # extract data for each iteration (30)
  
  START <- grep("Population Index", X) # line numbers
  END <- grep("Inbreeding", X) # line numbers
  
  X <- X[START:END]
  X
  
  # get pop codes into df
  
  pops <- X[2] %>%
    str_split(" ") %>%
    data.frame()
  
  colnames(pops) <- "pops"
  
  pops <- pops %>%
    separate(pops, c("code", "pop")) %>%
    filter(!is.na(pop))
  
  # get migration rates into df
  
  mig <- X[5:length(X) - 1] %>%
    data.frame()
  
  colnames(mig) <- "mig"
  
  mig <- mig %>%
    separate(mig, c("a", "b", "c", "d", "e", "f", "g"), sep = " m")
  
  # m[i][j] is the fraction of individuals in population i that are migrants derived from pop j (per generation)
  # pop1 == i & pop2 == j
  
  df <-
    data.frame(col = c(mig$b, mig$c, mig$d, mig$e, mig$f, mig$g)) %>%
    separate(col, c("pops", "mig_rates"), ": ") %>%
    separate(pops, c("pop1", "pop2"), "\\]\\[") %>%
    mutate(pop1 = gsub("\\[", "", pop1)) %>%
    mutate(pop2 = gsub("\\]", "", pop2)) %>%
    mutate(mig_rates = gsub("\\)", "", mig_rates)) %>%
    separate(mig_rates, c("mig_rate", "sd"), "\\(") %>%
    left_join(pops, by = c("pop1" = "code")) %>%
    left_join(pops, by = c("pop2" = "code")) %>%
    dplyr::select(!c(pop1, pop2)) %>%
    mutate(mig_rate = as.numeric(mig_rate)) %>%
    mutate(sd = as.numeric(sd))
  
  return(df)
}


# Extract migration rates

files <- list.files("data/ba3/alf_mig_rates/")
files <- paste0("data/ba3/alf_mig_rates/", files)

migs <- map(files, get_ba3_out) %>%
  bind_rows(.id = "column_label") %>%
  group_by(pop.x, pop.y) %>%
  summarise(mean_mig_rate = mean(mig_rate),
            mean_sd = mean(sd)) %>%
  mutate(CI_upper = mean_mig_rate + 1.96 * mean_sd,
         CI_lower = mean_mig_rate - 1.96 * mean_sd) %>%
  mutate(pop.x = case_when(pop.x == "Seychelles" ~"SEY",
                           pop.x == "Maldives" ~ "MAL",
                           pop.x == "Chagos" ~ "CHAG",
                           pop.x == "Australia" ~ "AP",
                           pop.x == "Hawaii" ~ "HAW",
                           pop.x == "Fiji" ~ "FIJI")) %>%
  mutate(pop.y = case_when(pop.y == "Seychelles" ~"SEY",
                           pop.y == "Maldives" ~ "MAL",
                           pop.y == "Chagos" ~ "CHAG",
                           pop.y == "Australia" ~ "AP",
                           pop.y == "Hawaii" ~ "HAW",
                           pop.y == "Fiji" ~ "FIJI")) %>%
  mutate(scale = case_when(pop.x == "MAL" & pop.y == "AUS" |
                             pop.x == "MAL" & pop.y == "FIJI" |
                             pop.x == "MAL" & pop.y == "HAW" |
                             pop.x == "SEY" & pop.y == "AUS"|
                             pop.x == "SEY" & pop.y == "FIJI"|
                             pop.x == "SEY" & pop.y == "HAW"|
                             pop.x == "CHAG" & pop.y == "AUS"|
                             pop.x == "CHAG" & pop.y == "FIJI"|
                             pop.x == "CHAG" & pop.y == "HAW"|
                             pop.x == "AUS" & pop.y == "MAL"|
                             pop.x == "AUS" & pop.y == "SEY"|
                             pop.x == "AUS" & pop.y == "CHAG"|
                             pop.x == "FIJI" & pop.y == "MAL"|
                             pop.x == "FIJI" & pop.y == "SEY"|
                             pop.x == "FIJI" & pop.y == "CHAG"|
                             pop.x == "HAW" & pop.y == "MAL"|
                             pop.x == "HAW" & pop.y == "SEY"|
                             pop.x == "HAW" & pop.y == "CHAG" ~
                             "Between", TRUE ~ "Within"))

mig_df_long <- migs %>%
  filter(pop.x != pop.y)

# Some numbers

summary(mig_df_long$mean_mig_rate)

mig_df_long %>%
  ungroup() %>%
  filter(scale == "Between") %>%
  summarise(mean = mean(mean_mig_rate),
            min = min(mean_mig_rate),
            max = max(mean_mig_rate))

mig_df_long %>%
  ungroup() %>%
  mutate(mig_rate = mean_mig_rate * 100) %>%
  filter(scale == "Between") %>%
  summarise(mean = mean(mig_rate),
            min = min(mig_rate),
            max = max(mig_rate))

mig_df_long %>%
  ungroup() %>%
  mutate(mig_rate = mean_mig_rate * 100) %>%
  filter(scale == "Within") %>%
  summarise(mean = mean(mig_rate),
            min = min(mig_rate),
            max = max(mig_rate))

# into
mig_df_long %>%
  mutate(mig_rate = mean_mig_rate * 100) %>%
  group_by(pop.x) %>%
  summarise(sum = sum(mig_rate),
            min = min(mig_rate),
            max = max(mig_rate)) %>%
  arrange(sum)

# average proportion of migrant individuals in a population

mig_df_long %>%
  mutate(mig_rate = mean_mig_rate * 100) %>%
  group_by(pop.x) %>%
  summarise(sum = sum(mig_rate),
            min = min(mig_rate),
            max = max(mig_rate)) %>%
  ungroup() %>%
  summarise(mean = mean(sum),
            min = min(sum),
            max = max(sum))

mig_df_long %>%
  mutate(mig_rate = mean_mig_rate * 100) %>%
  filter(scale == "Between") %>%
  group_by(pop.x) %>%
  summarise(sum = sum(mig_rate),
            min = min(mig_rate),
            max = max(mig_rate)) %>%
  ungroup() %>%
  summarise(mean = mean(sum),
            min = min(sum),
            max = max(sum))

# out of

mig_df_long %>%
  mutate(mig_rate = mean_mig_rate * 100) %>%
  group_by(pop.y) %>%
  summarise(sum = sum(mig_rate),
            min = min(mig_rate),
            max = max(mig_rate)) %>%
  arrange(sum)

#~~ Create plot

plot_df <- migs %>%
  select(mean_mig_rate, mean_sd, pop.x, pop.y) %>%
  pivot_wider(names_from = pop.x, values_from = c(mean_mig_rate, mean_sd))

# Convert to matrix
mig_rate_mat <- as.matrix(plot_df[, 2:7])
mig_rate_mat

sd_mat <- as.matrix(plot_df[, 8:13])
sd_mat

dimnames(mig_rate_mat) <-
  list(source = plot_df$pop.y, sink = plot_df$pop.y)

mig_rate_mat

# Create circular plot

cols <- c("#9986A5",
          "#9A8822",
          "#35274A",
          "#3B9AB2",
          "#0B775E",
          "#E2D200")

# Reset the circular layout parameters
circos.clear()

# Parameters for the circular layout
circos.par(start.degree = 90, gap.degree = 6)

# Plot chord diagram
chordDiagram(
  x = mig_rate_mat,
  grid.col = cols,
  grid.border = "black",
  transparency = 0.25,
  order = plot_df$pop.y ,
  directional = T,
  direction.type = "arrows",
  self.link = 1,
  preAllocateTracks = list(track.height = 0.1),
  annotationTrack = "grid",
  annotationTrackHeight = c(0.1, 0.1),
  link.border = "NA",
  link.sort = T,
  link.decreasing = T,
  link.arr.length = 0.15,
  link.arr.lty = 3,
  link.arr.col = "#252525",
  link.largest.ontop = F
)

# Add labels to chord diagram
circos.trackPlotRegion(
  track.index = 1,
  bg.border = NA,
  panel.fun = function(x, y) {
    xlim = get.cell.meta.data("xlim")
    sector.index = get.cell.meta.data("sector.index")
    # Text direction
    theta = circlize(mean(xlim), 1)[1, 1] %% 360
    dd = ifelse(theta < 180 ||
                  theta > 360, "bending.inside", "bending.outside")
    circos.text(
      x = mean(xlim),
      y = 0.8,
      labels = sector.index,
      facing = dd,
      niceFacing = TRUE,
      cex = 2,
      font = 1
    )
  }
)

# Add title
mtext("E", outer = FALSE, cex = 2.5, font = 1, 
      line = -1, adj = 0.1, col = "#333333")

dev.copy2pdf(file = "figs/ba3_alfredi.pdf", height = 10, width = 10)

circ_alf <- recordPlot()
circ_alf
saveRDS(circ_alf, "figs/circos_alf.RDS")



