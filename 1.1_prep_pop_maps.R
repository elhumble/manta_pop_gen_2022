# Script to prepare popmap files for denovomap optimisation
# And for full run of denovo map

# read process_radtag log files and determine depth of coverage per individual
# 280 individuals including outgroup

library(purrr)
library(dplyr)
library(ggplot2)
library(tidyr)
source("scripts/theme_emily.R")

# Convert processradtags output into csv files

lib_1 <- readLines("data/process_radtags_logs/process_radtags.lib_1.log") %>%
  .[c(14:109)] %>%
  write.table("data/process_radtags_logs/lib_1.csv", sep = "\t", 
              row.names = F, col.names = F,
              quote = F)

lib_2 <- readLines("data/process_radtags_logs/process_radtags.lib_2.log") %>%
  .[c(14:104)] %>%
  write.table("data/process_radtags_logs/lib_2.csv", sep = "\t", 
              row.names = F, col.names = F,
              quote = F)

lib_3 <- readLines("data/process_radtags_logs/process_radtags.lib_3.log") %>%
  .[c(14:109)] %>%
  write.table("data/process_radtags_logs/lib_3.csv", sep = "\t", 
              row.names = F, col.names = F,
              quote = F)

# Read in csv files

csv_files <- list.files(path = "data/process_radtags_logs/", pattern="*.csv")
csv_files <- paste0("data/process_radtags_logs/", csv_files)

# process rad tags out

prt <- csv_files %>% 
  map_dfr(read.csv, sep = "\t", .id = "source") %>%
  mutate(percent.retained = Retained/Total)

summary(prt$Retained)

# Combine with metadata

master <- read.csv("data/meta/sample_info_mastersheet_reformatted.csv")


df <- prt %>%
  separate(Filename, c("Filename", "Dup"), sep = "_") %>%
  left_join(master, by = c("Filename" = "Sample_code")) %>%
  group_by(Filename) %>%
  mutate(dup = n())


# Get some numbers

df %>% filter(Species == "Mobula alfredi") %>%
  distinct(Filename, .keep_all = T) %>%
  group_by(Location) %>%
  summarise(n = n())

df %>% filter(Species == "Mobula birostris") %>%
  distinct(Filename, .keep_all = T) %>%
  filter(Filename != 1055) %>%
  group_by(Location) %>%
  summarise(n = n())

df %>% filter(Species == "Mobula alfredi") %>%
  distinct(Filename, .keep_all = T) %>%
  group_by(Species) %>%
  summarise(n = n())

df %>% filter(Species == "Mobula birostris") %>%
  distinct(Filename, .keep_all = T) %>%
  filter(Filename != 1055) %>%
  group_by(Species) %>%
  summarise(n = n())

# Plot of percent reads retained

ggplot(prt, aes(percent.retained)) +
  geom_histogram() +
  facet_wrap(~source)

ggplot(prt, aes(percent.retained)) +
  geom_histogram()

# Filter individuals with low read numbers 

ggplot(prt, aes(Retained)) +
  geom_histogram()

prt <- prt %>%
  arrange(Retained) %>%
  mutate(rank = 1:nrow(.))

ggplot(prt, aes(Retained, rank)) +
  geom_point()

# Useful plot (percent retained but ranked by number of reads)

ggplot(prt, aes(percent.retained, rank)) +
  geom_point() +
  geom_vline(xintercept = 0.20)


prt_filter <- ggplot(prt, aes(percent.retained, rank)) +
  geom_point(size = 2, alpha = 0.3, shape = 21, fill = "black") +
  geom_vline(xintercept = 0.20, col = "red", linetype = 2) +
  theme_emily() +
  xlab("Percentage of reads retained") + 
  ylab(paste0("Samples ranked by", "\n", "percentage of reads retained"))

prt_filter

ggsave("figs/prt_filter.png", prt_filter, height = 4, width = 5)

nrow(filter(prt, percent.retained > 0.2))
nrow(filter(prt, percent.retained > 0.5))
nrow(filter(prt, percent.retained > 0.75))
nrow(filter(prt, percent.retained > 0.9))
nrow(filter(prt, percent.retained > 0))

# Filter at <0.20 percent reads retained

#low_qual <- filter(prt, percent.retained < 0.20) # Using this threshold
#low_qual_a <- filter(prt, Retained < 799753)

# using working set df from 0_prep to compare filtered inds to janes final list
#low_qual_a$Filename %in% working_set$Sample_code
#low_qual$Filename %in% low_qual_a$Filename


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#        Create population map files       #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Separate duplicate IDs into column

prt <- prt %>% 
  filter(percent.retained >= 0.20) %>% # 268 individuals
  separate(Filename, c("Filename", "Dup"), sep = "_")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#      Pop map for parameter optimisation       #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Read in master list

master <- read.csv("data/meta/sample_info_mastersheet_reformatted.csv")

# All individuals, no subsetting

# Combine with metadata

df <- prt %>%
 left_join(master, by = c("Filename" = "Sample_code")) %>%
 group_by(Filename) %>%
 mutate(dup = n())

# Remove possible duplicate sample 1055

df <- df %>%
  filter(Filename != 1055)

#~~~~~~~~~~~~~~~~~~~#
#   Across species  #
#~~~~~~~~~~~~~~~~~~~#

pop_map <- df %>%
  select(Filename, Dup, Species, Location, source, percent.retained)

pop_map %>%
  group_by(Species, Location) %>%
  summarise(n = n())

pop_map %>%
  group_by(Species) %>%
  summarise(n = n())

# Remove other species and fix filenames
# Write popmap for all individuals

full_dataset <- pop_map %>%
  mutate(Species = as.character(Species)) %>%
  dplyr::filter(Species != "Aetobatus narinari") %>%
  filter(Species != "Rhinoptera bonasus") %>%
  # Recode sample name
  unite("Filename", Filename, Dup) %>%
  mutate(Filename = gsub("_NA", "", Filename)) %>%
  unite("Pop", Species:Location) %>%
  mutate(Pop = gsub(" ", "_", Pop)) %>%
  mutate(opt_pop = "Manta")

write.table(full_dataset[c(1,5)], "data/file_lists/pop_map_across_sp_param_optimisation_full_dataset.txt",
            col.names = F, row.names = F, quote = F, sep = "\t")

# scp data/pop_map_across_sp_param_optimisation_full_dataset.txt ehumble@eddie.ecdf.ed.ac.uk:/exports/cmvm/eddie/eb/groups/ogden_grp/emily/manta_pop_gen_2020/file_lists/pop_map_across_sp_param_optimisation_full_dataset.txt


#~~ Full dataset with species information for running ustacks on all individuals (minus those with v low reads and dup)

full_dataset_sp_info <- pop_map %>%
  mutate(Species = as.character(Species)) %>%
  dplyr::filter(Species != "Aetobatus narinari") %>%
  filter(Species != "Rhinoptera bonasus") %>%
  # Recode sample name
  unite("Filename", Filename, Dup) %>%
  mutate(Filename = gsub("_NA", "", Filename)) %>%
  unite("Pop", Species:Location) %>%
  mutate(Pop = gsub(" ", "_", Pop))

write.table(full_dataset_sp_info[c(1,2)], "data/file_lists/pop_map_across_sp_full_dataset.txt",
            col.names = F, row.names = F, quote = F, sep = "\t")

# scp data/pop_map_across_sp_full_dataset.txt ehumble@eddie.ecdf.ed.ac.uk:/exports/cmvm/eddie/eb/groups/ogden_grp/emily/manta_pop_gen_2020/file_lists/pop_map_across_sp_full_dataset.txt

#~~ Full filtered alfredi

full_dataset_alfredi <- full_dataset_sp_info %>%
  filter(grepl("alfredi", Pop))


write.table(full_dataset_alfredi[c(1,2)], "data/file_lists/pop_map_alfredi_full_dataset.txt",
            col.names = F, row.names = F, quote = F, sep = "\t")

scp data/pop_map_alfredi_full_dataset.txt ehumble@eddie.ecdf.ed.ac.uk:/exports/cmvm/eddie/eb/groups/ogden_grp/emily/manta_pop_gen_2020/file_lists/pop_map_alfredi_full_dataset.txt

#~~ Full filtered birostris


full_dataset_birostris <- full_dataset_sp_info %>%
  filter(grepl("birostris", Pop))


write.table(full_dataset_birostris[c(1,2)], "data/file_lists/pop_map_birostris_full_dataset.txt",
            col.names = F, row.names = F, quote = F, sep = "\t")

scp data/pop_map_birostris_full_dataset.txt ehumble@eddie.ecdf.ed.ac.uk:/exports/cmvm/eddie/eb/groups/ogden_grp/emily/manta_pop_gen_2020/file_lists/pop_map_birostris_full_dataset.txt


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Subset of individuals with highest reads #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Across species

pop_map <- df %>%
  select(Filename, Dup, Species, Location, source, percent.retained)

across_sp_opt <- pop_map %>%
  # Remove other species
  filter(Species != "Rhinoptera bonasus") %>%
  filter(Species != "Aetobatus narinari") %>%
  # Group by ID
  group_by(Filename) %>%
  # Select duplicate individual with highest number of reads retained
  filter(percent.retained == max(percent.retained)) %>%
  # Select two individuals from each species and location
  group_by(Species, Location) %>%
  arrange(desc(percent.retained)) %>%
  slice(1:2) %>%
  # Recode sample name
  unite("Filename", Filename, Dup) %>%
  mutate(Filename = gsub("_NA", "", Filename)) %>%
  unite("Pop", Species:Location) %>%
  mutate(Pop = gsub(" ", "_", Pop)) %>%
  mutate(opt_pop = "Manta") # optimisation should assume one pop

write.table(across_sp_opt[c(1,5)], "data/file_lists/pop_map_across_sp_param_optimisation.txt",
            col.names = F, row.names = F, quote = F, sep = "\t")

#~~~~~~~~~~~~~~~~~~~#
#   Within species  #
#~~~~~~~~~~~~~~~~~~~#

#~~ Alfredi

alfredi_opt <- pop_map %>%
  filter(Species == "Mobula alfredi") %>%
  group_by(Filename,) %>%
  filter(percent.retained == max(percent.retained)) %>%
  # Select two individuals from each species and location
  group_by(Species, Location) %>%
  arrange(desc(percent.retained)) %>%
  slice(1:3) %>% # taking top 3 animals here
  # Recode sample name
  unite("Filename", Filename, Dup) %>%
  mutate(Filename = gsub("_NA", "", Filename)) %>%
  unite("Pop", Species:Location) %>%
  mutate(Pop = gsub(" ", "_", Pop)) %>%
  mutate(opt_pop = "Manta") # optimisation should assume one pop

write.table(alfredi_opt[c(1,5)], "data/file_lists/pop_map_alfredi_param_optimisation.txt",
            col.names = F, row.names = F, quote = F, sep = "\t")

#~~ Birostris

birostris_opt <- pop_map %>%
  filter(Species == "Mobula birostris") %>%
  group_by(Filename,) %>%
  filter(percent.retained == max(percent.retained)) %>%
  # Select two individuals from each species and location
  group_by(Species, Location) %>%
  arrange(desc(percent.retained)) %>%
  slice(1:3) %>%
  # Recode sample name
  unite("Filename", Filename, Dup) %>%
  mutate(Filename = gsub("_NA", "", Filename)) %>%
  unite("Pop", Species:Location) %>%
  mutate(Pop = gsub(" ", "_", Pop)) %>%
  mutate(opt_pop = "Manta") # optimisation should assume one pop

write.table(birostris_opt[c(1,5)], "data/file_lists/pop_map_birostris_param_optimisation.txt",
            col.names = F, row.names = F, quote = F, sep = "\t")

#~~ Upload to EDDIE

scp data/pop_map_across_sp_param_optimisation.txt ehumble@eddie.ecdf.ed.ac.uk:/exports/cmvm/eddie/eb/groups/ogden_grp/emily/manta_pop_gen_2020/file_lists/pop_map_across_sp_param_optimisation.txt
scp data/pop_map_alfredi_param_optimisation.txt ehumble@eddie.ecdf.ed.ac.uk:/exports/cmvm/eddie/eb/groups/ogden_grp/emily/manta_pop_gen_2020/file_lists/pop_map_alfredi_param_optimisation.txt
scp data/pop_map_birostris_param_optimisation.txt ehumble@eddie.ecdf.ed.ac.uk:/exports/cmvm/eddie/eb/groups/ogden_grp/emily/manta_pop_gen_2020/file_lists/pop_map_birostris_param_optimisation.txt


