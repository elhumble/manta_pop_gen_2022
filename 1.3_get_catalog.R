# Checking ustacks output
# Create catalog of good quality samples

library(tidyverse)
library(data.table)
library(stringr)

# ustacks_log <- "data/ustacks_logs/e_files/ustacks.e10003842.1"

# Function to read in ustacks e_files / log files 

get_table <- function(ustacks_log) {
  file_temp <- fread(ustacks_log, fill = T) %>% as_tibble()
  index1 <- which(file_temp$V1 == "Final")
  index2 <- which(file_temp$V1 == "Input")
  df <- file_temp[c(index1, index2), ] %>% 
    select(V1:V6) %>%
    mutate(ID = .[2,3]) %>%
    filter(V1 != "Input") %>%
    pivot_longer(V3:V6) %>%
    separate(value, c("statistic", "coverage"), sep = "=") %>%
    mutate(coverage = gsub(";", "", coverage),
           V2 = gsub(":", "", V2)) %>%
    mutate(ID = gsub("'/exports/eddie/scratch/ehumble/reads/", "", ID$V3)) %>%
    mutate(ID = gsub(".fq.gz'", "", ID)) %>%
    mutate(coverage = gsub("\\(.+\\)", "", coverage)) %>%
    select(-c(V1, V2, name)) %>%
    mutate(coverage = as.numeric(coverage))
}


#~~ across species PE
path <- "data/ustacks_logs/"

files <- list.files(path = "data/ustacks_logs/", 
                    pattern="ustacks*",
                    recursive = T)

#~~ Read in data

data <- tibble(filename = files) %>% # create a data frame holding the file names
  mutate(file_contents = map(filename, ~ get_table(file.path(path, .)))) %>% # a new data column
  unnest(cols = c(file_contents)) %>% # unnest into a df
  separate(filename, c("run", "e_file"), sep = "\\/e_files\\/")

# mean locus coverage per individuals

mean_per_locus_cov <- data %>%
  filter(statistic == "mean")

hist(mean_per_locus_cov$coverage)
summary(mean_per_locus_cov$coverage)

#~~ Read in processradtags output and compare to ustacks logs

#csv_files <- list.files(path = "temp/", pattern="*.csv")
#csv_files <- paste0("temp/", csv_files)

csv_files <- list.files(path = "data/process_radtags_logs/", pattern="*.csv")
csv_files <- paste0("data/process_radtags_logs/", csv_files)

# process rad tags out
prt <- csv_files %>% 
  map_dfr(read.csv, sep = "\t", .id = "source") %>%
  mutate(percent.retained = Retained/Total)

mean_per_locus_cov <- mean_per_locus_cov %>%
  left_join(prt, by = c("ID" = "Filename"))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#        Select individuals for catalog        #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Combine with metadata

mean_per_locus_cov <- mean_per_locus_cov %>% 
  separate(ID, c("Filename", "Dup"), sep = "_")

# Read in master list

master <- read.csv("data/meta/sample_info_mastersheet_reformatted.csv")

# Combine

mean_per_locus_cov <- mean_per_locus_cov %>%
  left_join(master, by = c("Filename" = "Sample_code")) %>%
  group_by(Filename) %>%
  mutate(dup = n())

ggplot(mean_per_locus_cov, aes(coverage, Retained)) +
  geom_point(aes(col = Location, shape = Species)) +
  facet_wrap(run ~ Location)

ggplot(mean_per_locus_cov, aes(coverage, Retained)) +
  geom_point(aes(col = Location, shape = Species)) +
  facet_wrap(~ run)

# All groupings have individuals with coverage over 50


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#      Across species pop map     #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Select top 100 individuals

pop_map <- mean_per_locus_cov %>%
  select(run, Filename, Dup, Species, Location, coverage, Retained)

across_sp_catalog <- pop_map %>%
  filter(run == "across_sp") %>%
  # Group by ID
  group_by(Filename) %>%
  # Select duplicate individual with highest coverage
  filter(coverage == max(coverage)) %>%
  # Filter individuals with coverage < 25
  filter(coverage >= 25) %>%
  # Remove putative hybrid individual from catalog
  filter(Filename != "1328") %>%
  # Select up to 10 individuals from each species and location
  group_by(Species, Location) %>%
  arrange(desc(coverage)) %>%
  slice(1:10) %>%
  # Recode sample name
  unite("Filename", Filename, Dup) %>%
  mutate(Filename = gsub("_NA", "", Filename)) %>%
  unite("Pop", Species:Location) %>%
  mutate(Pop = gsub(" ", "_", Pop))

write.table(across_sp_catalog[c(2,3)], "data/file_lists/across_sp_catalog.txt",
            col.names = F, row.names = F, quote = F, sep = "\t")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#         Alfredi pop map         #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

alfredi_catalog <- pop_map %>%
  filter(run == "alfredi") %>%
  # Group by ID
  group_by(Filename) %>%
  # Select duplicate individual with highest coverage
  filter(coverage == max(coverage)) %>%
  # Filter individuals with coverage < 25
  filter(coverage >= 25) %>%
  # Remove putative hybrid individual from catalog
  filter(Filename != "1328") %>%
  # Select up to 10 individuals from each species and location
  group_by(Species, Location) %>%
  arrange(desc(coverage)) %>%
  slice(1:10) %>%
  # Recode sample name
  unite("Filename", Filename, Dup) %>%
  mutate(Filename = gsub("_NA", "", Filename)) %>%
  unite("Pop", Species:Location) %>%
  mutate(Pop = gsub(" ", "_", Pop))

write.table(alfredi_catalog[c(2,3)], "data/file_lists/alfredi_catalog.txt",
            col.names = F, row.names = F, quote = F, sep = "\t")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#        Birostris pop map         #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


birostris_catalog <- pop_map %>%
  filter(run == "birostris") %>%
  # Group by ID
  group_by(Filename) %>%
  # Select duplicate individual with highest coverage
  filter(coverage == max(coverage)) %>%
  # Filter individuals with coverage < 25
  filter(coverage >= 25) %>%
  # Remove putative hybrid individual from catalog
  filter(Filename != "1328") %>%
  # Select up to 10 individuals from each species and location
  group_by(Species, Location) %>%
  arrange(desc(coverage)) %>%
  slice(1:10) %>%
  # Recode sample name
  unite("Filename", Filename, Dup) %>%
  mutate(Filename = gsub("_NA", "", Filename)) %>%
  unite("Pop", Species:Location) %>%
  mutate(Pop = gsub(" ", "_", Pop))

write.table(birostris_catalog[c(2,3)], "data/file_lists/birostris_catalog.txt",
            col.names = F, row.names = F, quote = F, sep = "\t")


# Transfer to EDDIE
scp data/across_sp_catalog.txt ehumble@eddie.ecdf.ed.ac.uk:/\
exports/cmvm/eddie/eb/groups/ogden_grp/emily/manta_pop_gen_2020/file_lists/across_sp_catalog.txt

scp data/alfredi_catalog.txt ehumble@eddie.ecdf.ed.ac.uk:/\
exports/cmvm/eddie/eb/groups/ogden_grp/emily/manta_pop_gen_2020/file_lists/alfredi_catalog.txt

scp data/birostris_catalog.txt ehumble@eddie.ecdf.ed.ac.uk:/\
exports/cmvm/eddie/eb/groups/ogden_grp/emily/manta_pop_gen_2020/file_lists/birostris_catalog.txt


