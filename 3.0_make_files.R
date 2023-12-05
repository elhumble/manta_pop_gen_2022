# Create files for running software

library(purrr)
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
library(wesanderson)
library(patchwork)
library(readxl)

#~~~~~~~~~~~~~~~~~~~~~~~~~#
#        Structure        #
#~~~~~~~~~~~~~~~~~~~~~~~~~#

# Handmade structure file

# scp ehumble@eddie.ecdf.ed.ac.uk:/exports/cmvm/eddie/eb/groups/ogden_grp/emily/manta_pop_gen_2020/data/out/stacks_PE/birostris_geno4_dp6_ind_geno_depth_no_dups_maf_unrelated_ld.ped data/structure/
# scp ehumble@eddie.ecdf.ed.ac.uk:/exports/cmvm/eddie/eb/groups/ogden_grp/emily/manta_pop_gen_2020/data/out/stacks_PE/alfredi_geno4_dp6_ind_geno_depth_no_dups_maf_unrelated_ld.ped data/structure/
  

#~~ BIROSTRIS

# get genotypes

system("mkdir data/structure/temp")
system("cut -d ' ' -f 7- data/structure/birostris_geno4_dp6_ind_geno_depth_no_dups_maf_unrelated_ld.ped > data/structure/temp/A")

# convert values
system("sed 's/A/1/g' data/structure/temp/A > data/structure/temp/B")
system("sed 's/T/2/g' data/structure/temp/B > data/structure/temp/C")
system("sed 's/G/3/g' data/structure/temp/C > data/structure/temp/D")
system("sed 's/C/4/g' data/structure/temp/D > data/structure/temp/E")
system("sed 's/0/-9/g' data/structure/temp/E > data/structure/temp/F")

# get ids
system("awk '{print $1}' data/structure/birostris_geno4_dp6_ind_geno_depth_no_dups_maf_unrelated_ld.ped > data/structure/temp/ids")

# make structure input file
gen <- fread("data/structure/temp/F", header = F)
meta <- fread("data/meta/sample_info_mastersheet_reformatted.csv") %>%
  dplyr::select(c(Sample_code, Location))

ids <- fread("data/structure/temp/ids", header = F) %>%
  separate(V1, c("Sample", "dup")) %>%
  left_join(meta, by = c("Sample" = "Sample_code")) %>%
  unite("V1", Sample:dup) %>%
  mutate(V1 = gsub("_NA", "", V1)) %>%
  mutate(Location = gsub(" ", "_", Location)) %>%
  mutate(POP = case_when(Location == "Sri_Lanka" ~ 1,
                         Location == "Mexico_Pacific" ~ 2,
                         Location == "the_Philippines" ~ 3,
                         Location == "Mexico_Caribbean" ~ 4,
                         Location == "Peru" ~ 5,
                         Location == "South_Africa" ~ 6)) %>%
  mutate(popinfo = 1) %>%
  mutate(LocData = POP) %>%
  select(V1, POP, popinfo, LocData)

write.table(ids, "data/structure/temp/ids_pop.txt", col.names = F, row.names = F, quote = F)
system("paste -d ' ' data/structure/temp/ids_pop.txt data/structure/temp/F > data/structure/temp/birostris_handmade_rad.stru")

# sort structure files by population

system("sort -k 2 data/structure/temp/birostris_handmade_rad.stru > data/structure/birostris_handmade_rad_sort.stru")

system("rm -r data/structure/temp")

scp data/structure/*.stru ehumble@eddie.ecdf.ed.ac.uk:/exports/cmvm/eddie/eb/groups/ogden_grp/emily/manta_pop_gen_2020/data/out/stacks_PE/structure/

#~~ ALFREDI

# get genotypes

system("mkdir data/structure/temp")
system("cut -d ' ' -f 7- data/structure/alfredi_geno4_dp6_ind_geno_depth_no_dups_maf_unrelated_ld.ped > data/structure/temp/A")

# convert values
system("sed 's/A/1/g' data/structure/temp/A > data/structure/temp/B")
system("sed 's/T/2/g' data/structure/temp/B > data/structure/temp/C")
system("sed 's/G/3/g' data/structure/temp/C > data/structure/temp/D")
system("sed 's/C/4/g' data/structure/temp/D > data/structure/temp/E")
system("sed 's/0/-9/g' data/structure/temp/E > data/structure/temp/F")

# get ids
system("awk '{print $1}' data/structure/alfredi_geno4_dp6_ind_geno_depth_no_dups_maf_unrelated_ld.ped > data/structure/temp/ids")

# make structure input file
gen <- fread("data/structure/temp/F", header = F)
meta <- fread("data/meta/sample_info_mastersheet_reformatted.csv") %>%
  dplyr::select(c(Sample_code, Location))

ids <- fread("data/structure/temp/ids", header = F) %>%
  separate(V1, c("Sample", "dup")) %>%
  left_join(meta, by = c("Sample" = "Sample_code")) %>%
  unite("V1", Sample:dup) %>%
  mutate(V1 = gsub("_NA", "", V1)) %>%
  mutate(Location = gsub(" ", "_", Location)) %>%
  mutate(POP = case_when(Location == "Maldives" ~ 1,
                         Location == "Seychelles" ~ 2,
                         Location == "Chagos" ~ 3,
                         Location == "Fiji" ~ 4,
                         Location == "Australia_Pacific" ~ 5,
                         Location == "Hawaii" ~ 6)) %>%
  mutate(popinfo = 1) %>%
  mutate(LocData = POP) %>%
  select(V1, POP, popinfo, LocData)

write.table(ids, "data/structure/temp/ids_pop.txt", col.names = F, row.names = F, quote = F)
system("paste -d ' ' data/structure/temp/ids_pop.txt data/structure/temp/F > data/structure/temp/alfredi_handmade_rad.stru")

# sort structure files by population

system("sort -k 2 data/structure/temp/alfredi_handmade_rad.stru > data/structure/alfredi_handmade_rad_sort.stru")
system("rm -r data/structure/temp")


scp data/structure/*.stru ehumble@eddie.ecdf.ed.ac.uk:/exports/cmvm/eddie/eb/groups/ogden_grp/emily/manta_pop_gen_2020/data/out/stacks_PE/structure/
  
#~~~~~~~~~~~~~~~~~~~~~~~~~#
#        BayesAss3        #
#~~~~~~~~~~~~~~~~~~~~~~~~~#

#scp ehumble@eddie.ecdf.ed.ac.uk:/exports/cmvm/eddie/eb/groups/ogden_grp/emily/manta_pop_gen_2020/data/out/stacks_PE/bayesass/*bayesass3.ped data/ba3/
#scp ehumble@eddie.ecdf.ed.ac.uk:/exports/cmvm/eddie/eb/groups/ogden_grp/emily/manta_pop_gen_2020/data/out/stacks_PE/bayesass/*bayesass3.map data/ba3/

meta <- fread("data/meta/sample_info_mastersheet_reformatted.csv") %>%
  dplyr::select(c(Sample_code, Location))

# ALFREDI

alf_ped <- fread("data/ba3/alfredi_bayesass3.ped") %>%
  separate(V1, c("Sample", "dup")) %>%
  left_join(meta, by = c("Sample" = "Sample_code")) %>%
  unite("V1", Sample:dup) %>%
  mutate(V1 = gsub("_NA", "", V1)) %>%
  mutate(Location = gsub(" ", "_", Location))

alf_ped$V1 <- alf_ped$Location

alf_ped <- alf_ped %>%
  dplyr::select(-c(Location))

write.table(alf_ped, "data/ba3/alfredi_bayesass3.ped",
            quote = F, row.names = F, col.names = F)

# Convert to immanc using PGDspider

# BIROSTRIS

bir_ped <- fread("data/ba3/birostris_bayesass3.ped") %>%
  separate(V1, c("Sample", "dup")) %>%
  left_join(meta, by = c("Sample" = "Sample_code")) %>%
  unite("V1", Sample:dup) %>%
  mutate(V1 = gsub("_NA", "", V1)) %>%
  mutate(Location = gsub(" ", "_", Location))

bir_ped$V1 <- bir_ped$Location

bir_ped <- bir_ped %>%
  dplyr::select(-c(Location))

write.table(bir_ped, "data/ba3/birostris_bayesass3.ped",
            quote = F, row.names = F, col.names = F)

# Convert to immanc using PGDspider

# scp data/ba3/*inp ehumble@eddie.ecdf.ed.ac.uk:/exports/cmvm/eddie/eb/groups/ogden_grp/emily/manta_pop_gen_2020/data/out/stacks_PE/bayesass/
  

#~~~~~~~~~~~~~~~~~~~~~~#
#        Treemix       #
#~~~~~~~~~~~~~~~~~~~~~~#

meta <- fread("data/meta/sample_info_mastersheet_reformatted.csv") %>%
  select(c(Sample_code, Location))

# prepare file for freq calculation in plink

#scp ehumble@eddie.ecdf.ed.ac.uk:/exports/cmvm/eddie/eb/groups/ogden_grp/emily/manta_pop_gen_2020/file_lists/*treemix data/file_lists/

fread("data/file_lists/alfredi_treemix") %>%
  separate(V1, c("Sample", "dup")) %>%
  left_join(meta, by = c("Sample" = "Sample_code")) %>%
  unite("V1", Sample:dup) %>%
  mutate(V1 = gsub("_NA", "", V1)) %>%
  mutate(Location = gsub(" ", "_", Location)) %>%
  write.table("data/file_lists/alfredi_treemix",
              quote = F, col.names = F, row.names = F)

fread("data/file_lists/birostris_treemix") %>%
  separate(V1, c("Sample", "dup")) %>%
  left_join(meta, by = c("Sample" = "Sample_code")) %>%
  unite("V1", Sample:dup) %>%
  mutate(V1 = gsub("_NA", "", V1)) %>%
  mutate(Location = gsub(" ", "_", Location)) %>%
  write.table("data/file_lists/birostris_treemix",
              quote = F, col.names = F, row.names = F)

#scp data/file_lists/*treemix ehumble@eddie.ecdf.ed.ac.uk:/exports/cmvm/eddie/eb/groups/ogden_grp/emily/manta_pop_gen_2020/file_lists/
  
  
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#        POPULATIONS for VCFtool       #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~ Alfredi

alfredi <- fread("data/file_lists/alfredi_treemix") |>
  filter(V3 != "Sri_Lanka") %>%
  unite("V1", c(V1, V2), sep = "_")

alfredi_pops <- sort(unique(alfredi$V3))

list_alf <- alfredi |>
  select(V1, V3) |>
  group_split(V3)

names(list_alf) <- alfredi_pops
  
list_alf %>%
  purrr::iwalk(~write.table(.x, paste0("data/file_lists/alfredi_IDs_", .y, ".txt"), row.names = FALSE, col.names = F, quote = F))

#scp data/file_lists/*IDs* ehumble@eddie.ecdf.ed.ac.uk:/exports/cmvm/eddie/eb/groups/ogden_grp/emily/manta_pop_gen_2020/file_lists/


#~~~~~~~~~~~~~~~~~~~~~~~#
#        Genepop        #
#~~~~~~~~~~~~~~~~~~~~~~~#

meta <- fread("data/meta/sample_info_mastersheet_reformatted.csv") %>%
  dplyr::select(c(Sample_code, Location))

# ALFREDI

alf_ped <- fread("data/plink/across_sp/alfredi_geno4_dp6_ind_geno_depth_no_dups_maf_unrelated_ld.ped") %>%
  separate(V1, c("Sample", "dup")) %>%
  left_join(meta, by = c("Sample" = "Sample_code")) %>%
  unite("V1", Sample:dup) %>%
  mutate(V1 = gsub("_NA", "", V1)) %>%
  mutate(Location = gsub(" ", "_", Location))

alf_ped$V1 <- alf_ped$Location

alf_ped <- alf_ped %>%
  dplyr::select(-c(Location))

write.table(alf_ped, "data/plink/across_sp/alfredi_for_genepop.ped",
            quote = F, row.names = F, col.names = F)

# Convert to genepop using PGDspider

# BIROSTRIS

bir_ped <- fread("data/plink/across_sp/birostris_geno4_dp6_ind_geno_depth_no_dups_maf_unrelated_ld.ped") %>%
  separate(V1, c("Sample", "dup")) %>%
  left_join(meta, by = c("Sample" = "Sample_code")) %>%
  unite("V1", Sample:dup) %>%
  mutate(V1 = gsub("_NA", "", V1)) %>%
  mutate(Location = gsub(" ", "_", Location))

bir_ped$V1 <- bir_ped$Location

bir_ped <- bir_ped %>%
  dplyr::select(-c(Location))

write.table(bir_ped, "data/plink/across_sp/birostris_for_genepop.ped",
            quote = F, row.names = F, col.names = F)

# Convert to genepop using PGDspider

# SPECIES

#scp ehumble@eddie.ecdf.ed.ac.uk:/exports/cmvm/eddie/eb/groups/ogden_grp/emily/manta_pop_gen_2020/data/out/stacks_PE/across_sp_geno4_dp6_ind_geno_depth_no_dups_mac_geno.ped data/plink/across_sp/
#scp ehumble@eddie.ecdf.ed.ac.uk:/exports/cmvm/eddie/eb/groups/ogden_grp/emily/manta_pop_gen_2020/data/out/stacks_PE/across_sp_geno4_dp6_ind_geno_depth_no_dups_mac_geno.map data/plink/across_sp/

meta <- fread("data/meta/sample_info_mastersheet_reformatted.csv") %>%
  dplyr::select(c(Sample_code, Species))

sp_ped <- fread("data/plink/across_sp/across_sp_geno4_dp6_ind_geno_depth_no_dups_mac_geno.ped") %>%
  separate(V1, c("Sample", "dup")) %>%
  left_join(meta, by = c("Sample" = "Sample_code")) %>%
  unite("V1", Sample:dup) %>%
  mutate(V1 = gsub("_NA", "", V1)) %>%
  mutate(Species = gsub(" ", "_", Species))

sp_ped$V1 <- sp_ped$Species

sp_ped <- sp_ped %>%
  dplyr::select(-c(Species))

write.table(sp_ped, "data/plink/across_sp/across_sp_for_genepop.ped",
            quote = F, row.names = F, col.names = F)


#~~ Birostris