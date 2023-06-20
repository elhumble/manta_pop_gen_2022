library(dplyr)
library(tidyr)
library(data.table)


#~~ Write file with dup IDs for plink

dups <- fread("data/vcftools/across_sp/across_sp_geno4_dp5.imiss") %>%
  separate(INDV, c("ID", "dup"), sep = "_") %>%
  # Group by ID
  group_by(ID) %>%
  filter(n()>1) %>%
  unite("INDV", ID, dup) %>%
  mutate(INDV = gsub("_NA", "", INDV))

write.table(dups[c(1,1)], "data/file_lists/duplicates.txt",
            col.names = F, row.names = F, quote = F, sep = "\t")


#~~ Get duplicate individuals with highest genotyping rate for filtering

# Select duplicate individual with highest number of reads retained

highest_dup <- fread("data/vcftools/across_sp/across_sp_geno4_dp5.imiss") %>%
  separate(INDV, c("ID", "dup"), sep = "_") %>%
  # Group by ID
  group_by(ID) %>%
  filter(n()>1) %>%
  filter(F_MISS == min(F_MISS)) %>%
  unite("INDV", ID, dup) %>%
  mutate(INDV = gsub("_NA", "", INDV))

non_dups <- fread("data/vcftools/across_sp/across_sp_geno4_dp5.imiss") %>%
  separate(INDV, c("ID", "dup"), sep = "_") %>%
  # Group by ID
  group_by(ID) %>%
  filter(n()==1) %>%
  unite("INDV", ID, dup) %>%
  mutate(INDV = gsub("_NA", "", INDV))

highest_dup <- rbind(highest_dup, non_dups)

write.table(highest_dup[c(1,1)], "data/file_lists/highest_dups.txt",
            col.names = F, row.names = F, quote = F, sep = "\t")

scp data/file_lists/duplicates.txt ehumble@eddie.ecdf.ed.ac.uk:/exports/cmvm/eddie/eb/groups/ogden_grp/emily/manta_pop_gen_2020/file_lists/duplicates.txt

scp data/file/lists/highest_dups.txt ehumble@eddie.ecdf.ed.ac.uk:/exports/cmvm/eddie/eb/groups/ogden_grp/emily/manta_pop_gen_2020/file_lists/highest_dups.txt
