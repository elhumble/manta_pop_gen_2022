# Script to calculate error rate from duplicate individuals

library(dplyr)
library(tidyr)
library(data.table)

#~~~~~~~~~~~~~~~~~~~~~~~~~~#
#      Across species      #
#~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~ Read in plink output

# scp ehumble@eddie.ecdf.ed.ac.uk:/exports/cmvm/eddie/eb/groups/ogden_grp/emily/manta_pop_gen_2020/data/out/stacks_PE/error_rate_dups.genome data/plink/across_sp/ 
  
error <- fread("data/plink/across_sp/error_rate_dups.genome", header = T) %>%
  separate(FID1, c("ID1", "dup1"), sep = "_") %>%
  separate(FID2, c("ID2", "dup2"), sep = "_") %>%
  filter(ID1 == ID2)

hist(error$PI_HAT, breaks = 10)

meta <- read.csv("data/meta/sample_info_mastersheet_reformatted.csv") %>%
  select(Sample_code, Species, Location)

error <- left_join(error, meta, by = c("ID1" = "Sample_code")) %>%
  left_join(meta, by = c("ID2" = "Sample_code"))

filter(error, Z2 > 0.9)

1 - (error %>%
       summarise(mean(Z2)))

# 0.05625556
# Z2 score: the proportion of SNPs at which two individuals (replicates) share both alleles IBD


