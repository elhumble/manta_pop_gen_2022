# Depth and popgen summary stats of SNP datasets

library(data.table)
library(dplyr)
library(tidyr)
library(diveRsity)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#      Across species       #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# scp ehumble@eddie.ecdf.ed.ac.uk:/exports/cmvm/eddie/eb/groups/ogden_grp/emily/manta_pop_gen_2020/data/out/stacks_PE/across_sp.idepth data/vcftools/across_sp/ 
# scp ehumble@eddie.ecdf.ed.ac.uk:/exports/cmvm/eddie/eb/groups/ogden_grp/emily/manta_pop_gen_2020/data/out/stacks_PE/across_sp.ldepth.mean data/vcftools/across_sp/ 
# scp ehumble@eddie.ecdf.ed.ac.uk:/exports/cmvm/eddie/eb/groups/ogden_grp/emily/manta_pop_gen_2020/data/out/stacks_PE/across_sp.gdepth data/vcftools/across_sp/ 
  
#~~ Overall individual depth

idepth <- fread("data/vcftools/across_sp/across_sp.idepth", colClasses = "numeric")

summary(idepth$MEAN_DEPTH)
hist(idepth$MEAN_DEPTH)

#~~ Mean locus depth

ldepth <- fread("data/vcftools/across_sp/across_sp.ldepth.mean", colClasses = "numeric")
summary(ldepth$MEAN_DEPTH)

CI <- 0.95
CI_meanDP <- stats::quantile(ldepth$MEAN_DEPTH, c((1-CI)/2,1-(1-CI)/2), na.rm=TRUE)
CI_sumDP <- stats::quantile(ldepth$MEAN_DEPTH, c((1-CI)/2,1-(1-CI)/2), na.rm=TRUE)
CI_meanDP
CI_sumDP


#~~ Genotype depth

gdepth <- fread("data/vcftools/across_sp/across_sp.gdepth", colClasses = "numeric") %>%
  unite(SNP, CHROM, POS) # join snp, chrom and pos cols

View(gdepth[1:100,1:21])

# Select genotypes only
gdepth <- select(gdepth, -SNP)

# Replace -1 with 0
gdepth[gdepth == -1] <- 0

# Coverage stats
gdepth <- gdepth %>%
  mutate(mean = rowMeans(.), 
         sum = rowSums(.))

View(gdepth[1:100,1:21])

summary(gdepth$mean)
hist(gdepth$mean)

CI <- 0.90
CI_meanDP <- stats::quantile(gdepth$mean, c((1-CI)/2,1-(1-CI)/2), na.rm=TRUE)
CI_sumDP <- stats::quantile(gdepth$sum, c((1-CI)/2,1-(1-CI)/2), na.rm=TRUE)
CI_meanDP
CI_sumDP


#~~ IMISS

# scp ehumble@eddie.ecdf.ed.ac.uk:/exports/cmvm/eddie/eb/groups/ogden_grp/emily/manta_pop_gen_2020/data/out/stacks_PE/across_sp_geno4_dp6.imiss data/vcftools/across_sp/ 

imiss <- fread("data/vcftools/across_sp/across_sp_geno4_dp6.imiss")
hist(imiss$F_MISS, breaks = 20)

nrow(filter(imiss, F_MISS < 0.9))
nrow(filter(imiss, F_MISS < 0.6))
nrow(filter(imiss, F_MISS < 0.5))
nrow(filter(imiss, F_MISS < 0.45))
nrow(filter(imiss, F_MISS < 0.4))
nrow(filter(imiss, F_MISS < 0.3))
nrow(filter(imiss, F_MISS < 0.2))

summary(imiss$F_MISS)

#~~ Filtered site mean depth

# scp ehumble@eddie.ecdf.ed.ac.uk:/exports/cmvm/eddie/eb/groups/ogden_grp/emily/manta_pop_gen_2020/data/out/stacks_PE/across_sp_geno4_dp6_ind.ldepth.mean data/vcftools/across_sp/ 
# scp ehumble@eddie.ecdf.ed.ac.uk:/exports/cmvm/eddie/eb/groups/ogden_grp/emily/manta_pop_gen_2020/data/out/stacks_PE/across_sp_geno4_dp6_ind_geno_depth.ldepth.mean data/vcftools/across_sp/ 

ldepth <- fread("data/vcftools/across_sp/across_sp_geno4_dp6_ind.ldepth.mean", colClasses = "numeric")
ldepth <- fread("data/vcftools/across_sp/across_sp_geno4_dp6_ind_geno_depth.ldepth.mean", colClasses = "numeric")

summary(ldepth$MEAN_DEPTH)

CI <- 0.95
CI_meanDP <- stats::quantile(ldepth$MEAN_DEPTH, c((1-CI)/2,1-(1-CI)/2), na.rm=TRUE)
CI_sumDP <- stats::quantile(ldepth$MEAN_DEPTH, c((1-CI)/2,1-(1-CI)/2), na.rm=TRUE)
CI_meanDP
CI_sumDP

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#      Population genetic statistics     #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

infile <- "data/plink/across_sp/alfredi_geno4_dp6_ind_geno_depth_no_dups_maf_unrelated_ld_genepop.txt"

genepop_alf <- readGenepop(infile = infile, gp = 3, bootstrap = FALSE)

genepop_alf$npops

alf_snp_stats <- basicStats(infile = "data/plink/across_sp/alfredi_geno4_dp6_ind_geno_depth_no_dups_maf_unrelated_ld_genepop.txt", 
                           outfile = "basic_stat", 
                           fis_ci = TRUE, ar_ci = TRUE, fis_boots = 999, 
                           ar_boots = 999, mc_reps = 9999, 
                           rarefaction = FALSE, ar_alpha = 0.05, 
                           fis_alpha = 0.05)

rownames(alf_snp_stats$main_tab$`1376`)
alf_snp_stats$main_tab$`1277`$overall
alf_snp_stats$main_tab$`1336`$overall
alf_snp_stats$main_tab$`0133`$overall
alf_snp_stats$main_tab$`0689`$overall
alf_snp_stats$main_tab$`0141`$overall
alf_snp_stats$main_tab$`1376`$overall

alf_snp_stats_df <- data.frame(statistic = rownames(alf_snp_stats$main_tab$`1376`),
                               maldives = alf_snp_stats$main_tab$`1277`$overall,
                               hawaii = alf_snp_stats$main_tab$`1336`$overall,
                               seychelles = alf_snp_stats$main_tab$`0133`$overall,
                               chagos = alf_snp_stats$main_tab$`0689`$overall,
                               fiji = alf_snp_stats$main_tab$`0141`$overall,
                               aus_pac = alf_snp_stats$main_tab$`1376`$overall) %>%
  pivot_longer(cols = maldives:aus_pac) %>%
  pivot_wider(names_from = name, values_from = value)

#~~ birostris

bir_snp_stats <- basicStats(infile = "data/plink/across_sp/birostris_geno4_dp6_ind_geno_depth_no_dups_maf_unrelated_ld_genepop.txt", 
                              outfile = "basic_stat", 
                              fis_ci = TRUE, ar_ci = TRUE, fis_boots = 999, 
                              ar_boots = 999, mc_reps = 9999, 
                              rarefaction = FALSE, ar_alpha = 0.05, 
                              fis_alpha = 0.05)

rownames(bir_snp_stats$main_tab$`0732`)
bir_snp_stats$main_tab$`0732`$overall
bir_snp_stats$main_tab$`1152`$overall
bir_snp_stats$main_tab$`1104`$overall
bir_snp_stats$main_tab$`0985`$overall
bir_snp_stats$main_tab$`1062`$overall
bir_snp_stats$main_tab$`1183`$overall

bir_snp_stats_df <- data.frame(statistic = rownames(bir_snp_stats$main_tab$`0732`),
                               sri_lanka = bir_snp_stats$main_tab$`0732`$overall,
                               mex_pac = bir_snp_stats$main_tab$`1152`$overall,
                               phil = bir_snp_stats$main_tab$`1104`$overall,
                               mex_car = bir_snp_stats$main_tab$`0985`$overall,
                               peru = bir_snp_stats$main_tab$`1062`$overall,
                               sa = bir_snp_stats$main_tab$`1183`$overall) %>%
  pivot_longer(cols = sri_lanka:sa) %>%
  pivot_wider(names_from = name, values_from = value)


#~~ species

infile <- "data/plink/across_sp/across_sp_geno4_dp6_ind_geno_depth_no_dups_mac_geno_genepop.txt"
genepop_sp <- readGenepop(infile = infile, gp = 3, bootstrap = FALSE)
genepop_sp$npops

sp_snp_stats <- basicStats(infile = "data/plink/across_sp/across_sp_geno4_dp6_ind_geno_depth_no_dups_mac_geno_genepop.txt", 
                            outfile = "basic_stat", 
                            fis_ci = TRUE, ar_ci = TRUE, fis_boots = 999, 
                            ar_boots = 999, mc_reps = 9999, 
                            rarefaction = FALSE, ar_alpha = 0.05, 
                            fis_alpha = 0.05)

rownames(sp_snp_stats$main_tab$`1277`)
sp_snp_stats$main_tab$`1277`$overall
sp_snp_stats$main_tab$`0732`$overall

sp_snp_stats_df <- data.frame(statistic = rownames(sp_snp_stats$main_tab$`1277`),
                               alfredi = sp_snp_stats$main_tab$`1277`$overall,
                               birostris = sp_snp_stats$main_tab$`0732`$overall) %>%
  pivot_longer(cols = alfredi:birostris) %>%
  pivot_wider(names_from = name, values_from = value)





