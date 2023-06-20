# Summary stats of raw SNPs

library(data.table)
library(dplyr)
library(tidyr)

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

ldepth <- fread("data/vcftools/across_sp/across_sp_geno4_dp6_ind.ldepth.mean", colClasses = "numeric")
summary(ldepth$MEAN_DEPTH)

CI <- 0.95
CI_meanDP <- stats::quantile(ldepth$MEAN_DEPTH, c((1-CI)/2,1-(1-CI)/2), na.rm=TRUE)
CI_sumDP <- stats::quantile(ldepth$MEAN_DEPTH, c((1-CI)/2,1-(1-CI)/2), na.rm=TRUE)
CI_meanDP
CI_sumDP
