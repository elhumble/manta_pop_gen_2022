#~~ Write script to cat R1 and R2 primary and remainder files output from process radtags

library(dplyr)
library(glue)
library(purrr)

# file names txt files downloaded from eddie into data/

file_names <- list.files(path = "data/file_lists/", pattern="^file_names*")
file_names <- paste0("data/", file_names)

#~~~~~~~~~~~~~~~~~~~~~~#
#      All libs        # 
#~~~~~~~~~~~~~~~~~~~~~~#

# Prepare dataframe to create command

run <- file_names %>% 
  map_dfr(read.table, header = F, .id = "lib") %>%
  filter(!grepl("\\.2", V1)) %>%
  filter(!grepl("rem", V1)) %>%
  mutate(V2 = V1,
         rem_R1 = V1,
         rem_R2 = V1) %>%
  mutate(V2 = gsub("\\.1.", ".2.", V1),
         rem_R1 = gsub("\\.1", ".rem.1", rem_R1),
         rem_R2 = gsub ("\\.1", ".rem.2", rem_R2)) %>%
  mutate(out = gsub("\\.1.fq", ".fq", V1))

# Create command

command <- glue("cat /exports/eddie/scratch/ehumble/lib_{run$lib}/{run$V1} ",
                "/exports/eddie/scratch/ehumble/lib_{run$lib}/{run$V2} ",
                "/exports/eddie/scratch/ehumble/lib_{run$lib}/{run$rem_R1} ",
                "/exports/eddie/scratch/ehumble/lib_{run$lib}/{run$rem_R2} > ", 
                "/exports/eddie/scratch/ehumble/reads/{run$out}")

# Collapse into individual rows

command <- glue_collapse(command, sep = "\n")


# Use glue to stick everything together and create gsub command

rt <- "1:00:00"
vmem <- "1G"

geo <-
  glue(
    "#!/bin/sh \n# Grid Engine Options \n#$ -N cat \n#$ -cwd \n#$ -l h_rt={rt}",
    "\n#$ -l h_vmem={vmem} \n#$ -R y \n#$ -e e_files \n#$ -o o_files",
    "\n\n# Jobscript to cat output of process radtags",
    "\n. /etc/profile.d/modules.sh \n\n",
    "\nmkdir /exports/eddie/scratch/ehumble/reads")

write.table(
  glue("{geo} \n\n {command}"),
  "data/scripts/1.2_cat_demultiplexed_files.sh",
  quote = F,
  col.names = F,
  row.names = F
)


#~~ Upload to EDDIE

scp data/1.2_cat_demultiplexed_files.sh ehumble@eddie.ecdf.ed.ac.uk:/exports/cmvm/eddie/eb/groups/ogden_grp/emily/manta_pop_gen_2020/1.2_cat_demultiplexed_files.sh

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#           For BWA           # 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


bwa <- file_names %>% 
  map_dfr(read.table, header = F, .id = "lib") %>%
 # filter(!grepl("\\.2", V1)) %>%
  filter(!grepl("rem", V1)) %>%
  mutate(V2 = V1) %>%
  mutate(V2 = gsub("\\.1", ".rem.1", V2),
         V2 = gsub ("\\.2", ".rem.2", V2)) %>%
  mutate(out = V1)


command <- glue("cat /exports/eddie/scratch/ehumble/lib_all/{bwa$V1} ",
                "/exports/eddie/scratch/ehumble/lib_all/{bwa$V2} > ", 
                "/exports/eddie/scratch/ehumble/bwa/{bwa$out}")

# Collapse into individual rows

command <- glue_collapse(command, sep = "\n")

# Use glue to stick everything together and create gsub command

rt <- "1:00:00"
vmem <- "1G"

geo <-
  glue(
    "#!/bin/sh \n# Grid Engine Options \n#$ -N cat \n#$ -cwd \n#$ -l h_rt={rt}",
    "\n#$ -l h_vmem={vmem} \n#$ -R y \n#$ -e e_files \n#$ -o o_files",
    "\n\n# Jobscript to cat output of process radtags",
    "\n. /etc/profile.d/modules.sh \n\n",
    "\nmkdir /exports/eddie/scratch/ehumble/bwa")

write.table(
  glue("{geo} \n\n {command}"),
  "data/scripts/4.0_cat_demultiplexed_files_bwa.sh",
  quote = F,
  col.names = F,
  row.names = F
)



#~~ Upload to EDDIE
scp data/4.0_cat_demultiplexed_files_bwa.sh ehumble@eddie.ecdf.ed.ac.uk:/exports/cmvm/eddie/eb/groups/ogden_grp/emily/manta_pop_gen_2020/4.0_cat_demultiplexed_files_bwa.sh

