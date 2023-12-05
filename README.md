**Analysis code for:**
-------------
Humble E, Hosegood J, Carvalho G, de Bruyn M, Creer S, Stevens GMW, Armstrong A, Bonfil R, Deakos M, Fernando D, Froman N, Peel LR, Pollett S, Ponzo A, Stewart JD, Wintner S, Ogden R **Comparative population genomics of manta rays has global implications for management.** *Molecular Ecology* https://doi.org/10.1111/mec.17220  

**Summary**
-------------

This repository contains the scripts required for analysis of population structure, differentiation, gene flow, heterozygosity and historical relationships. The repository also contains scripts for producing all figures in the manuscript.

For the STACKS SNP calling and filtering pipeline refer to `workflows.`

**Code structure**
-------------

*Scripts for use in STACKS pipeline*  
`1.0_cat_demultiplexed_files.R`  
`1.1_prep_pop_maps.R`  
`1.2_param_optimisation.R`  
`1.3_get_catalog.R`  

*SNP QC and file preparation*   
`2.0_snp_sumstats.R`  
`2.1_prep_dup_files.R`  
`2.2_error_rate.R`  
`2.3_relatedness.R`  
`3.0_make_files.R`  

*Population structure and gene flow*    
`4.0_PCA_DAPC.R`  
`4.1_admixture_out.R`  
`4.2_map.R`  

*Genetic diversity and downstream treemix*  
`5.0_IBD_dbMEM_RDA.R`  
`5.1_FST_WC.R`  
`6.0_het.R`  
`7.0_ba3.R`  
`8.0_treemix.R`  

*Figure generation*   
`9.0_figures.R`  

**Data**
-------------
Sequencing reads are available at the European Nucleotide Archive (https://www.ebi.ac.uk/ena/browser) under study accession [PRJEB66437](https://www.ebi.ac.uk/ena/browser/view/PRJEB66437).

###### Please feel free to [get in touch](mailto:emily.humble@ed.ac.uk) if you have any questions.
