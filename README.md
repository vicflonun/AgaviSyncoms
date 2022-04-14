# AgaviSyncoms
This repository contains the code used to perform the upstream and downstream statistical analyses in Flores-Nuñez, et al. 2022

ampagave2_5.bash : shell pipeline that processes 16SrRNAV4 and ITS2 amplicon sequencing paired reads to generate OTUs, OTU table and Taxonomic classification. see citation on manuscript. 
getOTUdata_6.R : R script that generates filtered OTU, taxa and metadata tables for upstream analysis (full OTU tables not provided)

./data : includes 16S and ITS filtered OTUs, OTU table,  taxonomic classification and metadata

agave_barplot_3.0.R : Script that generates relative abundance barplots for Figure 1 and more
agave_venn2.0.R : Script that generates venn diagrams of OTUs between seasons for Figure 1
agave_season_erichments3.0.R: Script that performs an OTU enrichment analysis between seasons and triplots for figure 1

agave.permanova1.0.R : Script for PERMANOVA analysis 
agave.diversity.3.0.R : Script that performs a diversity analysis and plots for Figure 2 and 3

agave_treatment_erichments3.0.R : Script that performs an OTU enrichment analysis between treatments and plots
agave_strain_erichments2.0.R : Script that performs an OTU (inoculated strains) enrichment analysis between seasons and tratments 

get_hubs7.0_syncoms.R: Script that performs a correlation network analysis between treatments 
