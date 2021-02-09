library(phyloseq)
library(dplyr)
library(metagMisc)
library(PhyloMeasures)
library(Rmisc)
library(tidyr)
library(nlme)
library(car)
library(emmeans)
library(multcomp)
library(ggpubr)

##Load tree_added.Rdata

setwd("/Users/mono5/Dropbox (GaTech)/PhD/GA_TECH/Dissertation/Spartina/2018-2019_Summary/16S_sequences/")
setwd("~/Dropbox/PhD/GA_TECH/Dissertation/Spartina/2018-2019_Summary/16S_sequences/")

phy.tree.prev.filt2 <- readRDS("./Sequencing_December_2019/phy.tree.prev.filt2.rds")

phy.bulk_sediment <- subset_samples(phy.tree.prev.filt2, (Compartment == "Bulk_sediment" & Location == "Sapelo") | 
                                      (Compartment == "Bulk_sediment" & Location == "Skidaway"))

metadata_summary <- data.frame(sample_data(phy.bulk_sediment))
metadata_summary[metadata_summary$Transect == 6 & metadata_summary$Point == 2,]$Spartina <- "Medium"

metadata_summary$Crab.burrows <- metadata_summary$Crab.burrows*4
metadata_summary$Snails.density <- metadata_summary$Snails.density*4

metadata_summary$Density <- NULL
metadata_summary$CB_rate <- NULL
metadata_summary$XYL_rate <- NULL

write.table(metadata_summary, "../Paper_spartina_assembly/Submission_Microbiome/Spartina_GA_Core_Microbiome/metadata_per_sampling_point.tsv", 
    quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
