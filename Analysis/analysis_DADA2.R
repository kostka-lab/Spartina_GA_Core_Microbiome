library(dada2); packageVersion("dada2")
library(phyloseq)
library(microbiome)
library(metagenomeSeq)
library(dplyr)
library(Glimma)
library(metagMisc)
library(PhyloMeasures)
library(Rmisc)
library(DESeq2)
library(beepr)
library(tidyr)
library(nlme)
library(car)
library(emmeans)
library(multcomp)
library(Rmisc)
library(ggpubr)
library(phangorn)
library(ggplot2)
library(ampvis)
library(lme4)
library(vegan)
library(ggtree)
library(scales)
library(phytools)
library(picante)

setwd(".")

#WILL RUN DADA 2 IN TWO PARTS BECAUSE THIS PROJECT WAS SEQUNCED IN TWO INDEPENDENT RUNS

#RUN 1: ALL BULK SEDIMENT FROM SAPELO YEAR 2018

# Forward and reverse fastq filenames have format: SAMPLENAME_L001_R1_001.fastq and SAMPLENAME_L001_R2_001.fastq
fnFs = sort(list.files(pattern="R1_trimmed.fastq.gz", full.names = TRUE))
fnRs = sort(list.files(pattern="R2_trimmed.fastq.gz", full.names = TRUE))
sample.names = sapply(strsplit(basename(fnFs), "_R"), `[`, 1)

# Examine quality of F reads
plotQualityProfile(fnFs[1:12])
## Will truncate at 190bp

# Examine quality of R reads
plotQualityProfile(fnRs[1:12])
## Will truncate at 175bp

# Assign filenames for the filtered files
filt_path = file.path("filtered")
filtFs = file.path(filt_path, paste0(sample.names, "_F_filt.fastq"))
filtRs = file.path(filt_path, paste0(sample.names, "_R_filt.fastq"))

# Filter the F and R reads
out = filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(175,150),
                    maxN=0, maxEE=c(2,2), truncQ=10, rm.phix=TRUE,
                    compress=TRUE, multithread=FALSE) 

# Examine quality of F reads
plotQualityProfile(filtFs[1:12])

# Examine quality of R reads
plotQualityProfile(filtRs[1:12])


##### ERROR RATES #####

errF = learnErrors(filtFs, multithread=FALSE, verbose = TRUE)
errR = learnErrors(filtRs, multithread=FALSE, verbose = TRUE)


# Plot the estimated error rates
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)

# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names


#Infer the sequence variants in each sample

dadaFs <- dada(derepFs, err=errF, multithread=TRUE, pool = FALSE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE, pool = FALSE)

#Inspect the returning dada-class object

dadaFs[[1]]
dadaRs[[1]]

##MERGE PAIRED READS

mergers_1 <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE, maxMismatch = 1)

# Inspect the merger data.frame from the first sample
head(mergers_1[[2]])

saveRDS(mergers_1, "../mergers_1.rds") # CHANGE ME to where you want sequence table saved



###############DO SAME ANALYSIS FOR THE SECOND RUN. ALL REMAINING SEQUENCES


# Forward and reverse fastq filenames have format: SAMPLENAME_L001_R1_001.fastq and SAMPLENAME_L001_R2_001.fastq
fnFs = sort(list.files(pattern="R1_trimmed.fastq.gz", full.names = TRUE))
fnRs = sort(list.files(pattern="R2_trimmed.fastq.gz", full.names = TRUE))
sample.names = sapply(strsplit(basename(fnFs), "_R1"), `[`, 1)

# Examine quality of F reads
plotQualityProfile(fnFs[1:8])
## Truncate at 175bp

# Examine quality of R reads
plotQualityProfile(fnRs[1:8])
## Truncate at 150bp

# Assign filenames for the filtered files
filt_path = file.path("../filtered")
filtFs = file.path(filt_path, paste0(sample.names, "_F_filt.fastq"))
filtRs = file.path(filt_path, paste0(sample.names, "_R_filt.fastq"))

# Filter the F and R reads
out_2 = filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(175,150),
                    maxN=0, maxEE=c(2,2), truncQ=10, rm.phix=TRUE,
                    compress=TRUE, multithread=FALSE) 

out <- rbind(out, out_2)

# Examine quality of F reads
plotQualityProfile(filtFs[9:20])

# Examine quality of R reads
plotQualityProfile(filtRs[9:20])


##### ERROR RATES #####

errF = learnErrors(filtFs, multithread=FALSE, verbose = TRUE)
beep()
errR = learnErrors(filtRs, multithread=FALSE, verbose = TRUE)
beep()


# Plot the estimated error rates
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

filtFs1 <- filtFs[1:64]
filtFs2 <- filtFs[65:128]
filtFs3 <- filtFs[129:184]

filtRs1 <- filtRs[1:64]
filtRs2 <- filtRs[65:128]
filtRs3 <- filtRs[129:184]


derepFs1 <- derepFastq(filtFs1, verbose=TRUE)
derepRs1 <- derepFastq(filtRs1, verbose=TRUE)

getN = function(x) sum(getUniques(x))
dadaFs <- sapply(dadaFs, getN)
dadaRs <- sapply(dadaRs, getN)

rm(derepFs)
rm(derepRs)

# Name the derep-class objects by the sample names
names(derepFs1) <- sample.names[1:64]
names(derepRs1) <- sample.names[1:64]


#Infer the sequence variants in each sample

dadaFs1 <- dada(derepFs1, err=errF, multithread=FALSE, pool = FALSE)
dadaRs1 <- dada(derepRs1, err=errR, multithread=FALSE, pool = FALSE)


#Inspect the returning dada-class object

dadaFs1[[1]]
dadaRs1[[1]]

##MERGE PAIRED READS

mergers_2a <- mergePairs(dadaFs1, derepFs1, dadaRs1, derepRs1, verbose=TRUE, maxMismatch = 1)

# Inspect the merger data.frame from the first sample
head(mergers_2a[[64]])

saveRDS(mergers_2a, "../mergers_2a.rds") # CHANGE ME to where you want sequence table saved


################ Merge the datasets

st1 <- readRDS("../../mergers_1.rds")
st2 <- readRDS("../mergers_2a.rds")

seqtab1 <- makeSequenceTable(st1)
seqtab2 <- makeSequenceTable(st2)


# We now have a data.frame for each sample with the merged $sequence, its $abundance, and the indices of the merged 
# $forward and  $reverse denoised sequences. Paired reads that did not exactly overlap were removed by mergePairs

##### CONSTRUCT SEQUENCE TABLE #####

# The seq table is analagous to the OTU table 
seqtab = mergeSequenceTables(seqtab1, seqtab2) 


# Select only sequences that are within the expected read length
seqtab.filt = seqtab[,nchar(colnames(seqtab)) %in% seq(251,255)]


# Plot
plot(table(nchar(getSequences(seqtab.filt))))


##### REMOVE CHIMERAS #####

seqtab.nochim = removeBimeraDenovo(seqtab.filt, method="consensus", multithread=TRUE, verbose=TRUE)

# Check proportion of chimeric sequences
sum(seqtab.nochim)/sum(seqtab.filt)

##### TRACK SEQUENCES THROUGH PIPELINE #####

getN = function(x) sum(getUniques(x))

s1_getN <- data.frame(sapply(st1, getN))
s2_getN <- data.frame(sapply(st2, getN))
names(s1_getN) = names(s2_getN) = "merged"
st_all <- rbind(s1_getN, s2_getN)

dadaFs_all <- as.data.frame(c(dadaFs, dadaFs1))
names(dadaFs_all) <- "denoised"

track = cbind(out, dadaFs_all, st_all, rowSums(seqtab),
              rowSums(seqtab.filt), rowSums(seqtab.nochim))


# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) = c("input", "filtered_quality", "denoised", "merged", "tabled", "filtered_read_size", "nonchim")
track

summary(track)
apply(track, 2, median)


##### ASSIGN TAXONOMY #####

sq.nochim <- getSequences(seqtab.nochim)

#Perfom in chunks if memory is an issue:

t1 <- assignTaxonomy(sq.nochim[1:5000],
                     "./silva_nr_v132_train_set.fa.gz", multithread=TRUE)
saveRDS(t1, "../t1.rds") # CHANGE ME to where you want sequence table saved
beep()

t2 <- assignTaxonomy(sq.nochim[5001:10000],
                     "./silva_nr_v132_train_set.fa.gz", multithread=TRUE)
saveRDS(t2, "../t2.rds") # CHANGE ME to where you want sequence table saved

t3 <- assignTaxonomy(sq.nochim[10001:15000],
                     "./silva_nr_v132_train_set.fa.gz", multithread=TRUE)
saveRDS(t3, "../t3.rds") # CHANGE ME to where you want sequence table saved
beep()

t4 <- assignTaxonomy(sq.nochim[15001:20000],
                     "./silva_nr_v132_train_set.fa.gz", multithread=TRUE)
saveRDS(t4, "../t4.rds") # CHANGE ME to where you want sequence table saved

t5 <- assignTaxonomy(sq.nochim[20001:25000],
                     "./silva_nr_v132_train_set.fa.gz", multithread=TRUE)
saveRDS(t5, "../t5.rds") # CHANGE ME to where you want sequence table saved

t6 <- assignTaxonomy(sq.nochim[25001:30000],
                     "./silva_nr_v132_train_set.fa.gz", multithread=TRUE)
saveRDS(t6, "../t6.rds") # CHANGE ME to where you want sequence table saved
beep()

t7 <- assignTaxonomy(sq.nochim[30001:35000],
                     "./silva_nr_v132_train_set.fa.gz", multithread=TRUE)
saveRDS(t7, "../t7.rds") # CHANGE ME to where you want sequence table saved

t8 <- assignTaxonomy(sq.nochim[35001:40000],
                     "./silva_nr_v132_train_set.fa.gz", multithread=TRUE)
saveRDS(t8, "../t8.rds") # CHANGE ME to where you want sequence table saved

t9 <- assignTaxonomy(sq.nochim[40001:45000],
                     "./silva_nr_v132_train_set.fa.gz", multithread=TRUE)
saveRDS(t9, "../t9.rds") # CHANGE ME to where you want sequence table saved
beep()

t10 <- assignTaxonomy(sq.nochim[45001:50000],
                     "./silva_nr_v132_train_set.fa.gz", multithread=TRUE)
saveRDS(t10, "../t10.rds") # CHANGE ME to where you want sequence table saved

t11 <- assignTaxonomy(sq.nochim[50001:55000],
                     "./silva_nr_v132_train_set.fa.gz", multithread=TRUE)
saveRDS(t11, "../t11.rds") # CHANGE ME to where you want sequence table saved

t12 <- assignTaxonomy(sq.nochim[55001:60000],
                     "./silva_nr_v132_train_set.fa.gz", multithread=TRUE)
saveRDS(t12, "../t12.rds") # CHANGE ME to where you want sequence table saved
beep()


t13 <- assignTaxonomy(sq.nochim[60001:65000],
                      "./silva_nr_v132_train_set.fa.gz", multithread=TRUE)
saveRDS(t13, "../t13.rds") # CHANGE ME to where you want sequence table saved
beep()

t14 <- assignTaxonomy(sq.nochim[65001:70000],
                      "./silva_nr_v132_train_set.fa.gz", multithread=TRUE)
saveRDS(t14, "../t14.rds") # CHANGE ME to where you want sequence table saved

t15 <- assignTaxonomy(sq.nochim[70001:74376],
                      "./silva_nr_v132_train_set.fa.gz", multithread=TRUE)
saveRDS(t15, "../t15.rds") # CHANGE ME to where you want sequence table saved
beep()

t1 <- readRDS("../t1.rds")
t2 <- readRDS("../t2.rds")
t3 <- readRDS("../t3.rds")
t4 <- readRDS("../t4.rds")
t5 <- readRDS("../t5.rds")
t6 <- readRDS("../t6.rds")
t7 <- readRDS("../t7.rds")
t8 <- readRDS("../t8.rds")
t9 <- readRDS("../t9.rds")
t10 <- readRDS("../t10.rds")
t11 <- readRDS("../t11.rds")
t12 <- readRDS("../t12.rds")
t13 <- readRDS("../t13.rds")
t14 <- readRDS("../t14.rds")
t15 <- readRDS("../t15.rds")

taxa <- rbind(t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15)


##################PHYLOSEQ######################

samples <- unique(row.names(seqtab.nochim))
write.table(samples, file = "../../samples_names.txt", sep = "\t")

#Upload metadata
samples <- read.csv("./metadata_new.csv", header = TRUE)
samples$Spartina = factor(samples$Spartina,levels(samples$Spartina)[c(2,1,3)])
samples$Compartment = factor(samples$Compartment,levels(samples$Compartment)[c(2,1,3)])
rownames(samples) <- samples$Seq_ID
samples$Transect <- as.factor(samples$Transect)
samples$Point <- as.factor(samples$Point)
samples$Year <- as.factor(samples$Year)


phy <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samples), 
               tax_table(taxa))

phy
summary(sample_sums(phy))


#Remove Chloroplast, Mitochondria, Eukaryots, and Prokaryotes with unknown Phylum
table(tax_table(phy)[, "Kingdom"], exclude = NULL)

phy.raw = subset_taxa(phy, Order != "Chloroplast" | is.na(Order))
phy.raw = subset_taxa(phy.raw, Family != "Mitochondria" | is.na(Family))
phy.raw
table(tax_table(phy.raw)[, "Kingdom"], exclude = NULL)
summary(sample_sums(phy.raw))


phy.filt = subset_taxa(phy.raw, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized") & !Kingdom %in% c("", "Eukaryota"))
phy.filt
table(tax_table(phy.filt)[, "Kingdom"], exclude = NULL)
summary(sample_sums(phy.filt))


# Select ASVs that are present in at least 5% of samples OR have at least 10 reads in all samples
phy.prev.filt = phyloseq_filter_prevalence(phy.filt, prev.trh = 0.05, abund.trh = 10,
                                           threshold_condition = "OR", abund.type = "total")
phy.prev.filt
summary(sample_sums(phy.prev.filt))


#Export FASTA file to make phylogenetic tree
writeFasta(row.names(tax_table(phy.prev.filt)), "../../alignment_tree/ASVs.fasta")
#aligned ASVs in mothur, made tree in fasttree, rooted tree in figtree
#Upload tree
ASVs_tree <- read_tree("../../alignment_tree/Spartina.ASVs.rooted.tre")
#Need to change Tree name to fit taxa table
ASV_names <- read.dna(file = "../../alignment_tree/ASVs.microbe.fasta", 
                      format = 'fasta', as.character = TRUE)

ASV_names_df <- tibble(name = names(ASV_names), sequence = rep(NA, length(ASV_names)))
for(i in 1:length(ASV_names)){
  ASV_names_df[i,2] <- toupper(paste(ASV_names[[i]], collapse = ''))
}
ASVs_tree$tip.label <- ASV_names_df[[2]][match(ASVs_tree$tip.label, ASV_names_df[[1]])]


samples_prev.filt <- sample_data(phy.prev.filt)
otu_prev.filt <- otu_table(phy.prev.filt)
tax_prev.filt <- tax_table(phy.prev.filt)

phy.tree.prev.filt <- phyloseq(otu_table(otu_prev.filt, taxa_are_rows=FALSE), 
                                 sample_data(samples_prev.filt), 
                                 tax_table(tax_prev.filt),
                                phy_tree(ASVs_tree))

save.image("../tree_added.RData")





###################ALPHA AND BETA DIVERSITY ANALYSES - FIGURE 3#############################

abundance_all <- as.data.frame(sample_sums(phy.tree.prev.filt))
names(abundance_all) <- "Counts"

samples_all <- data.frame(sample_data(phy.tree.prev.filt))
samples_all$Seq_ID <- as.character(samples_all$Seq_ID)

abundance_all$Seq_ID <- samples_all$Seq_ID


alpha_div <- estimate_richness(phy.tree.prev.filt, split = TRUE,
                               measures =c("Observed", "Shannon", "Simpson")) 
alpha_div$Seq_ID <- samples_all$Seq_ID

alpha_div <- merge(samples_all, alpha_div)
alpha_div <- merge(alpha_div, abundance_all)

#Plot Shannon

theme_alpha <- theme_bw() + theme(axis.title.x = element_blank(),
                                  axis.title.y = element_text(size = 18),
                                  axis.text = element_text(size = 16),
                                  panel.grid = element_blank(),
                                  legend.position = "top", legend.title = element_blank(),
                                  legend.text = element_text(size = 16))

alpha_div$Compartment <- gsub("Bulk_sediment", "Bulk sediment", alpha_div$Compartment)
alpha_div$Compartment <- as.factor(alpha_div$Compartment)
alpha_div$Compartment <- factor(alpha_div$Compartment,levels(alpha_div$Compartment)[c(1,3,2)])
alpha_div$Spartina <- factor(alpha_div$Spartina,levels(alpha_div$Spartina)[c(2,1,3)])


fig_Shannon <- ggplot(alpha_div, aes (x = Compartment, y = Shannon, fill = Spartina)) +
  geom_boxplot() + theme_alpha + labs(y = "Shannon index", x = "") +
  scale_y_continuous(limits = c(3,7))


#Plot qPCR abundance

qpcr_data <- read.csv("qpcr_data.csv", header = TRUE)
qpcr_data$Compartment <- factor(qpcr_data$Compartment, levels = levels(qpcr_data$Compartment)[c(1,3,2)])

gg_color <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

cols = gg_color(3)
cols <- cols[c(1,3)]

fig_qPCR <- ggplot(qpcr_data, aes (x = Compartment, y = logcopies_g_mass, fill = Spartina)) +
  geom_boxplot() + theme_alpha + labs(y = expression(Log[copies~'16S'~rRNA~gene]~g^-1), x = "") +
  theme(axis.text.x = element_text(angle = 0), axis.title.x = element_text(size = 18)) +
  scale_fill_manual(values=cols)

#Plot Rank abundance

dsn <- transform_sample_counts(phy.tree.prev.filt, function(x) x / sum(x) * 100)
sample_data(dsn)$Compartment <- gsub("Bulk_sediment", "Bulk sediment", sample_data(dsn)$Compartment)
otu_table(dsn) = t(otu_table(dsn))

data <- list(abund = as.data.frame(otu_table(dsn)@.Data), 
             tax = data.frame(tax_table(dsn)@.Data, OTU = rownames(tax_table(dsn))), 
             sample = suppressWarnings(as.data.frame(as.matrix(sample_data(dsn)))))


abund <- data[["abund"]]
tax <- data[["tax"]]
tax$OTU <- paste(rep("ASV", length(tax$OTU)), 1:length(tax$OTU), sep = "_")
sample <- data[["sample"]]


abund3 <- cbind.data.frame(Display = tax[, "OTU"], 
                           abund) %>% melt(id.var = "Display", value.name = "Abundance", 
                                           variable.name = "Sample")

abund3 <- data.table(abund3)[, `:=`(Abundance, sum(Abundance)), 
                             by = list(Display, Sample)] %>% setkey(Display, Sample) %>% 
  unique() %>% as.data.frame()


grp <- data.frame(Sample = rownames(sample), Group = sample[, 
                                                            "Compartment"])

abund3$Group <- grp$Group[match(abund3$Sample, grp$Sample)]
abund5 <- abund3
abund5 <- abund5[order(abund5$Abundance, decreasing = T),]
abund5_ord <- abund5[order(abund5$Sample),]
abund5_ord$Rank <- rep(1:length(unique(abund5_ord$Display)), length(unique(abund5_ord$Sample)))

temp3 <- dplyr::group_by(abund5_ord, Group, Rank) %>% dplyr::summarise(Mean = mean(Abundance))
temp3 <- temp3[temp3$Mean > 0, ]

TotalCounts <- temp3[with(temp3, order(-Mean)), ] %>% 
  group_by(Group) %>% mutate(dummy = 1) %>% mutate(Cumsum = cumsum(Mean), 
                                                   Rank = cumsum(dummy)) %>% as.data.frame()

TotalCounts$Group <- as.factor(TotalCounts$Group) 
TotalCounts$Group <- factor(TotalCounts$Group, levels(TotalCounts$Group)[c(1,3,2)])


rank_abund_plot <- ggplot(data = TotalCounts, aes(x = Rank, y = Cumsum, 
                                                  color = Group)) + geom_point(size = 2) + geom_line(size = 1) + ylim(0, 100) + 
  xlab("Rank abundance") + ylab("Cummulative read abundance") + scale_x_log10() + 
  theme(legend.position = c(0.90,0.12), axis.title = element_text(size = 18), axis.text = element_text(size = 16),
        legend.title = element_blank(), legend.text = element_text(size = 16), panel.background = element_blank(),
        panel.grid = element_blank()) + labs(x = "ASV Rank", y = "Cumulative rel. abundance (%)")


#Plot NDMS

sample_data(phy.tree.prev.filt)$Compartment <- 
  as.factor(gsub("Bulk_sediment", "Bulk sediment", sample_data(phy.tree.prev.filt)$Compartment))

sample_data(phy.tree.prev.filt)$Compartment <-
  factor(sample_data(phy.tree.prev.filt)$Compartment,
         levels(sample_data(phy.tree.prev.filt)$Compartment)[c(1,3,2)])

sample_data(phy.tree.prev.filt)$Spartina <-
  factor(sample_data(phy.tree.prev.filt)$Spartina,
         levels(sample_data(phy.tree.prev.filt)$Spartina)[c(2,1,3)])

theme_ord <- theme(axis.title.x = element_text(size = 18),
                   axis.title.y = element_text(size = 18),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   axis.text = element_text(size = 16),
                   plot.title = element_blank(),
                   legend.position = "bottom",
                   legend.title = element_blank(),
                   legend.text = element_text(size = 16)) 

ord.nmds.bray <- ordinate(phy.tree.prev.filt, method="NMDS", distance="bray")

fig_nmds_bray_Spartina <- plot_ordination(phy.tree.prev.filt, ord.nmds.bray, color="Spartina",
                                          shape = "Spartina", title="Bray NMDS") + 
  geom_point(size =4) + coord_fixed(expand= TRUE) +
  theme_bw() + theme_ord

fig_nmds_bray_Compartment <- plot_ordination(phy.tree.prev.filt, ord.nmds.bray, color="Compartment",
                                             shape = "Compartment", title="Bray NMDS") + 
  geom_point(size =4) + coord_fixed(expand= TRUE) +
  theme_bw() + theme_ord

Fig3_final <- ggarrange(ggarrange(fig_Shannon, fig_qPCR, ncol = 2, labels = c("a", "b"), align = "hv",
                            font.label = list(size = 28)),
                  rank_abund_plot, nrow = 3, 
                  ggarrange(fig_nmds_bray_Compartment, fig_nmds_bray_Spartina, ncol = 2,
                            labels = c("d", "e"), font.label = list(size = 28)),
                  labels = c("", "c", ""), font.label = list(size = 28))

jpeg("../Paper_spartina_assembly/Figures/fig_3.jpeg", width = 3800, height = 4500, res = 300)
Fig3_final
dev.off()

tiff("../Paper_spartina_assembly/Figures/fig_3.tiff", width = 3800, height = 4500, res = 300)
Fig3_final
dev.off()


########PERMANOVA##########

distance.bray <- phyloseq::distance(phy.tree.prev.filt, method="bray")
sample_metadata <- data.frame(sample_data(phy.tree.prev.filt))

adonis(distance.bray ~ Compartment + Spartina + Location + Depth + Year, data = sample_metadata)

#Permanova bulk
phy.bulk <- subset_samples(phy.tree.prev.filt, Compartment == "Bulk sediment")

distance.bray.bulk <- phyloseq::distance(phy.bulk, method="bray")
sample_metadata_bulk <- data.frame(sample_data(phy.bulk))

adonis(distance.bray.bulk ~ Spartina + Location + Depth + Year, data = sample_metadata_bulk)

#Permanova rhizo
phy.rhizo <- subset_samples(phy.tree.prev.filt, Compartment == "Rhizosphere")

distance.bray.rhizo <- phyloseq::distance(phy.rhizo, method="bray")
sample_metadata_rhizo <- data.frame(sample_data(phy.rhizo))

adonis(distance.bray.rhizo ~ Spartina + Location + Depth + Year, data = sample_metadata_rhizo)

#Permanova endo
phy.endo <- subset_samples(phy.tree.prev.filt, Compartment == "Endosphere")

distance.bray.endo <- phyloseq::distance(phy.endo, method="bray")
sample_metadata_endo <- data.frame(sample_data(phy.endo))

adonis(distance.bray.endo ~ Spartina + Location + Depth + Year, data = sample_metadata_endo)





##################### RELATIVE ABUNDANCE ANALYSIS ######################
##S and Fe oxidizers, sulfate reducers, and nitrifiers relative abundance analysis

function_taxa <- read.csv("../taxonomy_to_function/taxonomy_function.csv", stringsAsFactors = FALSE)
Fe_Ox <- function_taxa$Fe_oxidizer
S_Ox <- function_taxa$Sulfur_oxidizer
S_Re <- function_taxa$Sulfur_reducers

theme_ra <- theme_bw() + theme(axis.title.x = element_blank(),
                               axis.title.y = element_text(size = 18),
                               axis.text = element_text(size = 16),
                               legend.position = "top",
                               legend.text = element_text(size = 16),
                               legend.title = element_blank(),
                               panel.grid.major = element_blank(),
                               panel.grid.minor = element_blank(),
                               strip.text.x = element_text(size = 14),
                               strip.background = element_rect(fill = "grey90"))


library(scales)
mysqrt_trans <- function() {
  trans_new("mysqrt", 
            transform = base::sqrt,
            inverse = function(x) ifelse(x<0, 0, x^2),
            domain = c(0, Inf))
}


phy.prev.filt.ra  = transform_sample_counts(phy.tree.prev.filt, function(x) x / sum(x))
Genus_ra <- tax_glom(phy.prev.filt.ra, taxrank="Genus", NArm=TRUE, bad_empty=c(NA, "", " ", "\t"))

sample_data(Genus_ra)$Compartment <- gsub("Bulk_sediment", "Bulk",  sample_data(Genus_ra)$Compartment)
sample_data(Genus_ra)$Compartment <- gsub("Rhizosphere", "Rhizo",  sample_data(Genus_ra)$Compartment)
sample_data(Genus_ra)$Compartment <- gsub("Endosphere", "Endo",  sample_data(Genus_ra)$Compartment)
sample_data(Genus_ra)$Compartment <- as.factor(sample_data(Genus_ra)$Compartment)

Genus_ra_melt  <- psmelt(Genus_ra)

####Nitrifiers per Genus

nitri_genus_phy <- subset_taxa(Genus_ra, grepl("Nitro", Genus))

#Select most abundant genus
nitri_genus_phy_filt <- phyloseq_filter_prevalence(nitri_genus_phy, prev.trh = 0, 
                                                   abund.trh = 0.10, abund.type = "total",
                                                   threshold_condition = "AND")

nitri_genus_melt <- psmelt(nitri_genus_phy_filt)
nitri_genus_melt$Genus <- paste(nitri_genus_melt$Order, nitri_genus_melt$Genus, sep = ", ")

nitri_genus_agg <- aggregate(nitri_genus_melt$Abundance, by=list(Seq_ID=nitri_genus_melt$Seq_ID, Spartina=nitri_genus_melt$Spartina,
                                                                 Compartment=nitri_genus_melt$Compartment, Location=nitri_genus_melt$Location,
                                                                 Depth=nitri_genus_melt$Depth, Genus=nitri_genus_melt$Genus), FUN=sum)

names(nitri_genus_agg)[length(nitri_genus_agg)] <- "rel_abund"

nitri_genus_agg$Compartment <- as.factor(nitri_genus_agg$Compartment)
nitri_genus_agg$Compartment <- factor(nitri_genus_agg$Compartment,
                                      levels(nitri_genus_agg$Compartment)[c(1,3,2)])

nitri_genus_agg$Spartina <- as.factor(nitri_genus_agg$Spartina)
nitri_genus_agg$Spartina <- factor(nitri_genus_agg$Spartina,
                                   levels(nitri_genus_agg$Spartina)[c(2,1,3)])


fig_nitri <- ggplot(nitri_genus_agg, aes (x = Compartment, y = rel_abund*100, fill = Spartina)) +
  geom_boxplot() + facet_wrap(~Genus, ncol = 1) + coord_flip() + theme_ra +
  scale_y_continuous(trans="mysqrt", limits=c(0, 4.5), breaks = c(0,0.4,1.7, 4)) + 
  theme(axis.text.x = element_text(angle = 0), axis.title.y = element_blank(),
        axis.title.x = element_text(size = 18), legend.position = "none") + labs(y = "Relative abundance (%)", x = "")


##Nitrifiers total

nitri_ra_sum <- aggregate(nitri_genus_melt$Abundance, by=list(nitri_genus_melt=nitri_genus_melt$Seq_ID, Spartina=nitri_genus_melt$Spartina,
                                                              Compartment=nitri_genus_melt$Compartment, Location=nitri_genus_melt$Location,
                                                              Depth=nitri_genus_melt$Depth), FUN=sum)

names(nitri_ra_sum)[6] <- "nitri_ra_sum"

nitri_ra_sum$nitri_ra_sum <- nitri_ra_sum$nitri_ra_sum*100


##Model and plotting

nitri_ra_sum$Compartment <- as.factor(nitri_ra_sum$Compartment)
nitri_ra_sum$Compartment <- factor(nitri_ra_sum$Compartment,
                                   levels(nitri_ra_sum$Compartment)[c(1,3,2)])

nitri_ra_sum$Spartina <- as.factor(nitri_ra_sum$Spartina)
nitri_ra_sum$Spartina <- factor(nitri_ra_sum$Spartina,
                                levels(nitri_ra_sum$Spartina)[c(2,1,3)])

g <- ggplot(data = nitri_ra_sum, aes(x = Compartment, y = nitri_ra_sum,
                                     fill = factor(Spartina))) 
fig_0 <- g + geom_boxplot() + 
  theme_ra + labs(y = "Nitrifiers (%)") + 
  scale_y_continuous(trans="mysqrt", limits=c(0, 5.5), breaks = c(0,0.3, 1.3, 3 ,5))



#######
Fe_Ox_genus_phy <- subset_taxa(Genus_ra, Genus %in% Fe_Ox)
Fe_Ox_genus_phy_filt <- phyloseq_filter_prevalence(Fe_Ox_genus_phy, prev.trh = 0, 
                                                   abund.trh = 0.10, abund.type = "total",
                                                   threshold_condition = "AND")

Fe_Ox_genus_melt <- psmelt(Fe_Ox_genus_phy_filt)
Fe_Ox_genus_melt$Genus <- paste(Fe_Ox_genus_melt$Order, Fe_Ox_genus_melt$Genus, sep = ", ")

Fe_Ox_genus_agg <- aggregate(Fe_Ox_genus_melt$Abundance, by=list(Seq_ID=Fe_Ox_genus_melt$Seq_ID, Spartina=Fe_Ox_genus_melt$Spartina,
                                                                 Compartment=Fe_Ox_genus_melt$Compartment, Location=Fe_Ox_genus_melt$Location,
                                                                 Depth=Fe_Ox_genus_melt$Depth, Genus=Fe_Ox_genus_melt$Genus), FUN=sum)

names(Fe_Ox_genus_agg)[length(Fe_Ox_genus_agg)] <- "rel_abund"

Fe_Ox_genus_agg$Compartment <- as.factor(Fe_Ox_genus_agg$Compartment)
Fe_Ox_genus_agg$Compartment <- factor(Fe_Ox_genus_agg$Compartment,
                                      levels(Fe_Ox_genus_agg$Compartment)[c(1,3,2)])

Fe_Ox_genus_agg$Spartina <- as.factor(Fe_Ox_genus_agg$Spartina)
Fe_Ox_genus_agg$Spartina <- factor(Fe_Ox_genus_agg$Spartina,
                                   levels(Fe_Ox_genus_agg$Spartina)[c(2,1,3)])


fig_iron <- ggplot(Fe_Ox_genus_agg, aes (x = Compartment, y = rel_abund*100, fill = Spartina)) +
  geom_boxplot() + facet_wrap(~Genus, ncol = 1) + coord_flip() + theme_ra +
  scale_y_continuous(trans="mysqrt", limits=c(0, 11), breaks = c(0,0.5,2.5, 6, 10.5)) + 
  theme(axis.text.x = element_text(angle = 0), axis.title.y = element_blank(),
        axis.title.x = element_text(size = 18), legend.position = "none") + labs(y = "Relative abundance (%)", x = "")


#FeOx
Fe_ox_ra <- Genus_ra_melt %>% filter(Genus %in% Fe_Ox)
Fe_ox_ra$Seq_ID <- as.factor(Fe_ox_ra$Seq_ID )

Fe_ox_ra_sum <- aggregate(Fe_ox_ra$Abundance, by=list(Seq_ID=Fe_ox_ra$Seq_ID, Spartina=Fe_ox_ra$Spartina,
                                                      Compartment=Fe_ox_ra$Compartment, Location=Fe_ox_ra$Location,
                                                      Depth=Fe_ox_ra$Depth), FUN=sum)

names(Fe_ox_ra_sum)[6] <- "Fe_ox_ra_sum"

Fe_ox_ra_sum$Fe_ox_ra_sum <- Fe_ox_ra_sum$Fe_ox_ra_sum*100


##Model and plotting

Fe_ox_ra_sum$Compartment <- as.factor(Fe_ox_ra_sum$Compartment)
Fe_ox_ra_sum$Compartment <- factor(Fe_ox_ra_sum$Compartment,
                                   levels(Fe_ox_ra_sum$Compartment)[c(1,3,2)])

Fe_ox_ra_sum$Spartina <- as.factor(Fe_ox_ra_sum$Spartina)
Fe_ox_ra_sum$Spartina <- factor(Fe_ox_ra_sum$Spartina,
                                levels(Fe_ox_ra_sum$Spartina)[c(2,1,3)])

g <- ggplot(data = Fe_ox_ra_sum, aes(x = Compartment, y = Fe_ox_ra_sum,
                                     fill = factor(Spartina))) 

fig_1 <- g + geom_boxplot() + 
  theme_ra + labs(y = "Fe oxidizers (%)") + 
  scale_y_continuous(trans="mysqrt", limits=c(0, 11), breaks = c(0,0.5,2.5,6,10.5))


#SOx

S_Ox_genus_phy <- subset_taxa(Genus_ra, Genus %in% S_Ox)
S_Ox_genus_phy_filt <- phyloseq_filter_prevalence(S_Ox_genus_phy, prev.trh = 0, 
                                                  abund.trh = 0.5, abund.type = "total",
                                                  threshold_condition = "AND")

S_Ox_genus_phy_filt <- subset_taxa(S_Ox_genus_phy, Genus == "endosymbionts" | Genus == "Candidatus_Thiodiazotropha" |
                                     Genus == "Desulfovibrio")

S_Ox_genus_melt <- psmelt(S_Ox_genus_phy_filt)
S_Ox_genus_melt$Genus <- paste(S_Ox_genus_melt$Order, S_Ox_genus_melt$Genus, sep = ", ")

S_Ox_genus_agg <- aggregate(S_Ox_genus_melt$Abundance, by=list(Seq_ID=S_Ox_genus_melt$Seq_ID, Spartina=S_Ox_genus_melt$Spartina,
                                                               Compartment=S_Ox_genus_melt$Compartment, Location=S_Ox_genus_melt$Location,
                                                               Depth=S_Ox_genus_melt$Depth, Genus=S_Ox_genus_melt$Genus), FUN=sum)

names(S_Ox_genus_agg)[length(S_Ox_genus_agg)] <- "rel_abund"

S_Ox_genus_agg$Compartment <- as.factor(S_Ox_genus_agg$Compartment)
S_Ox_genus_agg$Compartment <- factor(S_Ox_genus_agg$Compartment,
                                     levels(S_Ox_genus_agg$Compartment)[c(1,3,2)])

S_Ox_genus_agg$Spartina <- as.factor(S_Ox_genus_agg$Spartina)
S_Ox_genus_agg$Spartina <- factor(S_Ox_genus_agg$Spartina,
                                  levels(S_Ox_genus_agg$Spartina)[c(2,1,3)])

fig_S_ox <- ggplot(S_Ox_genus_agg, aes (x = Compartment, y = rel_abund*100, fill = Spartina)) +
  geom_boxplot() + facet_wrap(~Genus, ncol = 1) + coord_flip() + theme_ra +
  scale_y_continuous(trans="mysqrt", limits=c(0, 60), breaks = c(0, 2, 8, 20, 38, 55)) + 
  theme(axis.text.x = element_text(angle = 0), axis.title.y = element_blank(),
        axis.title.x = element_text(size = 18),  legend.position = "none") + labs(y = "Relative abundance (%)", x = "")


#all of them

S_ox_ra <- Genus_ra_melt %>% filter(Genus %in% S_Ox)
S_ox_ra$Seq_ID <- as.factor(S_ox_ra$Seq_ID )

S_ox_ra_sum <- aggregate(S_ox_ra$Abundance, by=list(Seq_ID=S_ox_ra$Seq_ID, Spartina=S_ox_ra$Spartina,
                                                    Compartment=S_ox_ra$Compartment, Location=S_ox_ra$Location,
                                                    Depth=S_ox_ra$Depth), FUN=sum)

names(S_ox_ra_sum)[6] <- "S_ox_ra_sum"

S_ox_ra_sum$S_ox_ra_sum <- S_ox_ra_sum$S_ox_ra_sum*100

##Model and plotting

#Calculating SE

S_ox_ra_sum$Compartment <- as.factor(S_ox_ra_sum$Compartment)
S_ox_ra_sum$Compartment <- factor(S_ox_ra_sum$Compartment,
                                  levels(S_ox_ra_sum$Compartment)[c(1,3,2)])

S_ox_ra_sum$Spartina <- as.factor(S_ox_ra_sum$Spartina)
S_ox_ra_sum$Spartina <- factor(S_ox_ra_sum$Spartina,
                               levels(S_ox_ra_sum$Spartina)[c(2,1,3)])

g <- ggplot(data = S_ox_ra_sum, aes(x = Compartment, y = S_ox_ra_sum, fill = Spartina))
fig_2 <- g + geom_boxplot() + 
  expand_limits(y = 0) + theme_ra + labs(y = "Sulfur oxidizers (%)")  + 
  scale_y_continuous(trans="mysqrt", limits=c(0, 65), breaks = c(0,2,8,20, 38, 60))



#Sulfate Reducers
S_Re_genus_phy <- subset_taxa(Genus_ra, Genus %in% S_Re)
S_Re_genus_phy_filt <- phyloseq_filter_prevalence(S_Re_genus_phy, prev.trh = 0, 
                                                  abund.trh = 1, abund.type = "total",
                                                  threshold_condition = "AND")

S_Re_genus_phy_filt <- subset_taxa(S_Re_genus_phy_filt, Genus == "Desulfatitalea" | Genus == "Desulfopila" |
                                     Genus == "Desulfosarcina")

S_Re_genus_melt <- psmelt(S_Re_genus_phy_filt)
S_Re_genus_melt$Genus <- gsub("Sva0081_sediment_group", "Sva0081", S_Re_genus_melt$Genus)
S_Re_genus_melt$Genus <- paste(S_Re_genus_melt$Order, S_Re_genus_melt$Genus, sep = ", ")


S_Re_genus_agg <- aggregate(S_Re_genus_melt$Abundance, by=list(Seq_ID=S_Re_genus_melt$Seq_ID, Spartina=S_Re_genus_melt$Spartina,
                                                               Compartment=S_Re_genus_melt$Compartment, Location=S_Re_genus_melt$Location,
                                                               Depth=S_Re_genus_melt$Depth, Genus=S_Re_genus_melt$Genus), FUN=sum)

names(S_Re_genus_agg)[length(S_Re_genus_agg)] <- "rel_abund"

S_Re_genus_agg$Compartment <- as.factor(S_Re_genus_agg$Compartment)
S_Re_genus_agg$Compartment <- factor(S_Re_genus_agg$Compartment,
                                     levels(S_Re_genus_agg$Compartment)[c(1,3,2)])

S_Re_genus_agg$Spartina <- as.factor(S_Re_genus_agg$Spartina)
S_Re_genus_agg$Spartina <- factor(S_Re_genus_agg$Spartina,
                                  levels(S_Re_genus_agg$Spartina)[c(2,1,3)])

fig_S_Re <- ggplot(S_Re_genus_agg, aes (x = Compartment, y = rel_abund*100, fill = Spartina)) +
  geom_boxplot() + facet_wrap(~Genus, ncol = 1) + coord_flip() + theme_ra +
  scale_y_continuous(trans="mysqrt", limits=c(0, 43), breaks = c(0, 2.5, 9, 22, 40)) + 
  theme(axis.text.x = element_text(angle = 0), axis.title.y = element_blank(),
        axis.title.x = element_text(size = 18),  legend.position = "none") + labs(y = "Relative abundance (%)", x = "")


S_re_ra <- Genus_ra_melt %>% filter(Genus %in% S_Re)
S_re_ra$Seq_ID <- as.factor(S_re_ra$Seq_ID )

S_re_ra_sum <- aggregate(S_re_ra$Abundance, by=list(Seq_ID=S_re_ra$Seq_ID, Spartina=S_re_ra$Spartina,
                                                    Compartment=S_re_ra$Compartment, Location=S_re_ra$Location,
                                                    Depth=S_re_ra$Depth), FUN=sum)

names(S_re_ra_sum)[6] <- "S_re_ra_sum"

S_re_ra_sum$S_re_ra_sum <- S_re_ra_sum$S_re_ra_sum*100

##Model and plotting

S_re_ra_sum$Compartment <- as.factor(S_re_ra_sum$Compartment)
S_re_ra_sum$Compartment <- factor(S_re_ra_sum$Compartment,
                                  levels(S_re_ra_sum$Compartment)[c(1,3,2)])

S_re_ra_sum$Spartina <- as.factor(S_re_ra_sum$Spartina)
S_re_ra_sum$Spartina <- factor(S_re_ra_sum$Spartina,
                               levels(S_re_ra_sum$Spartina)[c(2,1,3)])

g <- ggplot(data = S_re_ra_sum, aes(x = Compartment, y = S_re_ra_sum, fill = Spartina))
fig_3 <- g + geom_boxplot() + 
  expand_limits(y = 0) + theme_ra + labs(y = "Sulfur reducers (%)")  + 
  scale_y_continuous(trans="mysqrt", limits=c(0, 55), breaks = c(0,3,12,28, 50))


Fig4_final <- ggarrange(fig_0, fig_1, fig_2, fig_3, fig_nitri, fig_iron, fig_S_ox, fig_S_Re,  labels = c("a", "b", "c", "d", "", "", "", ""),
                          nrow = 2, ncol = 4, heights = c(1,2), align = "v", font.label = c(size = 28))

jpeg("Fig4.jpeg", res = 300, width = 6000, height = 4000)
Fig4_final
dev.off()


tiff("Fig4.tiff", res = 300, width = 6000, height = 4000)
Fig4_final
dev.off()



###########CORE MICROBIOME ANALYSIS - FIG 5 AND S8###############

#Subset prevalence in 10% increments. Will calculate what is the relative abundance
#and richness retained at each level -compared to non-subsetted dataset-.
#This analysis was used to justify a threshold of 60%

phy_endo <- subset_samples(phy.tree.prev.filt, Compartment == "Endosphere")

phy_endo_initial <- phyloseq_filter_prevalence(phy_endo, prev.trh = 0, abund.trh = 1,
                                               threshold_condition = "AND", abund.type = "total")
phy_endo_core_100 = phyloseq_filter_prevalence(phy_endo, prev.trh = 1, abund.trh = 1,
                                               threshold_condition = "AND", abund.type = "total")
phy_endo_core_90 = phyloseq_filter_prevalence(phy_endo, prev.trh = 0.9, abund.trh = 1,
                                              threshold_condition = "AND", abund.type = "total")
phy_endo_core_80 = phyloseq_filter_prevalence(phy_endo, prev.trh = 0.8, abund.trh = 1,
                                              threshold_condition = "AND", abund.type = "total")
phy_endo_core_70 = phyloseq_filter_prevalence(phy_endo, prev.trh = 0.7, abund.trh = 1,
                                              threshold_condition = "AND", abund.type = "total")
phy_endo_core_60 = phyloseq_filter_prevalence(phy_endo, prev.trh = 0.6, abund.trh = 1,
                                              threshold_condition = "AND", abund.type = "total")
phy_endo_core_50 = phyloseq_filter_prevalence(phy_endo, prev.trh = 0.5, abund.trh = 1,
                                              threshold_condition = "AND", abund.type = "total")
phy_endo_core_40 = phyloseq_filter_prevalence(phy_endo, prev.trh = 0.4, abund.trh = 1,
                                              threshold_condition = "AND", abund.type = "total")
phy_endo_core_30 = phyloseq_filter_prevalence(phy_endo, prev.trh = 0.3, abund.trh = 1,
                                              threshold_condition = "AND", abund.type = "total")
phy_endo_core_20 = phyloseq_filter_prevalence(phy_endo, prev.trh = 0.2, abund.trh = 1,
                                              threshold_condition = "AND", abund.type = "total")
phy_endo_core_10 = phyloseq_filter_prevalence(phy_endo, prev.trh = 0.1, abund.trh = 1,
                                              threshold_condition = "AND", abund.type = "total")


#Relative abundance analysis

phy_endo_core_0_ra <- sample_sums(phy_endo_initial)/sample_sums(phy_endo_initial)
phy_endo_core_10_ra <- sample_sums(phy_endo_core_10)/sample_sums(phy_endo_initial)
phy_endo_core_20_ra <- sample_sums(phy_endo_core_20)/sample_sums(phy_endo_initial)
phy_endo_core_30_ra <- sample_sums(phy_endo_core_30)/sample_sums(phy_endo_initial)
phy_endo_core_40_ra <- sample_sums(phy_endo_core_40)/sample_sums(phy_endo_initial)
phy_endo_core_50_ra <- sample_sums(phy_endo_core_50)/sample_sums(phy_endo_initial)
phy_endo_core_60_ra <- sample_sums(phy_endo_core_60)/sample_sums(phy_endo_initial)
phy_endo_core_70_ra <- sample_sums(phy_endo_core_70)/sample_sums(phy_endo_initial)
phy_endo_core_80_ra <- sample_sums(phy_endo_core_80)/sample_sums(phy_endo_initial)
phy_endo_core_90_ra <- sample_sums(phy_endo_core_90)/sample_sums(phy_endo_initial)
phy_endo_core_100_ra <- sample_sums(phy_endo_core_100)/sample_sums(phy_endo_initial)


ra_core <- data_frame(sample = names(phy_endo_core_30_ra),
                      '0' = phy_endo_core_0_ra,
                      '10' = phy_endo_core_10_ra,
                      '20' = phy_endo_core_20_ra,
                      '30' = phy_endo_core_30_ra,
                      '40' = phy_endo_core_40_ra,
                      '50' = phy_endo_core_50_ra,
                      '60' = phy_endo_core_60_ra,
                      '70' = phy_endo_core_70_ra,
                      '80' = phy_endo_core_80_ra,
                      '90' = phy_endo_core_90_ra,
                      '100' = phy_endo_core_100_ra)

metadata <- sample_data(phy_endo_core_100)

ra_core <- cbind(metadata[,c(3:10)], ra_core)

ra_core_long_endo <- gather(ra_core, percent, relative_abundance, '0':'100')
ra_core_long_endo$percent <- as.numeric(ra_core_long_endo$percent)

ra_core_long_endo$Spartina = factor(ra_core_long_endo$Spartina,levels(ra_core_long_endo$Spartina)[c(2,1,3)])

fig1 <- ggplot(ra_core_long_endo, aes(x = percent, y = 100*relative_abundance)) + geom_boxplot(aes(group = percent)) +
  geom_point(size = 2, aes(color = Spartina)) + theme_bw() + 
  theme(legend.position = "top", axis.text = element_text(size =18), axis.title = element_text(size = 20),
        legend.text = element_text(size = 16), legend.title = element_blank()) +
  labs(x = "Species prevalence cutoff (%)",
       y = "Abundance maintained (%)")


#Richness analysis ra

richness_initial <- estimate_richness(phy_endo_initial, measures = "Observed")

phy_endo_core_0_richness_ra <- estimate_richness(phy_endo_initial, measures = "Observed")/richness_initial
phy_endo_core_10_richness_ra <- estimate_richness(phy_endo_core_10, measures = "Observed")/richness_initial
phy_endo_core_20_richness_ra <- estimate_richness(phy_endo_core_20, measures = "Observed")/richness_initial
phy_endo_core_30_richness_ra <- estimate_richness(phy_endo_core_30, measures = "Observed")/richness_initial
phy_endo_core_40_richness_ra <- estimate_richness(phy_endo_core_40, measures = "Observed")/richness_initial
phy_endo_core_50_richness_ra <- estimate_richness(phy_endo_core_50, measures = "Observed")/richness_initial
phy_endo_core_60_richness_ra <- estimate_richness(phy_endo_core_60, measures = "Observed")/richness_initial
phy_endo_core_70_richness_ra <- estimate_richness(phy_endo_core_70, measures = "Observed")/richness_initial
phy_endo_core_80_richness_ra <- estimate_richness(phy_endo_core_80, measures = "Observed")/richness_initial
phy_endo_core_90_richness_ra <- estimate_richness(phy_endo_core_90, measures = "Observed")/richness_initial
phy_endo_core_100_richness_ra <- estimate_richness(phy_endo_core_100, measures = "Observed")/richness_initial

richness_ra_core <- data_frame(sample = row.names(phy_endo_core_0_richness_ra),
                               '0' = phy_endo_core_0_richness_ra$Observed,
                               '10' = phy_endo_core_10_richness_ra$Observed,
                               '20' = phy_endo_core_20_richness_ra$Observed,
                               '30' = phy_endo_core_30_richness_ra$Observed,
                               '40' = phy_endo_core_40_richness_ra$Observed,
                               '50' = phy_endo_core_50_richness_ra$Observed,
                               '60' = phy_endo_core_60_richness_ra$Observed,
                               '70' = phy_endo_core_70_richness_ra$Observed,
                               '80' = phy_endo_core_80_richness_ra$Observed,
                               '90' = phy_endo_core_90_richness_ra$Observed,
                               '100' = phy_endo_core_100_richness_ra$Observed)

richness_ra_core <- cbind(metadata[,c(3:10)], richness_ra_core)

richness_ra_core_long_endo <- gather(richness_ra_core, percent, relative_abundance, '0':'100')
richness_ra_core_long_endo$percent <- as.numeric(richness_ra_core_long_endo$percent)

richness_ra_core_long_endo$Spartina = factor(richness_ra_core_long_endo$Spartina,levels(richness_ra_core_long_endo$Spartina)[c(2,1,3)])

fig2 <- ggplot(richness_ra_core_long_endo, aes(x = percent, y = 100*relative_abundance)) + geom_boxplot(aes(group = percent)) +
  geom_point(size = 2, aes(color = Spartina)) + theme_bw() + 
  theme(legend.position = "top", axis.text = element_text(size =18), axis.title = element_text(size = 20),
        legend.text = element_text(size = 16), legend.title = element_blank()) +
  labs(x = "Species prevalence cutoff (%)",
       y = "Richness maintained (%)")


####RHIZOSPHERE

phy_rhizo <- subset_samples(phy.tree.prev.filt, Compartment == "Rhizosphere")

phy_rhizo_core_100 = phyloseq_filter_prevalence(phy_rhizo, prev.trh = 1, abund.trh = 1,
                                                threshold_condition = "AND", abund.type = "total")
phy_rhizo_core_90 = phyloseq_filter_prevalence(phy_rhizo, prev.trh = 0.9, abund.trh = 1,
                                               threshold_condition = "AND", abund.type = "total")
phy_rhizo_core_80 = phyloseq_filter_prevalence(phy_rhizo, prev.trh = 0.8, abund.trh = 1,
                                               threshold_condition = "AND", abund.type = "total")
phy_rhizo_core_70 = phyloseq_filter_prevalence(phy_rhizo, prev.trh = 0.7, abund.trh = 1,
                                               threshold_condition = "AND", abund.type = "total")
phy_rhizo_core_60 = phyloseq_filter_prevalence(phy_rhizo, prev.trh = 0.6, abund.trh = 1,
                                               threshold_condition = "AND", abund.type = "total")
phy_rhizo_core_50 = phyloseq_filter_prevalence(phy_rhizo, prev.trh = 0.5, abund.trh = 1,
                                               threshold_condition = "AND", abund.type = "total")
phy_rhizo_core_40 = phyloseq_filter_prevalence(phy_rhizo, prev.trh = 0.4, abund.trh = 1,
                                               threshold_condition = "AND", abund.type = "total")
phy_rhizo_core_30 = phyloseq_filter_prevalence(phy_rhizo, prev.trh = 0.3, abund.trh = 1,
                                               threshold_condition = "AND", abund.type = "total")
phy_rhizo_core_20 = phyloseq_filter_prevalence(phy_rhizo, prev.trh = 0.2, abund.trh = 1,
                                               threshold_condition = "AND", abund.type = "total")
phy_rhizo_core_10 = phyloseq_filter_prevalence(phy_rhizo, prev.trh = 0.1, abund.trh = 1,
                                               threshold_condition = "AND", abund.type = "total")

#Relative abundance analysis

phy_rhizo_core_0_ra <- sample_sums(phy_rhizo_initial)/sample_sums(phy_rhizo_initial)
phy_rhizo_core_10_ra <- sample_sums(phy_rhizo_core_10)/sample_sums(phy_rhizo_initial)
phy_rhizo_core_20_ra <- sample_sums(phy_rhizo_core_20)/sample_sums(phy_rhizo_initial)
phy_rhizo_core_30_ra <- sample_sums(phy_rhizo_core_30)/sample_sums(phy_rhizo_initial)
phy_rhizo_core_40_ra <- sample_sums(phy_rhizo_core_40)/sample_sums(phy_rhizo_initial)
phy_rhizo_core_50_ra <- sample_sums(phy_rhizo_core_50)/sample_sums(phy_rhizo_initial)
phy_rhizo_core_60_ra <- sample_sums(phy_rhizo_core_60)/sample_sums(phy_rhizo_initial)
phy_rhizo_core_70_ra <- sample_sums(phy_rhizo_core_70)/sample_sums(phy_rhizo_initial)
phy_rhizo_core_80_ra <- sample_sums(phy_rhizo_core_80)/sample_sums(phy_rhizo_initial)
phy_rhizo_core_90_ra <- sample_sums(phy_rhizo_core_90)/sample_sums(phy_rhizo_initial)
phy_rhizo_core_100_ra <- sample_sums(phy_rhizo_core_100)/sample_sums(phy_rhizo_initial)


ra_core <- data_frame(sample = names(phy_rhizo_core_30_ra),
                      '0' = phy_rhizo_core_0_ra,
                      '10' = phy_rhizo_core_10_ra,
                      '20' = phy_rhizo_core_20_ra,
                      '30' = phy_rhizo_core_30_ra,
                      '40' = phy_rhizo_core_40_ra,
                      '50' = phy_rhizo_core_50_ra,
                      '60' = phy_rhizo_core_60_ra,
                      '70' = phy_rhizo_core_70_ra,
                      '80' = phy_rhizo_core_80_ra,
                      '90' = phy_rhizo_core_90_ra,
                      '100' = phy_rhizo_core_100_ra)

metadata <- sample_data(phy_rhizo_core_100)

ra_core <- cbind(metadata[,c(3:10)], ra_core)

ra_core_long_rhizo <- gather(ra_core, percent, relative_abundance, '0':'100')
ra_core_long_rhizo$percent <- as.numeric(ra_core_long_rhizo$percent)

ra_core_long_rhizo$Spartina = factor(ra_core_long_rhizo$Spartina,levels(ra_core_long_rhizo$Spartina)[c(2,1,3)])

fig3 <- ggplot(ra_core_long_rhizo, aes(x = percent, y = 100*relative_abundance)) + geom_boxplot(aes(group = percent)) +
  geom_point(size = 2, aes(color = Spartina)) + theme_bw() + 
  theme(legend.position = "top", axis.text = element_text(size =18), axis.title = element_text(size = 20),
        legend.text = element_text(size = 16), legend.title = element_blank()) +
  labs(x = "Species prevalence cutoff (%)",
       y = "Abundance maintained (%)")


#Richness analysis ra

richness_initial <- estimate_richness(phy_rhizo_initial, measures = "Observed")

phy_rhizo_core_0_richness_ra <- estimate_richness(phy_rhizo_initial, measures = "Observed")/richness_initial
phy_rhizo_core_10_richness_ra <- estimate_richness(phy_rhizo_core_10, measures = "Observed")/richness_initial
phy_rhizo_core_20_richness_ra <- estimate_richness(phy_rhizo_core_20, measures = "Observed")/richness_initial
phy_rhizo_core_30_richness_ra <- estimate_richness(phy_rhizo_core_30, measures = "Observed")/richness_initial
phy_rhizo_core_40_richness_ra <- estimate_richness(phy_rhizo_core_40, measures = "Observed")/richness_initial
phy_rhizo_core_50_richness_ra <- estimate_richness(phy_rhizo_core_50, measures = "Observed")/richness_initial
phy_rhizo_core_60_richness_ra <- estimate_richness(phy_rhizo_core_60, measures = "Observed")/richness_initial
phy_rhizo_core_70_richness_ra <- estimate_richness(phy_rhizo_core_70, measures = "Observed")/richness_initial
phy_rhizo_core_80_richness_ra <- estimate_richness(phy_rhizo_core_80, measures = "Observed")/richness_initial
phy_rhizo_core_90_richness_ra <- estimate_richness(phy_rhizo_core_90, measures = "Observed")/richness_initial
phy_rhizo_core_100_richness_ra <- estimate_richness(phy_rhizo_core_100, measures = "Observed")/richness_initial

richness_ra_core <- data_frame(sample = row.names(phy_rhizo_core_0_richness_ra),
                               '0' = phy_rhizo_core_0_richness_ra$Observed,
                               '10' = phy_rhizo_core_10_richness_ra$Observed,
                               '20' = phy_rhizo_core_20_richness_ra$Observed,
                               '30' = phy_rhizo_core_30_richness_ra$Observed,
                               '40' = phy_rhizo_core_40_richness_ra$Observed,
                               '50' = phy_rhizo_core_50_richness_ra$Observed,
                               '60' = phy_rhizo_core_60_richness_ra$Observed,
                               '70' = phy_rhizo_core_70_richness_ra$Observed,
                               '80' = phy_rhizo_core_80_richness_ra$Observed,
                               '90' = phy_rhizo_core_90_richness_ra$Observed,
                               '100' = phy_rhizo_core_100_richness_ra$Observed)

richness_ra_core <- cbind(metadata[,c(3:10)], richness_ra_core)

richness_ra_core_long_rhizo <- gather(richness_ra_core, percent, relative_abundance, '0':'100')
richness_ra_core_long_rhizo$percent <- as.numeric(richness_ra_core_long_rhizo$percent)

richness_ra_core_long_rhizo$Spartina = factor(richness_ra_core_long_rhizo$Spartina,levels(richness_ra_core_long_rhizo$Spartina)[c(2,1,3)])

fig4 <- ggplot(richness_ra_core_long_rhizo, aes(x = percent, y = 100*relative_abundance)) + geom_boxplot(aes(group = percent)) +
  geom_point(size = 2, aes(color = Spartina)) + theme_bw() + 
  theme(legend.position = "top", axis.text = element_text(size =18), axis.title = element_text(size = 20),
        legend.text = element_text(size = 16), legend.title = element_blank()) +
  labs(x = "Species prevalence cutoff (%)",
       y = "Richness maintained (%)")


FigS6 <- ggarrange(fig3, fig1, fig4, fig2,
                    ncol = 2, nrow = 2, label.x = 0.12,
                    font.label = list(size = 30))

jpeg("figS6.jpeg", res = 300, width = 5000, height = 3000)
FigS6
dev.off()



########STUDY THE TAXONOMY AND PHYLOGENY OF THE CORE MICROBIOME#########


########PHYLOGENETIC TREE###########
library(RColorBrewer)
n <- 24
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
col_vector_tree <- col_vector[c(1:11,13,14,17,22,27,56,53,60)]
pie(rep(1,74), col=(col_vector))

##PLOT TREE ENDOSPHERE
endo_ASVs <- colnames(otu_table(phy_endo_core_60))
rhizo_ASVs <- colnames(otu_table(phy_rhizo_core_60))

root_core_ASVs <- unique(c(endo_ASVs, rhizo_ASVs))

phy.tree.prev.filt_ra = transform_sample_counts(phy.tree.prev.filt, function(x)  100*x / sum(x) )

my_subset <- subset(t(otu_table(phy.tree.prev.filt_ra, taxa_are_rows = FALSE)),
                    rownames(t(otu_table(phy.tree.prev.filt_ra))) %in% root_core_ASVs)

new_physeq <- phyloseq(my_subset, tax_table(phy.tree.prev.filt_ra),
                       sample_data(phy.tree.prev.filt_ra), phy_tree(phy.tree.prev.filt_ra))

endo_rhizo_core <- subset_samples(new_physeq,
                                  Compartment == "Endosphere" | Compartment == "Rhizosphere")

endo_rhizo_core_plot = merge_samples(endo_rhizo_core, "Compartment")

endo_rhizo_core_plot = transform_sample_counts(endo_rhizo_core_plot, function(x)  x / 68 ) #To calculate average from 68 rhizosphere samples

sample_data(endo_rhizo_core_plot)$Compartment <- gsub(pattern = "1", replacement = "Endosphere",
                                                      sample_data(endo_rhizo_core_plot)$Compartment)
sample_data(endo_rhizo_core_plot)$Compartment <- gsub(pattern = "2", replacement = "Rhizosphere",
                                                      sample_data(endo_rhizo_core_plot)$Compartment)

otu_table_plot <- otu_table(endo_rhizo_core_plot)

log_endo <- (!(colnames(otu_table_plot) %in% endo_ASVs)) 
otu_table(endo_rhizo_core_plot)[1, log_endo] <- 0

log_rhizo <- (!(colnames(otu_table_plot) %in% rhizo_ASVs)) 
otu_table(endo_rhizo_core_plot)[2, log_rhizo] <- 0


write.csv(tax_table(endo_rhizo_core_plot), file = "Core_microbiome/tax_table_plot.csv")

#Changed ASVs taxonomy names manually to name them under their closest taxonomic identity (Ejm. Unk_Anaerolineaceae)

tax_table_plot <- read.csv(file = "Core_microbiome/tax_table_plot.csv", header = TRUE, row.names = 1)

tax_table(endo_rhizo_core_plot) <- as.matrix(tax_table_plot)

Fig5a <- ggtree(endo_rhizo_core_plot, layout = "circular") + 
  geom_text2(aes(subset=!isTip, label=label), hjust=-0.03, size=2) +
  geom_tiplab(aes(label=Genus), align = FALSE, offset = 0.07, size = 3.5) + 
  geom_point(aes(x=x+hjust, color=Class, size=Abundance, shape = Compartment), alpha = 0.75) +
  scale_color_manual(values = col_vector_tree) + xlim(0,1.7) + geom_treescale(width = 0.2, x = 0.4) +
  theme(legend.position="bottom") + guides(shape = guide_legend(override.aes = list(size = 3)), 
                                           color = guide_legend(override.aes = list(size = 5))) +
  labs(size="Relative abundance (%)")

jpeg(filename = "Fig5a.jpeg", width = 5000, height = 3000, res = 300)
Fig5a
dev.off()

tiff(filename = "Fig5a.tiff", width = 5000, height = 3000, res = 300)
Fig5a
dev.off()



######RELATIVE ABUNDANCE GROUPING ASVs PER GENUS

sample_data(endo_rhizo_core)$Treatment <- 
  paste(sample_data(endo_rhizo_core)$Spartina, sample_data(endo_rhizo_core)$Compartment, sep = "-")

endo_rhizo_core_barplot <- merge_samples(x = endo_rhizo_core, group = "Treatment", fun = mean)

otu_table(endo_rhizo_core_barplot)[c(1,2,5,6),] <- otu_table(endo_rhizo_core_barplot)[c(1,2,5,6),]/22
otu_table(endo_rhizo_core_barplot)[c(3,4),] <- otu_table(endo_rhizo_core_barplot)[c(3,4),]/24

otu_table(endo_rhizo_core_barplot)[c(1,3,5), log_endo] <- 0
otu_table(endo_rhizo_core_barplot)[c(2,4,6), log_rhizo] <- 0

#GLOM and MELT
endo_rhizo_core_barplot_genus <- tax_glom(endo_rhizo_core_barplot, taxrank="Genus", NArm=FALSE, bad_empty=c(NA, "", " ", "\t"))
endo_rhizo_core_barplot_genus_melt  <- psmelt(endo_rhizo_core_barplot_genus)

for(i in 1:length(endo_rhizo_core_barplot_genus_melt$Genus)){
  if(is.na(endo_rhizo_core_barplot_genus_melt$Genus[i]) == TRUE){
    endo_rhizo_core_barplot_genus_melt$Order[i] <- NA
  }
}

endo_rhizo_core_barplot_genus_melt$Genus <- 
  paste(endo_rhizo_core_barplot_genus_melt$Order, endo_rhizo_core_barplot_genus_melt$Genus, sep = ", ")

endo_rhizo_core_barplot_genus_melt$Genus <- gsub("Sva0081_sediment_group", "Sva0081",
                                                 endo_rhizo_core_barplot_genus_melt$Genus)

endo_rhizo_core_barplot_genus_melt$Genus <- gsub("NA, NA", "Unknown Genus",
                                                 endo_rhizo_core_barplot_genus_melt$Genus)

Genus_ra_rhizo_core_sum <- aggregate(endo_rhizo_core_barplot_genus_melt$Abundance, by=list(Seq_ID=endo_rhizo_core_barplot_genus_melt$Seq_ID, Spartina=endo_rhizo_core_barplot_genus_melt$Spartina,
                                                                                           Compartment=endo_rhizo_core_barplot_genus_melt$Compartment, Location=endo_rhizo_core_barplot_genus_melt$Location,
                                                                                           Depth=endo_rhizo_core_barplot_genus_melt$Depth, Genus = endo_rhizo_core_barplot_genus_melt$Genus), FUN=sum)

names(Genus_ra_rhizo_core_sum)[7] <- "relative_abundance"

Genus_ra_rhizo_core_plot <- summarySE(data = Genus_ra_rhizo_core_sum, measurevar = "relative_abundance", na.rm = TRUE,
                                      groupvars = c("Spartina", "Compartment", "Genus"))

library(RColorBrewer)
n <- 60
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

Genus_ra_rhizo_core_plot$Compartment <- gsub("1", "Endosphere", Genus_ra_rhizo_core_plot$Compartment)
Genus_ra_rhizo_core_plot$Compartment <- gsub("2", "Rhizosphere", Genus_ra_rhizo_core_plot$Compartment)

Genus_ra_rhizo_core_plot$Spartina <- gsub("1", "Medium", Genus_ra_rhizo_core_plot$Spartina)
Genus_ra_rhizo_core_plot$Spartina <- gsub("2", "Short", Genus_ra_rhizo_core_plot$Spartina)
Genus_ra_rhizo_core_plot$Spartina <- gsub("3", "Tall", Genus_ra_rhizo_core_plot$Spartina)

Genus_ra_rhizo_core_plot$Spartina <- as.factor(Genus_ra_rhizo_core_plot$Spartina)

Genus_ra_rhizo_core_plot$Spartina = factor(Genus_ra_rhizo_core_plot$Spartina,levels(Genus_ra_rhizo_core_plot$Spartina)[c(2,1,3)])

Fig5b <- ggplot(Genus_ra_rhizo_core_plot, aes(x = Spartina, y = relative_abundance, fill = Genus)) +
  geom_bar(stat = "identity")  + scale_fill_manual(values = col_vector) + facet_grid(~Compartment) +
  theme(axis.text = element_text(size = 14), axis.title = element_text(size = 18),
        legend.position = "bottom", legend.text = element_text(size = 10),
        legend.title = element_blank(), panel.background = element_blank(),
        panel.grid = element_blank(), strip.text = element_text(size = 20)) + 
  guides(fill=guide_legend(ncol=5)) + labs(y = "Relative abundance (%)", x = "")

jpeg("Fig5b.jpeg", width = 4250, height = 2250, res = 300)
Fig5b
dev.off()

tiff("Fig5b.tiff", width = 4250, height = 2250, res = 300)
Fig5b
dev.off()





###############PHYLOGENETIC ALFA AND BETA DIVERSITY, NTI AND BETA-NTI INDEXES#################

######CALCULATE NTI########

##CALCULATE EACH SAMPLING LOCATION AND YEAR INDEPENDENTLY

##SAPELO ISLAND 2018

## read phyloseq object in OTU table

phy_sapelo_2018 <- subset_samples(phy.tree.prev.filt, Location == "Sapelo" &
                                    Year == "2018")

#Filter out taxa with less than 5% prevalence
phy_sapelo_2018_prevfilt = phyloseq_filter_prevalence(phy_sapelo_2018, prev.trh = 0.05, abund.trh = 100,
                                                      threshold_condition = "OR", abund.type = "total")

phy_sapelo_2018_prevfilt_ra <- transform_sample_counts(phy_sapelo_2018_prevfilt, function(x) x / sum(x) )

taxa_names(phy_sapelo_2018_prevfilt_ra) <- paste("ASV", seq(ntaxa(phy_sapelo_2018_prevfilt_ra)), sep = "_")


otu = otu_table(phy_sapelo_2018_prevfilt_ra);
dim(otu); # this gives the dimensions
otu[1:5,1:5]; # this gives a look at the first 5 rows and columns

## read in the phylogeny

phylo = phy_tree(phy_sapelo_2018_prevfilt_ra);
phylo; # a summary of the phylogeny
## make sure the names on the phylogeny are ordered the same as the names in otu table

match.phylo.otu = match.phylo.data(phylo, t(otu));
str(match.phylo.otu);

## calculate NTI
asv_data <- as.matrix(as.data.frame(t(match.phylo.otu$data)))
NTI_Sapelo2018 <- ses.mntd(asv_data, cophenetic(match.phylo.otu$phy),
                         null.model="taxa.labels", abundance.weighted = T)

NTI_Sapelo2018$Seq_ID <- rownames(NTI_Sapelo2018)
phy_sapelo_2018_meta <- data.frame(sample_data(phy_sapelo_2019_prevfilt_ra))
NTI_Sapelo2018 <- merge(NTI_Sapelo2018, phy_sapelo_2018_meta)
NTI_Sapelo2018$NTI <- -NTI_Sapelo2018$mntd.obs.z


####SAPELO 2019

phy_sapelo_2019 <- subset_samples(phy.tree.prev.filt, Location == "Sapelo" &
                                    Year == "2019")

#Filter out taxa with less than 5% prevalence
phy_sapelo_2019_prevfilt = phyloseq_filter_prevalence(phy_sapelo_2019, prev.trh = 0.05, abund.trh = 100,
                                                      threshold_condition = "OR", abund.type = "total")

phy_sapelo_2019_prevfilt_ra <- transform_sample_counts(phy_sapelo_2019_prevfilt,
                                                       function(x) x / sum(x) )

taxa_names(phy_sapelo_2019_prevfilt_ra) <- paste("ASV", seq(ntaxa(phy_sapelo_2019_prevfilt_ra)), sep = "_")


otu = otu_table(phy_sapelo_2019_prevfilt_ra);
dim(otu); # this gives the dimensions
otu[1:5,1:5]; # this gives a look at the first 5 rows and columns

## read in the phylogeny

phylo = phy_tree(phy_sapelo_2019_prevfilt_ra);
phylo; # a summary of the phylogeny
## make sure the names on the phylogeny are ordered the same as the names in otu table

match.phylo.otu = match.phylo.data(phylo, t(otu));
str(match.phylo.otu);

## calculate NTI
asv_data <- as.matrix(as.data.frame(t(match.phylo.otu$data)))
NTI_Sapelo2019 <- ses.mntd(asv_data, cophenetic(match.phylo.otu$phy),
                         null.model="taxa.labels", abundance.weighted = T)

NTI_Sapelo2019$Seq_ID <- rownames(NTI_Sapelo2019)
phy_sapelo_2019_meta <- data.frame(sample_data(phy_sapelo_2019_prevfilt_ra))
NTI_Sapelo2019 <- merge(NTI_Sapelo2019, phy_sapelo_2019_meta)
NTI_Sapelo2019$NTI <- -NTI_Sapelo2019$mntd.obs.z

####SKIDAWAY 2019

phy_skidaway_2019 <- subset_samples(phy.tree.prev.filt, Location == "skidaway" &
                                    Year == "2019")

#Filter out taxa with less than 5% prevalence
phy_skidaway_2019_prevfilt = phyloseq_filter_prevalence(phy_skidaway_2019, prev.trh = 0.05, abund.trh = 100,
                                                      threshold_condition = "OR", abund.type = "total")

phy_skidaway_2019_prevfilt_ra <- transform_sample_counts(phy_skidaway_2019_prevfilt,
                                                       function(x) x / sum(x) )

taxa_names(phy_skidaway_2019_prevfilt_ra) <- paste("ASV", seq(ntaxa(phy_skidaway_2019_prevfilt_ra)), sep = "_")


otu = otu_table(phy_skidaway_2019_prevfilt_ra);
dim(otu); # this gives the dimensions
otu[1:5,1:5]; # this gives a look at the first 5 rows and columns

## read in the phylogeny

phylo = phy_tree(phy_skidaway_2019_prevfilt_ra);
phylo; # a summary of the phylogeny
## make sure the names on the phylogeny are ordered the same as the names in otu table

match.phylo.otu = match.phylo.data(phylo, t(otu));
str(match.phylo.otu);

## calculate NTI
asv_data <- as.matrix(as.data.frame(t(match.phylo.otu$data)))
NTI_Skidaway2019 <- ses.mntd(asv_data, cophenetic(match.phylo.otu$phy),
                         null.model="taxa.labels", abundance.weighted = T)

NTI_Skidaway2019$Seq_ID <- rownames(NTI_Skidaway2019)
phy_skidaway_2019_meta <- data.frame(sample_data(phy_skidaway_2019_prevfilt_ra))
NTI_Skidaway2019 <- merge(NTI_Skidaway2019, phy_skidaway_2019_meta)
NTI_Skidaway2019$NTI <- -NTI_Skidaway2019$mntd.obs.z


NTI_Sapelo2018$Location <- "Sapelo2018"
NTI_Sapelo2019$Location <- "Sapelo2019"
NTI_Skidaway2019$Location <- "Skidaway2019"

all_NTI <- rbind(NTI_Sapelo2018, NTI_Sapelo2019, NTI_Skidaway2019)

all_NTI$Compartment <- as.factor(all_NTI$Compartment)
all_NTI$Compartment <- factor(all_NTI$Compartment,
                              levels(all_NTI$Compartment)[c(1,3,2)])

all_NTI$Spartina <- as.factor(all_NTI$Spartina)
all_NTI$Spartina <- factor(all_NTI$Spartina,
                           levels(all_NTI$Spartina)[c(2,1,3)])

theme_NTI <- theme_bw() + theme(axis.title.x = element_blank(),
                                axis.title.y = element_text(size = 24),
                                axis.text = element_text(size = 20),
                                panel.grid.major = element_blank(),
                                panel.grid.minor = element_blank(),
                                strip.text.x = element_text(size = 20),
                                strip.background = element_rect(fill = "grey90"))

#Histogram
fig_NTI_all_hist <- ggplot(all_NTI, aes(x = NTI)) + 
  geom_histogram(aes(y = ..density..), bins = 18) + facet_wrap(~Compartment) +
  geom_vline(aes(xintercept = 2), linetype = "dashed") + 
  geom_density() + theme_NTI + labs(y = "Density", x = "NTI")


########CALCULATE BETA-NTI########

##Sapelo 2018

## read phylogeny

phylo = phy_tree(phy_sapelo_2018_prevfilt_ra);
phylo; # a summary of the phylogeny
## make sure the names on the phylogeny are ordered the same as the names in otu table

match.phylo.otu = match.phylo.data(phylo, t(otu));
str(match.phylo.otu);

## calculate empirical betaMNTD

beta.mntd.weighted = as.matrix(comdistnt(t(match.phylo.otu$data),
                                         cophenetic(match.phylo.otu$phy),abundance.weighted=T));
dim(beta.mntd.weighted);
beta.mntd.weighted[1:5,1:5];
write.csv(beta.mntd.weighted,'bMNTD/betaMNTD_weighted.csv',quote=F);

identical(colnames(match.phylo.otu$data),colnames(beta.mntd.weighted)); # just a check, should be TRUE
identical(colnames(match.phylo.otu$data),rownames(beta.mntd.weighted)); # just a check, should be TRUE

# calculate randomized betaMNTD

beta.reps = 999; # number of randomizations

rand.weighted.bMNTD.comp = array(c(-999),dim=c(ncol(match.phylo.otu$data),
                                               ncol(match.phylo.otu$data),beta.reps));
dim(rand.weighted.bMNTD.comp);

for (rep in 1:beta.reps) {
  
  rand.weighted.bMNTD.comp[,,rep] = as.matrix(comdistnt(t(match.phylo.otu$data),
                                                        taxaShuffle(cophenetic(match.phylo.otu$phy)),
                                                        abundance.weighted=T,
                                                        exclude.conspecifics = F));
  
  print(c(date(),rep));
  
}

bNTI = matrix(c(NA),nrow=ncol(match.phylo.otu$data),
                       ncol=ncol(match.phylo.otu$data));
dim(bNTI);


for (columns in 1:(ncol(match.phylo.otu$data)-1)) {
  for (rows in (columns+1):ncol(match.phylo.otu$data)) {
    
    rand.vals = rand.weighted.bMNTD.comp[rows,columns,];
    bNTI[rows,columns] = (beta.mntd.weighted[rows,columns] - mean(rand.vals)) / sd(rand.vals);
    rm("rand.vals");
    
  };
};

rownames(bNTI) = colnames(match.phylo.otu$data)
colnames(bNTI) = colnames(match.phylo.otu$data)

bNTI <- as.matrix(bNTI)
bNTI_long <-  melt(bNTI, varnames = c("Sample1", "Sample2"), value.name = "bNTI")
bNTI_long <- bNTI_long[!is.na(bNTI_long$bNTI),]

#Subset same compartment, compare across Spartina phenotype
#Bulk sediment
bNTI_long_Bulk <- bNTI_long[grepl("Bulk", bNTI_long$Sample1) & grepl("Bulk", bNTI_long$Sample2),]

bNTI_long_Bulk_tall <- bNTI_long_Bulk[grepl("P1", bNTI_long_Bulk$Sample1) & grepl("P1", bNTI_long_Bulk$Sample2),]
bNTI_long_Bulk_tall$Spartina <- "Tall"

bNTI_long_Bulk_medium <- bNTI_long_Bulk[grepl("P2", bNTI_long_Bulk$Sample1) & grepl("P2", bNTI_long_Bulk$Sample2),]
bNTI_long_Bulk_medium$Spartina <- "Medium"

bNTI_long_Bulk_short <- bNTI_long_Bulk[(grepl("P3", bNTI_long_Bulk$Sample1) | grepl("P4", bNTI_long_Bulk$Sample1) ) &
                                         (grepl("P3", bNTI_long_Bulk$Sample2) | grepl("P4", bNTI_long_Bulk$Sample2)),]
bNTI_long_Bulk_short$Spartina <- "Short"

Bulk_Spartina_bNTI <- rbind(bNTI_long_Bulk_tall[,c(3,4)], 
                            bNTI_long_Bulk_medium[,c(3,4)], bNTI_long_Bulk_short[,c(3,4)])
Bulk_Spartina_bNTI$Compartment <- "Bulk_sediment"

#Rhizosphere
bNTI_long_Rhizosphere <- bNTI_long[grepl("Rhizosphere", bNTI_long$Sample1) & grepl("Rhizosphere", bNTI_long$Sample2),]

bNTI_long_Rhizosphere_tall <- bNTI_long_Rhizosphere[grepl("P1", bNTI_long_Rhizosphere$Sample1) & grepl("P1", bNTI_long_Rhizosphere$Sample2),]
bNTI_long_Rhizosphere_tall$Spartina <- "Tall"

bNTI_long_Rhizosphere_medium <- bNTI_long_Rhizosphere[grepl("P2", bNTI_long_Rhizosphere$Sample1) & grepl("P2", bNTI_long_Rhizosphere$Sample2),]
bNTI_long_Rhizosphere_medium$Spartina <- "Medium"

bNTI_long_Rhizosphere_short <- bNTI_long_Rhizosphere[(grepl("P3", bNTI_long_Rhizosphere$Sample1) | grepl("P4", bNTI_long_Rhizosphere$Sample1) ) &
                                                       (grepl("P3", bNTI_long_Rhizosphere$Sample2) | grepl("P4", bNTI_long_Rhizosphere$Sample2)),]
bNTI_long_Rhizosphere_short$Spartina <- "Short"

Rhizosphere_Spartina_bNTI <- rbind(bNTI_long_Rhizosphere_tall[,c(3,4)], 
                                   bNTI_long_Rhizosphere_medium[,c(3,4)], bNTI_long_Rhizosphere_short[,c(3,4)])
Rhizosphere_Spartina_bNTI$Compartment <- "Rhizosphere"

#Endosphere

bNTI_long_Endosphere <- bNTI_long[grepl("Endosphere", bNTI_long$Sample1) & grepl("Endosphere", bNTI_long$Sample2),]

bNTI_long_Endosphere_tall <- bNTI_long_Endosphere[grepl("P1", bNTI_long_Endosphere$Sample1) & grepl("P1", bNTI_long_Endosphere$Sample2),]
bNTI_long_Endosphere_tall$Spartina <- "Tall"

bNTI_long_Endosphere_medium <- bNTI_long_Endosphere[grepl("P2", bNTI_long_Endosphere$Sample1) & grepl("P2", bNTI_long_Endosphere$Sample2),]
bNTI_long_Endosphere_medium$Spartina <- "Medium"

bNTI_long_Endosphere_short <- bNTI_long_Endosphere[(grepl("P3", bNTI_long_Endosphere$Sample1) | grepl("P4", bNTI_long_Endosphere$Sample1) ) &
                                                     (grepl("P3", bNTI_long_Endosphere$Sample2) | grepl("P4", bNTI_long_Endosphere$Sample2)),]
bNTI_long_Endosphere_short$Spartina <- "Short"

Endosphere_Spartina_bNTI <- rbind(bNTI_long_Endosphere_tall[,c(3,4)], 
                                  bNTI_long_Endosphere_medium[,c(3,4)], bNTI_long_Endosphere_short[,c(3,4)])
Endosphere_Spartina_bNTI$Compartment <- "Endosphere"

Compartment_same_Spartina_bNTI_Sapelo18 <- rbind(Bulk_Spartina_bNTI,Rhizosphere_Spartina_bNTI,Endosphere_Spartina_bNTI)

Compartment_same_Spartina_bNTI_Sapelo18$Location <- "Sapelo 2018"




#######Sapelo 2019

## read phylogeny

phylo = phy_tree(phy_sapelo_2019_prevfilt_ra);
phylo; # a summary of the phylogeny
## make sure the names on the phylogeny are ordered the same as the names in otu table

match.phylo.otu = match.phylo.data(phylo, t(otu));
str(match.phylo.otu);

## calculate empirical betaMNTD

beta.mntd.weighted = as.matrix(comdistnt(t(match.phylo.otu$data),
                                         cophenetic(match.phylo.otu$phy),abundance.weighted=T));
dim(beta.mntd.weighted);
beta.mntd.weighted[1:5,1:5];
write.csv(beta.mntd.weighted,'bMNTD/betaMNTD_weighted.csv',quote=F);

identical(colnames(match.phylo.otu$data),colnames(beta.mntd.weighted)); # just a check, should be TRUE
identical(colnames(match.phylo.otu$data),rownames(beta.mntd.weighted)); # just a check, should be TRUE

# calculate randomized betaMNTD

beta.reps = 999; # number of randomizations

rand.weighted.bMNTD.comp = array(c(-999),dim=c(ncol(match.phylo.otu$data),
                                               ncol(match.phylo.otu$data),beta.reps));
dim(rand.weighted.bMNTD.comp);

for (rep in 1:beta.reps) {
  
  rand.weighted.bMNTD.comp[,,rep] = as.matrix(comdistnt(t(match.phylo.otu$data),
                                                        taxaShuffle(cophenetic(match.phylo.otu$phy)),
                                                        abundance.weighted=T,
                                                        exclude.conspecifics = F));
  
  print(c(date(),rep));
  
}

bNTI = matrix(c(NA),nrow=ncol(match.phylo.otu$data),
              ncol=ncol(match.phylo.otu$data));
dim(bNTI);


for (columns in 1:(ncol(match.phylo.otu$data)-1)) {
  for (rows in (columns+1):ncol(match.phylo.otu$data)) {
    
    rand.vals = rand.weighted.bMNTD.comp[rows,columns,];
    bNTI[rows,columns] = (beta.mntd.weighted[rows,columns] - mean(rand.vals)) / sd(rand.vals);
    rm("rand.vals");
    
  };
};

rownames(bNTI) = colnames(match.phylo.otu$data)
colnames(bNTI) = colnames(match.phylo.otu$data)

bNTI <- as.matrix(bNTI)
bNTI_long <-  melt(bNTI, varnames = c("Sample1", "Sample2"), value.name = "bNTI")
bNTI_long <- bNTI_long[!is.na(bNTI_long$bNTI),]

#Subset same compartment, compare across Spartina phenotype
#Bulk sediment
bNTI_long_Bulk <- bNTI_long[grepl("Bulk", bNTI_long$Sample1) & grepl("Bulk", bNTI_long$Sample2),]

bNTI_long_Bulk_tall <- bNTI_long_Bulk[grepl("P1", bNTI_long_Bulk$Sample1) & grepl("P1", bNTI_long_Bulk$Sample2),]
bNTI_long_Bulk_tall$Spartina <- "Tall"

bNTI_long_Bulk_medium <- bNTI_long_Bulk[grepl("P2", bNTI_long_Bulk$Sample1) & grepl("P2", bNTI_long_Bulk$Sample2),]
bNTI_long_Bulk_medium$Spartina <- "Medium"

bNTI_long_Bulk_short <- bNTI_long_Bulk[(grepl("P3", bNTI_long_Bulk$Sample1) | grepl("P4", bNTI_long_Bulk$Sample1) ) &
                                         (grepl("P3", bNTI_long_Bulk$Sample2) | grepl("P4", bNTI_long_Bulk$Sample2)),]
bNTI_long_Bulk_short$Spartina <- "Short"

Bulk_Spartina_bNTI <- rbind(bNTI_long_Bulk_tall[,c(3,4)], 
                            bNTI_long_Bulk_medium[,c(3,4)], bNTI_long_Bulk_short[,c(3,4)])
Bulk_Spartina_bNTI$Compartment <- "Bulk_sediment"

#Rhizosphere
bNTI_long_Rhizosphere <- bNTI_long[grepl("Rhizosphere", bNTI_long$Sample1) & grepl("Rhizosphere", bNTI_long$Sample2),]

bNTI_long_Rhizosphere_tall <- bNTI_long_Rhizosphere[grepl("P1", bNTI_long_Rhizosphere$Sample1) & grepl("P1", bNTI_long_Rhizosphere$Sample2),]
bNTI_long_Rhizosphere_tall$Spartina <- "Tall"

bNTI_long_Rhizosphere_medium <- bNTI_long_Rhizosphere[grepl("P2", bNTI_long_Rhizosphere$Sample1) & grepl("P2", bNTI_long_Rhizosphere$Sample2),]
bNTI_long_Rhizosphere_medium$Spartina <- "Medium"

bNTI_long_Rhizosphere_short <- bNTI_long_Rhizosphere[(grepl("P3", bNTI_long_Rhizosphere$Sample1) | grepl("P4", bNTI_long_Rhizosphere$Sample1) ) &
                                                       (grepl("P3", bNTI_long_Rhizosphere$Sample2) | grepl("P4", bNTI_long_Rhizosphere$Sample2)),]
bNTI_long_Rhizosphere_short$Spartina <- "Short"

Rhizosphere_Spartina_bNTI <- rbind(bNTI_long_Rhizosphere_tall[,c(3,4)], 
                                   bNTI_long_Rhizosphere_medium[,c(3,4)], bNTI_long_Rhizosphere_short[,c(3,4)])
Rhizosphere_Spartina_bNTI$Compartment <- "Rhizosphere"

#Endosphere

bNTI_long_Endosphere <- bNTI_long[grepl("Endosphere", bNTI_long$Sample1) & grepl("Endosphere", bNTI_long$Sample2),]

bNTI_long_Endosphere_tall <- bNTI_long_Endosphere[grepl("P1", bNTI_long_Endosphere$Sample1) & grepl("P1", bNTI_long_Endosphere$Sample2),]
bNTI_long_Endosphere_tall$Spartina <- "Tall"

bNTI_long_Endosphere_medium <- bNTI_long_Endosphere[grepl("P2", bNTI_long_Endosphere$Sample1) & grepl("P2", bNTI_long_Endosphere$Sample2),]
bNTI_long_Endosphere_medium$Spartina <- "Medium"

bNTI_long_Endosphere_short <- bNTI_long_Endosphere[(grepl("P3", bNTI_long_Endosphere$Sample1) | grepl("P4", bNTI_long_Endosphere$Sample1) ) &
                                                     (grepl("P3", bNTI_long_Endosphere$Sample2) | grepl("P4", bNTI_long_Endosphere$Sample2)),]
bNTI_long_Endosphere_short$Spartina <- "Short"

Endosphere_Spartina_bNTI <- rbind(bNTI_long_Endosphere_tall[,c(3,4)], 
                                  bNTI_long_Endosphere_medium[,c(3,4)], bNTI_long_Endosphere_short[,c(3,4)])
Endosphere_Spartina_bNTI$Compartment <- "Endosphere"

Compartment_same_Spartina_bNTI_Sapelo19 <- rbind(Bulk_Spartina_bNTI,Rhizosphere_Spartina_bNTI,Endosphere_Spartina_bNTI)

Compartment_same_Spartina_bNTI_Sapelo19$Location <- "Sapelo 2019"





###########Skidaway 2019 

## read phylogeny

phylo = phy_tree(phy_skidaway_2019_prevfilt_ra);
phylo; # a summary of the phylogeny
## make sure the names on the phylogeny are ordered the same as the names in otu table

match.phylo.otu = match.phylo.data(phylo, t(otu));
str(match.phylo.otu);

## calculate empirical betaMNTD

beta.mntd.weighted = as.matrix(comdistnt(t(match.phylo.otu$data),
                                         cophenetic(match.phylo.otu$phy),abundance.weighted=T));
dim(beta.mntd.weighted);
beta.mntd.weighted[1:5,1:5];
write.csv(beta.mntd.weighted,'bMNTD/betaMNTD_weighted.csv',quote=F);

identical(colnames(match.phylo.otu$data),colnames(beta.mntd.weighted)); # just a check, should be TRUE
identical(colnames(match.phylo.otu$data),rownames(beta.mntd.weighted)); # just a check, should be TRUE

# calculate randomized betaMNTD

beta.reps = 999; # number of randomizations

rand.weighted.bMNTD.comp = array(c(-999),dim=c(ncol(match.phylo.otu$data),
                                               ncol(match.phylo.otu$data),beta.reps));
dim(rand.weighted.bMNTD.comp);

for (rep in 1:beta.reps) {
  
  rand.weighted.bMNTD.comp[,,rep] = as.matrix(comdistnt(t(match.phylo.otu$data),
                                                        taxaShuffle(cophenetic(match.phylo.otu$phy)),
                                                        abundance.weighted=T,
                                                        exclude.conspecifics = F));
  
  print(c(date(),rep));
  
}

bNTI = matrix(c(NA),nrow=ncol(match.phylo.otu$data),
              ncol=ncol(match.phylo.otu$data));
dim(bNTI);


for (columns in 1:(ncol(match.phylo.otu$data)-1)) {
  for (rows in (columns+1):ncol(match.phylo.otu$data)) {
    
    rand.vals = rand.weighted.bMNTD.comp[rows,columns,];
    bNTI[rows,columns] = (beta.mntd.weighted[rows,columns] - mean(rand.vals)) / sd(rand.vals);
    rm("rand.vals");
    
  };
};

rownames(bNTI) = colnames(match.phylo.otu$data)
colnames(bNTI) = colnames(match.phylo.otu$data)

bNTI <- as.matrix(bNTI)
bNTI_long <-  melt(bNTI, varnames = c("Sample1", "Sample2"), value.name = "bNTI")
bNTI_long <- bNTI_long[!is.na(bNTI_long$bNTI),]

#Subset same compartment, compare across Spartina phenotype
#Bulk sediment
bNTI_long_Bulk <- bNTI_long[grepl("Bulk", bNTI_long$Sample1) & grepl("Bulk", bNTI_long$Sample2),]

bNTI_long_Bulk_tall <- bNTI_long_Bulk[grepl("P1", bNTI_long_Bulk$Sample1) & grepl("P1", bNTI_long_Bulk$Sample2),]
bNTI_long_Bulk_tall$Spartina <- "Tall"

bNTI_long_Bulk_medium <- bNTI_long_Bulk[grepl("P2", bNTI_long_Bulk$Sample1) & grepl("P2", bNTI_long_Bulk$Sample2),]
bNTI_long_Bulk_medium$Spartina <- "Medium"

bNTI_long_Bulk_short <- bNTI_long_Bulk[(grepl("P3", bNTI_long_Bulk$Sample1) | grepl("P4", bNTI_long_Bulk$Sample1) ) &
                                         (grepl("P3", bNTI_long_Bulk$Sample2) | grepl("P4", bNTI_long_Bulk$Sample2)),]
bNTI_long_Bulk_short$Spartina <- "Short"

Bulk_Spartina_bNTI <- rbind(bNTI_long_Bulk_tall[,c(3,4)], 
                            bNTI_long_Bulk_medium[,c(3,4)], bNTI_long_Bulk_short[,c(3,4)])
Bulk_Spartina_bNTI$Compartment <- "Bulk_sediment"

#Rhizosphere
bNTI_long_Rhizosphere <- bNTI_long[grepl("Rhizosphere", bNTI_long$Sample1) & grepl("Rhizosphere", bNTI_long$Sample2),]

bNTI_long_Rhizosphere_tall <- bNTI_long_Rhizosphere[grepl("P1", bNTI_long_Rhizosphere$Sample1) & grepl("P1", bNTI_long_Rhizosphere$Sample2),]
bNTI_long_Rhizosphere_tall$Spartina <- "Tall"

bNTI_long_Rhizosphere_medium <- bNTI_long_Rhizosphere[grepl("P2", bNTI_long_Rhizosphere$Sample1) & grepl("P2", bNTI_long_Rhizosphere$Sample2),]
bNTI_long_Rhizosphere_medium$Spartina <- "Medium"

bNTI_long_Rhizosphere_short <- bNTI_long_Rhizosphere[(grepl("P3", bNTI_long_Rhizosphere$Sample1) | grepl("P4", bNTI_long_Rhizosphere$Sample1) ) &
                                                       (grepl("P3", bNTI_long_Rhizosphere$Sample2) | grepl("P4", bNTI_long_Rhizosphere$Sample2)),]
bNTI_long_Rhizosphere_short$Spartina <- "Short"

Rhizosphere_Spartina_bNTI <- rbind(bNTI_long_Rhizosphere_tall[,c(3,4)], 
                                   bNTI_long_Rhizosphere_medium[,c(3,4)], bNTI_long_Rhizosphere_short[,c(3,4)])
Rhizosphere_Spartina_bNTI$Compartment <- "Rhizosphere"

#Endosphere

bNTI_long_Endosphere <- bNTI_long[grepl("Endosphere", bNTI_long$Sample1) & grepl("Endosphere", bNTI_long$Sample2),]

bNTI_long_Endosphere_tall <- bNTI_long_Endosphere[grepl("P1", bNTI_long_Endosphere$Sample1) & grepl("P1", bNTI_long_Endosphere$Sample2),]
bNTI_long_Endosphere_tall$Spartina <- "Tall"

bNTI_long_Endosphere_medium <- bNTI_long_Endosphere[grepl("P2", bNTI_long_Endosphere$Sample1) & grepl("P2", bNTI_long_Endosphere$Sample2),]
bNTI_long_Endosphere_medium$Spartina <- "Medium"

bNTI_long_Endosphere_short <- bNTI_long_Endosphere[(grepl("P3", bNTI_long_Endosphere$Sample1) | grepl("P4", bNTI_long_Endosphere$Sample1) ) &
                                                     (grepl("P3", bNTI_long_Endosphere$Sample2) | grepl("P4", bNTI_long_Endosphere$Sample2)),]
bNTI_long_Endosphere_short$Spartina <- "Short"

Endosphere_Spartina_bNTI <- rbind(bNTI_long_Endosphere_tall[,c(3,4)], 
                                  bNTI_long_Endosphere_medium[,c(3,4)], bNTI_long_Endosphere_short[,c(3,4)])
Endosphere_Spartina_bNTI$Compartment <- "Endosphere"

Compartment_same_Spartina_bNTI_Skio19 <- rbind(Bulk_Spartina_bNTI,Rhizosphere_Spartina_bNTI,Endosphere_Spartina_bNTI)

Compartment_same_Spartina_bNTI_Skio19$Location <- "Skidaway 2019"








#Same Spartina bNTI

Compartment_same_Spartina_bNTI <- rbind(Compartment_same_Spartina_bNTI_Skio19, 
                                        Compartment_same_Spartina_bNTI_Sapelo19,
                                        Compartment_same_Spartina_bNTI_Sapelo18)

Compartment_same_Spartina_bNTI$Compartment <- as.factor(Compartment_same_Spartina_bNTI$Compartment)
Compartment_same_Spartina_bNTI$Compartment <- factor(Compartment_same_Spartina_bNTI$Compartment,
                                                     levels(Compartment_same_Spartina_bNTI$Compartment)[c(1,3,2)])

Compartment_same_Spartina_bNTI$Location <- as.factor(Compartment_same_Spartina_bNTI$Location)
Compartment_same_Spartina_bNTI$Location <- factor(Compartment_same_Spartina_bNTI$Location,
                                                  levels(Compartment_same_Spartina_bNTI$Location)[c(3,2,1)])


#PLOT
theme_bNTI <- theme_bw() + theme(axis.title.x = element_blank(),
                                 axis.title.y = element_text(size = 24),
                                 axis.text = element_text(size = 20),
                                 panel.grid.major = element_blank(),
                                 panel.grid.minor = element_blank(),
                                 strip.text.x = element_text(size = 20),
                                 strip.background = element_rect(fill = "grey90"))


#Histogram
fig_same_Spartina_hist <- ggplot(Compartment_same_Spartina_bNTI, aes(x = bNTI)) + 
  geom_histogram(aes(y = ..density..), bins = 18) + facet_wrap(~Location*Compartment) +
  geom_vline(aes(xintercept = 2), linetype = "dashed") + 
  geom_vline(aes(xintercept = -2), linetype = "dashed") + 
  geom_density() + theme_bNTI + labs(y = "Density", x = "bNTI")

jpeg("bNTI/Summary/same_Spartina_bNTI_hist.jpeg", res = 300, height = 3000, width = 4000)
fig_same_Spartina_hist
dev.off()



#########DeSEQ2 analysis













