#### Community Assembly and Functioning under Environmental Stress ####
#### Authors: Jessica R. Bernardin, Leonora S. Bittleston ####
#### last update : Feb 28, 2025 ####

#### Load Required Packages ####
packages_to_load <- c(
  "ggplot2", "vegan", "lme4", "tidyverse", "effects",
  "dplyr", "reshape2", "ape", "DiagrammeR",
  "tidybayes", "coefplot", "standardize", "bayesplot", "MCMCvis", "car",
  "patchwork", "ggpubr", "corrr", "ggcorrplot", "factoextra", "MASS",
  "pairwiseAdonis", "plotrix", "gridExtra", "multcompView", "ggeffects", "this.path", "brms", "ggmulti",
  "phyloseq", "qiime2R", "picante", "decontam", "performance", "janitor", "ANCOMBC", "pheatmap", "chron",
  "lubridate", "igraph", "Hmisc", "qgraph", "mia", "cowplot", "microbiome", "stringr", "ggtree", "sjPlot"
)

#install.packages("BiocManager")
library(BiocManager)

#install.packages("devtools")
library(devtools)

#install.packages("microbiome")
#BiocManager::install("mia")
#BiocManager::install("ggdraw")
#if (!requireNamespace("BiocManager", quietly = TRUE))
  #  install.packages("BiocManager")
#BiocManager::install("phyloseq")
#BiocManager::install("ggtree")
#if (!requireNamespace("devtools", quietly = TRUE)){install.packages("devtools")}
#devtools::install_github("jbisanz/qiime2R")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("decontam")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("ANCOMBC")


#install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")

# Load and install required packages
for (i in packages_to_load) { #Installs packages if not yet installed
  if (!require(i, character.only = TRUE)) install.packages(i)
}

setwd(this.path::here())

#### Community Structure ####
#read in 16S
asv16s <- read_tsv("exported-files/asv-table-dada2.txt")
asv16s <- asv16s %>% column_to_rownames(var="#OTU ID")

meta <- read_csv("data/Exp_2_metadata_tubes.csv")
meta <- meta %>% remove_rownames %>% column_to_rownames(var="ID")

#make a new column for the scaled data
meta$scaled_chit <- as.numeric(scale(meta$chitinase))
meta$scaled_prot <- as.numeric(scale(meta$protease))
meta$sample_id <- as.factor(meta$sample_id)

tax.16s <- read_tsv("exported-files/taxonomy.tsv")
tax.16s <- tax.16s %>% column_to_rownames(var="Feature ID")
tax.16s <- tax.16s[order(rownames(tax.16s)), ]
asv16s <- asv16s[order(rownames(asv16s)), ]

row.names(asv16s) == row.names(tax.16s) # sanity check
nrow(asv16s) #=3342 ASVs
nrow(tax.16s) #=3342 ASVs

asv16s.tree <- read.tree("exported-files/tree_sepp.nwk")
asv16s.tree <-root(asv16s.tree, "63ffc596781727668a90d6da135a33b4")## root with an archaeon

asv16s.exp2 <- asv16s %>%
  filter(rowSums(.) != 0)
dim(asv16s.exp2)#1585 asvs, 152 samples

tax.16s.exp2 <- subset(tax.16s, row.names(tax.16s) %in% rownames(asv16s.exp2)) 

tax.16s.exp2 <- tax.16s.exp2[order(rownames(tax.16s.exp2)), ]
asv16s.exp2 <- asv16s.exp2[order(rownames(asv16s.exp2)), ]

row.names(asv16s.exp2) == row.names(tax.16s.exp2) # sanity check
nrow(asv16s.exp2) #=1585 ASVs
nrow(tax.16s.exp2) #=1585 ASVs

#### asv16s.physeq3 = USING all data (raw) ####
asv16s.tax1raw <- data.frame(Feature.ID=row.names(tax.16s.exp2),Taxon=tax.16s.exp2[,1])
asv16s.tax2raw <- parse_taxonomy(asv16s.tax1raw)
asv16s.tax3raw <- cbind(ASVs=row.names(asv16s.tax2raw),asv16s.tax2raw)

pruned.tree<- ape::drop.tip(asv16s.tree,asv16s.tree$tip.label[-match(rownames(asv16s.exp2), asv16s.tree$tip.label)])
new_treeraw <- ape::multi2di(pruned.tree)

TAX.16sraw <- tax_table(as.matrix(asv16s.tax3raw))
OTU.16sraw <- otu_table(as.matrix(asv16s.exp2), taxa_are_rows = TRUE)
SAM.16sraw <- sample_data(meta)
Exp2.physeq1_raw <-merge_phyloseq(phyloseq(OTU.16sraw),SAM.16sraw,TAX.16sraw,new_treeraw)
Exp2.physeq1_raw #1585 taxa and 152 samples
saveRDS(Exp2.physeq1_raw, "RDS/Exp2.physeq1_raw.RDS")

#### DECONTAM PACKAGE for identifying contaminants####
#https://benjjneb.github.io/decontam/vignettes/decontam_intro.html
## Read in data and metadata as a phyloseq object (non-rarified data, non-filtered)
decon <- as.data.frame(sample_data(Exp2.physeq1_raw))
decon$LibrarySize <- sample_sums(Exp2.physeq1_raw)
decon <- decon[order(decon$LibrarySize),]
decon$Index <- seq(nrow(decon))
ggplot(data=decon, aes(x=Index, y=LibrarySize, color=day)) + geom_point()

# Step 1: Calculate the number of samples in which each ASV is present
present_in_samples <- as.data.frame(apply(asv16s.exp2 > 0, 1, sum))
present_in_samples$sample_number <- present_in_samples[[1]]
present_in_samples <- present_in_samples %>% dplyr::select(sample_number) %>% arrange(sample_number)
hist(present_in_samples)
dim(present_in_samples)#1585 ASV

counts_by_samples <- table(present_in_samples)
counts_by_samples
sorted_counts <- sort(counts_by_samples)
print(sorted_counts)#1989 ASVs are only present in one sample
hist(sorted_counts)

#Prevalance Method (this is raw data, no filtering)
sample_data(Exp2.physeq1_raw)$is.neg <- sample_data(Exp2.physeq1_raw)$container == "NEG"
contamdf.prev <- isContaminant(Exp2.physeq1_raw, method="prevalence", neg="is.neg")
contamdf.prev
table(contamdf.prev$contaminant) ### 37 ASVs are identified as contaminants, 1548 as not
head(which(contamdf.prev$contaminant), n=37)

ggplot(data=contamdf.prev, aes(x=p)) + 
  labs(x = 'decontam-prevalence Score', y='Number of ASVs') + 
  geom_histogram(binwidth=0.02)

###filter the df to samples id'd as contaminants and find out home many samples they are present in and the read abundance for these
cont.prev.true <-subset(contamdf.prev, contaminant=="TRUE")
abund <- asv16s.exp2
abund$reads <- rowSums(asv16s.exp2)
cont.prev.true.reads <- subset(abund, row.names(abund) %in% rownames(cont.prev.true)) 
cont.true.reads.only <- cont.prev.true.reads %>% dplyr::select(reads)  
cont.prev.true.reads <- cont.prev.true.reads[,-153]

#count up how many samples each ASV is present in
count_non_zero <- function(row) {
  sum(row != 0)
}
cont.prev.true.reads$NonZeroCount <- apply(cont.prev.true.reads[, -1], MARGIN = 1, count_non_zero)

#### Clean up data ####
#these 37 ASV have low frequency and low abundance and were closely associated with the negative controls, they will be removed from the dataset
cont.tax <- subset(tax.16s.exp2, row.names(tax.16s.exp2) %in% rownames(cont.prev.true.reads))

# remove chloroplasts
tax.16s.exp2.2 <- tax.16s.exp2[grep("Chloroplast",tax.16s.exp2$Taxon, invert = T),]#25 removed

# remove mitochondria
tax.16s.exp2.2 <- tax.16s.exp2.2[grep("Mitochondria",tax.16s.exp2.2$Taxon, invert = T),] #removed 8 ASVs

tax.16s.exp2.2 <- tax.16s.exp2.2[grep("Eukaryota",tax.16s.exp2.2$Taxon, invert = T),]#7 removed

tax.16s.exp2.2 <- tax.16s.exp2.2[grep("Archaea",tax.16s.exp2.2$Taxon, invert = T),]#2 removed

# remove unassigned
tax.16s.exp2.2 <- tax.16s.exp2.2[grep("Unassigned",tax.16s.exp2.2$Taxon, invert = T),]#18 removed

#remove contaminants identified above
cont.tax <- cont.tax %>% rownames_to_column(var = "ASV")
tax.16s.exp2.2 <- tax.16s.exp2.2 %>% rownames_to_column(var = "ASV")
tax.16s.exp2.2 <- anti_join(tax.16s.exp2.2, cont.tax, by="ASV") #removed 37 ASV
rownames(tax.16s.exp2.2) <- tax.16s.exp2.2$ASV
tax.16s.exp2.2$ASV<- NULL

#remove negative controls
asv16s.exp2_filt <- asv16s.exp2[,-c(149:152)]
meta.filt <- meta[-c(5:8),]

#remove the filtered taxa from the ASV table
asv16s.exp2_filt <- asv16s.exp2_filt %>%
  filter(rowSums(.) != 0)

dim(asv16s.exp2_filt)#1479  148

asv16s.exp2_filt <- subset(asv16s.exp2_filt, row.names(asv16s.exp2_filt) %in% row.names(tax.16s.exp2.2))
tax.16s.exp2.2 <- subset(tax.16s.exp2.2, row.names(tax.16s.exp2.2) %in% row.names(asv16s.exp2_filt))

rownames(asv16s.exp2_filt) == rownames(tax.16s.exp2.2)
dim(tax.16s.exp2.2) #1406 ASVs
dim(asv16s.exp2_filt) #1406 ASVs

#### Exp2.physeq2 = decontaminated data, filtered, negative controls removed ####
asv16s.tax1raw <- data.frame(Feature.ID=row.names(tax.16s.exp2.2),Taxon=tax.16s.exp2.2[,1])
asv16s.tax2raw <- parse_taxonomy(asv16s.tax1raw)
asv16s.tax3raw <- cbind(ASVs=row.names(asv16s.tax2raw),asv16s.tax2raw)

pruned.tree<- ape::drop.tip(asv16s.tree,asv16s.tree$tip.label[-match(rownames(asv16s.exp2_filt), asv16s.tree$tip.label)])
new_treeraw <- ape::multi2di(pruned.tree)
TAX.16sraw <- tax_table(as.matrix(asv16s.tax3raw))
OTU.16sraw <- otu_table(as.matrix(asv16s.exp2_filt), taxa_are_rows = TRUE)
SAM.16sraw <- sample_data(meta)
Exp2.physeq2 <-merge_phyloseq(phyloseq(OTU.16sraw),SAM.16sraw,TAX.16sraw,new_treeraw)
Exp2.physeq2 #1406 taxa and 148 samples
saveRDS(Exp2.physeq2, "RDS/Exp2.physeq2.RDS")
Exp2.physeq2 <- readRDS("RDS/Exp2.physeq2.RDS")

Exp2.physeq2.nomix <- subset_samples(Exp2.physeq2, day %in% c("1","8","15","22","29","36","43","50","57")) 
row_sums <- rowSums(otu_table(Exp2.physeq2.nomix))
nonzero_rows <- row_sums != 0
Exp2.physeq2.nomix <- prune_taxa(nonzero_rows, Exp2.physeq2.nomix)#1366 taxa and 144 samples

#### Rarify ####
asv16s.exp2_filt <- asv16s.exp2_filt[as.logical(rowSums(asv16s.exp2_filt != 0)), ] #1406 ASV left
tax.16s.exp2_filt <- subset(tax.16s.exp2.2, row.names(tax.16s.exp2.2) %in% row.names(asv16s.exp2_filt))

asv16s.exp2_filt.t <- t(asv16s.exp2_filt)
min(rowSums(asv16s.exp2_filt.t)) #min number of reads per sample = 6798
max(rowSums(asv16s.exp2_filt.t)) #max number of reads per sample = 112126
mean(rowSums(asv16s.exp2_filt.t))

set.seed(1115)
asv16s.rt <- rrarefy(asv16s.exp2_filt.t, 6798) ## rarefy at 6798, samples are rows
asv16s.rt <- asv16s.rt[,colSums(asv16s.rt) > 0] ## remove ASVs no longer present, 397 ASV removed
asv16s.rt <- asv16s.rt[order(row.names(asv16s.rt)),] # order samples alphabetically
asv16s.rt <- asv16s.rt[,order(colnames(asv16s.rt))] # order asvs alphabetically
dim(asv16s.rt) #148 samples, 1009 ASVs
summary(rowSums(asv16s.rt)) #=6798
summary(colSums(asv16s.rt)) #1

asv16s.r <- t(asv16s.rt) #rows are asvs and columns are samples

#filter taxonomy to match the 579 ASVs
tax.16s2.r <- subset(tax.16s.exp2_filt, row.names(tax.16s.exp2_filt) %in% row.names(asv16s.r)) #filter tax table to match asvs
nrow(asv16s.r)#1009
nrow(tax.16s2.r)#1009

meta.filt <- meta[order(row.names(meta)),] # order samples alphabetically
meta.r <- subset(meta.filt, row.names(meta.filt) %in% colnames(asv16s.r)) #filter meta to relevant samples

#now asv, tax, meta all the same size
row.names(asv16s.rt) == row.names(meta.r)
tax.16s2.r <- tax.16s2.r[order(row.names(tax.16s2.r)),] # order ASVs alphabetically
asv16s.r <- asv16s.r[order(row.names(asv16s.r)),] # order ASVs alphabetically
row.names(asv16s.r) == row.names(tax.16s2.r)

#### 16S alpha diversity ####
shannon.16s <- vegan::diversity(asv16s.rt, index="shannon")
ef.16s <- exp(shannon.16s)
ef.16s <- round(ef.16s, 3)
summary(ef.16s)
richness <- colSums(asv16s.r !=0)
md16s <- cbind(meta.r, richness)
md16s <- cbind(md16s, shannon.16s)
md16s <- cbind(md16s, ef.16s)
write.csv(md16s, "output_files/exp2_md16s.csv", row.names=TRUE)
md16s <- md16s[order(row.names(md16s)),] # order samples alphabetically
md16s$ID <- rownames(md16s)
row.names(asv16s.rt) == row.names(md16s)

#visualize alpha diversity metrics
md16s_mix <- filter(md16s, container == "MIX")
summary(md16s_mix)

## Figure S1 ##
ggplot(data=md16s_mix, aes(x=ID, y=richness)) +
  geom_jitter() + ylab("Richness (ASVs)") + theme_classic()+
  ylim(c(0,50))

summary(c(47, 45, 46, 48))
SE = sd(c(47, 45, 46, 48)) / sqrt(length(c(47, 45, 46, 48)))
md16s$food <- as.factor(md16s$food)
md16s$ph <- as.factor(md16s$ph)
md16s$temperature <- as.factor(md16s$temperature)

detach("package:plyr", unload=TRUE)

df.summaryrich <- md16s %>%
  group_by(treatment_combo,day, food, temperature, ph) %>%
  summarise(
    sd = sd(richness, na.rm = TRUE),
    richness = mean(richness, na.rm = TRUE))
df.summaryrich <- na.omit(df.summaryrich)

df.summaryrich$ph <- factor(df.summaryrich$ph, levels=c("3", "4", "5.6"))

## Figure 2A ##
ggplot(df.summaryrich, aes(day, richness, color = temperature, shape = food)) +
  geom_jitter(size = 3, alpha = 0.9, position = position_dodge(0.3)) +
  facet_wrap(~ph, nrow = 1) +
  geom_line(aes(group = treatment_combo), position = position_dodge(0.3)) +
  geom_errorbar(aes(ymin = richness - sd, ymax = richness + sd), position = position_dodge(0.3), width = 0.2) +
  theme_bw() +scale_shape_manual(values=c(16,2)) +
  scale_color_manual(values = c("#1f5776", "#c83126")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none") +
  scale_y_continuous(limits = c(20, 80), breaks = c(20, 40, 60, 80))


#### GLM ASV Richness ####
md16s$food <- relevel(md16s$food, ref = "3")
md16s$temperature <- relevel(md16s$temperature, ref = "22")
md16s$ph <- relevel(md16s$ph, ref = "5.6")
md16s$day <- as.integer(md16s$day)

#check variance to see if we use a negative binomial or poisson model distribution
mean_richness <- mean(md16s$richness)
var_richness <- var(md16s$richness)
var_richness / mean_richness #2.730272

#greater than 1 so overdispersed, use negbin
mrich <- brm(richness ~ ph*food*temperature+day,data=md16s, family=negbinomial(),iter = 10000, chains = 4, cores = 4)#,control = list(adapt_delta = 0.999, stepsize = 0.001, max_treedepth=20)
saveRDS(mrich, "RDS/Exp2_richness_tubes.RDS")
mrich <- readRDS("RDS/Exp2_richness_tubes.RDS")
summary(mrich)
mcmc_plot(mrich)
rich.pred <- predict_response(mrich, terms = c("food", "temperature", "ph"))
rich.pred$facet <- factor(rich.pred$facet, levels = c("3", "4", "5.6"))

## Figure 2B ##
ggplot(rich.pred, aes(x = x, y = predicted, shape = x, color = group)) + 
  facet_wrap(~facet) +
  scale_color_manual(values = c("#1f5776", "#c83126")) +
  theme_bw() +scale_shape_manual(values=c(16,2)) +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0, size=1,position = position_dodge(width = .7))+
  theme(legend.position="none")+ ylab("Predicted ASV Richness")+
  xlab("Food (g/L)")+ geom_point(size = 3, position = position_dodge(width = .7)) +theme(
    text = element_text(color = "black", size = 10 ))+
  scale_y_continuous(limits = c(20, 80), breaks = c(20, 40, 60, 80))

tab_model(mrich)

tab_model(mrich, show.est = TRUE, show.ci = TRUE, show.se = TRUE, 
          file = "output_files/brms_richness_table.doc")

posteriormrich <- mcmc_intervals_data(mrich, 
                                     prob_outer=0.95,
                                     prob=0.95)

posteriormrich$nonzero <- NA
posteriormrich$nonzero[posteriormrich$ll>0 & posteriormrich$hh>0] <- "nonzero"
posteriormrich$nonzero[posteriormrich$ll<0 & posteriormrich$hh<0] <- "nonzero"
posteriormrich$nonzero[is.na(posteriormrich$nonzero)] <- "zero"
posteriormrich<- posteriormrich[1:13,]
posteriormrich <- posteriormrich %>%
  mutate(parameter = str_remove(parameter, "^b_") %>%  # Remove "b_" prefix
           str_replace_all(":", ": ") %>%
           str_replace_all("ph3", "pH 3") %>%
           str_replace_all("ph4", "pH 4") %>%
           str_replace_all("food6", "high food") %>%
           str_replace_all("temperature37", "high temp"))

posteriormrich <- posteriormrich %>%
  mutate(parameter = factor(parameter, 
                            levels = c("Intercept", 
                                       "day",
                                       "high temp", 
                                       "high food", 
                                       "pH 3", 
                                       "pH 4", 
                                       "pH 3: high food", 
                                       "pH 4: high food", 
                                       "pH 3: high temp", 
                                       "pH 4: high temp", 
                                       "high food: high temp", 
                                       "pH 3: high food: high temp", 
                                       "pH 4: high food: high temp")))

## Figure S3 ##
ggplot(posteriormrich, aes(x = parameter, shape=nonzero)) +
  geom_hline(yintercept = 0, linetype = 2, size=.35, color = "black") +
  geom_pointrange(aes(ymin = ll, ymax = hh, y = m), position= position_dodge(width=0.75), size = 1, linewidth = 1) +
  scale_shape_manual(values=c(17, 19), labels=c("95% CI does\nnot contain zero", "95% CI\ncontains zero"))+
  coord_flip() +xlab(NULL) + ylab("Estimated effect on ASV Richness")+ theme(text = element_text(color = "black", size = 16))+
  theme(panel.background = element_rect(fill = "white", colour = "black"))+theme(legend.position="none")

#### MAKE PHYLOSEQ OBJECTS ####
#Exp2.physeq3.RDS = RARIFIED
row.names(asv16s.r) == row.names(tax.16s2.r)
r.asv16s.tax1 <- data.frame(Feature.ID=row.names(tax.16s2.r),Taxon=tax.16s2.r[,1])
r.asv16s.tax2 <- parse_taxonomy(r.asv16s.tax1)
r.asv16s.tax3 <- cbind(ASVs=row.names(r.asv16s.tax2),r.asv16s.tax2)

r.tree.16s <- picante::prune.sample(asv16s.rt, asv16s.tree)# Subset tree to relevant samples
r.new_tree <- ape::multi2di(r.tree.16s)

r.TAX.16s <- tax_table(as.matrix(r.asv16s.tax3))
r.OTU.16s <- otu_table(as.matrix(asv16s.r), taxa_are_rows = TRUE)
r.SAM.16s <- sample_data(meta.r)
Exp2.physeq3 <-merge_phyloseq(phyloseq(r.OTU.16s),r.SAM.16s,r.TAX.16s,r.new_tree)
Exp2.physeq3
#1009 taxa and 148 samples
saveRDS(Exp2.physeq3, "RDS/Exp2.physeq3.RDS")
Exp2.physeq3 <- readRDS("RDS/Exp2.physeq3.RDS")


Exp2.tubes.mix.physeq3 = subset_samples(Exp2.physeq3, day %in% c("0", "1")) 
row_sums <- rowSums(otu_table(Exp2.tubes.mix.physeq3))
nonzero_rows <- row_sums != 0
Exp2.tubes.mix.physeq3 <- prune_taxa(nonzero_rows, Exp2.tubes.mix.physeq3)#84 taxa and 52 samples


Exp2.tubes.mix2.physeq3 = subset_samples(Exp2.physeq3, day %in% c("0"))
row_sums <- rowSums(otu_table(Exp2.tubes.mix2.physeq3))
nonzero_rows <- row_sums != 0
Exp2.tubes.mix2.physeq3 <- prune_taxa(nonzero_rows, Exp2.tubes.mix2.physeq3)#51 taxa and 4 samples
mixonly.meta <- data.frame(sample_data(Exp2.tubes.mix2.physeq3))

Exp2.tubes.nomix.physeq3 = subset_samples(Exp2.physeq3, day %in% c("1", "15", "57"))
row_sums <- rowSums(otu_table(Exp2.tubes.nomix.physeq3))
nonzero_rows <- row_sums != 0
Exp2.tubes.nomix.physeq3 <- prune_taxa(nonzero_rows, Exp2.tubes.nomix.physeq3)#51 taxa and 4 samples

####  dbRDA BETA DIVERSITY ####
## Calculate weighted Unifrac distance
wu.dist.16s.tubes <- phyloseq::distance(Exp2.tubes.nomix.physeq3, method="wunifrac")
tubes.nomix.meta <- data.frame(sample_data(Exp2.tubes.nomix.physeq3))
tubes.nomix.meta$sample_id <- as.numeric(tubes.nomix.meta$sample_id)
tubes.nomix.meta$ph <- factor(tubes.nomix.meta$ph, levels = c(5.6, 4, 3))
tubes.nomix.meta$food <- factor(tubes.nomix.meta$food, levels = c(3, 6))
tubes.nomix.meta$temperature <- factor(tubes.nomix.meta$temperature, levels = c(22,37))
tubes.nomix.meta$day <- as.numeric(tubes.nomix.meta$day)

dbrda_beta <- dbrda(as.dist(wu.dist.16s.tubes) ~ temperature*ph*food+ day + Condition(sample_id),
                   data = tubes.nomix.meta,
                   distance = "NULL")

dbrda_beta #constrained prop = 0.36369
beta_sum_dbrda <- summary(dbrda_beta) #dbRDA1 and 2 Proportion Explained  0.5557 0.2380
beta_dbrda_anova <- anova.cca(dbrda_beta, by = "onedf", perm = 999)
beta_dbrda_anova
write.csv(beta_dbrda_anova, "output_files/beta_16s_dbrda_cca_aov.csv", row.names=TRUE)
significant_factors <- c("temperature37","food6", "ph3","day","temperature37:food6")

beta_est <- as.data.frame(cbind(x1 = dbrda_beta$CCA$biplot[,1], y1 = dbrda_beta$CCA$biplot[,2]))
beta_est$factor <- row.names(dbrda_beta$CCA$biplot)
beta_est$significant <- beta_est$factor %in% significant_factors
significant_beta_est <- beta_est[beta_est$significant, ]
sites_beta <- as.data.frame(beta_sum_dbrda$sites) 
sites_beta_meta <- cbind(sites_beta,tubes.nomix.meta )

sites_beta_meta$day <- as.factor(sites_beta_meta$day)
sites_beta_meta$temperature <- as.factor(sites_beta_meta$temperature)
sites_beta_meta$ph <- factor(sites_beta_meta$ph, levels=c(3, 4, 5.6))
sites_beta_meta$food <- factor(sites_beta_meta$food, levels=c(3, 6))

## Figure 2C ##
ggplot(data = sites_beta_meta, aes(x = dbRDA1, y = dbRDA2, color = temperature, alpha=day, shape=day)) +
  geom_point(size = 2) +scale_shape_manual(values = c(16, 17, 15)) +
  scale_color_manual(values = c("#1f5776", "#c83126")) +
  geom_segment(data = significant_beta_est, aes(x = 0, y = 0, xend = x1, yend = y1),
               arrow = arrow(length = unit(.2, "cm")), color = "black", size = 1, inherit.aes = FALSE) +
  geom_text(data = significant_beta_est, aes(x = x1, y = y1, label = factor), 
            hjust = 0, vjust = 0, inherit.aes = FALSE)+ theme_classic()+scale_alpha_discrete(range = c(1,.7))+
  theme(legend.title=element_blank())

## Figure 2D ##
ggplot(data = sites_beta_meta, aes(x = dbRDA1, y = dbRDA2, color = ph, alpha=day, shape=day)) +
  geom_point(size = 2) +scale_shape_manual(values = c(16, 17, 15)) +
  scale_color_manual(values = c("#CA64A3","#9A4EAE", "#301934")) +
  geom_segment(data = significant_beta_est, aes(x = 0, y = 0, xend = x1, yend = y1),
               arrow = arrow(length = unit(.2, "cm")), color = "black", size = 1, inherit.aes = FALSE) +
  geom_text(data = significant_beta_est, aes(x = x1, y = y1, label = factor), 
            hjust = 0, vjust = 0, inherit.aes = FALSE)+ theme_classic()+scale_alpha_discrete(range = c(1,.7))+
  theme(legend.title=element_blank())

## Figure 2E ##
ggplot(data = sites_beta_meta, aes(x = dbRDA1, y = dbRDA2, color = food, alpha=day, shape=day)) +
  geom_point(size = 2) +scale_shape_manual(values = c(16, 17, 15)) +
  scale_color_manual(values = c("#90ee90","#228822")) +
  geom_segment(data = significant_beta_est, aes(x = 0, y = 0, xend = x1, yend = y1),
               arrow = arrow(length = unit(.2, "cm")), color = "black", size = 1, inherit.aes = FALSE) +
  geom_text(data = significant_beta_est, aes(x = x1, y = y1, label = factor), 
            hjust = 0, vjust = 0, inherit.aes = FALSE)+ theme_classic()+scale_alpha_discrete(range = c(1,.7))+
  theme(legend.title=element_blank())

#### Mantel Test ####
Exp2.tubes.physeq3_noNA = subset_samples(Exp2.tubes.nomix.physeq3, chitinase != "NA")
Exp2.tubes.physeq3_noNA = subset_samples(Exp2.tubes.physeq3_noNA, protease != "NA")
row_sums <- rowSums(otu_table(Exp2.tubes.physeq3_noNA))
nonzero_rows <- row_sums != 0
Exp2.tubes.physeq3_noNA <- prune_taxa(nonzero_rows, Exp2.tubes.physeq3_noNA)
wu.dist.16s.tubes_noNA <- phyloseq::distance(Exp2.tubes.physeq3_noNA,"wUniFrac")

meta_mantel <- meta.r %>% drop_na(chitinase)
meta_mantel <- meta_mantel %>% drop_na(protease)

chit.dist <- vegdist(meta_mantel$chitinase,method="euclidean", na.rm=TRUE)
prot.dist <- vegdist(meta_mantel$protease,method="euclidean", na.rm=TRUE)

chit.wu.man <- mantel(wu.dist.16s.tubes_noNA,chit.dist, method = "spearman", permutations=999)
chit.wu.man
#Mantel statistic r: 0.2459  -1 strong negative, +1 strong positive, 0 none
#Significance: 0.001

prot.wu.man <- mantel(wu.dist.16s.tubes_noNA,prot.dist, method = "spearman", permutations=999)
prot.wu.man
#Mantel statistic r: 0.1225
#Significance: 0.006

#divide up the chit and prot mantel to extreme and normal like for ecoplate
meta_en_ht<- subset(meta_mantel, temperature %in% c("37"))
meta_en_lt<- subset(meta_mantel, temperature %in% c("22"))
meta_en_ph3<- subset(meta_mantel, ph %in% c("3"))
meta_en_ph4<- subset(meta_mantel, ph %in% c("4"))
meta_en_ph5.6<- subset(meta_mantel, ph %in% c("5.6"))
meta_en_food3<- subset(meta_mantel, food %in% c("3"))
meta_en_food6<- subset(meta_mantel, food %in% c("6"))

chit.distht <- vegdist(meta_en_ht$chitinase,method="bray", na.rm=TRUE)
chit.distlt <- vegdist(meta_en_lt$chitinase,method="bray", na.rm=TRUE)
chit.distph3 <- vegdist(meta_en_ph3$chitinase,method="bray", na.rm=TRUE)
chit.distph4 <- vegdist(meta_en_ph4$chitinase,method="bray", na.rm=TRUE)
chit.distph5.6 <- vegdist(meta_en_ph5.6$chitinase,method="bray", na.rm=TRUE)
chit.distfood3 <- vegdist(meta_en_food3$chitinase,method="bray", na.rm=TRUE)
chit.distfood6 <- vegdist(meta_en_food6$chitinase,method="bray", na.rm=TRUE)

prot.distht <- vegdist(meta_en_ht$protease,method="bray", na.rm=TRUE)
prot.distlt <- vegdist(meta_en_lt$protease,method="bray", na.rm=TRUE)
prot.distph3 <- vegdist(meta_en_ph3$protease,method="bray", na.rm=TRUE)
prot.distph4 <- vegdist(meta_en_ph4$protease,method="bray", na.rm=TRUE)
prot.distph5.6 <- vegdist(meta_en_ph5.6$protease,method="bray", na.rm=TRUE)
prot.distfood3 <- vegdist(meta_en_food3$protease,method="bray", na.rm=TRUE)
prot.distfood6 <- vegdist(meta_en_food6$protease,method="bray", na.rm=TRUE)

asven_ht <- subset(asv16s.rt, row.names(asv16s.rt) %in% rownames(meta_en_ht))
asven_lt <- subset(asv16s.rt, row.names(asv16s.rt) %in% rownames(meta_en_lt))
asven_ph3 <- subset(asv16s.rt, row.names(asv16s.rt) %in% rownames(meta_en_ph3))
asven_ph4 <- subset(asv16s.rt, row.names(asv16s.rt) %in% rownames(meta_en_ph4))
asven_ph5.6 <- subset(asv16s.rt, row.names(asv16s.rt) %in% rownames(meta_en_ph5.6))
asven_food3 <- subset(asv16s.rt, row.names(asv16s.rt) %in% rownames(meta_en_food3))
asven_food6 <- subset(asv16s.rt, row.names(asv16s.rt) %in% rownames(meta_en_food6))

rownames_1 <- rownames(meta_en_ht)
asven_ht <- asven_ht[rownames_1, ]
rownames(asven_ht) == rownames(meta_en_ht)
asven_ht <- asven_ht[,colSums(asven_ht) >0]

rownames_2 <- rownames(meta_en_lt)
asven_lt <- asven_lt[rownames_2, ]
rownames(asven_lt) == rownames(meta_en_lt)
asven_lt <- asven_lt[,colSums(asven_lt) >0]

rownames_3 <- rownames(meta_en_ph3)
asven_ph3 <- asven_ph3[rownames_3, ]
rownames(asven_ph3) == rownames(meta_en_ph3)
asven_ph3 <- asven_ph3[,colSums(asven_ph3) >0]

rownames_4 <- rownames(meta_en_ph4)
asven_ph4 <- asven_ph4[rownames_4, ]
rownames(asven_ph4) == rownames(meta_en_ph4)
asven_ph4 <- asven_ph4[,colSums(asven_ph4) >0]

rownames_5 <- rownames(meta_en_ph5.6)
asven_ph5.6 <- asven_ph5.6[rownames_5, ]
rownames(asven_ph5.6) == rownames(meta_en_ph5.6)
asven_ph5.6 <- asven_ph5.6[,colSums(asven_ph5.6) >0]

rownames_6 <- rownames(meta_en_food3)
asven_food3 <- asven_food3[rownames_6, ]
rownames(asven_food3) == rownames(meta_en_food3)
asven_food3 <- asven_food3[,colSums(asven_food3) >0]

rownames_7 <- rownames(meta_en_food6)
asven_food6 <- asven_food6[rownames_7, ]
rownames(asven_food6) == rownames(meta_en_food6)
asven_food6 <- asven_food6[,colSums(asven_food6) >0]

mantel_asv_test <- as.matrix(wu.dist.16s.tubes)

asv_dist_ht <- as.dist(mantel_asv_test[rownames_1, rownames_1])
asv_dist_lt <- as.dist(mantel_asv_test[rownames_2, rownames_2])
asv_dist_ph3 <- as.dist(mantel_asv_test[rownames_3, rownames_3])
asv_dist_ph4 <- as.dist(mantel_asv_test[rownames_4, rownames_4])
asv_dist_ph5.6 <- as.dist(mantel_asv_test[rownames_5, rownames_5])
asv_dist_food3 <- as.dist(mantel_asv_test[rownames_6, rownames_6])
asv_dist_food6 <- as.dist(mantel_asv_test[rownames_7, rownames_7])

mantel_ht <- mantel(asv_dist_ht,chit.distht, method = "spearman", permutations=999)
mantel_ht
#Mantel statistic r: 0.5821  -1 strong negative, +1 strong positive, 0 none
#Significance: 0.001

mantel_lt <- mantel(asv_dist_lt,chit.distlt, method = "spearman", permutations=999)
mantel_lt
#Mantel statistic r: 0.1767  -1 strong negative, +1 strong positive, 0 none
#Significance: 0.001

mantel_ph3 <- mantel(asv_dist_ph3,chit.distph3, method = "spearman", permutations=999)
mantel_ph3
#Mantel statistic r: 0.3676  -1 strong negative, +1 strong positive, 0 none
#Significance: 0.001

mantel_ph4 <- mantel(asv_dist_ph4,chit.distph4, method = "spearman", permutations=999)
mantel_ph4
#Mantel statistic r: 0.394  -1 strong negative, +1 strong positive, 0 none
#Significance: 0.001

mantel_ph5.6 <- mantel(asv_dist_ph5.6,chit.distph5.6, method = "spearman", permutations=999)
mantel_ph5.6
#Mantel statistic r: 0.3339  -1 strong negative, +1 strong positive, 0 none
#Significance: 0.001

mantel_food3 <- mantel(asv_dist_food3,chit.distfood3, method = "spearman", permutations=999)
mantel_food3
#Mantel statistic r: 0.3178  -1 strong negative, +1 strong positive, 0 none
#Significance: 0.003

mantel_food6 <- mantel(asv_dist_food6,chit.distfood6, method = "spearman", permutations=999)
mantel_food6
#Mantel statistic r: 0.4403  -1 strong negative, +1 strong positive, 0 none
#Significance: 0.001

########mantel prot
mantel_ht <- mantel(asv_dist_ht,prot.distht, method = "spearman", permutations=999)
mantel_ht
#Mantel statistic r: 0.3512  -1 strong negative, +1 strong positive, 0 none
#Significance: 0.001

mantel_lt <- mantel(asv_dist_lt,prot.distlt, method = "spearman", permutations=999)
mantel_lt
#Mantel statistic r: 0.1205  -1 strong negative, +1 strong positive, 0 none
#Significance: 0.001

mantel_ph3 <- mantel(asv_dist_ph3,prot.distph3, method = "spearman", permutations=999)
mantel_ph3
#Mantel statistic r: 0.0835  -1 strong negative, +1 strong positive, 0 none
#Significance: 0.075

mantel_ph4 <- mantel(asv_dist_ph4,prot.distph4, method = "spearman", permutations=999)
mantel_ph4
#Mantel statistic r: 0.03348  -1 strong negative, +1 strong positive, 0 none
#Significance: 0.277

mantel_ph5.6 <- mantel(asv_dist_ph5.6,prot.distph5.6, method = "spearman", permutations=999)
mantel_ph5.6
#Mantel statistic r: 0.1124  -1 strong negative, +1 strong positive, 0 none
#Significance: 0.037

mantel_food3 <- mantel(asv_dist_food3,prot.distfood3, method = "spearman", permutations=999)
mantel_food3
#Mantel statistic r: 0.1516  -1 strong negative, +1 strong positive, 0 none
#Significance: 0.015

mantel_food6 <- mantel(asv_dist_food6,prot.distfood6, method = "spearman", permutations=999)
mantel_food6
#Mantel statistic r: 0.1724  -1 strong negative, +1 strong positive, 0 none
#Significance: 0.005

#### 16S BETA DISPERSION ####
wu.bd.temperature <- betadisper(wu.dist.16s.tubes, tubes.nomix.meta$temperature)
wu.bd.temperature
boxplot(wu.bd.temperature)
anova(wu.bd.temperature)# p=0.2252 , no significant differences in the dispersion between temperatures

wu.bd.food <- betadisper(wu.dist.16s.tubes, tubes.nomix.meta$food)
wu.bd.food
boxplot(wu.bd.food)
anova(wu.bd.food)# p=0.2505, no significant differences in the dispersion between food

wu.bd.ph <- betadisper(wu.dist.16s.tubes, tubes.nomix.meta$ph)
wu.bd.ph
boxplot(wu.bd.ph)
anova(wu.bd.ph)# p=0.01143, significant differences in the dispersion between pH

wu.bd.day <- betadisper(wu.dist.16s.tubes, tubes.nomix.meta$day)
wu.bd.day
boxplot(wu.bd.day)
anova(wu.bd.day)# p=0.0001878 ***, significant differences in the dispersion between time


#### PHYLOGENETIC ANALYSES ####
# Phylogenetic tree of original inoculating community
Exp2.physeqmix <- Exp2.physeq3
sample_data_df <- data.frame(sample_data(Exp2.physeqmix))
mix_indices <- which(sample_data_df$container == "MIX")
new_names <- paste0("mix", seq_along(mix_indices))
sample_data_df$container[mix_indices] <- new_names
sample_data(Exp2.physeqmix) <- sample_data(sample_data_df)

Exp2.physeqmix = subset_samples(Exp2.physeqmix, container %in% c("mix1"))
row_sums <- rowSums(otu_table(Exp2.physeqmix))
nonzero_rows <- row_sums != 0
Exp2.physeqmix <- prune_taxa(nonzero_rows, Exp2.physeqmix)
# 46 taxa and 1 sample

treemix <- phy_tree(Exp2.physeqmix)

tax_table_df <- data.frame(tax_table(Exp2.physeqmix))
length(unique(tax_table_df$Genus))
length(unique(tax_table_df$Family))
length(unique(tax_table_df$Phylum))

tree_data <- fortify(treemix)
tree_data <- merge(tree_data, tax_table_df, by.x = "label", by.y = "row.names", all.x = TRUE)
newasvtaxnames <- read.csv("exported-files/taxonomy_NEW_ASV.csv", header=TRUE)
colnames(tree_data)[1] <- "Feature.ID"
tree_data_asvnames <- dplyr::left_join(tree_data, newasvtaxnames, by="Feature.ID")

## Figure 5A ##
# Family
ggtree(tree_data_asvnames)+
  geom_tiplab(aes(label=ASV,color = Family)) +
  theme_tree2() + scale_color_discrete(name = "Family")+
  geom_cladelab(data = tree_data_asvnames,mapping = aes(node = node, label = Family,color = Family), 
 offset = 1, offset.text=0.5)+theme(legend.position = "none")+
  geom_hilight(data=tree_data_asvnames[1:45,], aes(node=node, fill=Family),
               type = "roundrect")

## Figure 5A ##
# Class
ggtree(tree_data_asvnames)+
  geom_tiplab(aes(label=ASV,color = Class)) +
  theme_tree2() +
  geom_cladelab(data = tree_data_asvnames,mapping = aes(node = node, label = Class,color = Class), 
                offset = 1, offset.text=0.5)+theme(legend.position = "none")+
  geom_hilight(data=tree_data_asvnames[1:45,], aes(node=node, fill=Class),
               type = "roundrect")

##prune to day 57
Exp2.physeq357 = subset_samples(Exp2.physeq3, day %in% c("1"))
row_sums <- rowSums(otu_table(Exp2.physeq357))
nonzero_rows <- row_sums != 0
Exp2.physeq357 <- prune_taxa(nonzero_rows, Exp2.physeq357)

com <- t(as.data.frame(otu_table(Exp2.physeq3)))
tree <- phy_tree(Exp2.physeq3)
Ntip(tree)
meta <- as.data.frame(sample_data(Exp2.physeq3))
meta$ID <- rownames(meta)
all(colnames(com) %in% tree$tip.label)

#### SES PD ####
#SESPD <- ses.pd(com, tree, null.model = "taxa.labels",runs = 999)
#saveRDS(SESPD, "RDS/SES_PD_all.RDS")
SESPD <- readRDS("RDS/ses.pd.resultall.RDS")
SESPD$ID <- rownames(SESPD)
SESPD <- dplyr::left_join(SESPD, meta, by="ID")
SESPD$food <- as.factor(SESPD$food)
SESPD$temperature <- as.factor(SESPD$temperature)
SESPD$ph <- as.factor(SESPD$ph)
SESPD <- SESPD %>%
  mutate(ph = replace(ph, is.na(ph), c(5.6)))

SESPD <- SESPD %>%
  mutate(food = replace(food, is.na(food), c(3)))

SESPD <- SESPD %>%
  mutate(temperature = replace(temperature, is.na(temperature), c(22)))
SESPD$day <- as.numeric(SESPD$day)
ggplot(SESPD, aes(x=ntaxa, y=pd.obs.z, color=temperature, shape=food, group=sample_id))+geom_point(size=3)+
         theme_classic()+scale_color_manual(values = c("#1f5776", "#c83126"))+
  geom_hline(yintercept=0, linetype="dashed", color="black")+
  facet_wrap(~day, nrow=1)+geom_line()+ylim(c(-4, 4))+
  geom_hline(yintercept = -1.96, linetype="dashed", color="gray")+geom_hline(yintercept = 1.96, linetype="dashed", color="gray")


SESPDdaymix <- subset(SESPD, sample_id =="MIX")
mean(SESPDdaymix$pd.obs.z)#-2.015556
sd(SESPDdaymix$pd.obs.z)/sqrt(length(SESPDdaymix$pd.obs.z))# 0.1819155

SESPDday57 <- subset(SESPD, day =="57")

count_greater_than_1.98 <- sum(abs(SESPDday57$pd.obs.z) > 1.98)
total_numbers <- length(SESPDday57$pd.obs.z)
(count_greater_than_1.98 / total_numbers) * 100
#29.16667

stat_pd <- SESPDday57 %>% group_by(treatment_combo) %>% 
  mutate(mean_pd =mean(pd.obs.z))

#mean overall is -1.855

SESPDday57$ph <- factor(SESPDday57$ph, levels=c("5.6", "3", "4"))

### GLM
pd_model <- brm(pd.obs.z ~ temperature * food * ph ,data=SESPDday57, iter = 10000, chains = 4, cores = 4)
saveRDS(pd_model, "RDS/brms_pd_model.RDS")
pd_model <- readRDS("RDS/brms_pd_model.RDS")
pd_model

pdsum<- summary(pd_model)
summary_df <- as.data.frame(pdsum$fixed)
#temp
((summary_df[2,1]) / abs(summary_df[1,1]))*100

#food
((summary_df[3,1]) / abs(summary_df[1,1]))*100

#temp*food
((summary_df[6,1]) / abs(summary_df[1,1]))*100




tab_model(pd_model)
pd.pred <- ggpredict(pd_model, c("food", "temperature", "ph"))
pd.pred$facet <- factor(pd.pred$facet, levels=c("3", "4", "5.6"))
mcmc_plot(pd_model)

## Figure S5A ##
ggplot(pd.pred, aes(x = x, y = predicted, shape = x, color = group)) +
  facet_wrap(~facet) +
  scale_color_manual(values = c("#1f5776", "#c83126")) +
  theme_bw() +scale_shape_manual(values=c(16,2))+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0, size=1,position = position_dodge(width = .7))+
  theme(legend.position="none")+ ylab("SESPD")+
  xlab("Food (g/L)")+ geom_point(size = 3, position = position_dodge(width = .7)) +theme(
    text = element_text(color = "black", size = 10 )
  )+ylim(c(-4, 4))+
  geom_hline(yintercept = 0, linetype="dashed")+
  geom_hline(yintercept = -1.96, linetype="dashed", color="gray")+geom_hline(yintercept = 1.96, linetype="dashed", color="gray")

posteriorpd <- mcmc_intervals_data(pd_model, 
                                     prob_outer=0.95,
                                     prob=0.5)

posteriorpd$nonzero <- NA
posteriorpd$nonzero[posteriorpd$ll>0 & posteriorpd$hh>0] <- "nonzero"
posteriorpd$nonzero[posteriorpd$ll<0 & posteriorpd$hh<0] <- "nonzero"
posteriorpd$nonzero[is.na(posteriorpd$nonzero)] <- "zero"
posteriorpd<- posteriorpd[1:12,]

posteriorpd <- posteriorpd %>%
  mutate(parameter = str_remove(parameter, "^b_") %>%  # Remove "b_" prefix
           str_replace_all(":", ": ") %>%
           str_replace_all("ph3", "pH 3") %>%
           str_replace_all("ph4", "pH 4") %>%
           str_replace_all("food6", "high food") %>%
           str_replace_all("temperature37", "high temp"))

posteriorpd <- posteriorpd %>%
  mutate(parameter = factor(parameter, 
                            levels = c("Intercept", 
                                       "high temp", 
                                       "high food", 
                                       "pH 3", 
                                       "pH 4", 
                                       "high food: pH 3", 
                                       "high food: pH 4", 
                                       "high temp: pH 3", 
                                       "high temp: pH 4", 
                                       "high temp: high food", 
                                       "high temp: high food: pH 3", 
                                       "high temp: high food: pH 4")))

ggplot(posteriorpd, aes(x = parameter, shape=nonzero)) +
  geom_hline(yintercept = 0, linetype = 2, size=.35, color = "black") +
  geom_pointrange(aes(ymin = ll, ymax = hh, y = m), position= position_dodge(width=0.75), size = 1, linewidth = 1) +
  scale_shape_manual(values=c(17, 19), labels=c("95% CI does\nnot contain zero", "95% CI\ncontains zero"))+
  coord_flip() +xlab(NULL) + ylab("Estimated effect on SESPD")+ theme(text = element_text(color = "black", size = 16))+
  theme(panel.background = element_rect(fill = "white", colour = "black"))+theme(legend.position="none")

#### SES MPD ####
phydist <- cophenetic(tree)
#SESMPD <- ses.mpd(com, phydist, null.model = "taxa.labels",
                          abundance.weighted = TRUE, runs = 999)
#saveRDS(SESMPD, "RDS/ses.mpd.resultall.RDS")
SESMPD <- readRDS("RDS/ses.mpd.resultall.RDS")

SESMPD$ID <- rownames(SESMPD)
SESMPD <- dplyr::left_join(SESMPD, meta, by="ID")
SESMPD$food <- as.factor(SESMPD$food)
SESMPD$temperature <- as.factor(SESMPD$temperature)
SESMPD$ph <- as.factor(SESMPD$ph)


SESMPD <- SESMPD %>%
  mutate(ph = replace(ph, is.na(ph), c(5.6)))

SESMPD <- SESMPD %>%
  mutate(food = replace(food, is.na(food), c(3)))

SESMPD <- SESMPD %>%
  mutate(temperature = replace(temperature, is.na(temperature), c(22)))
SESMPD$day <- as.numeric(SESMPD$day)

SESMPDdaymix <- subset(SESMPD, sample_id =="MIX")
mean(SESMPDdaymix$mpd.obs.z)#0.3774597
sd(SESMPDdaymix$mpd.obs.z)/sqrt(length(SESMPDdaymix$mpd.obs.z))# 0.05695824

SESMPDday57 <- subset(SESMPD, day=="57")
SESMPDday57$ph <- factor(SESMPDday57$ph, levels=c("5.6", "3", "4"))

count_greater_than_1.98 <- sum(abs(SESMPDday57$mpd.obs.z) > 1.98)
total_numbers <- length(SESMPDday57$mpd.obs.z)
(count_greater_than_1.98 / total_numbers) * 100

#12.5 percent

stat_mpd <- SESMPDday57 %>% group_by(treatment_combo) %>% 
  mutate(mean_mpd =mean(mpd.obs.z))

#mean overall is -0.80479

#mpd_model <- brm(mpd.obs.z ~ temperature * food * ph ,data=SESMPDday57, iter = 10000, chains = 4, cores = 4)
#saveRDS(mpd_model, "RDS/brms_mpd_model.RDS")
mpd_model <- readRDS("RDS/brms_mpd_model.RDS")
tab_model(mpd_model)
mpdsum<- summary(mpd_model)
mpdsummary_df <- as.data.frame(mpdsum$fixed)
#temp
((mpdsummary_df[2,1]) / abs(mpdsummary_df[1,1]))*100

#food
((mpdsummary_df[3,1]) / abs(mpdsummary_df[1,1]))*100


#ph3
((mpdsummary_df[4,1]) / abs(mpdsummary_df[1,1]))*100


#temp*food
((mpdsummary_df[6,1]) / abs(mpdsummary_df[1,1]))*100
((mpdsummary_df[6,1]+ mpdsummary_df[2,1]+mpdsummary_df[3,1]) / abs(mpdsummary_df[1,1]))*100






mpd.pred <- ggpredict(mpd_model, c("food", "temperature", "ph"))
mpd.pred$facet <- factor(mpd.pred$facet, levels=c("3", "4", "5.6"))

## Figure S5B ##
ggplot(mpd.pred, aes(x = x, y = predicted, shape = x, color = group)) +
  facet_wrap(~facet) +
  scale_color_manual(values = c("#1f5776", "#c83126")) +
  theme_bw() +scale_shape_manual(values=c(16,2))+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0, size=1,position = position_dodge(width = .7))+
  theme(legend.position="none")+ ylab("SES MPD")+
  xlab("Food (g/L)")+ geom_point(size = 3, position = position_dodge(width = .7)) +theme(
    text = element_text(color = "black", size = 10 )
  )+ylim(c(-4, 4))+
  geom_hline(yintercept = 0, linetype="dashed")+
  geom_hline(yintercept = -1.96, linetype="dashed", color="gray")+geom_hline(yintercept = 1.96, linetype="dashed", color="gray")

posteriormpd <- mcmc_intervals_data(mpd_model, 
                                   prob_outer=0.95,
                                   prob=0.5)

posteriormpd$nonzero <- NA
posteriormpd$nonzero[posteriormpd$ll>0 & posteriormpd$hh>0] <- "nonzero"
posteriormpd$nonzero[posteriormpd$ll<0 & posteriormpd$hh<0] <- "nonzero"
posteriormpd$nonzero[is.na(posteriormpd$nonzero)] <- "zero"
posteriormpd<- posteriormpd[1:12,]

posteriormpd <- posteriormpd %>%
  mutate(parameter = str_remove(parameter, "^b_") %>%  # Remove "b_" prefix
           str_replace_all(":", ": ") %>%
           str_replace_all("ph3", "pH 3") %>%
           str_replace_all("ph4", "pH 4") %>%
           str_replace_all("food6", "high food") %>%
           str_replace_all("temperature37", "high temp"))

posteriormpd <- posteriormpd %>%
  mutate(parameter = factor(parameter, 
                            levels = c("Intercept", 
                                       "high temp", 
                                       "high food", 
                                       "pH 3", 
                                       "pH 4", 
                                       "high food: pH 3", 
                                       "high food: pH 4", 
                                       "high temp: pH 3", 
                                       "high temp: pH 4", 
                                       "high temp: high food", 
                                       "high temp: high food: pH 3", 
                                       "high temp: high food: pH 4")))

ggplot(posteriormpd, aes(x = parameter, shape=nonzero)) +
  geom_hline(yintercept = 0, linetype = 2, size=.35, color = "black") +
  geom_pointrange(aes(ymin = ll, ymax = hh, y = m), position= position_dodge(width=0.75), size = 1, linewidth = 1) +
  scale_shape_manual(values=c(17, 19), labels=c("95% CI does\nnot contain zero", "95% CI\ncontains zero"))+
  coord_flip() +xlab(NULL) + ylab("Estimated effect on SESMPD")+ theme(text = element_text(color = "black", size = 16))+
  theme(panel.background = element_rect(fill = "white", colour = "black"))+theme(legend.position="none")

#### SES MNTD ####

#SESMNTD <- ses.mntd(com, phydist, null.model = "taxa.labels", abundance.weighted = TRUE, runs = 999)
#saveRDS(SESMNTD, "RDS/ses.mntd.resultall.RDS")
SESMNTD <- readRDS("RDS/ses.mntd.resultall.RDS")

SESMNTD$ID <- rownames(SESMNTD)
SESMNTD <- dplyr::left_join(SESMNTD, meta, by="ID")
SESMNTD$food <- as.factor(SESMNTD$food)
SESMNTD$temperature <- as.factor(SESMNTD$temperature)
SESMNTD$ph <- as.factor(SESMNTD$ph)

SESMNTD <- SESMNTD %>%
  mutate(ph = replace(ph, is.na(ph), c(5.6)))

SESMNTD <- SESMNTD %>%
  mutate(food = replace(food, is.na(food), c(3)))

SESMNTD <- SESMNTD %>%
  mutate(temperature = replace(temperature, is.na(temperature), c(22)))
SESMNTD$day <- as.numeric(SESMNTD$day)

SESMNTDdaymix <- subset(SESMNTD, sample_id =="MIX")
mean(SESMNTDdaymix$mntd.obs.z)#-1.347052
sd(SESMNTDdaymix$mntd.obs.z)/sqrt(length(SESMNTDdaymix$mntd.obs.z))# 0.0851035

SESMNTDday57 <- subset(SESMNTD, day=="57")

SESMNTDday57$ph <- factor(SESMNTDday57$ph, levels=c("5.6", "3", "4"))

stat_mntd <- SESMNTDday57 %>% group_by(treatment_combo) %>% 
  mutate(mean_mntd =mean(mntd.obs.z))

#mean overall is -0.9528

count_greater_than_1.98 <- sum(abs(SESMNTDday57$mntd.obs.z) > 1.98)
total_numbers <- length(SESMNTDday57$mntd.obs.z)
(count_greater_than_1.98 / total_numbers) * 100

#0


#mntd_model <- brm(mntd.obs.z ~ temperature * food * ph,data=SESMNTDday57, iter = 10000, chains = 4, cores = 4)
#saveRDS(mntd_model, "RDS/brms_mntd_model.RDS")
mntd_model <- readRDS("RDS/brms_mntd_model.RDS")
tab_model(mntd_model)
mntdsum<- summary(mntd_model)
mntdsummary_df <- as.data.frame(mntdsum$fixed)
#temp
((mntdsummary_df[2,1]) / abs(mntdsummary_df[1,1]))*100

#ph3
((mntdsummary_df[4,1]) / abs(mntdsummary_df[1,1]))*100

#ph4
((mntdsummary_df[5,1]) / abs(mntdsummary_df[1,1]))*100

#temp*ph3
((mntdsummary_df[7,1]) / abs(mntdsummary_df[1,1]))*100
((mntdsummary_df[7,1]+ mntdsummary_df[2,1]+ mntdsummary_df[4,1]) / abs(mntdsummary_df[1,1]))*100

#temp*food*ph3
((mntdsummary_df[11,1]) / abs(mntdsummary_df[1,1]))*100

((mntdsummary_df[11,1]+ mntdsummary_df[2,1]+mntdsummary_df[3,1] + mntdsummary_df[4,1]) / abs(mntdsummary_df[1,1]))*100




plot(conditional_effects(mntd_model))
plot_model(mntd_model,
           type = "int")




mntd.pred <- predict_response(mntd_model, c("food", "temperature", "ph"))
mntd.pred$facet <- factor(mntd.pred$facet, levels=c("3", "4", "5.6"))

## Figure S5C ##
ggplot(mntd.pred, aes(x = x, y = predicted, shape = x, color = group)) +
  facet_wrap(~facet) +
  scale_color_manual(values = c("#1f5776", "#c83126")) +
  theme_bw() +scale_shape_manual(values=c(16,2))+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0, size=1,position = position_dodge(width = .7))+
  theme(legend.position="none")+ ylab("SES MNTD")+
  xlab("Food (g/L)")+ geom_point(size = 3, position = position_dodge(width = .7)) +theme(
    text = element_text(color = "black", size = 10 )
  )+ylim(c(-4, 4))+
  geom_hline(yintercept = 0, linetype="dashed")+
  geom_hline(yintercept = -1.96, linetype="dashed", color="gray")+geom_hline(yintercept = 1.96, linetype="dashed", color="gray")

posteriormntd <- mcmc_intervals_data(mntd_model, 
                                    prob_outer=0.95,
                                    prob=0.5)

posteriormntd$nonzero <- NA
posteriormntd$nonzero[posteriormntd$ll>0 & posteriormntd$hh>0] <- "nonzero"
posteriormntd$nonzero[posteriormntd$ll<0 & posteriormntd$hh<0] <- "nonzero"
posteriormntd$nonzero[is.na(posteriormntd$nonzero)] <- "zero"
posteriormntd<- posteriormntd[1:12,]

posteriormntd <- posteriormntd %>%
  mutate(parameter = str_remove(parameter, "^b_") %>%  # Remove "b_" prefix
           str_replace_all(":", ": ") %>%
           str_replace_all("ph3", "pH 3") %>%
           str_replace_all("ph4", "pH 4") %>%
           str_replace_all("food6", "high food") %>%
           str_replace_all("temperature37", "high temp"))

posteriormntd <- posteriormntd %>%
  mutate(parameter = factor(parameter, 
                            levels = c("Intercept", 
                                       "high temp", 
                                       "high food", 
                                       "pH 3", 
                                       "pH 4", 
                                       "high food: pH 3", 
                                       "high food: pH 4", 
                                       "high temp: pH 3", 
                                       "high temp: pH 4", 
                                       "high temp: high food", 
                                       "high temp: high food: pH 3", 
                                       "high temp: high food: pH 4")))

ggplot(posteriormntd, aes(x = parameter, shape=nonzero)) +
  geom_hline(yintercept = 0, linetype = 2, size=.35, color = "black") +
  geom_pointrange(aes(ymin = ll, ymax = hh, y = m), position= position_dodge(width=0.75), size = 1, linewidth = 1) +
  scale_shape_manual(values=c(17, 19), labels=c("95% CI does\nnot contain zero", "95% CI\ncontains zero"))+
  coord_flip() +xlab(NULL) + ylab("Estimated effect on SESMNTD")+ theme(text = element_text(color = "black", size = 16))+
  theme(panel.background = element_rect(fill = "white", colour = "black"))+theme(legend.position="none")

#Negative values of SESmpd and SESmntd indicate lower phylogenetic diversity than expected 
# under the assumption of the null model, whereas values greater than zero indicate higher 
# phylogenetic diversity than predicted by the null model.


posterior_m1a <- as.data.frame(mntd_model)
posterior_m1a_melt <- posterior_m1a[,1:12]
posterior_m1a_melt <- reshape2::melt(posterior_m1a_melt)
posterior_m1a_melt$model <- "mntd"

posterior_m2a <- as.data.frame(pd_model)
posterior_m2a_melt <- posterior_m2a[,1:12]
posterior_m2a_melt <- reshape2::melt(posterior_m2a_melt)
posterior_m2a_melt$model <- "pd"

posterior_m3a <- as.data.frame(mpd_model)
posterior_m3a_melt <- posterior_m3a[,1:12]
posterior_m3a_melt <- reshape2::melt(posterior_m3a_melt)
posterior_m3a_melt$model <- "mpd"

posterior_all <- rbind(posterior_m2a_melt, posterior_m3a_melt, posterior_m1a_melt)

ggplot(posterior_m2a_melt, aes(x = value,y=variable, fill = model)) +
  stat_halfeye(slab_alpha=0.50, .width=c(.95, 0.95)) +
  geom_vline(aes(xintercept=0), 
             color="black", size=1, linetype="dashed")+
  ylab("Probability density") +
  theme_classic() + guides(fill="none")+ 
  theme(text = element_text(color = "black", size = 16 ))+
  xlab("Predicted Effect")

ggplot(posterior_m3a_melt, aes(x = value,y=variable, fill = model)) +
  stat_halfeye(slab_alpha=0.50, .width=c(.95, 0.95)) +
  geom_vline(aes(xintercept=0), 
             color="black", size=1, linetype="dashed")+
  ylab("Probability density") +
  theme_classic() + guides(fill="none")+ 
  theme(text = element_text(color = "black", size = 16 ))+
  xlab("Predicted Effect")

ggplot(posterior_m1a_melt, aes(x = value,y=variable, fill = model)) +
  stat_halfeye(slab_alpha=0.50, .width=c(.95, 0.95)) +
  geom_vline(aes(xintercept=0), 
             color="black", size=1, linetype="dashed")+
  ylab("Probability density") +
  theme_classic() + guides(fill="none")+ 
  theme(text = element_text(color = "black", size = 16 ))+
  xlab("Predicted Effect")

posterior_all$model <- factor(posterior_all$model, levels=c("pd", "mpd", "mntd"))
posterior_all_filt <- posterior_all
posterior_all_filt <- posterior_all_filt %>% filter(variable != "b_Intercept")

posterior_all_filt <- posterior_all_filt %>%
  mutate(variable = str_remove(variable, "^b_") %>%  # Remove "b_" prefix
           str_replace_all(":", ": ") %>%
           str_replace_all("ph3", "pH 3") %>%
           str_replace_all("ph4", "pH 4") %>%
           str_replace_all("food6", "high food") %>%
           str_replace_all("temperature37", "high temp"))

posterior_all_filt <- posterior_all_filt %>%
  mutate(variable = factor(variable, 
                            levels = c("high temp", 
                                       "high food", 
                                       "pH 3", 
                                       "pH 4", 
                                       "high temp: pH 3", 
                                       "high temp: pH 4", 
                                       "high food: pH 3", 
                                       "high food: pH 4", 
                                       "high temp: high food", 
                                       "high temp: high food: pH 3", 
                                       "high temp: high food: pH 4")))



ggplot(posterior_all_filt, aes(x = value,y=variable, fill = model)) +
  stat_halfeye(slab_alpha=0.70, .width=c(.95, 0.95), position = "dodge") +
  geom_vline(aes(xintercept=0), 
             color="black", size=1, linetype="dashed")+
  ylab("Probability density") +
  theme_classic() + guides(fill="none")+ 
  theme(text = element_text(color = "black"))+
  xlab("Predicted Effect")+facet_wrap(~model, scales="free_x")+
scale_fill_brewer(palette = "Set2") +
  scale_color_brewer(palette = "Dark2")

posterior_m2a_melt2 <- posterior_m2a_melt %>% mutate(variable = paste0(variable, "A"))
posterior_m3a_melt2 <- posterior_m3a_melt %>% mutate(variable = paste0(variable, "B"))
posterior_m1a_melt2 <- posterior_m1a_melt %>% mutate(variable = paste0(variable, "C"))

## Figure 5B ##
posterior_all2 <- rbind(posterior_m2a_melt2, posterior_m3a_melt2, posterior_m1a_melt2)
ggplot(posterior_all2, aes(x = value,y=variable, fill = model)) +
  stat_halfeye(slab_alpha=0.50, .width=c(.95, 0.95)) +
  geom_vline(aes(xintercept=0), 
             color="black", size=1, linetype="dashed")+
  ylab("Probability density") +
  theme_classic() + guides(fill="none")+ 
  theme(text = element_text(color = "black", size = 16 ))+
  xlab("Predicted Effect")



#### Beta Diversity of Original Inoculating Communities ####
## Calculate weighted Unifrac distance and run NMDS
wu.dist.16s.mix <- phyloseq::distance(Exp2.tubes.mix.physeq3, method="wunifrac")
set.seed(123)
#wu.nmds.16s.mix <- metaMDS(wu.dist.16s.mix,
#                             k = 2, 
#                             trymax = 1000,
#                             wascores = TRUE)#0.09442731

#saveRDS(wu.nmds.16s.mix, file = "RDS/Exp2.wu.nmds.16s.mix.RDS")
wu.nmds.16s.mix <- readRDS("RDS/Exp2.wu.nmds.16s.mix.RDS")

## Plot NMDS
data.scores2 <- as.data.frame(scores(wu.nmds.16s.mix$points[,1:2]))
data.scores2$ID <- rownames(data.scores2)
meta_mix <- subset(md16s, day %in% c("0", "1"))
data_merge2 <- merge(data.scores2, meta_mix, by = c("ID"))

data_merge2$day <- as.factor(data_merge2$day)

## Figure S1B ##
ggplot(data_merge2, aes(x = MDS1, y = MDS2, color=day)) + 
  geom_point(size = 4)+
  labs(x = "NMDS1", y = "NMDS2")+scale_color_manual(values=c("black", "gray"))+
  theme_bw()

#### ANCOMBC DIFFERENTIAL ABUNDANCE #### 
#use non-rarified phyloseq object
#TUBES ONLY
tse = mia::convertFromPhyloseq(Exp2.physeq2.nomix)
tse$treatment_combo = factor(tse$treatment_combo, levels = c("temp22_ph5.6_food3","temp22_ph5.6_food6", "temp22_ph4_food3",
                                                             "temp22_ph4_food6", "temp22_ph3_food3", "temp22_ph3_food6", "temp37_ph5.6_food3",
                                                             "temp37_ph5.6_food6", "temp37_ph4_food3", "temp37_ph4_food6", "temp37_ph3_food3", "temp37_ph3_food6"))
tse$ph <- factor(tse$ph, levels=c("5.6", "4", "3"))
tse$temperature <- factor(tse$temperature, levels=c("22", "37"))
tse$food <- factor(tse$food, levels=c("3", "6"))
tse$sample_id <- as.factor(tse$sample_id)

#read in new asv names
newasv <- read.csv("Data/ANCOMBC_NEWNAMES.csv", header=TRUE)
newasv$newname2 <- paste(newasv$new_name, " - ", newasv$Genus)

#### differentially abundant taxa across chitinase activity
#output_chit = ancombc2(data = tse, assay_name = "counts", tax_level = NULL,
#                       fix_formula = "scaled_chit+scaled_prot",
#                       rand_formula = "(1 | day)",#cant use sample_id as random effect bc not enough observations per each factor (48)
#                       p_adj_method = "fdr", pseudo_sens = TRUE,
#                       prv_cut = 0.02,#4 tubes in 12 treatment groups at three time points, 144 tubes, so set prvcut to 2% meaning in at least 4 tubes
#                       alpha = 0.05, n_cl = 3, verbose = TRUE,
#                       global = FALSE, pairwise = FALSE)
#saveRDS(output_chit, "RDS/ancombc_chit.RDS")
output_chit <- readRDS("RDS/ancombc_chit.RDS")

res_prim_chit <- output_chit$res
chit_res_sig <- subset(res_prim_chit, diff_scaled_chit == "TRUE")
prot_res_sig <- subset(res_prim_chit, diff_scaled_prot == "TRUE")
dim(res_prim_chit) #69 taxa present in at least 2% of the samples
dim(chit_res_sig)#10 taxa differentially abundant
dim(prot_res_sig)#8 22

chit_res_sig_filt <- chit_res_sig[,c(1,3,5)]
colnames(chit_res_sig_filt) <- c("ASV", "lfc", "se")

chit_res_sig_filt_melt_names <- dplyr::left_join(chit_res_sig_filt, newasv, by="ASV")

chit_res_sig_filt_melt_names <- chit_res_sig_filt_melt_names %>%
  arrange(lfc)
taxon_levels <- unique(chit_res_sig_filt_melt_names$newname2)
chit_res_sig_filt_melt_names$newname2 <- factor(chit_res_sig_filt_melt_names$newname2, levels = taxon_levels)

prot_res_sig_filt <- prot_res_sig[,c(1,3,5)]
colnames(prot_res_sig_filt) <- c("ASV", "lfc", "se")

prot_res_sig_filt_melt_names <- dplyr::left_join(prot_res_sig_filt, newasv, by="ASV")

prot_res_sig_filt_melt_names <- prot_res_sig_filt_melt_names %>%
  arrange(lfc)
taxon_levels <- unique(prot_res_sig_filt_melt_names$newname2)
prot_res_sig_filt_melt_names$newname2 <- factor(prot_res_sig_filt_melt_names$newname2, levels = taxon_levels)

#combine them
chit_res_sig_filt_melt_names$df <- "chitinase"
prot_res_sig_filt_melt_names$df <- "protease"
df <- rbind(chit_res_sig_filt_melt_names, prot_res_sig_filt_melt_names)

ch_anc <- chit_res_sig_filt_melt_names
ch_anc$enzyme <- "chitinase"
pr_anc <- prot_res_sig_filt_melt_names
pr_anc$enzyme <- "protease"
enz_anc <- rbind(ch_anc, pr_anc)

## Figure 6E ##
ggplot(enz_anc, aes(x = enzyme, y = newname2, fill = lfc)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "#88B04B",mid= "#FFF8DC",high = "#FFC107",limits = c(-2, 2))+
  theme_minimal() +
  theme(text = element_text(color = "black")) + xlab("") + ylab("")

write.csv(df, "output_files/ANCOMBC_enzyme.csv", row.names=TRUE)

#### differentially abundant taxa between treatments
#output_temp = ancombc2(data = tse, assay_name = "counts", tax_level = NULL,
#                       fix_formula = "temperature + food + ph",
#                       rand_formula = "(1 | day)",
#                       p_adj_method = "fdr", pseudo_sens = TRUE,
#                       prv_cut = 0.02,#4 tubes in 12 treatment groups at three time points, 144 tubes, so set prvcut to 2% meaning in at least 4 tubes
#                      group = "temperature",
#                       alpha = 0.05, n_cl = 3, verbose = TRUE,
#                       global = TRUE, pairwise = TRUE)

#saveRDS(output_temp, "RDS/ancombc.RDS")
output_temp <- readRDS("RDS/ancombc.RDS")

food_prim = output_temp$res
all_res_sig <- food_prim[apply(food_prim[, 27:31], 1, function(x) any(grepl("TRUE", x))), ]
food_res_sig <- subset(food_prim, diff_food6 == "TRUE")
ph_res_sig <- food_prim[apply(food_prim[, 29:30], 1, function(x) any(grepl("TRUE", x))), ]
temp_res_sig <- subset(food_prim, diff_temperature37 == "TRUE")


dim(food_prim) #69 taxa present in 2% samples
dim(all_res_sig)#30 differentially abundant taxa across all treatments
dim(food_res_sig) #3 da taxa food
dim(ph_res_sig) #11 da taxa ph
dim(temp_res_sig) #23 da taxa temp


##all
all.filt <- as.data.frame(all_res_sig[,1:6])
colnames(all.filt)[1] <- "ASV"
all.filt.name <- dplyr::left_join(all.filt, newasv, by="ASV")
all.filt.name.filt <- all.filt.name[,c(3,4,5,6,9)]

all.filt.melt <- reshape2::melt(all.filt.name.filt, id.vars = "newname2")

##temp
temp.filt <- as.data.frame(temp_res_sig[,c(1,6)])
colnames(temp.filt)[1] <- "ASV"
temp.filt.name <- dplyr::left_join(temp.filt, newasv, by="ASV")
temp.filt.name.filt <- temp.filt.name[,c(2,5)]

temp.filt.melt <- reshape2::melt(temp.filt.name.filt, id.vars = "newname2")
temp.filt.melt$newname2 <- factor(temp.filt.melt$newname2, levels = unique(temp.filt.melt$newname2[order(temp.filt.melt$value, decreasing = TRUE)]))

## Figure 4F ##
ggplot(temp.filt.melt, aes(x = variable, y = newname2, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "#063648", mid="white",high = "#A33232")  +
  theme_minimal()+theme(text = element_text(color = "black", size = 10 ))+
  xlab("")+ ylab("")

##pH
ph.filt <- as.data.frame(ph_res_sig[,c(1,5)])
colnames(ph.filt)[1] <- "ASV"
ph.filt.name <- dplyr::left_join(ph.filt, newasv, by="ASV")
ph.filt.name.filt <- ph.filt.name[,c(2,5)]

ph.filt.melt <- reshape2::melt(ph.filt.name.filt, id.vars = "newname2")
ph.filt.melt$newname2 <- factor(ph.filt.melt$newname2, levels = unique(ph.filt.melt$newname2[order(ph.filt.melt$value, decreasing = TRUE)]))

## Figure 4F ##
ggplot(ph.filt.melt, aes(x = variable, y = newname2, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(high = "#CA64A3", mid="white",low = "#301934") +
  labs(title = "Heatmap of Taxa vs Variables", x = "Variables", y = "Taxa") +
  theme_minimal()+theme(text = element_text(color = "black", size = 10 ))+
  xlab("")+ ylab("")

##food
food.filt <- as.data.frame(food_res_sig[,c(1,3)])
colnames(food.filt)[1] <- "ASV"
food.filt.name <- dplyr::left_join(food.filt, newasv, by="ASV")
food.filt.name.filt <- food.filt.name[,c(2,5)]

food.filt.melt <- reshape2::melt(food.filt.name.filt, id.vars = "newname2")
food.filt.melt$newname2 <- factor(food.filt.melt$newname2, levels = unique(food.filt.melt$newname2[order(food.filt.melt$value, decreasing = TRUE)]))

## Figure 4F ##
ggplot(food.filt.melt, aes(x = variable, y = newname2, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(high = "#228822", mid="white",low = "#90ee90") +
  labs(title = "Heatmap of Taxa vs Variables", x = "Variables", y = "Taxa") +
  theme_minimal()+theme(text = element_text(color = "black", size = 16 ))+
  xlab("")+ ylab("")

#### RELATIVE ABUNDANCE PLOTS TUBES ####
newasvRA<- read.csv("exported-files/taxonomy_NEW_ASV.csv", header=TRUE)
newasv <- newasvRA$ASV
newasvRA <- parse_taxonomy(newasvRA)
newasvRA <- cbind(newasvRA,newasv )
newasvRA$newname2 <- paste(newasvRA$newasv, " - ", newasvRA$Genus)
otu <- data.frame(otu_table(Exp2.tubes.nomix.physeq3))
rownames(otu)==rownames(newasvRA)
#991 vs 3342
newasvRA_subset <- newasvRA[rownames(otu), , drop = FALSE]
identical(rownames(otu), rownames(newasvRA_subset))
rownames(otu)==rownames(newasvRA_subset)
otu$ASV <- rownames(otu)
newasvRA_subset$ASV <- rownames(newasvRA_subset)

otu.g <- dplyr::left_join(otu, newasvRA_subset, by="ASV")
otu.ga <- aggregate(. ~ otu.g$newname2, otu.g[,1:144], sum)
row.names(otu.ga) <- otu.ga[,1]
otu.ga <- otu.ga[,1:144]
otu.ga <- otu.ga[,-1]
otu.gao <- as.matrix(otu.ga[order(rowSums(otu.ga),decreasing = T),])
otu.gaop <- otu.gao[c(1:30),]#top 30 ASVs 
other <- colSums(otu.gao[c(31:nrow(otu.gao)),])#put the all the other genera in a "other" category
otu.gaopo <- rbind(otu.gaop, other)
otu.gaopo_gg <- as.data.frame(otu.gaopo)
otu.gaopo_ggp <- adorn_percentages(otu.gaopo_gg,, 1:ncol(otu.gaopo_gg), denominator="col")#number of columns (samples), calculate RA
other <- otu.gaopo_ggp[31,]
otu.gaopo_ggp <- otu.gaopo_ggp[-31,]
otu.gaopo_ggp <- otu.gaopo_ggp[order(rowSums(otu.gaopo_ggp),decreasing = F),]
otu.gaopo_ggp <- rbind(other, otu.gaopo_ggp)
otu.gaopo_ggp[ "ASV" ] <- rownames(otu.gaopo_ggp)
otu.gaopo_ggp$ASV <- factor(otu.gaopo_ggp$ASV, levels = otu.gaopo_ggp$ASV)
otu.gaopo_ggp.rmelt<- reshape2::melt(otu.gaopo_ggp, id.vars="ASV", value.name="Relative_Abundance", variable.name="Sample")

samp <- data.frame(sample_data(Exp2.tubes.nomix.physeq3))
samp.new <- samp
samp.new$Sample <- rownames(samp.new)
otu.gaopo_ggp.rmelt.meta <- dplyr::left_join(otu.gaopo_ggp.rmelt, samp.new, by="Sample")

otu.gaopo_ggp.rmelt.meta$day <- as.factor(otu.gaopo_ggp.rmelt.meta$day)
otu.gaopo_ggp.rmelt.meta <- otu.gaopo_ggp.rmelt.meta %>%
  dplyr::arrange(treatment_combo, day)


sample_order <- unique(otu.gaopo_ggp.rmelt.meta$Sample[order(otu.gaopo_ggp.rmelt.meta$day)])
otu.gaopo_ggp.rmelt.meta$Sample <- factor(otu.gaopo_ggp.rmelt.meta$Sample, levels = sample_order)
customcol28 <- colors <- c("black","royalblue3","darkblue","#DDDD77","#117744",
                           "#AAAA44","darkred","palegoldenrod","purple","mediumblue",
                           "lightgoldenrod","#77CCCC","yellow","purple4","darkgreen","#AA4455",
                           "lightsalmon","yellow3","purple2","lightblue","firebrick",
                           "navy","red4","skyblue","darkmagenta","mediumvioletred",
                           "#771155", "#44AAAA", "#CC99BB", "#114477",  "#777711")

## Figure S4 ##
ggplot(otu.gaopo_ggp.rmelt.meta, aes(x=Sample, y=Relative_Abundance, fill = ASV)) + 
  geom_bar( position = "fill", stat = "identity", width=1) + 
  theme_classic()+ theme(axis.text.x = element_text(angle = 90, size=5))+
  guides(fill = guide_legend(ncol = 5)) + theme(legend.position="bottom") +
  facet_wrap(~treatment_combo, scales = "free_x", nrow = 2)+ scale_fill_manual(values = customcol28)

 #### Top 16S ASVs, with taxonomy #####
asv16s.top <- data.frame(sort(colSums(asv16s.rt),decreasing = T))
names(asv16s.top) <- "Sequences"
asv16s.top.tax <- subset(tax.16s2.r, row.names(tax.16s2.r) %in% row.names(asv16s.top))
asv16s.top.tax <- asv16s.top.tax[match(rownames(asv16s.top), rownames(asv16s.top.tax)), ]
asv16s.top <- data.frame(asv16s.top,asv16s.top.tax)

Feature.ID <- rownames(asv16s.top)
asv16s.top.parse <- cbind(asv16s.top, Feature.ID)
asv16s.top.parse <- parse_taxonomy(asv16s.top.parse)

length(unique(asv16s.top.parse$Genus))
length(unique(asv16s.top.parse$Family))
length(unique(asv16s.top.parse$Phylum))

####look at top 20 asv abundance taxonomy in 16S
sum(asv16s.top$Sequences)
asv16s.top.20 <- asv16s.top[1:20,]
sum(asv16s.top.20$Sequences)

Feature.ID <- rownames(asv16s.top.20)
asv16s.top.20 <- cbind(asv16s.top.20, Feature.ID)
asv16s.top.20.tax <- parse_taxonomy(asv16s.top.20)
rownames(asv16s.top.20.tax) == rownames(asv16s.top.20)
asv16s.top.20.tax.abund <- cbind(asv16s.top.20.tax, asv16s.top.20)

#look at top 16S genera
asv.16S.mv <- subset(asv16s.top.20.tax.abund, Genus == "Delftia")
((sum(asv.16S.mv$Sequences)/1121010)*100)#11.51132
((sum(asv.16S.mv$Sequences)/1383702)*100)#9.325924

asv.16S.st <- subset(asv16s.top.20.tax.abund, Genus == "Elizabethkingia")
((sum(asv.16S.st$Sequences)/1121010)*100)#8.390915
((sum(asv.16S.st$Sequences)/1383702)*100)#6.797923

asv.16S.pd <- subset(asv16s.top.20.tax.abund, Genus == "Undibacterium")
((sum(asv.16S.pd$Sequences)/287729)*100)
((sum(asv.16S.pd$Sequences)/434238)*100)

### COMMUNITY FUNCTIONAL ANALYSIS ####

#### ECOPLATE ANALYSIS ####
#read in the data
final_raw <- read.csv("data/2022_Ecoplates_EXP2_end.csv", header = TRUE, row.names = 1, check.names = FALSE)
start_raw <- read.csv("data/2022_Ecoplates_Exp2_start.csv", header = TRUE, row.names = 1, check.names = FALSE)

#meta data
meta_wk8 <- read.csv("data/2022_Ecoplates_EXP2_end_META.csv", header = TRUE, check.names = FALSE)
meta_wk1 <- read.csv("data/2022_Ecoplates_EXP2_start_META.csv", header = TRUE, check.names = FALSE)

#transpose dataframe
final <-  t(final_raw)
start <- t(start_raw)

#convert all negative numbers to zeros
final[final < 0] <- 0
start[start < 0] <- 0

#combine with metadata then subset
meta_wk8$sample_id <- as.character(meta_wk8$sample_id)
final.df <- as.data.frame(final)
final.df$sample_id <- rownames(final.df)
final.df <- dplyr::left_join(final.df, meta_wk8, by= "sample_id")

meta_wk1$sample_id <- as.character(meta_wk1$sample_id)
start.df <- as.data.frame(start)
start.df$sample_id <- rownames(start.df)
start.df <- dplyr::left_join(start.df, meta_wk1, by= "sample_id")

#combine start and final
eco_all <- rbind(start.df, final.df)
eco.tubes.meta <- eco_all[,32:37]

#remove metadata
eco.all.tubes.nometa <- subset(eco_all, select = -c(sample_id, ph, temperature, food, container, wk))
eco.all.tubes.melt <- melt(eco_all, id.vars=c("sample_id", "ph", "temperature", "food", "container", "wk"))
eco.all.tubes.melt$value <- ifelse(eco.all.tubes.melt$value <= 0, 0.001, eco.all.tubes.melt$value)

eco.tubes.meta$sample_id <- as.numeric(eco_all$sample_id)
eco.tubes.meta$ph <- factor(eco.tubes.meta$ph, levels = c(5.6, 4, 3))
eco.tubes.meta$food <- factor(eco.tubes.meta$food, levels = c(3, 6))
eco.tubes.meta$temperature <- factor(eco.tubes.meta$temperature, levels = c(22,37))
eco.tubes.meta$wk <- factor(eco.tubes.meta$wk, levels = c(1,8))

dbrda_eco <- dbrda(eco.all.tubes.nometa ~ temperature*ph*food+ wk, Condition=sample_id,
                  data = eco.tubes.meta,
                  distance = "bray")

eco_sum_dbrda <- summary(dbrda_eco)
plot(dbrda_eco,cex=1)
#Partitioning of squared Bray distance:
#Inertia Proportion
#Total           4.269     1.0000
#Constrained     2.458     0.5757
#Unconstrained   1.811     0.4243

eco_dbrda_anova <- anova.cca(dbrda_eco, by = "onedf", perm = 999)
eco_dbrda_anova
#Df SumOfSqs       F Pr(>F)    
#temperature37            1  0.13102  6.0045  0.001 ***
#ph4                      1  0.18930  8.6753  0.001 ***
#ph3                      1  0.55226 25.3089  0.001 ***
#food6                    1  0.04634  2.1235  0.087 .  
#wk8                      1  1.29816 59.4922  0.001 ***
#temperature37:ph4        1  0.03994  1.8303  0.130    
#temperature37:ph3        1  0.11000  5.0409  0.002 ** 
#temperature37:food6      1  0.03229  1.4798  0.213    
#ph4:food6                1  0.03557  1.6301  0.168    
#ph3:food6                1  0.02775  1.2716  0.256    
#temperature37:ph4:food6  1 -0.00101 -0.0464  1.000    
#temperature37:ph3:food6  1 -0.00388 -0.1776  1.000    
#Residual                83  1.81112  

write.csv(eco_dbrda_anova, "output_files/ecoplate_tubes_dbrda_cca_aov.csv", row.names=TRUE)
significant_factors <- c("temperature37","ph3", "ph4","wk8", "temperature37:ph3")

eco_est <- as.data.frame(cbind(x1 = dbrda_eco$CCA$biplot[,1], y1 = dbrda_eco$CCA$biplot[,2]))
eco_est$factor <- row.names(dbrda_eco$CCA$biplot)
eco_est$significant <- eco_est$factor %in% significant_factors
significant_eco_est <- eco_est[eco_est$significant, ]

sites_eco <- as.data.frame(eco_sum_dbrda$sites) 
sites_eco_meta <- cbind(sites_eco,eco.tubes.meta )


#install.packages("RVAideMemoire")
library(RVAideMemoire)
pairwise_results <-pairwise.perm.manova(sites_eco_meta[,1:2],eco.tubes.meta$ph,nperm=999)
pairwise_results
#  5.6      4     
#4 0.4580.  -     
#3 0.0015   0.0015


sites_eco_meta$wk <- as.factor(sites_eco_meta$wk)
sites_eco_meta$temperature <- as.factor(sites_eco_meta$temperature)
sites_eco_meta$ph <- factor(sites_eco_meta$ph, levels=c(3, 4, 5.6))

## Figure 3A ##
ggplot(data = sites_eco_meta, aes(x = dbRDA1, y = dbRDA2, color = temperature, alpha=wk, shape=wk)) +
  geom_point(size = 3) +scale_shape_manual(values = c(19, 15)) +
  scale_color_manual(values = c("#1f5776", "#c83126")) +
  geom_segment(data = significant_eco_est, aes(x = 0, y = 0, xend = x1, yend = y1),
               arrow = arrow(length = unit(.2, "cm")), color = "black", size = 1, inherit.aes = FALSE) +
  geom_text(data = significant_eco_est, aes(x = x1, y = y1, label = factor), 
            hjust = 0, vjust = 0, inherit.aes = FALSE)+ theme_classic()+scale_alpha_discrete(range = c(1,.7))+
  theme(legend.title=element_blank())

## Figure 3B ##
ggplot(data = sites_eco_meta, aes(x = dbRDA1, y = dbRDA2, color = ph, alpha=wk,shape=wk)) +
  geom_point(size = 3) +scale_shape_manual(values = c(19, 15)) +
  scale_color_manual(values = c("#CA64A3","#9A4EAE", "#301934")) +
  geom_segment(data = significant_eco_est, aes(x = 0, y = 0, xend = x1, yend = y1),
               arrow = arrow(length = unit(.2, "cm")), color = "black", size = 1, inherit.aes = FALSE) +
  geom_text(data = significant_eco_est, aes(x = x1, y = y1, label = factor), 
            hjust = 0, vjust = 0, inherit.aes = FALSE)+ theme_classic()+scale_alpha_discrete(range = c(1,.7)) +
  theme(legend.title=element_blank())

#### PROCRUSTES ####
dim(eco_all)#96
sample_ids <- unique(eco_all$sample_id)
length(sample_ids) #48

meta_pro <- subset(meta)
meta_pro <- subset(meta_pro, day %in% c("1", "57"))

meta_pro_filt <- meta_pro
meta_pro_filt$id <- rownames(meta_pro_filt)
meta_pro_filt$wk <- NA
meta_pro_filt$wk[meta_pro_filt$day == 1] <- 1  
meta_pro_filt$wk[meta_pro_filt$day == 57] <- 8 
meta_pro_filt <- meta_pro_filt[,c(2, 19, 20)]
eco_pro_all <- dplyr::left_join(eco_all, meta_pro_filt, by=c("sample_id", "wk"))

rownames(eco_pro_all) <- eco_pro_all$id
eco_pro_filt <- eco_pro_all[,c(1:31)]

#format ASV data
asv_pro <- asv16s.rt
asv_pro_filt <- subset(asv_pro, row.names(asv_pro) %in% rownames(meta_pro))
dim(asv_pro_filt) #96
rownames(asv_pro_filt) == rownames(meta_pro) #not in order

#fix order
rownames_df1 <- rownames(meta_pro)
asv_pro_filt <- asv_pro_filt[rownames_df1, ]
rownames(asv_pro_filt) == rownames(meta_pro)#now in order

eco_pro_filt <- eco_pro_filt[rownames_df1, ]
rownames(asv_pro_filt) == rownames(eco_pro_filt)#now in order
#asv, meta, and eco all in same order

mds.eco <- metaMDS(eco_pro_filt)
mds.asv <- metaMDS(asv_pro_filt)

procrustes_result2 <- protest(mds.eco,mds.asv, scale=TRUE)
procrustes_result2

pro_result2 <- procrustes(mds.asv,mds.eco, symmetric = TRUE)
res <- residuals(pro_result2)

#Procrustes Sum of Squares (m12 squared):        0.9119 
#Correlation in a symmetric Procrustes rotation: 0.2968 
#Significance:  0.001 
plot(procrustes_result2)
plot(procrustes_result2, kind = 2)


Yscore <- procrustes_result2$Y
Xscore <- procrustes_result2$X


pro_final_merge <- as.data.frame(cbind(Yscore, Xscore))

pro_final_merge$sample_id <- rownames(pro_final_merge)
meta_pro$sample_id <- rownames(meta_pro)
pro_final_merge_meta <- dplyr::left_join(pro_final_merge, meta_pro, by="sample_id")

pro_final_merge_meta$temperature <- as.factor(pro_final_merge_meta$temperature)
pro_final_merge_meta$ph <- as.factor(pro_final_merge_meta$ph)
pro_final_merge_meta$food <- as.factor(pro_final_merge_meta$food)
pro_final_merge_meta$treatment_combo <- as.factor(pro_final_merge_meta$treatment_combo)

pro_a <- ggplot(pro_final_merge_meta, aes(group = sample_id, color=temperature)) +
  geom_point(aes(x = V1, y = V2)) + 
  geom_point(aes(x = NMDS1, y = NMDS2))+
  geom_segment(aes(x = V1, y = V2, xend = NMDS1, yend = NMDS2, group = sample_id))+
  theme_classic()+scale_color_manual(values=c("#1f5776", "#c83126"))+theme(legend.position = "none")

pro_b <-ggplot(pro_final_merge_meta, aes(group = sample_id, color=ph)) +
  geom_point(aes(x = V1, y = V2)) + 
  geom_point(aes(x = NMDS1, y = NMDS2))+
  geom_segment(aes(x = V1, y = V2, xend = NMDS1, yend = NMDS2, group = sample_id))+
  theme_classic()+scale_color_manual(values=c("#CA64A3","#9A4EAE", "#301934"))+theme(legend.position = "none")

pro_c <-ggplot(pro_final_merge_meta, aes(group = sample_id, color=food)) +
  geom_point(aes(x = V1, y = V2, alpha=day)) + 
  geom_point(aes(x = NMDS1, y = NMDS2))+
  geom_segment(aes(x = V1, y = V2, xend = NMDS1, yend = NMDS2, group = sample_id))+
  theme_classic()+scale_color_manual(values=c("#90ee90", "#228822"))+theme(legend.position = "none")

ggarrange(pro_a, pro_b, pro_c, nrow=1)




#repeat but just for the day 57
dim(eco_all)#96
sample_ids <- unique(eco_all$sample_id)
length(sample_ids) #48

meta_pro <- meta
meta_pro <- subset(meta_pro, day %in% c("57"))

meta_pro_filt <- meta_pro
meta_pro_filt$id <- rownames(meta_pro_filt)
meta_pro_filt$wk <- NA
meta_pro_filt$wk[meta_pro_filt$day == 1] <- 1  
meta_pro_filt$wk[meta_pro_filt$day == 57] <- 8 


meta_pro_filt <- meta_pro_filt[,c(2, 19, 20)]

eco_pro_all <- dplyr::left_join(eco_all, meta_pro_filt, by=c("sample_id", "wk"))
#eco_pro_all <- eco_pro_all[1:48,]
eco_pro_all <- eco_pro_all[49:96,]
rownames(eco_pro_all) <- eco_pro_all$id
eco_pro_filt <- eco_pro_all[,c(1:31)]

#format ASV data
asv_pro <- asv16s.rt
#rownames(meta_pro) <- meta_pro[,1]
asv_pro_filt <- subset(asv_pro, row.names(asv_pro) %in% rownames(meta_pro))
dim(asv_pro_filt) #48

rownames(asv_pro_filt) == rownames(meta_pro) #not in order

#fix order
rownames_df1 <- rownames(meta_pro)
asv_pro_filt <- asv_pro_filt[rownames_df1, ]
rownames(asv_pro_filt) == rownames(meta_pro)#now in order

eco_pro_filt <- eco_pro_filt[rownames_df1, ]
rownames(asv_pro_filt) == rownames(eco_pro_filt)#now in order
#asv, meta, and eco all in same order

asv_pro_filt <- asv_pro_filt[,colSums(asv_pro_filt) >0]

#procrustes analysis
mds.eco2 <- metaMDS(eco_pro_filt)
mds.asv2 <- metaMDS(asv_pro_filt)

procrustes_result3 <- protest(mds.eco2,mds.asv2, scale=TRUE)
procrustes_result3
plot(procrustes_result3, kind=2)
#day57 only
#Procrustes Sum of Squares (m12 squared):        0.9433 
#Correlation in a symmetric Procrustes rotation: 0.238
#Significance:  0.113


#day1 only
#Procrustes Sum of Squares (m12 squared):        0.5914 
#Correlation in a symmetric Procrustes rotation: 0.6392
#Significance:  0.001
plot(procrustes_result3)

Yrotscore <- as.data.frame(procrustes_result3$Yrot)
Xscore <- as.data.frame(procrustes_result3$X)

pro_final_merge <- cbind(Yrotscore, Xscore)

pro_final_merge$sample_id <- rownames(pro_final_merge)
meta_pro$sample_id <- rownames(meta_pro)
pro_final_merge_meta <- dplyr::left_join(pro_final_merge, meta_pro, by="sample_id")

pro_final_merge_meta$temperature <- as.factor(pro_final_merge_meta$temperature)
pro_final_merge_meta$ph <- as.factor(pro_final_merge_meta$ph)
pro_final_merge_meta$food <- as.factor(pro_final_merge_meta$food)
pro_final_merge_meta$treatment_combo <- as.factor(pro_final_merge_meta$treatment_combo)

ggplot(pro_final_merge_meta, aes(group = sample_id, shape=food)) +
  geom_point(aes(x = V1, y = V2, size=2, color=temperature)) + 
  geom_point(aes(x = NMDS1, y = NMDS2, size=2, color=temperature))+
  geom_segment(aes(x = V1, y = V2, xend = NMDS1, yend = NMDS2,group = sample_id, color="gray"),
               arrow = arrow(type = "closed", length = unit(0.1, "inches")))+
  facet_wrap(~ph)+scale_color_manual(values=c("#1f5776", "#c83126", "gray"))+
  theme_classic()+theme(text = element_text(color = "black", size = 16 ))+
  xlab("")+ ylab("")
  
pro_a <- ggplot(pro_final_merge_meta, aes(group = sample_id, color=temperature)) +
  geom_point(aes(x = V1, y = V2)) + 
  geom_point(aes(x = NMDS1, y = NMDS2))+
  geom_segment(aes(x = V1, y = V2, xend = NMDS1, yend = NMDS2, group = sample_id))+
  theme_classic()+scale_color_manual(values=c("#1f5776", "#c83126"))+theme(legend.position = "none")

pro_b <-ggplot(pro_final_merge_meta, aes(group = sample_id, color=ph)) +
  geom_point(aes(x = V1, y = V2)) + 
  geom_point(aes(x = NMDS1, y = NMDS2))+
  geom_segment(aes(x = V1, y = V2, xend = NMDS1, yend = NMDS2, group = sample_id))+
  theme_classic()+scale_color_manual(values=c("#CA64A3","#9A4EAE", "#301934"))+theme(legend.position = "none")

pro_c <-ggplot(pro_final_merge_meta, aes(group = sample_id, color=food)) +
  geom_point(aes(x = V1, y = V2)) + 
  geom_point(aes(x = NMDS1, y = NMDS2))+
  geom_segment(aes(x = V1, y = V2, xend = NMDS1, yend = NMDS2, group = sample_id))+
  theme_classic()+scale_color_manual(values=c("#90ee90", "#228822"))+theme(legend.position = "none")

ggarrange(pro_a, pro_b, pro_c, nrow=1)


#### ECOPLATE MANTEL TESTS ####
mantel_meta <- meta.filt
mantel_meta$id <- rownames(mantel_meta)
mantel_meta$wk <- NA
mantel_meta$wk[mantel_meta$day == 1] <- 1  
mantel_meta$wk[mantel_meta$day == 57] <- 8 

mantel_meta <- subset(mantel_meta, wk %in% c("1", "8"))

meta_ht<- subset(mantel_meta, temperature %in% c("37"))
meta_lt<- subset(mantel_meta, temperature %in% c("22"))
meta_ph3<- subset(mantel_meta, ph %in% c("3"))
meta_ph4<- subset(mantel_meta, ph %in% c("4"))
meta_ph5.6<- subset(mantel_meta, ph %in% c("5.6"))
meta_food3<- subset(mantel_meta, food %in% c("3"))
meta_food6<- subset(mantel_meta, food %in% c("6"))

eco_mantel <- dplyr::left_join(mantel_meta, eco_all, by=c("sample_id", "wk"))
rownames(eco_mantel) <- eco_mantel$id
eco_mantel <- eco_mantel[,20:50]

#format ASV data
asv_pro <- asv16s.rt

eco_ht <- subset(eco_mantel, row.names(eco_mantel) %in% rownames(meta_ht))
eco_lt <- subset(eco_mantel, row.names(eco_mantel) %in% rownames(meta_lt))
eco_ph3 <- subset(eco_mantel, row.names(eco_mantel) %in% rownames(meta_ph3))
eco_ph4 <- subset(eco_mantel, row.names(eco_mantel) %in% rownames(meta_ph4))
eco_ph5.6 <- subset(eco_mantel, row.names(eco_mantel) %in% rownames(meta_ph5.6))
eco_food3 <- subset(eco_mantel, row.names(eco_mantel) %in% rownames(meta_food3))
eco_food6 <- subset(eco_mantel, row.names(eco_mantel) %in% rownames(meta_food6))

asv_ht <- subset(asv16s.rt, row.names(asv16s.rt) %in% rownames(meta_ht))
asv_lt <- subset(asv16s.rt, row.names(asv16s.rt) %in% rownames(meta_lt))
asv_ph3 <- subset(asv16s.rt, row.names(asv16s.rt) %in% rownames(meta_ph3))
asv_ph4 <- subset(asv16s.rt, row.names(asv16s.rt) %in% rownames(meta_ph4))
asv_ph5.6 <- subset(asv16s.rt, row.names(asv16s.rt) %in% rownames(meta_ph5.6))
asv_food3 <- subset(asv16s.rt, row.names(asv16s.rt) %in% rownames(meta_food3))
asv_food6 <- subset(asv16s.rt, row.names(asv16s.rt) %in% rownames(meta_food6))

rownames_1 <- rownames(eco_ht)
asv_ht <- asv_ht[rownames_1, ]
rownames(asv_ht) == rownames(eco_ht)
asv_ht <- asv_ht[,colSums(asv_ht) >0]

rownames_2 <- rownames(eco_lt)
asv_lt <- asv_lt[rownames_2, ]
rownames(asv_lt) == rownames(eco_lt)
asv_lt <- asv_lt[,colSums(asv_lt) >0]

rownames_3 <- rownames(eco_ph3)
asv_ph3 <- asv_ph3[rownames_3, ]
rownames(asv_ph3) == rownames(eco_ph3)
asv_ph3 <- asv_ph3[,colSums(asv_ph3) >0]

rownames_4 <- rownames(eco_ph4)
asv_ph4 <- asv_ph4[rownames_4, ]
rownames(asv_ph4) == rownames(eco_ph4)
asv_ph4 <- asv_ph4[,colSums(asv_ph4) >0]

rownames_5 <- rownames(eco_ph5.6)
asv_ph5.6 <- asv_ph5.6[rownames_5, ]
rownames(asv_ph5.6) == rownames(eco_ph5.6)
asv_ph5.6 <- asv_ph5.6[,colSums(asv_ph5.6) >0]

rownames_6 <- rownames(eco_food3)
asv_food3 <- asv_food3[rownames_6, ]
rownames(asv_food3) == rownames(eco_food3)
asv_food3 <- asv_food3[,colSums(asv_food3) >0]

rownames_7 <- rownames(eco_food6)
asv_food6 <- asv_food6[rownames_7, ]
rownames(asv_food6) == rownames(eco_food6)
asv_food6 <- asv_food6[,colSums(asv_food6) >0]

mantel_asv_test <- as.matrix(wu.dist.16s.tubes)

asv_dist_ht <- as.dist(mantel_asv_test[rownames_1, rownames_1])
asv_dist_lt <- as.dist(mantel_asv_test[rownames_2, rownames_2])
asv_dist_ph3 <- as.dist(mantel_asv_test[rownames_3, rownames_3])
asv_dist_ph4 <- as.dist(mantel_asv_test[rownames_4, rownames_4])
asv_dist_ph5.6 <- as.dist(mantel_asv_test[rownames_5, rownames_5])
asv_dist_food3 <- as.dist(mantel_asv_test[rownames_6, rownames_6])
asv_dist_food6 <- as.dist(mantel_asv_test[rownames_7, rownames_7])

eco_dist_ht <- vegdist(eco_ht, method="bray")
eco_dist_lt <- vegdist(eco_lt, method="bray")
eco_dist_ph3 <- vegdist(eco_ph3, method="bray")
eco_dist_ph4 <- vegdist(eco_ph4, method="bray")
eco_dist_ph5.6 <- vegdist(eco_ph5.6, method="bray")
eco_dist_food3 <- vegdist(eco_food3, method="bray")
eco_dist_food6 <- vegdist(eco_food6, method="bray")

mantel_ht <- mantel(asv_dist_ht,eco_dist_ht, method = "spearman", permutations=999)
mantel_ht
#Mantel statistic r: 0.4507  -1 strong negative, +1 strong positive, 0 none
#Significance: 0.001
plot(asv_dist_ht,eco_dist_ht)

asv_dist_vectorht <- as.vector(asv_dist_ht)
eco_dist_vectorht <- as.vector(eco_dist_ht)
data_lt <- data.frame(asv_dist_vectorht = asv_dist_vectorht, eco_dist_vectorht = eco_dist_vectorht)

## Figure 6A ##
ggplot(data_lt, aes(x = asv_dist_vectorht, y = eco_dist_vectorht)) +
  geom_point(col="#c83126") +xlim(c(0,.82))+ylim(c(0,.82))+theme_classic()+
  geom_smooth(method = "lm", col = "black") +
  labs( x = "High Temperature Weighted UniFrac Dissimilarity\n(16S Community Composition)",
       y = "High Temperature Bray-Curtis Dissimilarity\n(EcoPlate Carbon Substrate Use)")

mantel_lt <- mantel(asv_dist_lt,eco_dist_lt, method = "spearman", permutations=999)
mantel_lt
#Mantel statistic r: 0.3023  -1 strong negative, +1 strong positive, 0 none
#Significance: 0.001
plot(asv_dist_lt,eco_dist_lt)

asv_dist_vectorlt <- as.vector(asv_dist_lt)
eco_dist_vectorlt <- as.vector(eco_dist_lt)
data_lt <- data.frame(asv_dist_vectorlt = asv_dist_vectorlt, eco_dist_vectorlt = eco_dist_vectorlt)

## Figure 6B ##
ggplot(data_lt, aes(x = asv_dist_vectorlt, y = eco_dist_vectorlt)) +
  geom_point(col="#1f5776") +xlim(c(0,.82))+ylim(c(0,.82))+theme_classic()+
  geom_smooth(method = "lm", col = "black") +
  labs( x = "Normal Temperature Weighted UniFrac Dissimilarity\n(16S Community Composition)",
        y = "Normal Temperature Bray-Curtis Dissimilarity\n(EcoPlate Carbon Substrate Use)")


mantel_ph3 <- mantel(asv_dist_ph3,eco_dist_ph3, method = "spearman", permutations=999)
mantel_ph3
#Mantel statistic r: 0.3603  -1 strong negative, +1 strong positive, 0 none
#Significance: 0.001
plot(asv_dist_ph3,eco_dist_ph3)

mantel_ph4 <- mantel(asv_dist_ph4,eco_dist_ph4, method = "spearman", permutations=999)
mantel_ph4
#Mantel statistic r: 0.4423  -1 strong negative, +1 strong positive, 0 none
#Significance: 0.001
plot(asv_dist_ph4,eco_dist_ph4)

mantel_ph5.6 <- mantel(asv_dist_ph5.6,eco_dist_ph5.6, method = "spearman", permutations=999)
mantel_ph5.6
#Mantel statistic r: 0.4311  -1 strong negative, +1 strong positive, 0 none
#Significance: 0.001
plot(asv_dist_ph5.6,eco_dist_ph5.6)

mantel_food3 <- mantel(asv_dist_food3,eco_dist_food3, method = "spearman", permutations=999)
mantel_food3
#Mantel statistic r: 0.2626  -1 strong negative, +1 strong positive, 0 none
#Significance: 0.003
plot(asv_dist_food3,eco_dist_food3)

mantel_food6 <- mantel(asv_dist_food6,eco_dist_food6, method = "spearman", permutations=999)
mantel_food6
#Mantel statistic r: 0.3882  -1 strong negative, +1 strong positive, 0 none
#Significance: 0.001
plot(asv_dist_food6,eco_dist_food6)



#procrustes analysis
mds.ecoht <- metaMDS(eco_ht)
mds.asvht <- metaMDS(asv_ht)

procrustes_resulttemmp37 <- protest(mds.asvht,mds.ecoht, scale=TRUE)
procrustes_resulttemmp37
plot(procrustes_resulttemmp37)
#high temp
#Procrustes Sum of Squares (m12 squared):        0.6454
#Correlation in a symmetric Procrustes rotation: 0.5955
#Significance:  0.001

Yrotscoreht <- as.data.frame(procrustes_resulttemmp37$Yrot)
Xscoreht <- as.data.frame(procrustes_resulttemmp37$X)
pro_ht <- cbind(Xscoreht, Yrotscoreht)
pro_ht$sample_id <- rownames(pro_ht)
meta_ht$sample_id <- rownames(meta_ht)
pro_meta_ht <- dplyr::left_join(pro_ht, meta_ht, by="sample_id")
pro_meta_ht$temperature <- as.factor(pro_meta_ht$temperature)

## Figure 6C ##
ggplot(pro_meta_ht) +
  geom_segment(aes(x = NMDS1, y = NMDS2, xend = V1, yend = V2, color="gray"), color="gray",
               arrow = arrow(type = "open", length = unit(0.1, "inches")))+
  geom_point(aes(x = V1, y = V2), size=3, col="#c83126",shape=1) + 
  geom_point(aes(x = NMDS1, y = NMDS2), size=3, col="#c83126",shape=19)+theme_classic()+
  xlim(c(-.26,.42))+ylim(c(-.1,.25))


max(pro_meta_ht$NMDS1)#0.1643132
min(pro_meta_ht$NMDS2)#-0.1167475
max(pro_meta_ht$V1)#0.245463
min(pro_meta_ht$V2)#-0.05705124





mds.ecolt <- metaMDS(eco_lt)
mds.asvlt <- metaMDS(asv_lt)

procrustes_resulttemmp22 <- protest(mds.asvlt,mds.ecolt, scale=TRUE)
procrustes_resulttemmp22
plot(procrustes_resulttemmp22)
#low temp
#Procrustes Sum of Squares (m12 squared):        0.7987
#Correlation in a symmetric Procrustes rotation: 0.4486
#Significance:  0.001

Yrotscorelt <- as.data.frame(procrustes_resulttemmp22$Yrot)
Xscorelt <- as.data.frame(procrustes_resulttemmp22$X)
pro_lt <- cbind(Xscorelt, Yrotscorelt)
pro_lt$sample_id <- rownames(pro_lt)
meta_lt$sample_id <- rownames(meta_lt)
pro_meta_lt <- left_join(pro_lt, meta_lt, by="sample_id")
pro_meta_lt$temperature <- as.factor(pro_meta_lt$temperature)

## Figure 6D ##
ggplot(pro_meta_lt) +
  geom_segment(aes(x = NMDS1, y = NMDS2, xend = V1, yend = V2, color="gray"), color="gray",
               arrow = arrow(type = "open", length = unit(0.1, "inches")))+
  geom_point(aes(x = V1, y = V2), size=3, col="#1f5776",shape=1) + 
  geom_point(aes(x = NMDS1, y = NMDS2), size=3, col="#1f5776",shape=19)+theme_classic()+
  xlim(c(-.26,.42))+ylim(c(-.1,.25))


max(pro_meta_lt$NMDS1)#0.4004384
min(pro_meta_lt$NMDS2)#-0.2561536
max(pro_meta_lt$V1)#0.08987522
min(pro_meta_lt$V2)#-0.05467201


procrustes_resultph4 <- protest(mds.eco2,mds.asv2, scale=TRUE)
procrustes_resultph4
plot(procrustes_resultph4)
#ph4
#Procrustes Sum of Squares (m12 squared):        0.6354
#Correlation in a symmetric Procrustes rotation: 0.6038
#Significance:  0.001


procrustes_resultfood3 <- protest(mds.eco2,mds.asv2, scale=TRUE)
procrustes_resultfood3
plot(procrustes_resultfood3)
#low food
#Procrustes Sum of Squares (m12 squared):        0.9246
#Correlation in a symmetric Procrustes rotation: 0.2746
#Significance:  0.047

procrustes_resultfood6 <- protest(mds.eco2,mds.asv2, scale=TRUE)
procrustes_resultfood6
plot(procrustes_resultfood6)
#highfood
#Procrustes Sum of Squares (m12 squared):        0.8056
#Correlation in a symmetric Procrustes rotation: 0.4409
#Significance:  0.002


procrustes_resulttemmp37 <- protest(mds.eco2,mds.asv2, scale=TRUE)
procrustes_resulttemmp37
plot(procrustes_resulttemmp37)
#high temp
#Procrustes Sum of Squares (m12 squared):        0.6481
#Correlation in a symmetric Procrustes rotation: 0.5932
#Significance:  0.001

procrustes_resulttemmp22 <- protest(mds.eco2,mds.asv2, scale=TRUE)
procrustes_resulttemmp22
plot(procrustes_resulttemmp22)
#low temp
#Procrustes Sum of Squares (m12 squared):        0.7929
#Correlation in a symmetric Procrustes rotation: 0.4551
#Significance:  0.001


procrustes_resultph3 <- protest(mds.eco2,mds.asv2, scale=TRUE)
procrustes_resultph3
plot(procrustes_resultph3)
#low ph
#Procrustes Sum of Squares (m12 squared):        0.8707
#Correlation in a symmetric Procrustes rotation: 0.3596
#Significance:  0.036

procrustes_resultph5.6 <- protest(mds.eco2,mds.asv2, scale=TRUE)
procrustes_resultph5.6
plot(procrustes_resultph5.6)
#highph
#Procrustes Sum of Squares (m12 squared):        0.7847
#Correlation in a symmetric Procrustes rotation: 0.464
#Significance:  0.003






meta_pro <- subset(meta, food %in% c("6"))
meta_pro$food <- as.factor(meta_pro$food)
Yrotscore <- as.data.frame(procrustes_resultfood6$Yrot)
Xscore <- as.data.frame(procrustes_resultfood6$X)
pro_final_merge <- cbind(Yrotscore, Xscore)
pro_final_merge$sample_id <- rownames(pro_final_merge)
meta_pro$sample_id <- rownames(meta_pro)
pro_final_merge_meta <- left_join(pro_final_merge, meta_pro, by="sample_id")
pro_final_merge_meta$temperature <- as.factor(pro_final_merge_meta$temperature)
pro_final_merge_meta$ph <- as.factor(pro_final_merge_meta$ph)
ggplot(pro_final_merge_meta) +
  geom_segment(aes(x = V1, y = V2, xend = NMDS1, yend = NMDS2, color="gray"), color="gray",
               arrow = arrow(type = "open", length = unit(0.1, "inches")))+
  geom_point(aes(x = V1, y = V2, color=food,shape=day), size=3) + 
  geom_point(aes(x = NMDS1, y = NMDS2, color=food,shape=day), size=3)+
  theme_classic()+scale_color_manual(values=c("#228822"))



#### CHITINASE ####
#filter metadata
meta<- read.csv("data/Exp_2_metadata_tubes.csv", header=TRUE)
meta_all_tubes <- subset(meta, container =="tube")
meta_all_tubes$food <- as.factor(meta_all_tubes$food)
meta_all_tubes$ph <- as.factor(meta_all_tubes$ph)
meta_all_tubes$temperature <- as.factor(meta_all_tubes$temperature)
meta_all_tubes$treatment_combo <- as.factor(meta_all_tubes$treatment_combo)
meta_all_tubes$day <- as.integer(meta_all_tubes$day)

detach("package:plyr", unload=TRUE)

df.summarychit <- meta_all_tubes %>%
  group_by(treatment_combo,day, food, temperature, ph) %>%
  summarise(
    sd = sd(chitinase, na.rm = TRUE),
    chitinase = mean(chitinase, na.rm = TRUE))
df.summarychit <- na.omit(df.summarychit)

## Figure 2A ##
ggplot(df.summarychit, aes(day, chitinase, color = temperature, shape = food)) +
  geom_jitter(size = 3, alpha = 0.9, position = position_dodge(0.3)) +
  facet_wrap(~ph, nrow = 1) +
  geom_line(aes(group = treatment_combo), position = position_dodge(0.3)) +
  geom_errorbar(aes(ymin = chitinase - sd, ymax = chitinase + sd), position = position_dodge(0.3), width = 0.2) +
  theme_bw() +scale_shape_manual(values=c(16,2)) +
  scale_color_manual(values = c("#1f5776", "#c83126")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none") 

#### PROTEASE ####
df.summaryprot <- meta_all_tubes %>%
  group_by(treatment_combo,day, food, temperature, ph) %>%
  summarise(
    sd = sd(protease, na.rm = TRUE),
    protease = mean(protease, na.rm = TRUE))
df.summaryprot <- na.omit(df.summaryprot)

## Figure 2C ##
ggplot(df.summaryprot, aes(day, protease, color = temperature, shape = food)) +
  geom_jitter(size = 3, alpha = 0.9, position = position_dodge(0.3)) +
  facet_wrap(~ph, nrow = 1) +
  geom_line(aes(group = treatment_combo), position = position_dodge(0.3)) +
  geom_errorbar(aes(ymin = protease - sd, ymax = protease + sd), position = position_dodge(0.3), width = 0.2) +
  theme_bw() +scale_shape_manual(values=c(16,2)) +
  scale_color_manual(values = c("#1f5776", "#c83126")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none") 

#### RESPIRATION ####
meta_all_tubes.resp <- meta_all_tubes[complete.cases(meta_all_tubes$respiration_co2_ppm_hr),]

df.summaryresp <- meta_all_tubes.resp %>%
  group_by(treatment_combo,day, food, temperature, ph) %>%
  summarise(
    sd = sd(respiration_co2_ppm_hr, na.rm = TRUE),
    respiration_co2_ppm_hr = mean(respiration_co2_ppm_hr, na.rm = TRUE))
df.summaryresp <- na.omit(df.summaryresp)

## Figure 2E ##
ggplot(df.summaryresp, aes(day, respiration_co2_ppm_hr, color = temperature, shape = food)) +
  geom_jitter(size = 3, alpha = 0.9, position = position_dodge(0.3)) +
  facet_wrap(~ph, nrow = 1) +
  geom_line(aes(group = treatment_combo), position = position_dodge(0.3)) +
  geom_errorbar(aes(ymin = respiration_co2_ppm_hr - sd, ymax = respiration_co2_ppm_hr + sd), position = position_dodge(0.3), width = 0.2) +
  theme_bw() +scale_shape_manual(values=c(16,2)) +
  scale_color_manual(values = c("#1f5776", "#c83126")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none") +scale_x_continuous(breaks = c(0, 20, 40))

#### GLMS ####
#response = chitinase, protease, respiration
#predictor = ph, food, temp, day
#repeated measure = sample_id = removed becasue not enough data/replication to control for repeated measures with all the interactive effects

#relevel the predictor variables to baselines
meta_all_tubes$food <- as.factor(meta_all_tubes$food)
meta_all_tubes$temperature <- as.factor(meta_all_tubes$temperature)
meta_all_tubes$ph <- as.factor(meta_all_tubes$ph)
meta_all_tubes$food <- relevel(meta_all_tubes$food, ref = "3")
meta_all_tubes$temperature <- relevel(meta_all_tubes$temperature, ref = "22")
meta_all_tubes$ph <- relevel(meta_all_tubes$ph, ref = "5.6")

informpriors_chit <- c(prior(normal(-5.241044, 0.1770824), 
                        class = "Intercept"))

informpriors_prot <- c(prior(normal(6.928552, 0.2161154), 
                             class = "Intercept"))

informpriors_resp <- c(prior(normal(7.850019, 1.511284), 
                             class = "Intercept"))

mchit <- brm(chitinase ~ ph*food*temperature + day,data=meta_all_tubes,
             family = Gamma(link = "log"), iter = 10000, chains = 4, cores = 4,control = list(adapt_delta = 0.99, stepsize = 0.001, max_treedepth=15),
             prior=informpriors_chit)

meta_all_tubes$protease <- ifelse(meta_all_tubes$protease <= 0, 0.001, meta_all_tubes$protease)
mprot <- brm(protease ~ ph*food*temperature + day,data=meta_all_tubes,
             family = Gamma(link = "log"), iter = 10000, chains = 4, cores = 4,control = list(adapt_delta = 0.999, stepsize = 0.001, max_treedepth=20),
             prior=informpriors_prot)


mresp <- brm(respiration_co2_ppm_hr ~ ph*food*temperature+ day,data=meta_all_tubes,
             family = Gamma(link = "log"), iter = 10000, chains = 4, cores = 4,control = list(adapt_delta = 0.999, stepsize = 0.001, max_treedepth=20),
             prior=informpriors_resp)

summary(mchit)
summary(mprot)
summary(mresp)

saveRDS(mchit, file = "RDS/brms_exp2_chit.RDS")
saveRDS(mprot, file = "RDS/brms_exp2_prot.RDS")
saveRDS(mresp, file = "RDS/brms_exp2_resp.RDS")

mchit <- readRDS("RDS/brms_exp2_chit.RDS")
mprot <- readRDS("RDS/brms_exp2_prot.RDS")
mresp <- readRDS("RDS/brms_exp2_resp.RDS")

pp_check(mchit)
pp_check(mprot)
pp_check(mresp)

#### GLM INTERACTION PLOTS ####
mchit
tab_model(mchit)
tab_model(mchit, show.est = TRUE, show.ci = TRUE, ci.lvl = 0.95, show.se = TRUE, 
          file = "output_files/brms_chit_table.doc")



chit.pred <- predict_response(mchit, terms = c("food", "temperature", "ph"))
chit.pred$facet <- factor(chit.pred$facet, levels = c("3", "4", "5.6"))

## Figure 2B ##
ggplot(chit.pred, aes(x = x, y = predicted, shape = x, color = group)) + 
  facet_wrap(~facet) +
  scale_color_manual(values = c("#1f5776", "#c83126")) +
  theme_bw() +scale_shape_manual(values=c(16,2)) +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0, size=1,position = position_dodge(width = .7))+
  theme(legend.position="none")+ ylab("Predicted Chitinase Activity")+
  xlab("Food (g/L)")+ geom_point(size = 3, position = position_dodge(width = .7)) +theme(
    text = element_text(color = "black", size = 10 ))


tab_model(mprot, show.est = TRUE, show.ci = TRUE, show.se = TRUE, 
          file = "output_files/brms_prot_table.doc")
prot.pred <- ggpredict(mprot, terms = c("food", "temperature", "ph"))
prot.pred$facet <- factor(prot.pred$facet, levels = c("3", "4", "5.6"))

## Figure 2D ##
ggplot(prot.pred, aes(x = x, y = predicted, shape = x, color = group)) +
  facet_wrap(~facet) +
  scale_color_manual(values = c("#1f5776", "#c83126")) +
  theme_bw() +scale_shape_manual(values=c(16,2))+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0, size=1,position = position_dodge(width = .7))+
  theme(legend.position="none")+ ylab("Predicted Protease Activity")+
  xlab("Food (g/L)")+ geom_point(size = 3, position = position_dodge(width = .7)) +theme(
    text = element_text(color = "black", size = 10 )
  )

tab_model(mresp, show.est = TRUE, show.ci = TRUE, show.se = TRUE, 
          file = "output_files/brms_resp_table.doc")
resp.pred <- ggpredict(mresp, terms = c("food", "temperature", "ph"))
resp.pred$facet <- factor(resp.pred$facet, levels = c("3", "4", "5.6"))

## Figure 2F ##
ggplot(resp.pred, aes(x = x, y = predicted, shape = x, color = group)) +
  facet_wrap(~facet) +
  scale_color_manual(values = c("#1f5776", "#c83126")) +
  theme_bw() +scale_shape_manual(values=c(16,2))+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0, size=1,position = position_dodge(width = .7))+
  theme(legend.position="none")+ ylab("Predicted Respir")+
  xlab("Food (g/L)")+ geom_point(size = 3, position = position_dodge(width = .7)) +theme(
    text = element_text(color = "black", size = 10 ))+ylim(c(0,35000))

#Directional Effect of food on respiration
posterior_mresp <- as.data.frame(mresp)
highfood <- posterior_mresp %>% filter(b_food6 >0)
nrow(highfood)/nrow(posterior_mresp)*100

# POSTERIOR PROBABILITY GRAPHS
posteriorchit <- mcmc_intervals_data(mchit, 
                                  prob_outer=0.95,
                                  prob=0.5)
posteriorchit$nonzero <- NA
posteriorchit$nonzero[posteriorchit$ll>0 & posteriorchit$hh>0] <- "nonzero"
posteriorchit$nonzero[posteriorchit$ll<0 & posteriorchit$hh<0] <- "nonzero"
posteriorchit$nonzero[is.na(posteriorchit$nonzero)] <- "zero"
posteriorchit<- posteriorchit[1:13,]

posteriorchit <- posteriorchit %>%
  mutate(parameter = str_remove(parameter, "^b_") %>%  # Remove "b_" prefix
           str_replace_all(":", ": ") %>%
           str_replace_all("ph3", "pH 3") %>%
           str_replace_all("ph4", "pH 4") %>%
           str_replace_all("food6", "high food") %>%
           str_replace_all("temperature37", "high temp"))

posteriorchit <- posteriorchit %>%
  mutate(parameter = factor(parameter, 
                            levels = c("Intercept", 
                                       "day",
                                       "high temp", 
                                       "high food", 
                                       "pH 3", 
                                       "pH 4", 
                                       "pH 3: high food", 
                                       "pH 4: high food", 
                                       "pH 3: high temp", 
                                       "pH 4: high temp", 
                                       "high food: high temp", 
                                       "pH 3: high food: high temp", 
                                       "pH 4: high food: high temp")))

## Figure S2 ##
ggplot(posteriorchit, aes(x = parameter, shape=nonzero)) +
  geom_hline(yintercept = 0, linetype = 2, size=.35, color = "black") +
  geom_pointrange(aes(ymin = ll, ymax = hh, y = m), position= position_dodge(width=0.75), size = 1, linewidth = 1) +
  scale_shape_manual(values=c(17, 19), labels=c("95% CI does\nnot contain zero", "95% CI\ncontains zero"))+
  coord_flip() +xlab(NULL) + ylab("Estimated effect on Chitinase Activity")+ theme(text = element_text(color = "black", size = 16))+
  theme(panel.background = element_rect(fill = "white", colour = "black"))+theme(legend.position="none")

#PROTEASE MODEL
posteriorprot <- mcmc_intervals_data(mprot, 
                                     prob_outer=0.95,
                                     prob=0.5)

posteriorprot$nonzero <- NA
posteriorprot$nonzero[posteriorprot$ll>0 & posteriorprot$hh>0] <- "nonzero"
posteriorprot$nonzero[posteriorprot$ll<0 & posteriorprot$hh<0] <- "nonzero"
posteriorprot$nonzero[is.na(posteriorprot$nonzero)] <- "zero"
posteriorprot<- posteriorprot[1:13,]

posteriorprot <- posteriorprot %>%
  mutate(parameter = str_remove(parameter, "^b_") %>%  # Remove "b_" prefix
           str_replace_all(":", ": ") %>%
           str_replace_all("ph3", "pH 3") %>%
           str_replace_all("ph4", "pH 4") %>%
           str_replace_all("food6", "high food") %>%
           str_replace_all("temperature37", "high temp"))

posteriorprot <- posteriorprot %>%
  mutate(parameter = factor(parameter, 
                            levels = c("Intercept", 
                                       "day",
                                       "high temp", 
                                       "high food", 
                                       "pH 3", 
                                       "pH 4", 
                                       "pH 3: high food", 
                                       "pH 4: high food", 
                                       "pH 3: high temp", 
                                       "pH 4: high temp", 
                                       "high food: high temp", 
                                       "pH 3: high food: high temp", 
                                       "pH 4: high food: high temp")))

## Figure S2 ##
ggplot(posteriorprot, aes(x = parameter, shape=nonzero)) +
  geom_hline(yintercept = 0, linetype = 2, size=.35, color = "black") +
  geom_pointrange(aes(ymin = ll, ymax = hh, y = m), position= position_dodge(width=0.75), size = 1, linewidth = 1) +
  scale_shape_manual(values=c(17, 19), labels=c("95% CI does\nnot contain zero", "95% CI\ncontains zero"))+
  coord_flip() +xlab(NULL) + ylab("Estimated effect on Protease Activity")+ theme(text = element_text(color = "black", size = 16))+
  theme(panel.background = element_rect(fill = "white", colour = "black"))+theme(legend.position="none")



#RESPIRATION MODEL
posteriorresp <- mcmc_intervals_data(mresp, 
                                     prob_outer=0.95,
                                     prob=0.5)

posteriorresp$nonzero <- NA
posteriorresp$nonzero[posteriorresp$ll>0 & posteriorresp$hh>0] <- "nonzero"
posteriorresp$nonzero[posteriorresp$ll<0 & posteriorresp$hh<0] <- "nonzero"
posteriorresp$nonzero[is.na(posteriorresp$nonzero)] <- "zero"
posteriorresp<- posteriorresp[1:13,]

posteriorresp <- posteriorresp %>%
  mutate(parameter = str_remove(parameter, "^b_") %>%  # Remove "b_" prefix
           str_replace_all(":", ": ") %>%
           str_replace_all("ph3", "pH 3") %>%
           str_replace_all("ph4", "pH 4") %>%
           str_replace_all("food6", "high food") %>%
           str_replace_all("temperature37", "high temp"))

posteriorresp <- posteriorresp %>%
  mutate(parameter = factor(parameter, 
                            levels = c("Intercept", 
                                       "day",
                                       "high temp", 
                                       "high food", 
                                       "pH 3", 
                                       "pH 4", 
                                       "pH 3: high food", 
                                       "pH 4: high food", 
                                       "pH 3: high temp", 
                                       "pH 4: high temp", 
                                       "high food: high temp", 
                                       "pH 3: high food: high temp", 
                                       "pH 4: high food: high temp")))

## Figure S2 ##
ggplot(posteriorresp, aes(x = parameter, shape=nonzero)) +
  geom_hline(yintercept = 0, linetype = 2, size=.35, color = "black") +
  geom_pointrange(aes(ymin = ll, ymax = hh, y = m), position= position_dodge(width=0.75), size = 1, linewidth = 1) +
  scale_shape_manual(values=c(17, 19), labels=c("95% CI does\nnot contain zero", "95% CI\ncontains zero"))+
  coord_flip() +xlab(NULL) + ylab("Estimated effect on Respiration")+ theme(text = element_text(color = "black", size = 16))+
  theme(panel.background = element_rect(fill = "white", colour = "black"))+theme(legend.position="none")
