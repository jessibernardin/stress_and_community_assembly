#### Community Assembly and Functioning under Environmental Stress ####
#### Authors: Jessica R. Bernardin, Leonora S. Bittleston ####
#### last update : December 19th, 2024 ####
#### Physiological Functions and Diversity Analysis

#### Load Required Packages ####
packages_to_load <- c(
  "ggplot2", "vegan", "lme4", "tidyverse", "effects",
  "dplyr", "reshape2", "ape", "DiagrammeR",
  "tidybayes", "coefplot", "standardize", "bayesplot", "MCMCvis", "car",
  "patchwork", "ggpubr", "corrr", "ggcorrplot", "factoextra", "MASS",
  "pairwiseAdonis", "plotrix", "gridExtra", "multcompView", "ggeffects", "this.path", "brms", "ggmulti",
  "phyloseq", "qiime2R", "picante", "decontam", "performance", "janitor", "ANCOMBC", "pheatmap", "chron",
  "lubridate", "igraph", "Hmisc", "qgraph", "mia", "cowplot", "microbiome"
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


#### Community Structure
#read in 16S
asv16s <- read_tsv("exported-files/asv-table-dada2.txt")
asv16s <- asv16s %>% column_to_rownames(var="#OTU ID")

meta <- read_csv("Exp_2_metadata_tubes.csv")
meta <- meta %>% remove_rownames %>% column_to_rownames(var="ID")

#make a new column for the scaled data
meta$scaled_chit <- as.numeric(scale(meta$chitinase))
meta$scaled_prot <- as.numeric(scale(meta$protease))
meta$sample_id <- as.factor(meta$sample_id)
str(meta$sample_id)

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
#what percentage of the asv are only present in one sample before prev filt
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

Exp2.physeq2.nomix <- subset_samples(Exp2.physeq2, day %in% c("1","8","15","22","29","36","43","50","57")) 
row_sums <- rowSums(otu_table(Exp2.physeq2.nomix))
nonzero_rows <- row_sums != 0
Exp2.physeq2.nomix <- prune_taxa(nonzero_rows, Exp2.physeq2.nomix)#1366 taxa and 144 samples



#### Rarify ####
#take out asv less than 10
#asv16s.exp2_filt[asv16s.exp2_filt < 10] <- 0
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
ggplot(data=md16s, aes(x=day, y=ef.16s)) +
  geom_jitter() + ylab("Effective species ASVs") + theme_classic() +
  facet_wrap(~container, nrow=2)

md16s_mix <- filter(md16s, container == "MIX")
summary(md16s_mix)
ggplot(data=md16s_mix, aes(x=ID, y=richness)) +
  geom_jitter() + ylab("Richness (ASVs)") + theme_classic()+
  ylim(c(0,50))

summary(c(47, 45, 46, 48))
SE = sd(c(47, 45, 46, 48)) / sqrt(length(c(47, 45, 46, 48)))
md16s$food <- as.factor(md16s$food)
md16s$ph <- as.factor(md16s$ph)
md16s$temperature <- as.factor(md16s$temperature)

ggplot(data=md16s, aes(x=day, y=ef.16s, color=temperature, shape=food)) +
  geom_jitter(size=3, alpha=.75) + ylab("Effective species ASVs") + theme_classic() +
  facet_wrap(~ph, nrow=1)+scale_color_manual(values=c("#1f5776", "#c83126"))

ggplot(data=md16s, aes(x=day, y=richness, group=treatment_combo, color=temperature, shape=food)) +
  geom_jitter(size=3, alpha=.75) + geom_smooth()+ylab("Richness") + theme_classic() +
  facet_wrap(~ph, nrow=1)+scale_color_manual(values=c("#1f5776", "#c83126"))

detach("package:plyr", unload=TRUE)

df.summaryrich <- md16s %>%
  group_by(treatment_combo,day, food, temperature, ph) %>%
  summarise(
    sd = sd(richness, na.rm = TRUE),
    richness = mean(richness, na.rm = TRUE))
df.summaryrich <- na.omit(df.summaryrich)

df.summaryrich$ph <- factor(df.summaryrich$ph, levels=c("3", "4", "5.6"))
ggplot(df.summaryrich, aes(day, richness, color = temperature, shape = food)) +
  geom_jitter(size = 3, alpha = 0.9, position = position_dodge(0.3)) +
  facet_wrap(~ph, nrow = 1) +
  geom_line(aes(group = treatment_combo), position = position_dodge(0.3)) +
  geom_errorbar(aes(ymin = richness - sd, ymax = richness + sd), position = position_dodge(0.3), width = 0.2) +
  theme_bw() +
  scale_color_manual(values = c("#1f5776", "#c83126")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none") +
  ylim(c(0, 100))


#model for alpha diversity
md16s$food <- relevel(md16s$food, ref = "3")
md16s$temperature <- relevel(md16s$temperature, ref = "22")
md16s$ph <- relevel(md16s$ph, ref = "5.6")
md16s$day <- as.numeric(as.character(md16s$day))

mens <- brm(ef.16s ~ ph*food*temperature+day,data=md16s, family=Gamma(link = "log"),iter = 5000, chains = 4, cores = 4,control = list(adapt_delta = 0.999, stepsize = 0.001, max_treedepth=20))
saveRDS(mens, "RDS/Exp2_ENS_tubes.RDS")
#mens <- readRDS("RDS/Exp2_ENS_tubes.RDS")
mcmc_plot(mens, regex_pars="b_",  
          prob_outer=0.95,
         prob=0.95)+theme_classic()

posteriormens <- mcmc_intervals_data(mens, 
                                     prob_outer=0.95,
                                     prob=0.5) 

posteriormens$nonzero <- NA
posteriormens$nonzero[posteriormens$ll>0 & posteriormens$hh>0] <- "nonzero"
posteriormens$nonzero[posteriormens$ll<0 & posteriormens$hh<0] <- "nonzero"
posteriormens$nonzero[is.na(posteriormens$nonzero)] <- "zero"
posteriormens<- posteriormens[1:13,]

ggplot(posteriormens, aes(x = parameter,
                          shape=nonzero)) +
  geom_hline(yintercept = 0, linetype = 3, 
             size=1, color = "#b0b5b3") +
  geom_pointrange(aes(ymin = ll, ymax = hh, y = m),
                  position= position_dodge(width=0.75),
                  size = 3/4) +
  scale_shape_manual(values=c(17, 19), 
                     labels=c("95% CI does\nnot contain zero", 
                              "95% CI\ncontains zero"))+
  coord_flip() +
  theme(axis.text.y = element_text( size=7), 
        axis.text.x=element_text(size=7),
        axis.title = element_text(size=7), 
        legend.text = element_text(size=7)) +
  xlab(NULL) +
  ylab("Estimated effect on Effective Number of ASVs")+
  theme(panel.background = element_rect(fill = "white", colour = "black"))

mens.pred <- ggpredict(mens, terms = c("day", "ph", "temperature", "food"))

mens.pred$group <- factor(mens.pred$group, levels = c("3", "4", "5.6"))
ggplot(mens.pred, aes(x=x, y=predicted, color=facet, shape=panel))+geom_point()+
  facet_wrap(~group)+theme_classic()+geom_line()+
  scale_color_manual(values=c("#1f5776", "#c83126"))



###richness
mrich2 <- brm(richness ~ ph*food*temperature+day + (1|sample_id),data=md16s, family=Gamma(link = "log"),iter = 1000, chains = 4, cores = 4)#,control = list(adapt_delta = 0.999, stepsize = 0.001, max_treedepth=20)
#saveRDS(mrich2, "RDS/Exp2_richness_tubes.RDS")
mrich <- readRDS("RDS/Exp2_richness_tubes.RDS")
summary(mrich2)
mcmc_plot(mrich2)
rich.pred <- ggpredict(mrich2, terms = c("food", "temperature", "ph", "day"))
rich.pred$facet <- factor(rich.pred$facet, levels = c("3", "4", "5.6"))
plot(rich.pred, facets = TRUE, line.size=2, dot.size=4, dodge=1) +
  scale_color_manual(values=c("#1f5776", "#c83126"))

ggplot(rich.pred, aes(x = x, y = predicted, shape = x, color = group)) +
  geom_point(size = 5, position = position_dodge(width = .7)) + # Dodging the points
  facet_wrap(~facet) +
  scale_color_manual(values = c("#1f5776", "#c83126")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0, size=2,position = position_dodge(width = .7))+
  theme(legend.position="none")+ylim(c(0, 75))


posteriormrich <- mcmc_intervals_data(mrich2, 
                                     prob_outer=0.95,
                                     prob=0.95)

posteriormrich$nonzero <- NA
posteriormrich$nonzero[posteriormrich$ll>0 & posteriormrich$hh>0] <- "nonzero"
posteriormrich$nonzero[posteriormrich$ll<0 & posteriormrich$hh<0] <- "nonzero"
posteriormrich$nonzero[is.na(posteriormrich$nonzero)] <- "zero"
posteriormrich<- posteriormrich[1:13,]

ggplot(posteriormrich, aes(x = parameter,
                          shape=nonzero)) +
  geom_hline(yintercept = 0, linetype = 3, 
             size=1, color = "#b0b5b3") +
  geom_pointrange(aes(ymin = ll, ymax = hh, y = m),
                  position= position_dodge(width=0.75),
                  size = 3/4) +
  scale_shape_manual(values=c(17, 19), 
                     labels=c("95% CI does\nnot contain zero", 
                              "95% CI\ncontains zero"))+
  coord_flip() +
  theme(axis.text.y = element_text( size=14), 
        axis.text.x=element_text(size=14),
        axis.title = element_text(size=14), 
        legend.text = element_text(size=14)) +
  xlab(NULL) +
  ylab("Estimated effect on Bacterial Richness")+
  theme(panel.background = element_rect(fill = "white", colour = "black"))






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

Exp2.tubes.nomix.physeq3 = subset_samples(Exp2.physeq3, day %in% c("1", "15", "57"))
row_sums <- rowSums(otu_table(Exp2.tubes.nomix.physeq3))
nonzero_rows <- row_sums != 0
Exp2.tubes.nomix.physeq3 <- prune_taxa(nonzero_rows, Exp2.tubes.nomix.physeq3)#51 taxa and 4 samples


#### BETA DIVERSITY ####
## Calculate weighted Unifrac distance and run NMDS
#wu.dist.16s.tubes <- distance(Exp2.tubes.nomix.physeq3,"wUniFrac") not working for some reason
wu.dist.16s.tubes <- phyloseq::distance(Exp2.tubes.nomix.physeq3, method="wunifrac")
tubes.nomix.meta <- data.frame(sample_data(Exp2.tubes.nomix.physeq3))
set.seed(123)
wu.nmds.16s.tubes <- metaMDS(wu.dist.16s.tubes,
                             k = 2, 
                             trymax = 1000,
                             wascores = TRUE)

saveRDS(wu.nmds.16s.tubes, file = "RDS/Exp2.wu.nmds.16s.tubes.RDS")
#wu.nmds.16s.tubes <- readRDS("Exp2.wu.nmds.16s.tubes.RDS")

## Plot NMDS
data.scores2 <- as.data.frame(scores(wu.nmds.16s.tubes$points[,1:2]))
data.scores2$ID <- rownames(data.scores2)
meta.r$ID <- rownames(meta.r)
data_merge2 <- merge(data.scores2, meta.r, by = c("ID"))
data_merge2$food <- ifelse(is.na(data_merge2$food), "inoculant", data_merge2$food)
data_merge2$ph <- ifelse(is.na(data_merge2$ph), "inoculant", data_merge2$ph)
data_merge2$temperature <- ifelse(is.na(data_merge2$temperature), "inoculant", data_merge2$temperature)
data_merge2$day <- as.factor(data_merge2$day)
data_merge2$ph <- as.factor(data_merge2$ph)
data_merge2$food <- as.factor(data_merge2$food)
data_merge2$temperature <- as.factor(data_merge2$temperature)

plot(1, type = "n", xlab="NMDS Axis 1", ylab="NMDS Axis 2", xlim=range(wu.nmds.16s.tubes$points[,1]), ylim=range(wu.nmds.16s.tubes$points[,2]), main="Exp2 16S")
ordiarrows(wu.nmds.16s.tubes, groups = data_merge2$sample_id, startmark = 1, lwd = 1, col="gray")
points(wu.nmds.16s.tubes$points[,1:2], col= c("#CA64A3","#9A4EAE", "#301934", "slategray")[data_merge2$ph], pch= 19, cex= 1)
legend("bottomleft", 
       legend=c("3.0","4.0", "5.6"),
       col= c("#CA64A3","#9A4EAE", "#301934", "slategray"),
       pch= 19,
       cex=1,
       title="pH",
       bty = "n")


plot(1, type = "n", xlab="NMDS Axis 1", ylab="NMDS Axis 2", xlim=range(wu.nmds.16s.tubes$points[,1]), ylim=range(wu.nmds.16s.tubes$points[,2]), main="Exp2 16S")
ordiarrows(wu.nmds.16s.tubes, groups = data_merge2$sample_id, startmark = 1, lwd = 1, col="gray")
points(wu.nmds.16s.tubes$points[,1:2], col= c("#90ee90", "#228822")[data_merge2$food], pch= 19, cex= 1)
legend("bottomleft", 
       legend=c("3g/L","6g/L"),
       col= c("#90ee90", "#228822"),
       title="Food",
       pch= 19,
       cex=1,
       bty = "n")


plot(1, type = "n", xlab="NMDS Axis 1", ylab="NMDS Axis 2", xlim=range(wu.nmds.16s.tubes$points[,1]), ylim=range(wu.nmds.16s.tubes$points[,2]), main="Exp2 16S")
ordiarrows(wu.nmds.16s.tubes, groups = data_merge2$sample_id, startmark = 1, lwd = 1, col="gray")
points(wu.nmds.16s.tubes$points[,1:2], col= c("#1f5776", "#c83126")[data_merge2$temperature], pch= 19, cex= 1)
legend("bottomleft", 
       legend=c("22°C","37°C"),
       col= c("#1f5776", "#c83126"),
       pch= 19,
       cex=1,
       title="Temperature",
       bty = "n")

data_merge2$day<- as.factor(data_merge2$day)
ggplot(data_merge2, aes(x = MDS1, y = MDS2, color = temperature, shape=food)) +
  geom_point(size = 4, aes(alpha = day)) + 
  scale_color_manual(values = c("#1f5776", "#c83126", "black")) + 
  scale_alpha_discrete(range = c(1,0.1)) +  # Adjust range as needed
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
        legend.text = element_text(size = 12, face = "bold", colour = "black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 14), 
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
        legend.title = element_text(size = 14, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2), 
        legend.key = element_blank()) + 
  labs(x = "NMDS1", y = "NMDS2", color = "Temperature", shape="Food", alpha = "Day")


ggplot(data_merge2, aes(x = MDS1, y = MDS2, color = ph, shape=food)) +
  geom_point(size = 4, aes(alpha = day)) + 
  scale_color_manual(values = c("#CA64A3","#9A4EAE", "#301934", "black")) + 
  scale_alpha_discrete(range = c(1,0.1)) +  # Adjust range as needed
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
        legend.text = element_text(size = 12, face = "bold", colour = "black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 14), 
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
        legend.title = element_text(size = 14, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2), 
        legend.key = element_blank()) + 
  labs(x = "NMDS1", y = "NMDS2", color = "pH", shape="Food", alpha = "Day")


#### dbRDA beta diversity ####
#wu.dist.16s.tubes
#tubes.nomix.meta
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
plot(dbrda_beta)

beta_dbrda_anova <- anova.cca(dbrda_beta, by = "onedf", perm = 999)
beta_dbrda_anova
write.csv(beta_dbrda_anova, "beta_16s_dbrda_cca_aov.csv", row.names=TRUE)
#Model: dbrda(formula = as.dist(wu.dist.16s.tubes) ~ temperature * ph * food + day + Condition(sample_id), data = tubes.nomix.meta, distance = "NULL")
#Df SumOfSqs       F Pr(>F)    
#temperature37             1   1.4385 19.7392  0.001 ***
#ph4                       1   0.0331  0.4547  0.824    
#ph3                       1   0.5707  7.8311  0.001 ***
#food6                     1   0.2706  3.7130  0.007 ** 
#day                       1   2.9019 39.8200  0.001 ***
#temperature37:ph4         1   0.0366  0.5017  0.821    
#temperature37:ph3         1   0.1362  1.8684  0.097 .  
#temperature37:food6       1   0.2198  3.0165  0.014 *  
#ph4:food6                 1   0.0216  0.2960  0.939    
#ph3:food6                 1   0.0367  0.5039  0.800    
#temperature37:ph4:food6   1   0.0255  0.3501  0.920    
#temperature37:ph3:food6   1   0.0232  0.3179  0.934    
#Residual                130   9.4739 
significant_factors <- c("temperature37","food6", "ph3","day","temperature37:food6")

beta_est <- as.data.frame(cbind(x1 = dbrda_beta$CCA$biplot[,1], y1 = dbrda_beta$CCA$biplot[,2]))
beta_est$factor <- row.names(dbrda_beta$CCA$biplot)
beta_est$significant <- beta_est$factor %in% significant_factors
significant_beta_est <- beta_est[beta_est$significant, ]

sites_beta <- as.data.frame(beta_sum_dbrda$sites) 
sites_beta_meta <- cbind(sites_beta,tubes.nomix.meta )
sites_beta_meta$temperature <- as.factor(sites_beta_meta$temperature)

temperature_palette <- colorRampPalette(c("#1f5776", "#c83126"))(length(unique(sites_beta_meta$temperature)))
temperature_colors <- temperature_palette[as.numeric(factor(sites_beta_meta$temperature))]
par(mar = c(5, 5, 3, 2)) # setting figure parameters
plot(dbRDA2 ~ dbRDA1, data = sites_beta_meta, pch = 20, type = "p", cex = 1.5, col = temperature_colors)
arrows(x0 = rep(0, nrow(significant_beta_est)), y0 = rep(0, nrow(significant_beta_est)), 
       x1 = significant_beta_est$x1, y1 = significant_beta_est$y1, lwd = 1, col = adjustcolor("black"))
text(x = significant_beta_est$x1 , y = significant_beta_est$y1 , labels = significant_beta_est$factor)
legend("topright", legend = levels(factor(sites_beta_meta$temperature)), fill = temperature_palette, 
       title = "Temperature", cex = 0.8)

sites_beta_meta$day <- as.factor(sites_beta_meta$day)
ggplot(data = sites_beta_meta, aes(x = dbRDA1, y = dbRDA2, color = temperature, alpha=day)) +
  geom_point(size = 3) +
  scale_color_manual(values = c("#1f5776", "#c83126")) +
  geom_segment(data = significant_beta_est, aes(x = 0, y = 0, xend = x1, yend = y1),
               arrow = arrow(length = unit(0.2, "cm")), color = "black", size = 1, inherit.aes = FALSE) +
  geom_text(data = significant_beta_est, aes(x = x1 * 1.2, y = y1 * 1.2, label = factor), 
            hjust = 0, vjust = 0, inherit.aes = FALSE)+ theme_classic()+scale_alpha_discrete(range = c(1,.2))



temperature_palette <- colorRampPalette(c("gray", "slategray", "black"))(length(unique(sites_beta_meta$day)))
temperature_colors <- temperature_palette[as.numeric(factor(sites_beta_meta$day))]
par(mar = c(5, 5, 3, 2)) # setting figure parameters
plot(dbRDA2 ~ dbRDA1, data = sites_beta_meta, pch = 20, type = "p", cex = 1.5, col = temperature_colors)
arrows(x0 = rep(0, nrow(significant_beta_est)), y0 = rep(0, nrow(significant_beta_est)), 
       x1 = significant_beta_est$x1, y1 = significant_beta_est$y1, lwd = 1, col = adjustcolor("black"))
text(x = significant_beta_est$x1 * 1.2, y = significant_beta_est$y1 * 1.2, labels = significant_beta_est$factor)
legend("topright", legend = levels(factor(sites_beta_meta$day)), fill = temperature_palette, 
       title = "Day", cex = 0.8)


sites_beta_meta$day <- as.factor(sites_beta_meta$day)
sites_beta_meta$temperature <- as.factor(sites_beta_meta$temperature)
sites_beta_meta$ph <- factor(sites_beta_meta$ph, levels=c(3, 4, 5.6))
sites_beta_meta$food <- factor(sites_beta_meta$food, levels=c(3, 6))

ggplot(data = sites_beta_meta, aes(x = dbRDA1, y = dbRDA2, color = temperature, alpha=day, shape=day)) +
  geom_point(size = 2) +scale_shape_manual(values = c(16, 17, 15)) +
  scale_color_manual(values = c("#1f5776", "#c83126")) +
  geom_segment(data = significant_beta_est, aes(x = 0, y = 0, xend = x1, yend = y1),
               arrow = arrow(length = unit(.2, "cm")), color = "black", size = 1, inherit.aes = FALSE) +
  geom_text(data = significant_beta_est, aes(x = x1, y = y1, label = factor), 
            hjust = 0, vjust = 0, inherit.aes = FALSE)+ theme_classic()+scale_alpha_discrete(range = c(1,.7))+
  theme(legend.title=element_blank())

ggplot(data = sites_beta_meta, aes(x = dbRDA1, y = dbRDA2, color = ph, alpha=day, shape=day)) +
  geom_point(size = 2) +scale_shape_manual(values = c(16, 17, 15)) +
  scale_color_manual(values = c("#CA64A3","#9A4EAE", "#301934")) +
  geom_segment(data = significant_beta_est, aes(x = 0, y = 0, xend = x1, yend = y1),
               arrow = arrow(length = unit(.2, "cm")), color = "black", size = 1, inherit.aes = FALSE) +
  geom_text(data = significant_beta_est, aes(x = x1, y = y1, label = factor), 
            hjust = 0, vjust = 0, inherit.aes = FALSE)+ theme_classic()+scale_alpha_discrete(range = c(1,.7))+
  theme(legend.title=element_blank())

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
wu.dist.16s.tubes_noNA <- distance(Exp2.tubes.physeq3_noNA,"wUniFrac")

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












# beta disper to check dispersion between treatments
wu.bd.temperature <- betadisper(wu.dist.16s.tubes, data_merge2$temperature)
wu.bd.temperature
boxplot(wu.bd.temperature)
anova(wu.bd.temperature)# p=0.4644 , no significant differences in the dispersion between temperatures

wu.bd.food <- betadisper(wu.dist.16s.tubes, data_merge2$food)
wu.bd.food
boxplot(wu.bd.food)
anova(wu.bd.food)# p=0.2952, no significant differences in the dispersion between food

wu.bd.ph <- betadisper(wu.dist.16s.tubes, data_merge2$ph)
wu.bd.ph
boxplot(wu.bd.ph)
anova(wu.bd.ph)# p=0.006923, significant differences in the dispersion between pH

wu.bd.day <- betadisper(wu.dist.16s.tubes, data_merge2$day)
wu.bd.day
boxplot(wu.bd.day)
anova(wu.bd.day)# p=0.0003871 ***, significant differences in the dispersion between time

#### PERMANOVA for categorical variables (factors) ####
# set number of permutations
ad.16s.treat2 <- adonis2(wu.dist.16s.tubes ~ temperature*ph*food + day, data=data_merge2, by="terms")
ad.16s.treat2

write.csv(ad.16s.treat2, "output_files/exp2_beta_inter.csv")  


#### PHYLOGENETIC ANALYSES ####
##MIX
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
#68 taxa and 1 sample

treemix <- phy_tree(Exp2.physeqmix)


# Plot the phylogenetic tree using ggtree
ggtree(treemix) +
  geom_tiplab() +
  theme_tree2()

tax_table_df <- data.frame(tax_table(Exp2.physeqmix))
length(unique(tax_table_df$Genus))
length(unique(tax_table_df$Family))
length(unique(tax_table_df$Phylum))


tree_data <- fortify(treemix)
tree_data <- merge(tree_data, tax_table_df, by.x = "label", by.y = "row.names", all.x = TRUE)
newasvtaxnames <- read.csv("exported-files/taxonomy_NEW_ASV.csv", header=TRUE)
colnames(tree_data)[1] <- "Feature.ID"
tree_data_asvnames <- left_join(tree_data, newasvtaxnames, by="Feature.ID")


ggtree(tree_data_asvnames)+
  geom_tiplab(aes(label=ASV,color = Family)) +
  theme_tree2() + scale_color_discrete(name = "Family")+
  geom_cladelab(data = tree_data_asvnames,mapping = aes(node = node, label = Family,color = Family), 
 offset = 1, offset.text=0.5)+theme(legend.position = "none")+
  geom_hilight(data=tree_data_asvnames[1:45,], aes(node=node, fill=Family),
               type = "roundrect")

ggtree(tree_data_asvnames, branch.length = "none")+
  geom_tiplab(aes(label=ASV,color = Family), align=TRUE) +
  theme_tree2() + scale_color_discrete(name = "Family")+
  geom_cladelab(data = tree_data_asvnames,mapping = aes(node = node, label = Family,color = Family), 
                offset = 1, offset.text=0.5, align=TRUE)+theme(legend.position = "none")+
  geom_treescale(width = 0.25, offset = 0.5)


ggtree(tree_data_asvnames)+
  geom_tiplab(aes(label=ASV,color = Order)) +
  theme_tree2() +
  geom_cladelab(data = tree_data_asvnames,mapping = aes(node = node, label = Order,color = Order), 
                offset = 1, offset.text=0.5)+theme(legend.position = "none")+
  geom_hilight(data=tree_data_asvnames[1:45,], aes(node=node, fill=Order),
               type = "roundrect")


ggtree(tree_data_asvnames)+
  geom_tiplab(aes(label=ASV,color = Class)) +
  theme_tree2() +
  geom_cladelab(data = tree_data_asvnames,mapping = aes(node = node, label = Class,color = Class), 
                offset = 1, offset.text=0.5)+theme(legend.position = "none")+
  geom_hilight(data=tree_data_asvnames[1:45,], aes(node=node, fill=Class),
               type = "roundrect")

ggtree(tree_data_asvnames, layout="circular")+
  geom_tiplab(aes(label=ASV,color = Family)) +
  theme_tree2() + scale_color_discrete(name = "Family")+
  geom_hilight(data=tree_data_asvnames[1:45,], aes(node=node, fill=Family),
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


## SES PD
SESPD <- ses.pd(com, tree, null.model = "taxa.labels",runs = 999)
saveRDS(SESPD, "RDS/SES_PD_all.RDS")
SESPD <- readRDS("RDS/SES_PD_all.RDS")
SESPD$ID <- rownames(SESPD)
SESPD <- left_join(SESPD, meta, by="ID")
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

ggplot(stat_pd, aes(x=food, y=mean_pd, color=temperature, shape=food))+geom_point(size=3)+
  facet_wrap(~ph)+theme_classic()+ ylim(c(-3, 3))+scale_color_manual(values = c("#1f5776", "#c83126"))


ggplot(SESPDday57, aes(x=ntaxa, y=pd.obs.z, color=temperature, shape=food, group=sample_id))+geom_point(size=3)+
  theme_classic()+scale_color_manual(values = c("#1f5776", "#c83126"))+
  geom_hline(yintercept=0, linetype="dashed", color="black")+
  facet_wrap(~day, nrow=1)+geom_line()+ylim(c(-4, 4))+
  geom_hline(yintercept = -1.96, linetype="dashed", color="gray")+geom_hline(yintercept = 1.96, linetype="dashed", color="gray")

SESPDday57$ph <- factor(SESPDday57$ph, levels=c("5.6", "3", "4"))

pd_model <- brm(pd.obs.z ~ temperature * food * ph ,data=SESPDday57, iter = 10000, chains = 4, cores = 4)
saveRDS(pd_model, "RDS/pd_model.RDS")

mcmc_plot(pd_model, regex_pars="b_",
          prob_outer=0.95,
          prob=0.95)+theme_classic()+scale_color_manual(values="black")

pd.pred <- ggpredict(pd_model, c("food", "temperature", "ph"))
pd.pred$facet <- factor(pd.pred$facet, levels=c("3", "4", "5.6"))

ggplot(pd.pred, aes(x = x, y = predicted, shape = x, color = group)) +
  geom_point(size = 5, position = position_dodge(width = .7)) + # Dodging the points
  facet_wrap(~facet) +
  scale_color_manual(values = c("#1f5776", "#c83126")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0, size=2,position = position_dodge(width = .7))+
  theme(legend.position="none")+ ylab("SES PD")+
  xlab("Food (g/L)")+ theme(
    text = element_text(color = "black", size = 16 ))+ylim(c(-4, 4))+
  geom_hline(yintercept = 0, linetype="dashed")+
  geom_hline(yintercept = -1.96, linetype="dashed", color="gray")+geom_hline(yintercept = 1.96, linetype="dashed", color="gray")


pd.pred1 <- ggpredict(pd_model, "food")
pd.pred1$treatment <- "food"

pd.pred2 <- ggpredict(pd_model, "ph")
pd.pred2$treatment <- "ph"
pd.pred2$x <- as.factor(pd.pred2$x)
pd.pred2 <- pd.pred2 %>%
  mutate(x = as.factor(ifelse(x == 3, "3.0", ifelse(x == 4, "4.0", as.character(x)))))
pd.pred2<- pd.pred2[1:3,]

pd.pred3 <- ggpredict(pd_model, "temperature")
pd.pred3$treatment <- "temperature"

pd.pred.all <- rbind(pd.pred1, pd.pred2, pd.pred3)
pd.pred.all$x <- as.factor(pd.pred.all$x)

pd_filt <- SESPD[,c(6, 15, 16, 17)]
pd_filt$food <- factor(pd_filt$food, levels = c("3", "6"))
pd_filt$ph <- factor(pd_filt$ph, levels = c("3", "4", "5.6"))
pd_filt <- pd_filt %>%
  mutate(ph = as.factor(ifelse(ph == 3, "3.0", ifelse(ph == 4, "4.0", as.character(ph)))))

pd_filt$ph <- factor(pd_filt$ph, levels = c("3.0", "4.0", "5.6"))

pd_melted <- melt(pd_filt, id.vars = "pd.obs.z", 
                  measure.vars = c("food", "ph", "temperature"),
                  variable.name = "treatment", 
                  value.name = "x")

names(pd_melted)[names(pd_melted) == "pd.obs.z"] <- "predicted"

pd_melted$x <- factor(pd_melted$x, levels = c("3", "6", "3.0", "4.0", "5.6", "22", "37"))

pd.pred.all$x <- factor(pd.pred.all$x, levels = c("3", "6", "3.0", "4.0", "5.6", "22", "37"))

ggplot(pd_melted, mapping = aes(x = x, y = predicted, group=treatment, color=x))+
  scale_color_manual(values=c( "#90ee90","#228822", "#CA64A3", "#9A4EAE","#301934", "#1f5776", "#c83126"))+
  geom_jitter(size=3) + geom_point(data=pd.pred.all, aes(x=x, y=predicted, group=treatment), color="black", size=5) +
  geom_linerange(data=pd.pred.all,aes(ymin=conf.low, ymax=conf.high), size=2, color="black",
                 position=position_dodge(width = 0.5)) +facet_wrap(~treatment, scales="free_x")+
  theme_classic()+ylab("SES PD")+theme(legend.position="none")+ylim(c(-3.5,3.5))+geom_hline(yintercept = 0, linetype = "dashed", color = "black")


posteriorpd <- mcmc_intervals_data(pd_model, 
                                     prob_outer=0.95,
                                     prob=0.5)

posteriorpd$nonzero <- NA
posteriorpd$nonzero[posteriorpd$ll>0 & posteriorpd$hh>0] <- "nonzero"
posteriorpd$nonzero[posteriorpd$ll<0 & posteriorpd$hh<0] <- "nonzero"
posteriorpd$nonzero[is.na(posteriorpd$nonzero)] <- "zero"
posteriorpd<- posteriorpd[1:12,]

ggplot(posteriorpd, aes(x = parameter,
                          shape=nonzero)) +
  geom_hline(yintercept = 0, linetype = 3, 
             size=1, color = "#b0b5b3") +
  geom_pointrange(aes(ymin = ll, ymax = hh, y = m),
                  position= position_dodge(width=0.75),
                  size = 3/4) +
  scale_shape_manual(values=c(17, 19), 
                     labels=c("95% CI does\nnot contain zero", 
                              "95% CI\ncontains zero"))+
  coord_flip() +
  xlab(NULL) +
  ylab("Estimated effect on SESPD")+
  theme(
    text = element_text(color = "black", size = 16 ))+
  theme(panel.background = element_rect(fill = "white", colour = "black"))+
  theme(legend.position="none")















## SES MPD
phydist <- cophenetic(tree)
SESMPD <- ses.mpd(com, phydist, null.model = "taxa.labels",
                          abundance.weighted = TRUE, runs = 999)
saveRDS(SESMPD, "RDS/ses.mpd.resultall.RDS")
#saveRDS(SESMPD, "RDS/ses.mpd.result.RDS")
SESMPD <- readRDS("RDS/ses.mpd.resultall.RDS")

SESMPD$ID <- rownames(SESMPD)
SESMPD <- left_join(SESMPD, meta, by="ID")
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
ggplot(SESMPD, aes(x=day, y=mpd.obs.z, color=temperature, shape=food, group=sample_id))+geom_point()+
  theme_classic()+scale_color_manual(values = c("#1f5776", "#c83126"))+
  geom_hline(yintercept=0, linetype="dashed", color="black")+
  facet_wrap(~ph)+geom_line()+ylim(c(-4, 4))+
  geom_hline(yintercept = -1.96, linetype="dashed", color="gray")+geom_hline(yintercept = 1.96, linetype="dashed", color="gray")



ggplot(SESMPD, aes(x=ntaxa, y=mpd.obs.z, color=temperature, shape=food))+geom_jitter(size=3)+
  facet_wrap(~day, nrow=1)+theme_classic()+
  scale_color_manual(values = c("#1f5776", "#c83126"))+
  ylim(c(-4,4))+geom_hline(yintercept = 0, linetype="dashed")+
  geom_hline(yintercept = -1.96, linetype="dashed", color="gray")+geom_hline(yintercept = 1.96, linetype="dashed", color="gray")

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

ggplot(stat_mpd, aes(x=food, y=mean_mpd, color=temperature, shape=food))+geom_point(size=3)+
  facet_wrap(~ph)+theme_classic()+ ylim(c(-3, 3))+scale_color_manual(values = c("#1f5776", "#c83126"))





mpd_model <- brm(mpd.obs.z ~ temperature * food * ph ,data=SESMPDday57, iter = 10000, chains = 4, cores = 4)
saveRDS(mpd_model, "mpd_model.RDS")

mcmc_plot(mpd_model, regex_pars="b_",
          prob_outer=0.95,
          prob=0.95)+theme_classic()


mpd.pred <- ggpredict(mpd_model, c("food", "temperature", "ph"))

ggplot(mpd.pred, aes(x = x, y = predicted, shape = x, color = group)) +
  geom_point(size = 5, position = position_dodge(width = .7)) + # Dodging the points
  facet_wrap(~facet) +
  scale_color_manual(values = c("#1f5776", "#c83126")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0, size=2,position = position_dodge(width = .7))+
  theme(legend.position="none")+ ylab("SES MPD")+
  xlab("Food (g/L)")+ theme(
    text = element_text(color = "black", size = 16 ))+ylim(c(-4, 4))+
  geom_hline(yintercept = 0, linetype="dashed")+
  geom_hline(yintercept = -1.96, linetype="dashed", color="gray")+geom_hline(yintercept = 1.96, linetype="dashed", color="gray")







mpd.pred1 <- ggpredict(mpd_model, "food")
mpd.pred1$treatment <- "food"

mpd.pred2 <- ggpredict(mpd_model, "ph")
mpd.pred2$treatment <- "ph"
mpd.pred2$x <- as.factor(mpd.pred2$x)
mpd.pred2 <- mpd.pred2 %>%
  mutate(x = as.factor(ifelse(x == 3, "3.0", ifelse(x == 4, "4.0", as.character(x)))))
mpd.pred2<- mpd.pred2[1:3,]

mpd.pred3 <- ggpredict(mpd_model, "temperature")
mpd.pred3$treatment <- "temperature"

mpd.pred.all <- rbind(mpd.pred1, mpd.pred2, mpd.pred3)
mpd.pred.all$x <- as.factor(mpd.pred.all$x)

mpd_filt <- SESMPD[,c(6, 15, 16, 17)]
mpd_filt$food <- factor(mpd_filt$food, levels = c("3", "6"))
mpd_filt$ph <- factor(mpd_filt$ph, levels = c("3", "4", "5.6"))
mpd_filt <- mpd_filt %>%
  mutate(ph = as.factor(ifelse(ph == 3, "3.0", ifelse(ph == 4, "4.0", as.character(ph)))))

mpd_filt$ph <- factor(mpd_filt$ph, levels = c("3.0", "4.0", "5.6"))

mpd_melted <- melt(mpd_filt, id.vars = "mpd.obs.z", 
                  measure.vars = c("food", "ph", "temperature"),
                  variable.name = "treatment", 
                  value.name = "x")

names(mpd_melted)[names(mpd_melted) == "mpd.obs.z"] <- "predicted"

mpd_melted$x <- factor(mpd_melted$x, levels = c("3", "6", "3.0", "4.0", "5.6", "22", "37"))

mpd.pred.all$x <- factor(mpd.pred.all$x, levels = c("3", "6", "3.0", "4.0", "5.6", "22", "37"))

ggplot(mpd_melted, mapping = aes(x = x, y = predicted, group=treatment, color=x))+
  scale_color_manual(values=c( "#90ee90","#228822", "#CA64A3", "#9A4EAE","#301934", "#1f5776", "#c83126"))+
  geom_jitter(size=3) + geom_point(data=mpd.pred.all, aes(x=x, y=predicted, group=treatment), color="black", size=5) +
  geom_linerange(data=mpd.pred.all,aes(ymin=conf.low, ymax=conf.high), size=2, color="black",
                 position=position_dodge(width = 0.5)) +facet_wrap(~treatment, scales="free_x")+
  theme_classic()+ylab("SES MPD")+theme(legend.position="none")+ylim(c(-3.5,3.5))+geom_hline(yintercept = 0, linetype = "dashed", color = "black")

posteriormpd <- mcmc_intervals_data(mpd_model, 
                                   prob_outer=0.95,
                                   prob=0.5)

posteriormpd$nonzero <- NA
posteriormpd$nonzero[posteriormpd$ll>0 & posteriormpd$hh>0] <- "nonzero"
posteriormpd$nonzero[posteriormpd$ll<0 & posteriormpd$hh<0] <- "nonzero"
posteriormpd$nonzero[is.na(posteriormpd$nonzero)] <- "zero"
posteriormpd<- posteriormpd[1:12,]

ggplot(posteriormpd, aes(x = parameter,
                        shape=nonzero)) +
  geom_hline(yintercept = 0, linetype = 3, 
             size=1, color = "#b0b5b3") +
  geom_pointrange(aes(ymin = ll, ymax = hh, y = m),
                  position= position_dodge(width=0.75),
                  size = 3/4) +
  scale_shape_manual(values=c(17, 19), 
                     labels=c("95% CI does\nnot contain zero", 
                              "95% CI\ncontains zero"))+
  coord_flip() +
  xlab(NULL) +
  ylab("Estimated effect on SESMPD")+
  theme(
    text = element_text(color = "black", size = 16 ))+
  theme(panel.background = element_rect(fill = "white", colour = "black"))+
  theme(legend.position="none")

#### SES MNTD ####

SESMNTD <- ses.mntd(com, phydist, null.model = "taxa.labels",
                           abundance.weighted = TRUE, runs = 999)
saveRDS(SESMNTD, "RDS/ses.mntd.resultall.RDS")

#saveRDS(SESMNTD, "RDS/ses.mntd.result.RDS")
SESMNTD <- readRDS("RDS/ses.mntd.resultall.RDS")

SESMNTD$ID <- rownames(SESMNTD)
SESMNTD <- left_join(SESMNTD, meta, by="ID")
SESMNTD$food <- as.factor(SESMNTD$food)
SESMNTD$temperature <- as.factor(SESMNTD$temperature)
SESMNTD$ph <- as.factor(SESMNTD$ph)

ggplot(SESMNTD, aes(x=food, y=mntd.obs.z, color=temperature, shape=food))+geom_jitter()+
  facet_wrap(~ph)+theme_classic()+scale_color_manual(values = c("#1f5776", "#c83126"))



SESMNTD <- SESMNTD %>%
  mutate(ph = replace(ph, is.na(ph), c(5.6)))

SESMNTD <- SESMNTD %>%
  mutate(food = replace(food, is.na(food), c(3)))

SESMNTD <- SESMNTD %>%
  mutate(temperature = replace(temperature, is.na(temperature), c(22)))
SESMNTD$day <- as.numeric(SESMNTD$day)


ggplot(SESMNTD, aes(x=day, y=mntd.obs.z, color=temperature, shape=food, group=sample_id))+geom_point()+
  theme_classic()+scale_color_manual(values = c("#1f5776", "#c83126"))+
  geom_hline(yintercept=0, linetype="dashed", color="black")+
  facet_wrap(temperature~ph)+geom_line()+ylim(c(-3.5, 3.5))



ggplot(SESMNTD, aes(x=ntaxa, y=mntd.obs.z, color=temperature, shape=food))+geom_jitter(size=3)+
  facet_wrap(~day, nrow=1)+theme_classic()+
  scale_color_manual(values = c("#1f5776", "#c83126"))+
  ylim(c(-4,4))+geom_hline(yintercept = 0, linetype="dashed")+
  geom_hline(yintercept = -1.96, linetype="dashed", color="gray")+geom_hline(yintercept = 1.96, linetype="dashed", color="gray")


SESMNTDdaymix <- subset(SESMNTD, sample_id =="MIX")
mean(SESMNTDdaymix$mntd.obs.z)#-1.347052
sd(SESMNTDdaymix$mntd.obs.z)/sqrt(length(SESMNTDdaymix$mntd.obs.z))# 0.0851035



SESMNTDday57 <- subset(SESMNTD, day=="57")

SESMNTDday57$ph <- factor(SESMNTDday57$ph, levels=c("5.6", "3", "4"))




stat_mntd <- SESMNTDday57 %>% group_by(treatment_combo) %>% 
  mutate(mean_mntd =mean(mntd.obs.z))

#mean overall is -0.9528

ggplot(stat_mntd, aes(x=food, y=mean_mntd, color=temperature, shape=food))+geom_point(size=3)+
  facet_wrap(~ph)+theme_classic()+ ylim(c(-3, 3))+scale_color_manual(values = c("#1f5776", "#c83126"))


count_greater_than_1.98 <- sum(abs(SESMNTDday57$mntd.obs.z) > 1.98)
total_numbers <- length(SESMNTDday57$mntd.obs.z)
(count_greater_than_1.98 / total_numbers) * 100

#0



mntd_model <- brm(mntd.obs.z ~ temperature * food * ph,data=SESMNTDday57, iter = 10000, chains = 4, cores = 4)
saveRDS(mntd_model, "mntd_model.RDS")
mcmc_plot(mntd_model, regex_pars="b_",
          prob_outer=0.95,
          prob=0.95)+theme_classic()


mntd.pred <- ggpredict(mntd_model, c("food", "temperature", "ph"))

ggplot(mntd.pred, aes(x = x, y = predicted, shape = x, color = group)) +
  geom_point(size = 5, position = position_dodge(width = .7)) + # Dodging the points
  facet_wrap(~facet) +
  scale_color_manual(values = c("#1f5776", "#c83126")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0, size=2,position = position_dodge(width = .7))+
  theme(legend.position="none")+ ylab("SES MNTD")+
  xlab("Food (g/L)")+ theme(
    text = element_text(color = "black", size = 16 ))+ylim(c(-4, 4))+
  geom_hline(yintercept = 0, linetype="dashed")+
  geom_hline(yintercept = -1.96, linetype="dashed", color="gray")+geom_hline(yintercept = 1.96, linetype="dashed", color="gray")















mntd.pred1 <- ggpredict(mntd_model, "food")
mntd.pred1$treatment <- "food"

mntd.pred2 <- ggpredict(mntd_model, "ph")
mntd.pred2$treatment <- "ph"
mntd.pred2$x <- as.factor(mntd.pred2$x)
mntd.pred2 <- mntd.pred2 %>%
  mutate(x = as.factor(ifelse(x == 3, "3.0", ifelse(x == 4, "4.0", as.character(x)))))
mntd.pred2<- mntd.pred2[1:3,]

mntd.pred3 <- ggpredict(mntd_model, "temperature")
mntd.pred3$treatment <- "temperature"

mntd.pred.all <- rbind(mntd.pred1, mntd.pred2, mntd.pred3)
mntd.pred.all$x <- as.factor(mntd.pred.all$x)

mntd_filt <- SESMNTD[,c(6, 15, 16, 17)]
mntd_filt$food <- factor(mntd_filt$food, levels = c("3", "6"))
mntd_filt$ph <- factor(mntd_filt$ph, levels = c("3", "4", "5.6"))
mntd_filt <- mntd_filt %>%
  mutate(ph = as.factor(ifelse(ph == 3, "3.0", ifelse(ph == 4, "4.0", as.character(ph)))))

mntd_filt$ph <- factor(mntd_filt$ph, levels = c("3.0", "4.0", "5.6"))

mntd_melted <- melt(mntd_filt, id.vars = "mntd.obs.z", 
                   measure.vars = c("food", "ph", "temperature"),
                   variable.name = "treatment", 
                   value.name = "x")

names(mntd_melted)[names(mntd_melted) == "mntd.obs.z"] <- "predicted"

mntd_melted$x <- factor(mntd_melted$x, levels = c("3", "6", "3.0", "4.0", "5.6", "22", "37"))

mntd.pred.all$x <- factor(mntd.pred.all$x, levels = c("3", "6", "3.0", "4.0", "5.6", "22", "37"))

ggplot(mntd_melted, mapping = aes(x = x, y = predicted, group=treatment, color=x))+
  scale_color_manual(values=c( "#90ee90","#228822", "#CA64A3", "#9A4EAE","#301934", "#1f5776", "#c83126"))+
  geom_jitter(size=3) + geom_point(data=mntd.pred.all, aes(x=x, y=predicted, group=treatment), color="black", size=5) +
  geom_linerange(data=mntd.pred.all,aes(ymin=conf.low, ymax=conf.high), size=2, color="black",
                 position=position_dodge(width = 0.5)) +facet_wrap(~treatment, scales="free_x")+
  theme_classic()+ylab("SES MNTD")+theme(legend.position="none")+ylim(c(-3.5,3.5))+geom_hline(yintercept = 0, linetype = "dashed", color = "black")

mcmc_plot(mntd_model)
posteriormntd <- mcmc_intervals_data(mntd_model, 
                                    prob_outer=0.95,
                                    prob=0.5)

posteriormntd$nonzero <- NA
posteriormntd$nonzero[posteriormntd$ll>0 & posteriormntd$hh>0] <- "nonzero"
posteriormntd$nonzero[posteriormntd$ll<0 & posteriormntd$hh<0] <- "nonzero"
posteriormntd$nonzero[is.na(posteriormntd$nonzero)] <- "zero"
posteriormntd<- posteriormntd[1:12,]

ggplot(posteriormntd, aes(x = parameter,
                         shape=nonzero)) +
  geom_hline(yintercept = 0, linetype = 3, 
             size=1, color = "#b0b5b3") +
  geom_pointrange(aes(ymin = ll, ymax = hh, y = m),
                  position= position_dodge(width=0.75),
                  size = 3/4) +
  scale_shape_manual(values=c(17, 19), 
                     labels=c("95% CI does\nnot contain zero", 
                              "95% CI\ncontains zero"))+
  coord_flip() +
  xlab(NULL) +
  ylab("Estimated effect on SESMNTD")+
  theme(
    text = element_text(color = "black", size = 16 ))+
  theme(panel.background = element_rect(fill = "white", colour = "black"))+
  theme(legend.position="none")








#Negative values of SESmpd and SESmntd indicate lower phylogenetic diversity than expected under the assumption of the null model, whereas values greater than zero indicate higher phylogenetic diversity than predicted by the null model.



### try nti again with relative abundance
cm_relative <- cm / rowSums(cm)



#### SES BetaMNTD ####
#Ran the betaMNTD code on the cluster, see github for scripts
#Read in model object
beta_nti<- readRDS("betamntd/beta-NTI-relative-abundance.RDS")
beta_nti <- as.data.frame(beta_nti) #148X148

#look at differences between pairwise comparisons
#DAY1
pairwise_dfday1 <- beta_nti
rows_with_day1 <- grepl("DAY1_", row_names)
cols_with_day1 <- grepl("DAY1_", col_names)
pairwise_dfday1<- pairwise_dfday1[rows_with_day1, cols_with_day1]#48X48

pairwise_dfday1$rowname <- rownames(pairwise_dfday1)
pairwise_dfday1_melt<- melt(pairwise_dfday1, id.vars = "rowname", value.name = "BetaNTI")
colnames(pairwise_dfday1_melt)[2] <- "columnname"
pairwise_dfday1_melt <- na.omit(pairwise_dfday1_melt)

pairwise_dfday1_meltR <- pairwise_dfday1_melt[,c(1,3)]
colnames(pairwise_dfday1_meltR)[1] <- "ID"
pairwise_dfday1_meltRmeta<- left_join(pairwise_dfday1_meltR, meta.r, by="ID" )
pairwise_dfday1_meltRmeta<- pairwise_dfday1_meltRmeta[,c(1,2,7:10)]

pairwise_dfday1_meltC <- pairwise_dfday1_melt[,c(2,3)]
colnames(pairwise_dfday1_meltC)[1] <- "ID"
pairwise_dfday1_meltCmeta<- left_join(pairwise_dfday1_meltC, meta.r, by="ID" )
pairwise_dfday1_meltCmeta<- pairwise_dfday1_meltCmeta[,c(1,2,7:10)]

pairwise_dfday1_merge <- as.data.frame(cbind(pairwise_dfday1_meltRmeta, pairwise_dfday1_meltCmeta))
colnames(pairwise_dfday1_merge)[7] <- "ID.1"
colnames(pairwise_dfday1_merge)[8] <- "BetaNTI.1"
colnames(pairwise_dfday1_merge)[9] <- "treatment_combo.1"
colnames(pairwise_dfday1_merge)[10] <- "food.1"
colnames(pairwise_dfday1_merge)[11] <- "ph.1"
colnames(pairwise_dfday1_merge)[12] <- "temperature.1"

pairwise_dfday1_merge$temperature <- as.numeric(pairwise_dfday1_merge$temperature)
pairwise_dfday1_merge$temperature.1 <- as.numeric(as.character(pairwise_dfday1_merge$temperature.1))
pairwise_dfday1_merge$ph <- as.numeric(pairwise_dfday1_merge$ph)
pairwise_dfday1_merge$ph.1 <- as.numeric(pairwise_dfday1_merge$ph.1)
pairwise_dfday1_merge$food <- as.numeric(pairwise_dfday1_merge$food)
pairwise_dfday1_merge$food.1 <- as.numeric(pairwise_dfday1_merge$food.1)

pairwise_dfday1_merge$temp_dif <- abs(pairwise_dfday1_merge$temperature - pairwise_dfday1_merge$temperature.1)
pairwise_dfday1_merge$ph_dif <- abs(pairwise_dfday1_merge$ph - pairwise_dfday1_merge$ph.1)
pairwise_dfday1_merge$food_dif <- abs(pairwise_dfday1_merge$food - pairwise_dfday1_merge$food.1)


ggplot(pairwise_dfday1_merge, aes(x=temp_dif, y=BetaNTI)) + geom_jitter()
ggplot(pairwise_dfday1_merge, aes(x=ph_dif, y=BetaNTI)) + geom_jitter()
ggplot(pairwise_dfday1_merge, aes(x=food_dif, y=BetaNTI, color=temperature, shape=as.factor(food))) + geom_jitter()



#DAY15
pairwise_dfday15 <- beta_nti
rows_with_day15 <- grepl("DAY15", row_names)
cols_with_day15 <- grepl("DAY15", col_names)
pairwise_dfday15<- pairwise_dfday15[rows_with_day15, cols_with_day15]#48x48

pairwise_dfday15$rowname <- rownames(pairwise_dfday1)
pairwise_dfday15_melt<- melt(pairwise_dfday15, id.vars = "rowname", value.name = "BetaNTI")
colnames(pairwise_dfday15_melt)[2] <- "columnname"
pairwise_dfday15_melt <- na.omit(pairwise_dfday15_melt)

pairwise_dfday15_meltR <- pairwise_dfday15_melt[,c(1,3)]
colnames(pairwise_dfday15_meltR)[1] <- "ID"
pairwise_dfday15_meltRmeta<- left_join(pairwise_dfday15_meltR, meta.r, by="ID" )
pairwise_dfday15_meltRmeta<- pairwise_dfday15_meltRmeta[,c(1,2,7:10)]

pairwise_dfday15_meltC <- pairwise_dfday15_melt[,c(2,3)]
colnames(pairwise_dfday15_meltC)[1] <- "ID"
pairwise_dfday15_meltCmeta<- left_join(pairwise_dfday15_meltC, meta.r, by="ID" )
pairwise_dfday15_meltCmeta<- pairwise_dfday15_meltCmeta[,c(1,2,7:10)]

pairwise_dfday15_merge <- as.data.frame(cbind(pairwise_dfday15_meltRmeta, pairwise_dfday15_meltCmeta))
colnames(pairwise_dfday15_merge)[7] <- "ID.1"
colnames(pairwise_dfday15_merge)[8] <- "BetaNTI.1"
colnames(pairwise_dfday15_merge)[9] <- "treatment_combo.1"
colnames(pairwise_dfday15_merge)[10] <- "food.1"
colnames(pairwise_dfday15_merge)[11] <- "ph.1"
colnames(pairwise_dfday15_merge)[12] <- "temperature.1"

pairwise_dfday15_merge$temperature <- as.numeric(pairwise_dfday15_merge$temperature)
pairwise_dfday15_merge$temperature.1 <- as.numeric(as.character(pairwise_dfday15_merge$temperature.1))
pairwise_dfday15_merge$ph <- as.numeric(pairwise_dfday15_merge$ph)
pairwise_dfday15_merge$ph.1 <- as.numeric(pairwise_dfday15_merge$ph.1)
pairwise_dfday15_merge$food <- as.numeric(pairwise_dfday15_merge$food)
pairwise_dfday15_merge$food.1 <- as.numeric(pairwise_dfday15_merge$food.1)

pairwise_dfday15_merge$temp_dif <- abs(pairwise_dfday15_merge$temperature - pairwise_dfday15_merge$temperature.1)
pairwise_dfday15_merge$ph_dif <- abs(pairwise_dfday15_merge$ph - pairwise_dfday15_merge$ph.1)
pairwise_dfday15_merge$food_dif <- abs(pairwise_dfday15_merge$food - pairwise_dfday15_merge$food.1)


ggplot(pairwise_dfday15_merge, aes(x=temp_dif, y=BetaNTI)) + geom_jitter()
ggplot(pairwise_dfday15_merge, aes(x=ph_dif, y=BetaNTI)) + geom_jitter()
ggplot(pairwise_dfday15_merge, aes(x=food_dif, y=BetaNTI)) + geom_jitter()




#DAY57
pairwise_dfday57 <- beta_nti
rows_with_day57 <- grepl("DAY57", row_names)
cols_with_day57 <- grepl("DAY57", col_names)
pairwise_dfday57<- pairwise_dfday57[rows_with_day57, cols_with_day57]#48x48

pairwise_dfday57$rowname <- rownames(pairwise_dfday57)
pairwise_dfday57_melt<- melt(pairwise_dfday57, id.vars = "rowname", value.name = "BetaNTI")
colnames(pairwise_dfday57_melt)[2] <- "columnname"
pairwise_dfday57_melt <- na.omit(pairwise_dfday57_melt)

pairwise_dfday57_meltR <- pairwise_dfday57_melt[,c(1,3)]
colnames(pairwise_dfday57_meltR)[1] <- "ID"
pairwise_dfday57_meltRmeta<- left_join(pairwise_dfday57_meltR, meta.r, by="ID" )
pairwise_dfday57_meltRmeta<- pairwise_dfday57_meltRmeta[,c(1,2,7:10)]

pairwise_dfday57_meltC <- pairwise_dfday57_melt[,c(2,3)]
colnames(pairwise_dfday57_meltC)[1] <- "ID"
pairwise_dfday57_meltCmeta<- left_join(pairwise_dfday57_meltC, meta.r, by="ID" )
pairwise_dfday57_meltCmeta<- pairwise_dfday57_meltCmeta[,c(1,2,7:10)]

pairwise_dfday57_merge <- as.data.frame(cbind(pairwise_dfday57_meltRmeta, pairwise_dfday57_meltCmeta))
colnames(pairwise_dfday57_merge)[7] <- "ID.1"
colnames(pairwise_dfday57_merge)[8] <- "BetaNTI.1"
colnames(pairwise_dfday57_merge)[9] <- "treatment_combo.1"
colnames(pairwise_dfday57_merge)[10] <- "food.1"
colnames(pairwise_dfday57_merge)[11] <- "ph.1"
colnames(pairwise_dfday57_merge)[12] <- "temperature.1"

pairwise_dfday57_merge$temperature <- as.numeric(pairwise_dfday57_merge$temperature)
pairwise_dfday57_merge$temperature.1 <- as.numeric(as.character(pairwise_dfday57_merge$temperature.1))
pairwise_dfday57_merge$ph <- as.numeric(pairwise_dfday57_merge$ph)
pairwise_dfday57_merge$ph.1 <- as.numeric(pairwise_dfday57_merge$ph.1)
pairwise_dfday57_merge$food <- as.numeric(pairwise_dfday57_merge$food)
pairwise_dfday57_merge$food.1 <- as.numeric(pairwise_dfday57_merge$food.1)

pairwise_dfday57_merge$temp_dif <- abs(pairwise_dfday57_merge$temperature - pairwise_dfday57_merge$temperature.1)
pairwise_dfday57_merge$ph_dif <- abs(pairwise_dfday57_merge$ph - pairwise_dfday57_merge$ph.1)
pairwise_dfday57_merge$food_dif <- abs(pairwise_dfday57_merge$food - pairwise_dfday57_merge$food.1)


ggplot(pairwise_dfday57_merge, aes(x=temp_dif, y=BetaNTI)) + geom_jitter()
ggplot(pairwise_dfday57_merge, aes(x=ph_dif, y=BetaNTI)) + geom_jitter()
ggplot(pairwise_dfday57_merge, aes(x=food_dif, y=BetaNTI)) + geom_jitter()



## COMPARED TO MIX (repeat with each day)
pairwise_dfMIX <- beta_nti
rows_with_mix <- grepl("GROUP", row_names)
cols_with_day57 <- grepl("DAY15_", col_names)
pairwise_dfMIX<- pairwise_dfMIX[rows_with_mix, cols_with_day57]#48x4

pairwise_dfMIX$rowname <- rownames(pairwise_dfMIX)
pairwise_dfMIX_melt<- melt(pairwise_dfMIX, id.vars = "rowname", value.name = "BetaNTI")
colnames(pairwise_dfMIX_melt)[2] <- "columnname"
pairwise_dfMIX_melt <- na.omit(pairwise_dfMIX_melt)

colnames(pairwise_dfMIX_melt)[2] <- "ID"
pairwise_dfMIX_meltmeta<- left_join(pairwise_dfMIX_melt, meta.r, by="ID" )

ggplot(pairwise_dfMIX_meltmeta, aes(x=temperature, y=BetaNTI)) + geom_jitter()
ggplot(pairwise_dfMIX_meltmeta, aes(x=ph, y=BetaNTI)) + geom_jitter()
ggplot(pairwise_dfMIX_meltmeta, aes(x=food, y=BetaNTI)) + geom_jitter()

ggplot(pairwise_dfMIX_meltmeta, aes(x=temperature, y=BetaNTI, color=as.factor(temperature), shape=as.factor(food)))+
  geom_jitter(size=2)+
  theme_classic()+ facet_wrap(~ph)+geom_hline(linetype="dashed", yintercept = 0, color="black")+ylim(c(-3,3))+
  scale_color_manual(values=c("#1f5776", "#c83126"))


## COMPARED TO MIX ONE GRAPH
pairwise_dfMIX <- beta_nti
rows_with_mix <- grepl("GROUP", row_names)
cols_with_day57 <- grepl("DAY", col_names)
pairwise_dfMIX<- pairwise_dfMIX[rows_with_mix, cols_with_day57]#4x144

pairwise_dfMIX$rowname <- rownames(pairwise_dfMIX)
pairwise_dfMIX_melt<- melt(pairwise_dfMIX, id.vars = "rowname", value.name = "BetaNTI")
colnames(pairwise_dfMIX_melt)[2] <- "columnname"
pairwise_dfMIX_melt <- na.omit(pairwise_dfMIX_melt)

colnames(pairwise_dfMIX_melt)[2] <- "ID"
pairwise_dfMIX_meltmeta<- left_join(pairwise_dfMIX_melt, meta.r, by="ID" )
pairwise_dfMIX_meltmeta$day <- as.numeric(pairwise_dfMIX_meltmeta$day)

ggplot(pairwise_dfMIX_meltmeta, aes(x=day, y=BetaNTI, group=temperature, color=as.factor(temperature), shape=as.factor(food)))+
  geom_jitter(size=2)+geom_line()+
  theme_classic()+ facet_wrap(~ph)+geom_hline(linetype="dashed", yintercept = 0, color="black")+ylim(c(-3,3))+
  scale_color_manual(values=c("#1f5776", "#c83126"))

ggplot(pairwise_dfMIX_meltmeta, aes(x = day, y = BetaNTI, group = sample_id, color = as.factor(temperature), shape = as.factor(food))) + 
  geom_jitter(size = 2, width = 0.2) +  # Added width to jitter to avoid overplotting
  geom_line() + 
  theme_classic() + 
  facet_grid(temperature ~ ph) +  # Facet by ph and temperature
  geom_hline(linetype = "dashed", yintercept = 0, color = "black") + 
  ylim(c(-3, 3)) + 
  scale_color_manual(values = c("#1f5776", "#c83126")) +
  labs(color = "Temperature", shape = "Food") 


ggplot(pairwise_dfMIX_meltmeta, aes(x = day, y = BetaNTI, group=day))+geom_boxplot()






#write.csv(pairwise_dfMIX_meltmeta, "pairwise_dfMIX_meltmeta.csv")
#manually matched up samples to correct inoculating mix
pairwise_dfMIX_matched <- read.csv("pairwise_dfMIX_meltmeta.csv", header=TRUE)
ggplot(pairwise_dfMIX_matched, aes(x = day, y = BetaNTI, group = sample_id, color = as.factor(temperature), shape = as.factor(food))) + 
  geom_jitter(size = 2, width = 0.2) +  # Added width to jitter to avoid overplotting
  geom_line() + 
  theme_classic() + 
  facet_grid(temperature ~ ph) +  # Facet by ph and temperature
  geom_hline(linetype = "dashed", yintercept = 0, color = "black") + 
  ylim(c(-2, 2)) + 
  scale_color_manual(values = c("#1f5776", "#c83126")) +
  labs(color = "Temperature", shape = "Food") 

pairwise_dfMIX_matched$ph  <- as.character(pairwise_dfMIX_matched$ph)
pairwise_dfMIX_matched$ph <- factor(pairwise_dfMIX_matched$ph, levels = c( "3", "4","5.6"))
pairwise_dfMIX_matched$food <- factor(pairwise_dfMIX_matched$food, levels = c("3", "6"))
pairwise_dfMIX_matched$temperature <- factor(pairwise_dfMIX_matched$temperature, levels = c("22", "37"))

ggplot(pairwise_dfMIX_matched, aes(x = day, y = BetaNTI, group=day))+geom_boxplot()+
  facet_wrap(~temperature)

ggplot(pairwise_dfMIX_matched, aes(x = day, y = BetaNTI, group=day))+geom_boxplot()+
  facet_wrap(temperature~ph)+geom_hline(yintercept = 0, linetype="dashed")+
  theme_classic()

ggplot(pairwise_dfMIX_matched, aes(x = day, y = BetaNTI, group=day))+geom_boxplot()+
  facet_wrap(~food)

pairwise_dfMIX_matched$day <- as.factor(pairwise_dfMIX_matched$day)
ggplot(pairwise_dfMIX_matched, aes(x = temperature, y = BetaNTI, fill=day)) + 
  geom_boxplot() + scale_fill_manual(values=c("gray", "slategray","gray48"))+ 
  geom_smooth(method = "lm", se = FALSE)+ facet_wrap(food~ph) + 
  scale_color_manual(values=c("gray", "slategray","gray48")) +
  theme_minimal() + 
  theme(axis.text.y = element_text(colour = "black", size = 16, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 16), 
        legend.text = element_text(size = 16, face ="bold", colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 16), 
        axis.title.x = element_text(face = "bold", size = 16, colour = "black"), 
        legend.title = element_text(size = 16, colour = "black", face = "bold"),
        legend.key=element_blank()) +
  ylab("Beta NTI") + xlab("")+
  geom_hline(yintercept = 0, linetype="dashed")+ylim(c(-2,2))

hist(pairwise_dfMIX_matched$BetaNTI, main = "Histogram of BetaNTI", xlab = "BetaNTI")
qqnorm(pairwise_dfMIX_matched$BetaNTI)
qqline(pairwise_dfMIX_matched$BetaNTI, col = "red")

# Density Plot
plot(density(pairwise_dfMIX_matched$BetaNTI), main = "Density Plot of BetaNTI")

#variance test
# Levene's Test (from the car package)
if (!require(car)) install.packages("car")
library(car)

leveneTest(BetaNTI ~ interaction(temperature, food, ph, day), data = pairwise_dfMIX_matched)
# Shapiro-Wilk Test
shapiro.test(pairwise_dfMIX_matched$BetaNTI)

aovbnti<- aov(BetaNTI ~ temperature*food*ph*day,data=pairwise_dfMIX_matched)
summary(aovbnti)
tukey_results <- TukeyHSD(aovbnti)
tukey_results
default_priors <- default_prior(BetaNTI ~ temperature * food * ph * day, data = pairwise_dfMIX_matched)
default_priors
priors <- c(
  set_prior("student_t(3, 0, 0.5)", class = "b"),       # Tighter Student-t priors for coefficients
  set_prior("normal(0, 1)", class = "Intercept"),       # Normal prior for the intercept
  set_prior("exponential(1)", class = "sigma")          # Exponential prior for the residual standard deviation
)

metanti_model <- brm(BetaNTI ~ temperature*food*ph*day,data=pairwise_dfMIX_matched, iter = 10000, chains = 4, cores = 4, prior=priors)

pp_check(metanti_model)
mcmc_plot(metanti_model, regex_pars="b_",
          prob_outer=0.95,
          prob=0.95)+theme_classic()

conditional_effects_interaction <- ggpredict(metanti_model, terms = c("day", "food", "ph", "temperature"))
plot(conditional_effects_interaction, facets = TRUE,  line_size=2, dot_size=4, dodge=7, shape="group")



ggplot(resp.pred, aes(x = x, y = predicted, shape = x, color = group)) +
  geom_point(size = 5, position = position_dodge(width = .7)) + # Dodging the points
  facet_wrap(~facet) +
  scale_color_manual(values = c("#1f5776", "#c83126")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0, size=2,position = position_dodge(width = .7))+
  theme(legend.position="none")+ ylab("Predicted Respiration")+
  xlab("Food (g/L)")+ theme(
    text = element_text(color = "black", size = 16 )
  )

if (!require(shinystan)) install.packages("shinystan")
library(shinystan)
launch_shinystan(metanti_model)











n <- nrow(beta_nti)  # Assuming it's a square matrix

# Fill in missing values in the upper triangle using values from the lower triangle
for (i in 1:n) {
  for (j in (i + 1):n) {  # Only loop through the upper triangle excluding the diagonal
    if (is.na(beta_nti[i, j])) {
      beta_nti[i, j] <- beta_nti[j, i]
    }
  }
}
beta_nti[is.na(beta_nti)] <- 0

row_names <- rownames(beta_nti)
col_names <- colnames(beta_nti)

# Identify rows and columns that contain "Day57"
rows_with_day57 <- grepl("DAY57", row_names)
cols_with_day57 <- grepl("DAY57", col_names)
day57_beta_nti<- beta_nti[rows_with_day57, cols_with_day57]


melt_beta<- melt(day57_beta_nti)
colnames(melt_beta)[1] <- "ID"
melt_beta_meta <- left_join(melt_beta, meta.r, by="ID")

melt_beta_meta$food <- as.factor(melt_beta_meta$food)
melt_beta_meta$temperature <- as.factor(melt_beta_meta$temperature)
melt_beta_meta$ph <- as.factor(melt_beta_meta$ph)

ggplot(melt_beta_meta, aes(x=food, y=value))+ geom_jitter()+
  geom_boxplot()+ theme_classic()

ggplot(melt_beta_meta, aes(x=temperature, y=value))+ geom_jitter()+
  geom_boxplot()+ theme_classic()

ggplot(melt_beta_meta, aes(x=ph, y=value))+ geom_jitter()+
  geom_boxplot()+ theme_classic()






































ggplot(melt_beta_meta, mapping = aes(x = x, y = predicted, group=treatment, color=x))+
  scale_color_manual(values=c( "#90ee90","#228822", "#CA64A3", "#9A4EAE","#301934", "#1f5776", "#c83126"))+
  geom_jitter(size=3) + geom_point(data=mntd.pred.all, aes(x=x, y=predicted, group=treatment), color="black", size=5) +
  geom_linerange(data=mntd.pred.all,aes(ymin=conf.low, ymax=conf.high), size=2, color="black",
                 position=position_dodge(width = 0.5)) +facet_wrap(~treatment, scales="free_x")+
  theme_classic()+ylab("SES MNTD")+theme(legend.position="none")+ylim(c(-3.5,3.5))+geom_hline(yintercept = 0, linetype = "dashed", color = "black")


# Identify rows and columns that contain "Day15"
rows_with_day15 <- grepl("DAY15", row_names)
cols_with_day15 <- grepl("DAY15", col_names)
day15_beta_nti<- beta_nti[rows_with_day15, cols_with_day15]

# Identify rows and columns that contain "Day1"
rows_with_day1 <- grepl("DAY1_", row_names)
cols_with_day1 <- grepl("DAY1_", col_names)
day1_beta_nti<- beta_nti[rows_with_day1, cols_with_day1]

# Identify and remove rows and columns that contain "MIX"
rows_with_mix <- grepl("MIX", row_names)
cols_with_mix <- grepl("MIX", col_names)
nomix_beta_nti <- beta_nti[!rows_with_mix, !cols_with_mix]

#### NMDS ####
beta_nti_nmds <- metaMDS(beta_nti,
                             k = 2, 
                             trymax = 1000,
                             wascores = TRUE)


## Plot NMDS
data.scores2 <- as.data.frame(scores(beta_nti_nmds$points[,1:2]))
data.scores2$ID <- rownames(data.scores2)
data_merge2 <- left_join(data.scores2, meta.r, by = c("ID"))
data_merge2$food <- ifelse(is.na(data_merge2$food), "inoculant", data_merge2$food)
data_merge2$ph <- ifelse(is.na(data_merge2$ph), "inoculant", data_merge2$ph)
data_merge2$temperature <- ifelse(is.na(data_merge2$temperature), "inoculant", data_merge2$temperature)
data_merge2$day <- as.factor(data_merge2$day)
data_merge2$ph <- factor(data_merge2$ph, levels=c("3", "4", "5.6"))
data_merge2$food <- factor(data_merge2$food, levels=c("3", "6"))
data_merge2$temperature <- factor(data_merge2$temperature, levels =c("22", "37"))

plot(1, type = "n", xlab="NMDS Axis 1", ylab="NMDS Axis 2", xlim=range(beta_nti_nmds$points[,1]), ylim=range(beta_nti_nmds$points[,2]), main="Beta NTI, temperature")
ordiarrows(beta_nti_nmds, groups = data_merge2$sample_id, startmark = 1, lwd = 1, col="gray")
points(beta_nti_nmds$points[,1:2], col= c("blue", "red")[data_merge2$temperature], pch= 19, cex= 1)

plot(1, type = "n", xlab="NMDS Axis 1", ylab="NMDS Axis 2", xlim=range(beta_nti_nmds$points[,1]), ylim=range(beta_nti_nmds$points[,2]), main="Beta NTI, food")
ordiarrows(beta_nti_nmds, groups = data_merge2$sample_id, startmark = 1, lwd = 1, col="gray")
points(beta_nti_nmds$points[,1:2], col= c("lightgreen","darkgreen")[data_merge2$food], pch= 19, cex= 1)

plot(1, type = "n", xlab="NMDS Axis 1", ylab="NMDS Axis 2", xlim=range(beta_nti_nmds$points[,1]), ylim=range(beta_nti_nmds$points[,2]), main="Beta NTI, ph")
ordiarrows(beta_nti_nmds, groups = data_merge2$sample_id, startmark = 1, lwd = 1, col="gray")
points(beta_nti_nmds$points[,1:2], col= c("pink","magenta", "purple")[data_merge2$ph], pch= 19, cex= 1)


## try bnti with dbrda
#remove MIX samples
data_merge2_nomix<- data_merge2[-(145:148),]
beta_nti_dist <- vegdist(nomix_beta_nti, method="bray")
beta_nti_dbrda <- dbrda(beta_nti_dist ~ ph+food+temperature,
                   data = data_merge2_nomix,
                   distance = "bray")
plot(beta_nti_dbrda)

#day1 samples
data_merge2_day1<- subset(data_merge2, day=="1")
beta_nti_dist <- vegdist(day1_beta_nti, method="bray")
beta_nti_dbrda <- dbrda(beta_nti_dist ~ ph+food+temperature,
                        data = data_merge2_day1,
                        distance = "bray")
plot(beta_nti_dbrda)

#day15 samples
data_merge2_day15<- subset(data_merge2, day=="15")
beta_nti_dist <- vegdist(day15_beta_nti, method="bray")
beta_nti_dbrda <- dbrda(beta_nti_dist ~ ph+food+temperature,
                        data = data_merge2_day15,
                        distance = "bray")
plot(beta_nti_dbrda)

#day57 samples
data_merge2_day57<- subset(data_merge2, day=="57")
beta_nti_dist <- vegdist(day57_beta_nti, method="bray")
beta_nti_dbrda <- dbrda(beta_nti_dist ~ ph+food+temperature,
                        data = data_merge2_day57,
                        distance = "bray")
plot(beta_nti_dbrda)



##anova.cca for the significance





#### beta diversity mix samples ####
## Calculate weighted Unifrac distance and run NMDS
wu.dist.16s.mix <- distance(Exp2.tubes.mix.physeq3,"wUniFrac")

set.seed(123)
wu.nmds.16s.mix <- metaMDS(wu.dist.16s.mix,
                             k = 2, 
                             trymax = 1000,
                             wascores = TRUE)#0.09442731

saveRDS(wu.nmds.16s.mix, file = "Exp2.wu.nmds.16s.mix.RDS")
wu.nmds.16s.mix <- readRDS("Exp2.wu.nmds.16s.mix.RDS")

## Plot NMDS
data.scores2 <- as.data.frame(scores(wu.nmds.16s.mix$points[,1:2]))
data.scores2$ID <- rownames(data.scores2)
meta_mix <- subset(md16s, day %in% c("0", "1"))
data_merge2 <- merge(data.scores2, meta_mix, by = c("ID"))

data_merge2$day <- as.factor(data_merge2$day)

ggplot(data_merge2, aes(x = MDS1, y = MDS2, color=day)) + 
  geom_point(size = 4)+
  labs(x = "NMDS1", y = "NMDS2")+scale_color_manual(values=c("black", "gray"))+
  theme_bw()

#### PERMANOVA for categorical variables (factors) ####
ad.16s.mix <- adonis2(wu.dist.16s.mix ~ X, data=data_merge2, by="margin")
ad.16s.mix



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
output_chit = ancombc2(data = tse, assay_name = "counts", tax_level = NULL,
                       fix_formula = "scaled_chit+scaled_prot",
                       rand_formula = "(1 | day)",#cant use sample_id as random effect bc not enough observations per each factor (48)
                       p_adj_method = "fdr", pseudo_sens = TRUE,
                       prv_cut = 0.02,#4 tubes in 12 treatment groups at three time points, 144 tubes, so set prvcut to 2% meaning in at least 4 tubes
                       alpha = 0.05, n_cl = 3, verbose = TRUE,
                       global = FALSE, pairwise = FALSE)
saveRDS(output_chit, "RDS/ancombc_chit.RDS")
#output_chit <- readRDS("RDS/ancombc_chit.RDS")

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


ggplot(chit_res_sig_filt_melt_names, aes(x = lfc, y = newname2)) +
  geom_point() +
  xlab("Log Fold Change") +
  geom_errorbarh(aes(xmin = lfc - se, xmax = lfc + se), height = 0.2) +
  theme_classic()+geom_vline(xintercept = 0, linetype="dashed")

prot_res_sig_filt <- prot_res_sig[,c(1,3,5)]
colnames(prot_res_sig_filt) <- c("ASV", "lfc", "se")

prot_res_sig_filt_melt_names <- dplyr::left_join(prot_res_sig_filt, newasv, by="ASV")

prot_res_sig_filt_melt_names <- prot_res_sig_filt_melt_names %>%
  arrange(lfc)
taxon_levels <- unique(prot_res_sig_filt_melt_names$newname2)
prot_res_sig_filt_melt_names$newname2 <- factor(prot_res_sig_filt_melt_names$newname2, levels = taxon_levels)

dim(prot_res_sig_filt_melt_names)
dim(chit_res_sig_filt_melt_names)

#combine them
chit_res_sig_filt_melt_names$df <- "chitinase"
prot_res_sig_filt_melt_names$df <- "protease"
df <- rbind(chit_res_sig_filt_melt_names, prot_res_sig_filt_melt_names)

 ggplot(df, aes(x = newname2, y = lfc, color = newname2)) +
   geom_point(shape=1) +
   ylab("Log Fold Change") +
   geom_errorbar(aes(ymin = lfc - se, ymax = lfc + se), width = 0.2) +
   theme_classic() +
   geom_hline(yintercept = 0, linetype = "dashed")+
   theme(axis.text.x = element_text(angle = 45, hjust=1))+
   theme(legend.position="none")
 
 ggplot(chit_res_sig_filt_melt_names, aes(x = newname2, y = lfc)) +
   geom_point(color = "blue") +
   geom_errorbar(aes(ymin = lfc - se, ymax = lfc + se), width = 0.2, color = "blue") +
   geom_hline(yintercept = 0, linetype = "dashed") +
   ylab("Primary Y-Axis (Log Fold Change)") +
   theme_classic()+
   theme(axis.text.x = element_text(angle = 45, hjust=1))+
   theme(legend.position="none")
 
ggplot(prot_res_sig_filt_melt_names, aes(x = newname2, y = lfc)) +
   geom_point(color = "red") +
   geom_errorbar(aes(ymin = lfc - se, ymax = lfc + se), width = 0.2, color = "red") +
   geom_hline(yintercept = 0, linetype = "dashed") +
   ylab("Secondary Y-Axis") +
   theme_classic() +
   theme(axis.title.y = element_text(color = "red"),
         axis.text.y = element_text(color = "red"))+
   theme(axis.text.x = element_text(angle = 45,hjust=1))+
   theme(legend.position="none")
 
ch_anc <- chit_res_sig_filt_melt_names
ch_anc$enzyme <- "chitinase"
pr_anc <- prot_res_sig_filt_melt_names
pr_anc$enzyme <- "protease"
enz_anc <- rbind(ch_anc, pr_anc)

ggplot(ch_anc, aes(x = enzyme, y = newname2, fill = lfc)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "#efe0a8", mid="#dfd8d2",high = "#5f555a")  +
  theme_minimal()+theme(text = element_text(color = "black", size = 10 ))+
  xlab("")+ ylab("")

ggplot(pr_anc, aes(x = enzyme, y = newname2, fill = lfc)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "#007FFF",high = "#FF8000")  +
  theme_minimal()+theme(text = element_text(color = "black", size = 10 ))+
  xlab("")+ ylab("")

ggplot(enz_anc, aes(x = enzyme, y = newname2, fill = lfc)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "#88B04B",mid= "#FFF8DC",high = "#FFC107",limits = c(-2, 2))+
  theme_minimal() +
  theme(text = element_text(color = "black")) + xlab("") + ylab("")

write.csv(df, "ANCOMBC_enzyme.csv", row.names=TRUE)

#### differentially abundant taxa between treatments
output_temp = ancombc2(data = tse, assay_name = "counts", tax_level = NULL,
                       fix_formula = "food + ph + temperature",
                       rand_formula = "(1 | day)",
                       p_adj_method = "fdr", pseudo_sens = TRUE,
                       prv_cut = 0.02,#4 tubes in 12 treatment groups at three time points, 144 tubes, so set prvcut to 2% meaning in at least 4 tubes
                       group = "temperature",
                       alpha = 0.05, n_cl = 3, verbose = TRUE,
                       global = TRUE, pairwise = TRUE)

saveRDS(output_temp, "RDS/ancombc.RDS")
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
ggplot(all.filt.melt, aes(x = variable, y = newname2, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low="yellow",mid = "white", high = "black") +
  labs(title = "Heatmap of Taxa vs Variables", x = "Variables", y = "Taxa") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

##temp
temp.filt <- as.data.frame(temp_res_sig[,c(1,6)])
colnames(temp.filt)[1] <- "ASV"
temp.filt.name <- dplyr::left_join(temp.filt, newasv, by="ASV")
temp.filt.name.filt <- temp.filt.name[,c(2,5)]

temp.filt.melt <- reshape2::melt(temp.filt.name.filt, id.vars = "newname2")
temp.filt.melt$newname2 <- factor(temp.filt.melt$newname2, levels = unique(temp.filt.melt$newname2[order(temp.filt.melt$value, decreasing = TRUE)]))

ggplot(temp.filt.melt, aes(x = variable, y = newname2, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "#063648", mid="white",high = "#A33232")  +
  theme_minimal()+theme(text = element_text(color = "black", size = 10 ))+
  xlab("")+ ylab("")

##pH

ph.filt <- as.data.frame(ph_res_sig[,c(1,5)])
colnames(ph.filt)[1] <- "ASV"
ph.filt.name <- left_join(ph.filt, newasv, by="ASV")
ph.filt.name.filt <- ph.filt.name[,c(2,5)]

ph.filt.melt <- reshape2::melt(ph.filt.name.filt, id.vars = "newname2")
ph.filt.melt$newname2 <- factor(ph.filt.melt$newname2, levels = unique(ph.filt.melt$newname2[order(ph.filt.melt$value, decreasing = TRUE)]))

ggplot(ph.filt.melt, aes(x = variable, y = newname2, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(high = "#CA64A3", mid="white",low = "#301934") +
  labs(title = "Heatmap of Taxa vs Variables", x = "Variables", y = "Taxa") +
  theme_minimal()+theme(text = element_text(color = "black", size = 10 ))+
  xlab("")+ ylab("")

##food
food.filt <- as.data.frame(food_res_sig[,c(1,3)])
colnames(food.filt)[1] <- "ASV"
food.filt.name <- left_join(food.filt, newasv, by="ASV")
food.filt.name.filt <- food.filt.name[,c(2,5)]

food.filt.melt <- reshape2::melt(food.filt.name.filt, id.vars = "newname2")
food.filt.melt$newname2 <- factor(food.filt.melt$newname2, levels = unique(food.filt.melt$newname2[order(food.filt.melt$value, decreasing = TRUE)]))

ggplot(food.filt.melt, aes(x = variable, y = newname2, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(high = "#228822", mid="white",low = "#90ee90") +
  labs(title = "Heatmap of Taxa vs Variables", x = "Variables", y = "Taxa") +
  theme_minimal()+theme(text = element_text(color = "black", size = 16 ))+
  xlab("")+ ylab("")


#taxonomy
tax <- as.data.frame(tax.16s2.r)
tax <- cbind(Feature.ID=rownames(tax),tax)
tax.p <- parse_taxonomy(tax)
tax.p$ASV <- row.names(tax.p)

#asvs
asv16s.anc <- asv16s.r
asv16s.anc <- as.data.frame(asv16s.anc)
asv16s.anc$ASV <- row.names(asv16s.anc)

#ancombc
colnames(all_res_sig)[1] <- "ASV"
ancom_tax <- merge(all_res_sig, tax.p, by="ASV")
ancom_tax_count <- merge(ancom_tax, asv16s.anc, by="ASV")

##########need to fix below
ancom_tax_count_filt <- ancom_tax_count[,c(1,44:249)]
ancom_tax_count_filt.melt <- ancom_tax_count_filt %>% melt(id=c( "ASV"))#,
tube_meta<- subset(meta, container == "tube")
tube_meta<- subset(tube_meta, day %in% c("1", "15", "57"))
tube_meta$variable <- rownames(tube_meta)

ancom_tax.filt.melt_meta <- left_join(ancom_tax_count_filt.melt, tube_meta, by="variable")
ancom_tax.filt.melt_meta <- subset (ancom_tax.filt.melt_meta, container == "tube")
ancom_tax.filt.melt_meta <- subset (ancom_tax.filt.melt_meta, day %in% c("1", "15", "57"))

ancom_tax.filt.melt_meta$temperature <- as.factor(ancom_tax.filt.melt_meta$temperature)
ancom_tax.filt.melt_meta$ph <- as.factor(ancom_tax.filt.melt_meta$ph)
ancom_tax.filt.melt_meta$food <- as.factor(ancom_tax.filt.melt_meta$food)
ancom_tax.filt.melt_meta$treatment_combo <- as.factor(ancom_tax.filt.melt_meta$treatment_combo)

##sort by temperature
ancom_tax.filt.melt_meta <- ancom_tax.filt.melt_meta %>%
  arrange(day, temperature, ph, food)

ancom_tax.filt.melt_meta$value <- log2(ancom_tax.filt.melt_meta$value + 0.5)

ancom_tax.filt.melt_meta2 <- ancom_tax.filt.melt_meta[,c(1,2,3,6,9,10,11)]


#get new names on heatmap


ancommerge <- left_join(ancom_tax.filt.melt_meta2, newasv, by = "ASV")
ancommerge.filt <- ancommerge[,c(2:7, 10)]

ancom_tax.filt.wide.meta <- ancommerge.filt %>%
  pivot_wider(
    id_cols = c(variable,day, temperature, ph, food), # Columns to keep as identifiers
    names_from = newname2, # Column to spread into wide format
    values_from = value # Column containing values corresponding to the ASVs
  )

ancom_all_meta <- ancom_tax.filt.wide.meta[,1:5]
ancom_all_meta <- as.data.frame(tube_meta)
rownames(ancom_all_meta) <- ancom_all_meta$variable
ancom_all_meta <- ancom_all_meta[,-1]
ancom_all_meta$ph <- as.factor(ancom_all_meta$ph)
ancom_all_meta$temperature <- as.factor(ancom_all_meta$temperature)
ancom_all_meta$food <- as.factor(ancom_all_meta$food)



ancom_all_matrix <- ancom_tax.filt.wide.meta[,c(1,6:35)]
ancom_all_matrix <- as.data.frame(ancom_all_matrix)
rownames(ancom_all_matrix) <- ancom_all_matrix$variable
ancom_all_matrix <- ancom_all_matrix[,-1]
ancom_all_matrix.t <- t(ancom_all_matrix)

pheatmap(ancom_all_matrix.t,
         cluster_cols = FALSE,cluster_rows = TRUE,
         col = colorRampPalette(c("white", "gray","black"))(40), border_color=NA,
         cex=1.5, fontsize=4, annotation_col = ancom_all_meta, gaps_col = c(48, 96))



#### RELATIVE ABUNDANCE PLOTS TUBES ####
tax <- data.frame(tax_table(Exp2.tubes.nomix.physeq3))

#filter to look at genus NA
tax_genus_na <- tax %>% filter(is.na(Genus))
library(Biostrings)
fasta_file <- readDNAStringSet("exported-files/dna-sequences.fasta")
asv_ids <- tax_genus_na$ASV
filtered_fasta <- fasta_file[names(fasta_file) %in% asv_ids]
writeXStringSet(filtered_fasta, filepath = "dna_seq_genus_na.fasta")

samp <- data.frame(sample_data(Exp2.tubes.nomix.physeq3))
otu <- data.frame(otu_table(Exp2.tubes.nomix.physeq3))
rownames(otu)==rownames(tax)

jewel_tone_palette <- c("#000000","#a9a19c","#493829","#117a65","#bdd09f", "#0c2461", "#d35400",
                        "#1f8a13","#8f3b1b","#855723","#e67e22","#740058",
                        "palegoldenrod","#3498db","#f39c12", "#28b463","#af7ac5",
                        "#1f618d","#d35400","#154360","darkmagenta")
# Plotting the relative abundance at Genus level
otu.g <- data.frame(Genus=tax$Genus,otu)


otu.g$Genus[is.na(otu.g$Genus)] <- "Unknown"
otu.ga <- aggregate(. ~ otu.g$Genus, otu.g[,2:ncol(otu.g)], sum)
row.names(otu.ga) <- otu.ga[,1]
otu.ga <- otu.ga[,2:ncol(otu.ga)]
otu.gao <- as.matrix(otu.ga[order(rowSums(otu.ga),decreasing = T),])

otu.gaop <- otu.gao[c(2:21),]#top 20 genera unknown is the top row, so I am skipping it and doing 2:21 and I'll put unknown in the other column
other <- colSums(otu.gao[c(1,22:nrow(otu.gao)),])#put the all the other genera in a "other" category
otu.gaopo <- rbind(otu.gaop, other)
otu.gaopo_gg <- as.data.frame(otu.gaopo)
otu.gaopo_ggp <- adorn_percentages(otu.gaopo_gg,, 1:ncol(otu.gaopo_gg), denominator="col")#number of columns (samples), calculate RA
other <- otu.gaopo_ggp[21,]
otu.gaopo_ggp <- otu.gaopo_ggp[-21,]
otu.gaopo_ggp <- otu.gaopo_ggp[order(rowSums(otu.gaopo_ggp),decreasing = F),]
otu.gaopo_ggp <- rbind(other, otu.gaopo_ggp)
otu.gaopo_ggp[ "Genus" ] <- rownames(otu.gaopo_ggp)
otu.gaopo_ggp$Genus <- factor(otu.gaopo_ggp$Genus, levels = otu.gaopo_ggp$Genus)
otu.gaopo_ggp.rmelt<- reshape2::melt(otu.gaopo_ggp, id.vars="Genus", value.name="Relative_Abundance", variable.name="Sample")

samp.new <- samp
samp.new$Sample <- rownames(samp.new)
otu.gaopo_ggp.rmelt.meta <- left_join(otu.gaopo_ggp.rmelt, samp.new, by="Sample")

otu.gaopo_ggp.rmelt.meta$day <- as.factor(otu.gaopo_ggp.rmelt.meta$day)
otu.gaopo_ggp.rmelt.meta <- otu.gaopo_ggp.rmelt.meta %>%
  dplyr::arrange(treatment_combo, day)

sample_order <- unique(otu.gaopo_ggp.rmelt.meta$Sample[order(otu.gaopo_ggp.rmelt.meta$day)])

# Convert Sample to a factor with the new order
otu.gaopo_ggp.rmelt.meta$Sample <- factor(otu.gaopo_ggp.rmelt.meta$Sample, levels = sample_order)

# Create the plot
ggplot(otu.gaopo_ggp.rmelt.meta, aes(x = Sample, y = Relative_Abundance, fill = Genus)) + 
  geom_bar(position = "fill", stat = "identity", width = 1) + 
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 90, size = 5, hjust = 1, vjust = 0.5)) +
  guides(fill = guide_legend(ncol = 7)) + 
  theme(legend.position = "bottom") +
  facet_wrap(~treatment_combo, scales = "free_x", nrow = 2)+
  scale_fill_manual(values = jewel_tone_palette)


# Plotting the relative abundance at ASV level
otu.g <- data.frame(ASV=tax$ASV,otu)
otu.g$ASV[is.na(otu.g$ASV)] <- "Unknown"
otu.ga <- aggregate(. ~ otu.g$ASV, otu.g[,2:ncol(otu.g)], sum)
row.names(otu.ga) <- otu.ga[,1]
otu.ga <- otu.ga[,2:ncol(otu.ga)]
otu.gao <- as.matrix(otu.ga[order(rowSums(otu.ga),decreasing = T),])

otu.gaop <- otu.gao[c(1:30),]#top 20 genera unknown is the top row, so I am skipping it and doing 2:21 and I'll put unknown in the other column
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


samp.new <- samp
samp.new$Sample <- rownames(samp.new)
otu.gaopo_ggp.rmelt.meta <- left_join(otu.gaopo_ggp.rmelt, samp.new, by="Sample")

otu.gaopo_ggp.rmelt.meta$day <- as.factor(otu.gaopo_ggp.rmelt.meta$day)
otu.gaopo_ggp.rmelt.meta <- otu.gaopo_ggp.rmelt.meta %>%
  dplyr::arrange(treatment_combo, day)


sample_order <- unique(otu.gaopo_ggp.rmelt.meta$Sample[order(otu.gaopo_ggp.rmelt.meta$day)])
otu.gaopo_ggp.rmelt.meta$Sample <- factor(otu.gaopo_ggp.rmelt.meta$Sample, levels = sample_order)
customcol28 <- colors <- c("black","royalblue3","darkblue","#DDDD77","#117744",
                           "#AAAA44","darkred","purple","mediumblue","palegoldenrod",
                           "lightgoldenrod","#77CCCC","yellow","purple4","darkgreen","#AA4455",
                           "lightsalmon","yellow3","purple2","lightblue","firebrick",
                           "navy","red4","skyblue","darkmagenta","mediumvioletred",
                           "#771155", "#44AAAA", "#CC99BB", "#114477",  "#777711")


ggplot(otu.gaopo_ggp.rmelt.meta, aes(x=Sample, y=Relative_Abundance, fill = ASV)) + 
  geom_bar( position = "fill", stat = "identity", width=1) + 
  theme_classic()+ theme(axis.text.x = element_text(angle = 90, size=5))+
  guides(fill = guide_legend(ncol = 5)) + theme(legend.position="bottom") +
  facet_wrap(~treatment_combo, scales = "free_x", nrow = 2)+ scale_fill_manual(values = customcol28)


#### make RA ASV graph again but with new ASV names ####
newasvRA<- read.csv("exported-files/taxonomy_NEW_ASV.csv", header=TRUE)
newasv <- newasvRA$ASV
newasvRA <- parse_taxonomy(newasvRA)
newasvRA <- cbind(newasvRA,newasv )
newasvRA$newname2 <- paste(newasvRA$newasv, " - ", newasvRA$Genus)
rownames(otu)==rownames(newasvRA)
#991 vs 3342
newasvRA_subset <- newasvRA[rownames(otu), , drop = FALSE]
identical(rownames(otu), rownames(newasvRA_subset))
rownames(otu)==rownames(newasvRA_subset)
otu$ASV <- rownames(otu)
newasvRA_subset$ASV <- rownames(newasvRA_subset)

otu.g <- left_join(otu, newasvRA_subset, by="ASV")
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


samp.new <- samp
samp.new$Sample <- rownames(samp.new)
otu.gaopo_ggp.rmelt.meta <- left_join(otu.gaopo_ggp.rmelt, samp.new, by="Sample")

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

### look at acinetobacter ###
acineto_tax <- tax.16s2.r[apply(tax.16s2.r, 1, function(x) any(grepl("Acinetobacter", x))), ]

acineto_tax$asv <- rownames(acineto_tax)
acineto_com <- as.data.frame(asv16s.r)
acineto_com$asv <- rownames(acineto_com)
acineto_combo <- left_join(acineto_tax, acineto_com, by="asv")
acineto_long <- acineto_combo %>% 
  pivot_longer(cols=-c(Taxon, Confidence, asv),
               names_to="sample_id", values_to="counts")
acineto_meta<- meta.r
acineto_meta$sample_id <- rownames(acineto_meta)

acineto_all <- left_join(acineto_long, acineto_meta, by="sample_id")
acineto_all <- acineto_all %>% filter(container == "tube")

ggplot(acineto_all, aes(x=day, y=counts, color=asv))+geom_jitter()+
  facet_wrap(temperature~ph)

ggplot(acineto_all, aes(x=day, y=log(0.1+counts), color=asv, label=sample_id))+geom_jitter()+
  facet_wrap(temperature~ph)+geom_text(size=2.5, color="black")

baff_all <- acineto_all[apply(acineto_all, 1, function(x) any(grepl("baff", x))), ]

baff_all$temperature <- as.factor(baff_all$temperature)
baff_all$food <- as.factor(baff_all$food)


ggplot(baff_all, aes(x=day, y=log(0.1+counts), color=temperature, shape=food, label=sample_id))+geom_jitter(size=3)+
  facet_wrap(~ph)+theme_classic()+geom_text(size=2.5, color="black")

ggplot(baff_all, aes(x=chitinase, y=log(0.1+counts), color=temperature, shape=food, label=sample_id))+geom_jitter(size=3)+
  facet_wrap(~ph)+theme_classic()+geom_text(size=2, color="black")

baff_day57 <- baff_all %>% filter(day == "57")
ggplot(baff_day57, aes(x=chitinase, y=log(0.1+counts), color=temperature, shape=food, label=sample_id))+geom_jitter(size=3)+
  facet_wrap(~ph)+theme_classic()+
  scale_color_manual(values = c("#1f5776", "#c83126"))+geom_text(size=2, color="black")



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





### by phylum
tax_olig_cop <- tax
copiotrophs <- c("Proteobacteria", "Firmicutes", "Actinobacteriota", "Bacteroidota", "Bdellovibrionota")
oligotrophs <- c("Microsporidia", "Cyanobacteria", "Planctomycetota", "Verrucomicrobiota", "Acidobacteriota", "Desulfobacterota", "Chloroflexi", "Euryarchaeota", "Crenarchaeota")

# Create new column based on Phylum categories
tax_olig_cop <- tax_olig_cop %>%
  mutate(c_o = if_else(Phylum %in% copiotrophs, "C", 
                       if_else(Phylum %in% oligotrophs, "O", NA_character_)))

otu.g <- data.frame(c_o=tax_olig_cop$c_o,otu)
otu.g$c_o[is.na(otu.g$c_o)] <- "Unknown"
otu.ga <- aggregate(. ~ otu.g$c_o, otu.g[,2:ncol(otu.g)], sum)
row.names(otu.ga) <- otu.ga[,1]
otu.ga <- otu.ga[,2:ncol(otu.ga)]
otu.gao <- as.matrix(otu.ga[order(rowSums(otu.ga),decreasing = T),])


otu.gaopo_gg <- as.data.frame(otu.gao)
otu.gaopo_ggp <- adorn_percentages(otu.gaopo_gg,, 1:ncol(otu.gaopo_gg), denominator="col")#number of columns (samples), calculate RA

otu.gaopo_ggp[ "c_o" ] <- rownames(otu.gaopo_ggp)
otu.gaopo_ggp$c_o <- factor(otu.gaopo_ggp$c_o, levels = otu.gaopo_ggp$c_o)
otu.gaopo_ggp.rmelt<- reshape2::melt(otu.gaopo_ggp, id.vars="c_o", value.name="Relative_Abundance", variable.name="Sample")


samp.new <- samp
samp.new$Sample <- rownames(samp.new)
otu.gaopo_ggp.rmelt.meta <- left_join(otu.gaopo_ggp.rmelt, samp.new, by="Sample")

otu.gaopo_ggp.rmelt.meta$day <- as.factor(otu.gaopo_ggp.rmelt.meta$day)
otu.gaopo_ggp.rmelt.meta <- otu.gaopo_ggp.rmelt.meta %>%
  dplyr::arrange(treatment_combo, day)


sample_order <- unique(otu.gaopo_ggp.rmelt.meta$Sample[order(otu.gaopo_ggp.rmelt.meta$day)])
otu.gaopo_ggp.rmelt.meta$Sample <- factor(otu.gaopo_ggp.rmelt.meta$Sample, levels = sample_order)

ggplot(otu.gaopo_ggp.rmelt.meta, aes(x=Sample, y=Relative_Abundance, fill = c_o)) + 
  geom_bar( position = "fill", stat = "identity", width=1) + 
  theme_classic()+ theme(axis.text.x = element_text(angle = 90, size=5))+
  scale_fill_manual(values = c( "navy","orange", "slategray")) +
  guides(fill = guide_legend(ncol = 7)) + theme(legend.position="bottom") +
  facet_wrap(~treatment_combo, scales = "free_x", nrow = 1)

#### TAX SUMMARY ####
abundance_data <- t(as.data.frame(otu_table(Exp2.physeq3)))

# Get the sample data
sample_data <- as.data.frame(sample_data(Exp2.physeq3))

# Combine abundance data with sample data
combined_data <- cbind(abundance_data, sample_data)

summarized_data <- combined_data %>%
  group_by(treatment_combo) %>%
  summarize(across(1:566, sum, na.rm = TRUE))
summarized_data <- as.data.frame(summarized_data)

rownames(summarized_data) <- summarized_data$treatment_combo

summarized_data <- summarized_data[,-1]
summarized_data.t <- t(summarized_data)
summarized_data.t <- as.data.frame(summarized_data.t)

summarized_data.t$ASV <- rownames(summarized_data.t)
tax.new<-tax 
tax.new$ASV <- rownames(tax.new)
tax.new <-as.data.frame(tax.new)
summarized_data.t.tax <- left_join(summarized_data.t, tax.new, by="ASV")                       
                                   
summarized_data.t.tax.filt <- summarized_data.t.tax[,c(1:12, 20)]                           
summarized_data.t.tax.filt <- summarized_data.t.tax.filt %>%
  group_by(Genus) %>%
  summarize(across(1:12, sum, na.rm = TRUE))               

melt.ra <- melt(summarized_data.t.tax.filt)

ggplot(melt.ra, aes(x=Genus, y=value)) + geom_boxplot()+
  facet_wrap(~variable)



### Community Function ####

#### ECOPLATE ANALYSIS ####
#read in the data
final_raw <- read.csv("Data/2022_Ecoplates_EXP2_end.csv", header = TRUE, row.names = 1, check.names = FALSE)
start_raw <- read.csv("Data/2022_Ecoplates_Exp2_start.csv", header = TRUE, row.names = 1, check.names = FALSE)

#meta data
meta_wk8 <- read.csv("Data/2022_Ecoplates_EXP2_end_META.csv", header = TRUE, check.names = FALSE)
meta_wk1 <- read.csv("Data/2022_Ecoplates_EXP2_start_META.csv", header = TRUE, check.names = FALSE)

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

#remove metadata
eco_wk0 <- subset(eco_all, wk == "1")
eco_wk0.nometa <- subset(eco_wk0, select = -c(sample_id, ph, temperature, food, container, wk))

eco_wk7 <- subset(eco_all, wk == "8")
eco_wk7.nometa <- subset(eco_wk7, select = -c(sample_id, ph, temperature, food, container, wk))

eco.all.tubes.nometa <- subset(eco_all, select = -c(sample_id, ph, temperature, food, container, wk))

eco.all.tubes.melt <- melt(eco_all, id.vars=c("sample_id", "ph", "temperature", "food", "container", "wk"))
eco.all.tubes.melt$value <- ifelse(eco.all.tubes.melt$value <= 0, 0.001, eco.all.tubes.melt$value)


#WEEK 0 AND 7 COMBINED
bc.nmds.eco.tubesall <- metaMDS(eco.all.tubes.nometa, k=2, trymax=100)### Bray-Curtis is the default metric, k = 2 dimensions
ordiplot(bc.nmds.eco.tubesall, type = "t", display = "sites", cex = 0.7)

data.scoresecoall <- as.data.frame(scores(bc.nmds.eco.tubesall$points[,1:2]))
ecomerge <- cbind(data.scoresecoall, eco_all)
ecomerge$ph <- as.factor(ecomerge$ph)
ecomerge$food <- as.factor(ecomerge$food)
ecomerge$temperature <- as.factor(ecomerge$temperature)

plot(bc.nmds.eco.tubesall$points[,1:2], type="n", xlab="NMDS Axis 1", ylab="NMDS Axis 2", xlim=range(bc.nmds.eco.tubesall$points[,1]), ylim=range(bc.nmds.eco.tubesall$points[,2]), main="EcoTubesAll")
ordiarrows(bc.nmds.eco.tubesall, groups = ecomerge$sample_id, startmark = 1, lwd = 1, col="gray")
points(bc.nmds.eco.tubesall$points[,1:2], col=c("#CA64A3", "#9A4EAE", "#301934")[ecomerge$ph], pch=19, cex=1)
legend("bottomright", legend=c("pH 3", "pH 4", "pH 5.6"), col=c("#CA64A3", "#9A4EAE", "#301934"), pch=19, cex=1, bty="n")

plot(bc.nmds.eco.tubesall$points[,1:2], type="n", xlab="NMDS Axis 1", ylab="NMDS Axis 2", xlim=range(bc.nmds.eco.tubesall$points[,1]), ylim=range(bc.nmds.eco.tubesall$points[,2]), main="EcoTubesAll")
ordiarrows(bc.nmds.eco.tubesall, groups = ecomerge$sample_id, startmark = 1, lwd = 1, col="gray")
points(bc.nmds.eco.tubesall$points[,1:2], col=c("#90ee90", "#228822")[ecomerge$food], pch=19, cex=1)
legend("bottomright", legend=c("3 g/L", "6 g/L"), col=c("#90ee90", "#228822"), pch=19, cex=1, bty="n")

plot(bc.nmds.eco.tubesall$points[,1:2], type="n", xlab="NMDS Axis 1", ylab="NMDS Axis 2", xlim=range(bc.nmds.eco.tubesall$points[,1]), ylim=range(bc.nmds.eco.tubesall$points[,2]), main="EcoTubesAll")
ordiarrows(bc.nmds.eco.tubesall, groups = ecomerge$sample_id, startmark = 1, lwd = 1, col="gray")
points(bc.nmds.eco.tubesall$points[,1:2], col=c("#1f5776", "#c83126")[ecomerge$temperature], pch=19, cex=1)
legend("bottomright", legend=c("22C", "37C"), col=c("#1f5776", "#c83126"), pch=19, cex=1, bty="n")



ecomerge$wk <- as.factor(ecomerge$wk)
day1 <- subset(ecomerge, wk == 1)
day57 <- subset(ecomerge, wk == 8)

# Merge on sample_id to get start and end points for arrows
arrows <- merge(day1, day57, by = "sample_id", suffixes = c(".day1", ".day57"))

ggplot() + geom_segment(data = arrows, aes(x = MDS1.day1, y = MDS2.day1, 
  xend = MDS1.day57, yend = MDS2.day57), arrow = arrow(length = unit(0.2, "cm")), size = .4, color="gray") +
  geom_point(data = ecomerge, aes(x = MDS1, y = MDS2, color = ph, shape=food, alpha = wk), size = 4)+
  scale_color_manual(values = c("#CA64A3", "#9A4EAE", "#301934")) +
  labs(x = "NMDS Axis 1", y = "NMDS Axis 2", title = "EcoTubesAll") +
  theme_classic() +theme(legend.position = "bottomright")+
  scale_alpha_discrete(range = c(1,.65))

ggplot() + geom_segment(data = arrows, aes(x = MDS1.day1, y = MDS2.day1, 
                       xend = MDS1.day57, yend = MDS2.day57), arrow = arrow(length = unit(0.2, "cm")),
                        size = .4, color="gray") +
  geom_point(data = ecomerge, aes(x = MDS1, y = MDS2, color = temperature,shape=food, alpha = wk), size = 4)+
  scale_color_manual(values = c("#1f5776", "#c83126")) +
  labs(x = "NMDS Axis 1", y = "NMDS Axis 2", title = "EcoTubesAll") +
  theme_classic() +theme(legend.position = "bottomright")+
  scale_alpha_discrete(range = c(1,.6))

ggplot() + geom_segment(data = arrows, aes(x = MDS1.day1, y = MDS2.day1, 
                         xend = MDS1.day57, yend = MDS2.day57), arrow = arrow(length = unit(0.2, "cm")),
                        size = .4, color="gray") +
  geom_point(data = ecomerge, aes(x = MDS1, y = MDS2, color = food, shape=food, alpha = wk), size = 4)+
  scale_color_manual(values = c("#90ee90", "#228822")) +
  labs(x = "NMDS Axis 1", y = "NMDS Axis 2", title = "EcoTubesAll") +
  theme_classic() +theme(legend.position = "bottomright")+
  scale_alpha_discrete(range = c(1,.6))



#WEEK 0 ECOPLATE
bc.nmds.eco.tubeswk0 <- metaMDS(eco_wk0.nometa, k=2, trymax=100)### Bray-Curtis is the default metric, k = 2 dimensions
ordiplot(bc.nmds.eco.tubeswk0, type = "t", display = "sites", cex = 0.7)
data.scoresecotubwk0 <- as.data.frame(scores(bc.nmds.eco.tubeswk0$points[,1:2]))
ecomergewk0 <- cbind(data.scoresecotubwk0, eco_wk0)
ecomergewk0$ph <- as.factor(ecomergewk0$ph)
plot(bc.nmds.eco.tubeswk0$points[,1:2], type="n", xlab="NMDS Axis 1", ylab="NMDS Axis 2", xlim=range(bc.nmds.eco.tubeswk0$points[,1]), ylim=range(bc.nmds.eco.tubeswk0$points[,2]), main="EcoTubeswk0")
points(bc.nmds.eco.tubeswk0$points[,1:2], col=c("#CA64A3", "#9A4EAE", "#301934")[ecomergewk0$ph], pch=19, cex=1)
legend("bottomright", legend=c("pH 3", "pH 4", "pH 5.6"), col=c("#CA64A3", "#9A4EAE", "#301934"), pch=19, cex=1, bty="n")

plot(bc.nmds.eco.tubeswk0$points[,1:2], type="n", xlab="NMDS Axis 1", ylab="NMDS Axis 2", xlim=range(bc.nmds.eco.tubeswk0$points[,1]), ylim=range(bc.nmds.eco.tubeswk0$points[,2]), main="EcoTubeswk0")
points(bc.nmds.eco.tubeswk0$points[,1:2], col=c("#90ee90", "#228822")[ecomergewk0$food], pch=19, cex=1)
legend("bottomright", legend=c("3 g/L", "6 g/L"), col=c("#90ee90", "#228822"), pch=19, cex=1, bty="n")

plot(bc.nmds.eco.tubeswk0$points[,1:2], type="n", xlab="NMDS Axis 1", ylab="NMDS Axis 2", xlim=range(bc.nmds.eco.tubeswk0$points[,1]), ylim=range(bc.nmds.eco.tubeswk0$points[,2]), main="EcoTubeswk0")
points(bc.nmds.eco.tubeswk0$points[,1:2], col=c("#1f5776", "#c83126")[ecomergewk0$temperature], pch=19, cex=1)
legend("bottomright", legend=c("22C", "37C"), col=c("#1f5776", "#c83126"), pch=19, cex=1, bty="n")

#WEEK 7 ECOPLATE
bc.nmds.eco.tubeswk7 <- metaMDS(eco_wk7.nometa, k=2, trymax=100)### Bray-Curtis is the default metric, k = 2 dimensions
ordiplot(bc.nmds.eco.tubeswk7, type = "t", display = "sites", cex = 0.7)

data.scoresecotubwk7 <- as.data.frame(scores(bc.nmds.eco.tubeswk7$points[,1:2]))
ecomergewk7 <- cbind(data.scoresecotubwk0, eco_wk7)
ecomergewk7$ph <- as.factor(ecomergewk7$ph)
ecomergewk7$food <- as.factor(ecomergewk7$food)
ecomergewk7$temperature <- as.factor(ecomergewk7$temperature)

plot(bc.nmds.eco.tubeswk7$points[,1:2], type="n", xlab="NMDS Axis 1", ylab="NMDS Axis 2", xlim=range(bc.nmds.eco.tubeswk7$points[,1]), ylim=range(bc.nmds.eco.tubeswk7$points[,2]), main="EcoTubeswk7")
points(bc.nmds.eco.tubeswk7$points[,1:2], col=c("#CA64A3", "#9A4EAE", "#301934")[ecomergewk7$ph], pch=19, cex=1)
legend("bottomright", legend=c("pH 3", "pH 4", "pH 5.6"), col=c("#CA64A3", "#9A4EAE", "#301934"), pch=19, cex=1, bty="n")

plot(bc.nmds.eco.tubeswk7$points[,1:2], type="n", xlab="NMDS Axis 1", ylab="NMDS Axis 2", xlim=range(bc.nmds.eco.tubeswk7$points[,1]), ylim=range(bc.nmds.eco.tubeswk7$points[,2]), main="EcoTubeswk7")
points(bc.nmds.eco.tubeswk7$points[,1:2], col=c("#90ee90", "#228822")[ecomergewk7$food], pch=19, cex=1)
legend("bottomright", legend=c("3 g/L", "6 g/L"), col=c("#90ee90", "#228822"), pch=19, cex=1, bty="n")

plot(bc.nmds.eco.tubeswk7$points[,1:2], type="n", xlab="NMDS Axis 1", ylab="NMDS Axis 2", xlim=range(bc.nmds.eco.tubeswk7$points[,1]), ylim=range(bc.nmds.eco.tubeswk7$points[,2]), main="EcoTubeswk7")
points(bc.nmds.eco.tubeswk7$points[,1:2], col=c("#1f5776", "#c83126")[ecomergewk7$temperature], pch=19, cex=1)
legend("bottomright", legend=c("22C", "37C"), col=c("#1f5776", "#c83126"), pch=19, cex=1, bty="n")


#### Statistical Analysis ####
### PERMANOVA for categorical variables (factors)
eco.tubes.meta <- eco_all

#remove metadata
eco.tubes.nometa <- subset(eco.tubes.meta, select = -c(sample_id, ph, temperature, food, container, wk))
tube_dist <- vegdist(eco.tubes.nometa, method="bray")

## Calculate multivariate dispersions
disp_tubes_wk <- betadisper(tube_dist, eco.tubes.meta$wk)
disp_tubes_wk
boxplot(disp_tubes_wk)
anova(disp_tubes_wk)#0.01482 *

disp_tubes_ph <- betadisper(tube_dist, eco.tubes.meta$ph)
disp_tubes_ph
boxplot(disp_tubes_ph)
anova(disp_tubes_ph)#1.191e-11 ***

disp_tubes_food <- betadisper(tube_dist, eco.tubes.meta$food)
disp_tubes_food
boxplot(disp_tubes_food)
anova(disp_tubes_food)#0.264

disp_tubes_temp <- betadisper(tube_dist, eco.tubes.meta$temperature)
disp_tubes_temp
boxplot(disp_tubes_temp)
anova(disp_tubes_temp)#0.02426 *

#EcoPlate PERMANOVA
ad_tubes2 <- adonis2(tube_dist ~ temperature*ph*food + wk, data=eco.tubes.meta, by="terms")
ad_tubes2 
write_csv(ad_tubes2, "ecoplate_tubes_adonis2.csv")

eco.tubes.pw.all <- pairwise.adonis(eco.all.tubes.nometa, eco.all.tubes$ph , p.adjust.m='holm')
eco.tubes.pw.all
write_csv(eco.tubes.pw.all, "ecoplate_tubes_pw_ph.csv")

#pairs Df  SumsOfSqs   F.Model        R2 p.value p.adjusted sig
#1 4 vs 5.6  1 0.04536638  2.205683 0.0343534   0.084      0.084    
#2   4 vs 3  1 0.51471228 10.710689 0.1473056   0.001      0.003   *
#  3 5.6 vs 3  1 0.55225694 12.229011 0.1647471   0.001      0.003   *

eco.tubes.meta$sample_id <- as.numeric(eco.tubes.meta$sample_id)
eco.tubes.meta$ph <- factor(eco.tubes.meta$ph, levels = c(5.6, 4, 3))
eco.tubes.meta$food <- factor(eco.tubes.meta$food, levels = c(3, 6))
eco.tubes.meta$temperature <- factor(eco.tubes.meta$temperature, levels = c(22,37))
eco.tubes.meta$wk <- factor(eco.tubes.meta$wk, levels = c(1,8))

dbrda_eco <- dbrda(eco.tubes.nometa ~ temperature*ph*food+ wk, Condition=sample_id,
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

write.csv(eco_dbrda_anova, "ecoplate_tubes_dbrda_cca_aov.csv", row.names=TRUE)
significant_factors <- c("temperature37","ph3", "ph4","wk8", "temperature37:ph3")

eco_est <- as.data.frame(cbind(x1 = dbrda_eco$CCA$biplot[,1], y1 = dbrda_eco$CCA$biplot[,2]))
eco_est$factor <- row.names(dbrda_eco$CCA$biplot)
eco_est$significant <- eco_est$factor %in% significant_factors
significant_eco_est <- eco_est[eco_est$significant, ]

sites_eco <- as.data.frame(eco_sum_dbrda$sites) 
sites_eco_meta <- cbind(sites_eco,eco.tubes.meta )


#install.packages("RVAideMemoire")
#library(RVAideMemoire)
pairwise_results <-pairwise.perm.manova(sites_eco_meta[,1:2],eco.tubes.meta$ph,nperm=999)
pairwise_results
#  5.6      4     
#4 0.4580.  -     
#3 0.0015   0.0015


sites_eco_meta$wk <- as.factor(sites_eco_meta$wk)
sites_eco_meta$temperature <- as.factor(sites_eco_meta$temperature)
sites_eco_meta$ph <- factor(sites_eco_meta$ph, levels=c(3, 4, 5.6))

ggplot(data = sites_eco_meta, aes(x = dbRDA1, y = dbRDA2, color = temperature, alpha=wk, shape=wk)) +
  geom_point(size = 3) +scale_shape_manual(values = c(19, 15)) +
  scale_color_manual(values = c("#1f5776", "#c83126")) +
  geom_segment(data = significant_eco_est, aes(x = 0, y = 0, xend = x1, yend = y1),
               arrow = arrow(length = unit(.2, "cm")), color = "black", size = 1, inherit.aes = FALSE) +
  geom_text(data = significant_eco_est, aes(x = x1, y = y1, label = factor), 
            hjust = 0, vjust = 0, inherit.aes = FALSE)+ theme_classic()+scale_alpha_discrete(range = c(1,.7))+
  theme(legend.title=element_blank())


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
meta_pro_filt <- meta_pro_filt[,c(2, 16, 17)]
eco_pro_all <- left_join(eco_all, meta_pro_filt, by=c("sample_id", "wk"))

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
pro_final_merge_meta <- left_join(pro_final_merge, meta_pro, by="sample_id")

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

meta_pro <- subset(meta)
meta_pro <- subset(meta_pro, day %in% c("57"))

meta_pro_filt <- meta_pro
meta_pro_filt$id <- rownames(meta_pro_filt)
meta_pro_filt$wk <- NA
meta_pro_filt$wk[meta_pro_filt$day == 1] <- 1  
meta_pro_filt$wk[meta_pro_filt$day == 57] <- 8 


meta_pro_filt <- meta_pro_filt[,c(2, 16, 17)]

eco_pro_all <- left_join(eco_all, meta_pro_filt, by=c("sample_id", "wk"))
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
pro_final_merge_meta <- left_join(pro_final_merge, meta_pro, by="sample_id")

pro_final_merge_meta$temperature <- as.factor(pro_final_merge_meta$temperature)
pro_final_merge_meta$ph <- as.factor(pro_final_merge_meta$ph)
pro_final_merge_meta$food <- as.factor(pro_final_merge_meta$food)
pro_final_merge_meta$treatment_combo <- as.factor(pro_final_merge_meta$treatment_combo)

ggplot(pro_final_merge_meta, aes(group = sample_id, shape=food)) +
  geom_point(aes(x = V1, y = V2, size=2, color=temperature)) + 
  geom_point(aes(x = MDS1, y = MDS2, size=2, color=temperature))+
  geom_segment(aes(x = V1, y = V2, xend = MDS1, yend = MDS2,group = sample_id, color="gray"),
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
mantel_meta <- meta
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

eco_mantel <- left_join(mantel_meta, eco_all, by=c("sample_id", "wk"))
rownames(eco_mantel) <- eco_mantel$id
eco_mantel <- eco_mantel[,18:48]

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
pro_meta_ht <- left_join(pro_ht, meta_ht, by="sample_id")
pro_meta_ht$temperature <- as.factor(pro_meta_ht$temperature)
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


# red "#c83126"
# blue "#1f5776"
# pink "#CA64A3"
#light purple "#9A4EAE"
#dark purple "#301934"
#light green "#90ee90"
#dark green "#228822"



#### CHITINASE ####
#filter metadata
meta<- read.csv("Exp_2_metadata_tubes.csv", header=TRUE)
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

ggplot(df.summarychit, aes(day, chitinase, color = temperature, shape = food)) +
  geom_jitter(size = 3, alpha = 0.9, position = position_dodge(0.3)) +
  facet_wrap(~ph, nrow = 1) +
  geom_line(aes(group = treatment_combo), position = position_dodge(0.3)) +
  geom_errorbar(aes(ymin = chitinase - sd, ymax = chitinase + sd), position = position_dodge(0.3), width = 0.2) +
  theme_bw() +scale_shape_manual(values=c(16,1)) +
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

ggplot(df.summaryprot, aes(day, protease, color = temperature, shape = food)) +
  geom_jitter(size = 3, alpha = 0.9, position = position_dodge(0.3)) +
  facet_wrap(~ph, nrow = 1) +
  geom_line(aes(group = treatment_combo), position = position_dodge(0.3)) +
  geom_errorbar(aes(ymin = protease - sd, ymax = protease + sd), position = position_dodge(0.3), width = 0.2) +
  theme_bw() +scale_shape_manual(values=c(16,1)) +
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


ggplot(df.summaryresp, aes(day, respiration_co2_ppm_hr, color = temperature, shape = food)) +
  geom_jitter(size = 3, alpha = 0.9, position = position_dodge(0.3)) +
  facet_wrap(~ph, nrow = 1) +
  geom_line(aes(group = treatment_combo), position = position_dodge(0.3)) +
  geom_errorbar(aes(ymin = respiration_co2_ppm_hr - sd, ymax = respiration_co2_ppm_hr + sd), position = position_dodge(0.3), width = 0.2) +
  theme_bw() +scale_shape_manual(values=c(16,1)) +
  scale_color_manual(values = c("#1f5776", "#c83126")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none") 

#### GLMMS ####
#response = chitinase, protease, respiration
#predictor = ph, food, temp
#repeated measure = day

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

mchit <- brm(chitinase ~ ph*food*temperature+day+ (1|sample_id),data=meta_all_tubes,
             family = Gamma(link = "log"), iter = 10000, chains = 4, cores = 4,control = list(adapt_delta = 0.999, stepsize = 0.001, max_treedepth=20),
             prior=informpriors_chit)
pp_check(mchit)

r2(mchit)
# Bayesian R2 with Compatibility Interval
#Conditional R2: 0.669 (95% CI [0.635, 0.695])
#Marginal R2: 0.390 (95% CI [0.187, 0.605])

meta_all_tubes$protease <- ifelse(meta_all_tubes$protease <= 0, 0.001, meta_all_tubes$protease)
mprot <- brm(protease ~ ph*food*temperature + (1|day),data=meta_all_tubes,
             family = Gamma(link = "log"), iter = 10000, chains = 4, cores = 4,control = list(adapt_delta = 0.999, stepsize = 0.001, max_treedepth=20),
             prior=informpriors_prot)


mresp <- brm(respiration_co2_ppm_hr ~ ph*food*temperature+ (1|day),data=meta_all_tubes,
             family = Gamma(link = "log"), iter = 10000, chains = 4, cores = 4,control = list(adapt_delta = 0.999, stepsize = 0.001, max_treedepth=20),
             prior=informpriors_resp)


summary(mchit)
summary(mprot)
summary(mresp)

saveRDS(mchit, file = "brms_exp2_chit.RDS")
saveRDS(mprot, file = "brms_exp2_prot.RDS")
saveRDS(mresp, file = "brms_exp2_resp.RDS")

mchit <- readRDS("brms_exp2_chit.RDS")
mprot <- readRDS("brms_exp2_prot.RDS")
mresp <- readRDS("brms_exp2_resp.RDS")

pp_check(mchit)
pp_check(mprot)
pp_check(mresp)

##INTERACTION PLOTS
chit.pred <- predict_response(mchit, terms = c("food", "temperature", "ph"))

ggplot(chit.pred, aes(x = x, y = predicted, shape = x, color = group)) +
  geom_point(size = 5, position = position_dodge(width = .7)) + # Dodging the points
  facet_wrap(~facet) +
  scale_color_manual(values = c("#1f5776", "#c83126")) +
  theme_bw() +scale_shape_manual(values=c(16,1)) +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0, size=2,position = position_dodge(width = .7))+
  theme(legend.position="none")+ ylab("Predicted Chitinase Activity")+
  xlab("Food (g/L)")+ theme(
    text = element_text(color = "black", size = 10 )
  )




prot.pred <- ggpredict(mprot, terms = c("food", "temperature", "ph"))
prot.pred$facet <- factor(prot.pred$facet, levels = c("3", "4", "5.6"))

ggplot(prot.pred, aes(x = x, y = predicted, shape = x, color = group)) +
  geom_point(size = 5, position = position_dodge(width = .7)) + # Dodging the points
  facet_wrap(~facet) +
  scale_color_manual(values = c("#1f5776", "#c83126")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0, size=2,position = position_dodge(width = .7))+
  theme(legend.position="none")+ ylab("Predicted Protease Activity")+
  xlab("Food (g/L)")+ theme(
    text = element_text(color = "black", size = 10 )
  )


resp.pred <- ggpredict(mresp, terms = c("food", "temperature", "ph"))
resp.pred$facet <- factor(resp.pred$facet, levels = c("3", "4", "5.6"))

ggplot(resp.pred, aes(x = x, y = predicted, shape = x, color = group)) +
  geom_point(size = 5, position = position_dodge(width = .7)) + # Dodging the points
  facet_wrap(~facet) +
  scale_color_manual(values = c("#1f5776", "#c83126")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0, size=2,position = position_dodge(width = .7))+
  theme(legend.position="none")+ ylab("Predicted Respiration")+
  xlab("Food (g/L)")+ theme(
    text = element_text(color = "black", size = 10 ))+ylim(c(0,35000))

#Directional Effect of food on respiration
posterior_mresp <- as.data.frame(mresp)
highfood <- posterior_mresp %>% filter(b_food6 >0)
nrow(highfood)/nrow(posterior_mresp)*100

##Cleaned up posterior distribution graphs

#CHITINASE MODEL
posteriorchit <- mcmc_intervals_data(mchit, 
                                  prob_outer=0.95,
                                  prob=0.5)

posteriorchit$nonzero <- NA
posteriorchit$nonzero[posteriorchit$ll>0 & posteriorchit$hh>0] <- "nonzero"
posteriorchit$nonzero[posteriorchit$ll<0 & posteriorchit$hh<0] <- "nonzero"
posteriorchit$nonzero[is.na(posteriorchit$nonzero)] <- "zero"
posteriorchit<- posteriorchit[1:13,]

ggplot(posteriorchit, aes(x = parameter,
                       shape=nonzero)) +
  geom_hline(yintercept = 0, linetype = 3, 
             size=1, color = "#b0b5b3") +
  geom_pointrange(aes(ymin = ll, ymax = hh, y = m),
                  position= position_dodge(width=0.75),
                  size = 3/4) +
  scale_shape_manual(values=c(17, 19), 
                     labels=c("95% CI does\nnot contain zero", 
                              "95% CI\ncontains zero"))+
  coord_flip() +
  xlab(NULL) +
  ylab("Estimated effect on Chitinase Activity")+
  theme(
    text = element_text(color = "black", size = 16 ))+
theme(panel.background = element_rect(fill = "white", colour = "black"))+
  theme(legend.position="none")

#PROTEASE MODEL
posteriorprot <- mcmc_intervals_data(mprot, 
                                     prob_outer=0.95,
                                     prob=0.5)

posteriorprot$nonzero <- NA
posteriorprot$nonzero[posteriorprot$ll>0 & posteriorprot$hh>0] <- "nonzero"
posteriorprot$nonzero[posteriorprot$ll<0 & posteriorprot$hh<0] <- "nonzero"
posteriorprot$nonzero[is.na(posteriorprot$nonzero)] <- "zero"
posteriorprot<- posteriorprot[1:12,]

ggplot(posteriorprot, aes(x = parameter,
                          shape=nonzero)) +
  geom_hline(yintercept = 0, linetype = 3, 
             size=1, color = "#b0b5b3") +
  geom_pointrange(aes(ymin = ll, ymax = hh, y = m),
                  position= position_dodge(width=0.75),
                  size = 3/4) +
  scale_shape_manual(values=c(17, 19), 
                     labels=c("95% CI does\nnot contain zero", 
                              "95% CI\ncontains zero"))+
  coord_flip() +
  xlab(NULL) +
  ylab("Estimated effect on Protease Activity")+
  theme(text = element_text(color = "black", size = 16 ))+
  theme(panel.background = element_rect(fill = "white", colour = "black"))+
  theme(legend.position="none")

#RESPIRATION MODEL
posteriorresp <- mcmc_intervals_data(mresp, 
                                     prob_outer=0.95,
                                     prob=0.5)

posteriorresp$nonzero <- NA
posteriorresp$nonzero[posteriorresp$ll>0 & posteriorresp$hh>0] <- "nonzero"
posteriorresp$nonzero[posteriorresp$ll<0 & posteriorresp$hh<0] <- "nonzero"
posteriorresp$nonzero[is.na(posteriorresp$nonzero)] <- "zero"
posteriorresp<- posteriorresp[1:12,]

ggplot(posteriorresp, aes(x = parameter, shape=nonzero)) +
  geom_hline(yintercept = 0, linetype = 3, size=1, color = "#b0b5b3") +
  geom_pointrange(aes(ymin = ll, ymax = hh, y = m), position= position_dodge(width=0.75), size = 3/4) +
  scale_shape_manual(values=c(17, 19), labels=c("95% CI does\nnot contain zero", "95% CI\ncontains zero"))+
  coord_flip() +xlab(NULL) + ylab("Estimated effect on Respiration")+ theme(text = element_text(color = "black", size = 16 ))+
  theme(panel.background = element_rect(fill = "white", colour = "black"))+theme(legend.position="none")


#EFFECTIVE NUMBER OF SPECIES MODEL
posteriormens <- mcmc_intervals_data(mens, 
                                     prob_outer=0.95,
                                     prob=0.5)

posteriormens$nonzero <- NA
posteriormens$nonzero[posteriormens$ll>0 & posteriormens$hh>0] <- "nonzero"
posteriormens$nonzero[posteriormens$ll<0 & posteriormens$hh<0] <- "nonzero"
posteriormens$nonzero[is.na(posteriormens$nonzero)] <- "zero"
posteriormens<- posteriormens[1:12,]

ggplot(posteriormens, aes(x = parameter,
                          shape=nonzero)) +
  geom_hline(yintercept = 0, linetype = 3, 
             size=1, color = "#b0b5b3") +
  geom_pointrange(aes(ymin = ll, ymax = hh, y = m),
                  position= position_dodge(width=0.75),
                  size = 3/4) +
  scale_shape_manual(values=c(17, 19), 
                     labels=c("95% CI does\nnot contain zero", 
                              "95% CI\ncontains zero"))+
  coord_flip() +
  theme(axis.text.y = element_text( size=7), 
        axis.text.x=element_text(size=7),
        axis.title = element_text(size=7), 
        legend.text = element_text(size=7)) +
  xlab(NULL) +
  ylab("Estimated effect on Effective Number of ASVs")+
  theme(panel.background = element_rect(fill = "white", colour = "black"))
