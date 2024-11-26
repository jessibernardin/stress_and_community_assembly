library(phylofactor)

#### practice example ####
data("FTmicrobiome")
OTUTable <- FTmicrobiome$OTUTable   #Our OTU table. Rows are OTUs and columns are samples
OTUTable[1:5,1:5]
body.site <- FTmicrobiome$X         #Our independent variable
str(body.site)
tree <- FTmicrobiome$tree           #Our phylogeny
taxonomy <- FTmicrobiome$taxonomy   #Our taxonomy
taxonomy[1:3,]
PF <- PhyloFactor(OTUTable,tree,body.site,nfactors=3)


gtree <- pf.tree(PF)
ggtree::rotate_tree(gtree$ggplot,-45)

library(phytools)
clr <- function(Matrix) apply(Matrix,MARGIN=2,FUN=function(x) log(x)-mean(log(x)))
par(mfrow=c(1,1))
phylo.heatmap(tree,clr(PF$Data))

#### Exp 2 Data ####
phylo_exp2 <- readRDS("~/Dropbox/2022_Comm_Exp2/stress_and_community_assembly/Exp2.physeq3.RDS")

Exp2.tubes.physeq3 = subset_samples(phylo_exp2, container == "tube") 
row_sums <- rowSums(otu_table(Exp2.tubes.physeq3))
nonzero_rows <- row_sums != 0
Exp2.tubes.physeq3 <- prune_taxa(nonzero_rows, Exp2.tubes.physeq3)#565 taxa and 144 samples
