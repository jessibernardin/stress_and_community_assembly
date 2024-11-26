packages_to_load <- c(
  "phyloseq", "dplyr",  "picante", "ape"
)

# Load and install required packages
for (i in packages_to_load) { #Installs packages if not yet installed
  if (!require(i, character.only = TRUE)) install.packages(i)
}

Exp2.physeq3 <- readRDS("Exp2.physeq3.RDS")

#make two seperate phyloseq objects, one for tubes and one for plants
Exp2.pitcher.physeq3 = subset_samples(Exp2.physeq3, container == "pitcher") #852 taxa and 62 samples

Exp2.tubes.physeq3 = subset_samples(Exp2.physeq3, container == "tube") #852 taxa and 144 samples
row_sums <- rowSums(otu_table(Exp2.tubes.physeq3))
nonzero_rows <- row_sums != 0
Exp2.tubes.physeq3 <- prune_taxa(nonzero_rows, Exp2.tubes.physeq3)#565 taxa and 144 samples

Exp2.tubesmix.physeq3 = subset_samples(Exp2.physeq3, container !="pitcher") #852 taxa and 144 samples
row_sums <- rowSums(otu_table(Exp2.tubesmix.physeq3))
nonzero_rows <- row_sums != 0
Exp2.tubesmix.physeq3 <- prune_taxa(nonzero_rows, Exp2.tubesmix.physeq3)#570 taxa and 148 samples



ps = Exp2.tubesmix.physeq3 %>% transform_sample_counts(fun = function(x) x/sum(x))

# make sure ps is in relative abundance!!
sample_sums(ps)

otu = otu_table(ps)
phylo = phy_tree(ps)

## make sure the names on the phylogeny are ordered the same as the names in otu table

match.phylo.otu = match.phylo.data(phylo, otu);
str(match.phylo.otu);

## calculate empirical betaMNTD

beta.mntd.weighted = as.matrix(comdistnt(t(match.phylo.otu$data),cophenetic(match.phylo.otu$phy),abundance.weighted=T));
dim(beta.mntd.weighted);
beta.mntd.weighted[1:5,1:5]

identical(colnames(match.phylo.otu$data),colnames(beta.mntd.weighted)); # just a check, should be TRUE
identical(colnames(match.phylo.otu$data),rownames(beta.mntd.weighted)); # just a check, should be TRUE

# calculate randomized betaMNTD

beta.reps = 999; # number of randomizations

rand.weighted.bMNTD.comp = array(c(-999),dim=c(ncol(match.phylo.otu$data),ncol(match.phylo.otu$data),beta.reps));
dim(rand.weighted.bMNTD.comp);

for (rep in 1:beta.reps) {
  
  rand.weighted.bMNTD.comp[,,rep] = as.matrix(comdistnt(t(match.phylo.otu$data),taxaShuffle(cophenetic(match.phylo.otu$phy)),abundance.weighted=T,exclude.conspecifics = F));
  
  print(c(date(),rep));
  
}
###
weighted.bNTI = matrix(c(NA),nrow=ncol(match.phylo.otu$data),ncol=ncol(match.phylo.otu$data));
dim(weighted.bNTI);

for (columns in 1:(ncol(match.phylo.otu$data)-1)) {
  for (rows in (columns+1):ncol(match.phylo.otu$data)) {
    
    rand.vals = rand.weighted.bMNTD.comp[rows,columns,];
    weighted.bNTI[rows,columns] = (beta.mntd.weighted[rows,columns] - mean(rand.vals)) / sd(rand.vals);
    rm("rand.vals");
    
  };
};

rownames(weighted.bNTI) = colnames(match.phylo.otu$data);
colnames(weighted.bNTI) = colnames(match.phylo.otu$data);
weighted.bNTI;
write.csv(weighted.bNTI,"weighted_bNTI.csv",quote=F);
saveRDS(weighted.bNTI, "beta-NTI-relative-abundance.RDS")













