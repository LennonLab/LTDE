rm(list = ls())
getwd()
setwd("~/GitHub/LTDE/")

#library('ape')
library('seqinr')
#library('phytools')
#source("https://bioconductor.org/biocLite.R")
#biocLite("ggtree")
#install.packages("BiocManager")
#BiocManager::install(c("Biostrings", "ggtree", "treeio"))
library('ggtree')
library('treeio')
library('latex2exp')
library('gridExtra')

ml.tree <- ggtree::read.raxml("data/tree/RAxML_bipartitionsBranchLabels.ltde")
df.species <- read.table("data/demography/weibull_results_clean_species.csv", 
                         header = TRUE, sep = ",", row.names = 1, stringsAsFactors = FALSE)
to.keep <- df.species$Species
#to.keep <- append(to.keep, "NC_005042.1.353331-354795")
to.keep <- to.keep[to.keep!= 'KBS0727']
pruned.tree <- treeio::drop.tip(ml.tree, ml.tree@phylo$tip.label[-match(to.keep, ml.tree@phylo$tip.label)])
## use  ggtree::drop.tip() to remove tips from ml.tree
pruned.tree@phylo$tip.label <- replace(pruned.tree@phylo$tip.label, pruned.tree@phylo$tip.label=="KBS0721", "Flavobacterium sp. KBS0721")
pruned.tree@phylo$tip.label <- replace(pruned.tree@phylo$tip.label, pruned.tree@phylo$tip.label=="KBS0707", "Pseudomonas sp. KBS0707")
pruned.tree@phylo$tip.label <- replace(pruned.tree@phylo$tip.label, pruned.tree@phylo$tip.label=="KBS0702", "Arthrobacter sp. KBS0702")
pruned.tree@phylo$tip.label <- replace(pruned.tree@phylo$tip.label, pruned.tree@phylo$tip.label=="ATCC13985", "Pseudomonas sp. ATCC13985")
pruned.tree@phylo$tip.label <- replace(pruned.tree@phylo$tip.label, pruned.tree@phylo$tip.label=="ATCC43928", "Pseudomonas sp. ATCC43928")
pruned.tree@phylo$tip.label <- replace(pruned.tree@phylo$tip.label, pruned.tree@phylo$tip.label=="KBS0701", "Pedobacter sp. KBS0701")
pruned.tree@phylo$tip.label <- replace(pruned.tree@phylo$tip.label, pruned.tree@phylo$tip.label=="KBS0703", "Arthrobacter sp. KBS0703")
pruned.tree@phylo$tip.label <- replace(pruned.tree@phylo$tip.label, pruned.tree@phylo$tip.label=="KBS0705", "Inquilinus sp. KBS0705")
pruned.tree@phylo$tip.label <- replace(pruned.tree@phylo$tip.label, pruned.tree@phylo$tip.label=="KBS0706", "Mycobacterium sp. KBS0706")
pruned.tree@phylo$tip.label <- replace(pruned.tree@phylo$tip.label, pruned.tree@phylo$tip.label=="KBS0710", "Pseudomonas sp. KBS0710")
pruned.tree@phylo$tip.label <- replace(pruned.tree@phylo$tip.label, pruned.tree@phylo$tip.label=="KBS0711", "Janthinobacterium sp. KBS0711")
pruned.tree@phylo$tip.label <- replace(pruned.tree@phylo$tip.label, pruned.tree@phylo$tip.label=="KBS0712", "Variovorax sp. KBS0712")
pruned.tree@phylo$tip.label <- replace(pruned.tree@phylo$tip.label, pruned.tree@phylo$tip.label=="KBS0713", "Yersinia sp. KBS0713")
pruned.tree@phylo$tip.label <- replace(pruned.tree@phylo$tip.label, pruned.tree@phylo$tip.label=="KBS0714", "Micrococcus sp. KBS0714")
pruned.tree@phylo$tip.label <- replace(pruned.tree@phylo$tip.label, pruned.tree@phylo$tip.label=="KBS0715", "Curtobacterium sp. KBS0715")
pruned.tree@phylo$tip.label <- replace(pruned.tree@phylo$tip.label, pruned.tree@phylo$tip.label=="KBS0721", "Flavobacterium sp. KBS0721")
pruned.tree@phylo$tip.label <- replace(pruned.tree@phylo$tip.label, pruned.tree@phylo$tip.label=="KBS0722", "Oerskovia sp. KBS0722")
pruned.tree@phylo$tip.label <- replace(pruned.tree@phylo$tip.label, pruned.tree@phylo$tip.label=="KBS0724", "Rhodococcus sp. KBS0724")
pruned.tree@phylo$tip.label <- replace(pruned.tree@phylo$tip.label, pruned.tree@phylo$tip.label=="KBS0725", "Bradyrhizobium sp. KBS0725")
pruned.tree@phylo$tip.label <- replace(pruned.tree@phylo$tip.label, pruned.tree@phylo$tip.label=="KBS0801", "Burkholderia sp. KBS0801")
pruned.tree@phylo$tip.label <- replace(pruned.tree@phylo$tip.label, pruned.tree@phylo$tip.label=="KBS0802", "Pseudomonas sp. KBS0802")
pruned.tree@phylo$tip.label <- replace(pruned.tree@phylo$tip.label, pruned.tree@phylo$tip.label=="KBS0812", "Bacillus sp. KBS0812")

tree.plot <- ggtree(pruned.tree, branch.length='none') +
  #geom_text(aes(label=bootstrap)) +
  geom_tiplab(size = 2.5) + 
  geom_label2(aes(label = bootstrap), size=2.5, label.size = 0.3) +
  geom_treescale(x = 0.2, y = -1.3) + xlim(0, 11)


ggsave(file="figs/tree.png", tree.plot, units='in', dpi=600)


