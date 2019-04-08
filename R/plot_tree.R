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



ml.tree <- ggtree::read.raxml("data/tree/RAxML_bipartitionsBranchLabels.ltde")
df.species <- read.table("data/demography/weibull_results_clean_species.csv", 
                         header = TRUE, sep = ",", row.names = 1, stringsAsFactors = FALSE)
to.keep <- df.species$Species
#to.keep <- append(to.keep, "NC_005042.1.353331-354795")
to.keep <- to.keep[to.keep!= 'KBS0727']
pruned.tree <- treeio::drop.tip(ml.tree, ml.tree@phylo$tip.label[-match(to.keep, ml.tree@phylo$tip.label)])
## use  ggtree::drop.tip() to remove tips from ml.tree


tree.plot <- ggtree(pruned.tree, branch.length='none') +
  #geom_text(aes(label=bootstrap)) +
  geom_tiplab(size = 3) + 
  geom_label2(aes(label = bootstrap), size=2.5, label.size = 0.3) +
  geom_treescale(x = 0, y = -0.5) + xlim(0, 10)


ggsave(file="figs/tree.png", tree.plot, units='in', dpi=600)


