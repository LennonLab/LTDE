rm(list = ls())
getwd()
setwd("~/GitHub/LTDE/")

library('phytools')
library('ape')

tree <- read.tree("data/tree/ribosomal_protein_fasttree.tre")
tree.root <- midpoint.root(tree)
#plot(tree.root)
identify.phylo(tree.root)
# print the outgroup
tree.root$tip.label[length(tree.root$tip.label)]
write.tree(tree.root, file = "data/tree/ribosomal_protein_fasttree_midpoint.tre")
