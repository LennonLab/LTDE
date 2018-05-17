rm(list = ls())
getwd()
setwd("~/GitHub/LTDE/")

library('phytools')
library('ape')
library('adephylo')

tree <- read.tree("data/tree/ribosomal_protein_fasttree.tre")
tree.root <- midpoint.root(tree)
identify.phylo(tree.root)
# print the outgroup
tree.root$tip.label[length(tree.root$tip.label)]
write.tree(tree.root, file = "data/tree/ribosomal_protein_fasttree_midpoint.tre")

length(tree.root$edge.length)
length(tree.root$node.label)

tree.root.rtt <- as.data.frame(distRoot(tree.root, method = 'patristic'))
colnames(tree.root.rtt) <- "distance"
tree.root.rtt[order(tree.root.rtt$distance), c(1)]
tree.root.rtt[order(distance),] 

order.taxa <- rev(rownames(tree.root.rtt)[order(tree.root.rtt$distance)])
order.taxa[1:5]
  
# longest branches on each side of the split
# Bacteria_Proteobacteria_Gammaproteobacteria_Enterobacteriales_Enterobacteriaceae_Candidatus_Riesia_pediculicola_USDA
# Archaea_Thaumarchaeota_Nitrosopumilales_Nitrosopumilaceae_Nitrosopumilus_sp__MY1_MY1
