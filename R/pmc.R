rm(list = ls())
getwd()
setwd("~/GitHub/LTDE/")

# unblock to install package
#library("devtools")
#install_github("cboettig/pmc")

library("pmc")
library("ape")

iter <- 1000

df <- read.table("data/demography/weibull_results_clean_species.csv", 
                 header = TRUE, sep = ",", row.names = 1, stringsAsFactors = FALSE)
df.no812 <- df[!(df$Species=="KBS0812" ),]
rownames(df) <- df$Species
rownames(df.no812) <- df.no812$Species

# Load ML tree
ml.tree <- read.tree("data/tree/RAxML_bipartitionsBranchLabels.ltde_no_contamination")
# Define the outgroup
outgroup <- match("NC_005042.1.353331-354795", ml.tree$tip.label)
# Create a rooted tree {ape} using the outgroup
ml.rooted <- root(ml.tree, outgroup, resolve.root = TRUE)
# Keep rooted but drop outgroup branch
ml.rooted <- drop.tip(ml.rooted, c("NC_005042.1.353331-354795"))
ml.rooted.no812 <- drop.tip(ml.rooted, c("KBS0812"))

is.ultrametric(ml.rooted.no812)
ml.rooted.no812.um  <- chronos(ml.rooted.no812)
is.ultrametric(ml.rooted.no812.um)

beta.log10.no812 <- df.no812$beta.log10
names(beta.log10.no812) <- df.no812$Species


beta.BM.lambda <- pmc(ml.rooted.no812, beta.log10.no812, "BM", "lambda", nboot = iter)
p_value.BM.lambda <- (length(beta.BM.lambda$null[beta.BM.lambda$null > beta.BM.lambda$lr] )+1) / (iter+1)

print(beta.BM.lambda$lr)
print(beta.BM.lambda$B$opt$lambda)
print(p_value.BM.lambda)



