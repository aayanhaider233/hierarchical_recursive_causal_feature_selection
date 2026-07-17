if (!requireNamespace("WGCNA", quietly = TRUE)) install.packages("WGCNA")
library(WGCNA)
options(stringsAsFactors = FALSE)
allowWGCNAThreads(nThreads=2)

# ------------------- Load Data -------------------
dataFile <- "C:/Users/USER/Documents/Thesis/Datasets/Final/train_gene_aggregate.csv"
data <- read.csv(dataFile, check.names = FALSE)

# ------------------- Separate genes and traits -------------------
n_traits <- 6
n_total <- ncol(data)

geneData <- data[, !(colnames(data) %in% c("age","sex","disease","hannum_EAA","horvath2013_EAA","grimage2_EAA"))]
traits <- data[, c("age","sex","disease","hannum_EAA","horvath2013_EAA","grimage2_EAA")]

# ------------------- Check genes and samples -------------------
gsg <- goodSamplesGenes(geneData, verbose = 3)
if (!gsg$allOK) {
  geneData <- geneData[gsg$goodSamples, gsg$goodGenes]
}

# ------------------- Choose soft-threshold power -------------------
powers <- c(1:10, seq(12, 20, 2))
sft <- pickSoftThreshold(geneData, powerVector = powers, verbose = 5)

# Pick the smallest power where scale-free topology fit R^2 >= 0.9
softPower <- sft$powerEstimate
if (is.na(softPower)) {
  softPower <- powers[which.max(sft$fitIndices[, 2])]
}
cat("Chosen soft-thresholding power:", softPower, "\n")

# ------------------- Construct network and detect modules -------------------
net <- blockwiseModules(
  geneData,
  power = softPower,
  TOMType = "unsigned",
  minModuleSize = 30,
  reassignThreshold = 0,
  mergeCutHeight = 0.25,
  numericLabels = TRUE,
  pamRespectsDendro = FALSE,
  saveTOMs = TRUE,
  verbose = 3
)

# ------------------- Module colors -------------------
moduleColors <- labels2colors(net$colors)
table(moduleColors)

# Save module assignments: gene -> module
geneModuleAssignment <- data.frame(
  Gene = colnames(geneData),
  ModuleColor = moduleColors,
  ModuleLabel = net$colors
)
write.csv(geneModuleAssignment,
          "C:/Users/USER/Documents/Thesis/Datasets/Final/gene_module_assignments.csv",
          row.names = FALSE)

# ------------------- Module eigengenes -------------------
MEs <- net$MEs
# Ensure column names match module colors
colnames(MEs) <- paste0("ME", substring(colnames(MEs), 3)) # e.g., "ME1" -> "MEblue"

write.csv(MEs, "C:/Users/USER/Documents/Thesis/Datasets/Final/module_eigengenes.csv", row.names = FALSE)

# ------------------- Module-trait correlations -------------------
moduleTraitCor <- cor(MEs, traits, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(geneData))

write.csv(moduleTraitCor,
          "C:/Users/USER/Documents/Thesis/Datasets/Final/module_trait_cor.csv",
          row.names = TRUE)
write.csv(moduleTraitPvalue,
          "C:/Users/USER/Documents/Thesis/Datasets/Final/module_trait_pvalues.csv",
          row.names = TRUE)

# ------------------- Gene-module membership (kME) -------------------
geneModuleMembership <- signedKME(geneData, MEs)
write.csv(geneModuleMembership,
          "C:/Users/USER/Documents/Thesis/Datasets/Final/gene_module_membership.csv",
          row.names = TRUE)

# ------------------- Extract genes per module -------------------
modules <- unique(moduleColors)
for (mod in modules) {
  genesInModule <- colnames(geneData)[moduleColors == mod]
  outFile <- paste0("C:/Users/USER/Documents/Thesis/Datasets/Final/genes_module_", mod, ".csv")
  write.csv(genesInModule, outFile, row.names = FALSE)
}

# ------------------- Done -------------------
print("WGCNA analysis complete! All outputs saved as CSVs.")
