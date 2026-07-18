if (!requireNamespace("WGCNA", quietly = TRUE)) install.packages("WGCNA")
library(WGCNA)
options(stringsAsFactors = FALSE)
allowWGCNAThreads(nThreads=2)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 4) {
  stop("Usage: Rscript wgcna.R <gene_matrix path> <module_assignment out path> <me_matrix out path> <module_membership out path>")
}

gene_matrix <- args[1]
assignment_path <- args[2]
me_matrix_path <- args[3]
membership_path <- args[4]

genes <- read.csv(gene_matrix, check.names = FALSE)

sample_ids <- genes$sample_id
gene_data <- genes[, !(colnames(genes) %in% c("sample_id"))]

gsg <- goodSamplesGenes(gene_data, verbose = 3)

if (!gsg$allOK) {
  gene_data <- gene_data[gsg$goodSamples, gsg$goodGenes]
  sample_ids <- sample_ids[gsg$goodSamples]
}

powers <- c(1:10, seq(12, 20, 2))
sft <- pickSoftThreshold(gene_data, powerVector = powers, verbose = 5)

soft_power <- sft$powerEstimate
if (is.na(soft_power)) {
  soft_power <- powers[which.max(sft$fitIndices[, 2])]
}

net <- blockwiseModules(
  gene_data,
  power = soft_power,
  TOMType = "unsigned",
  minModuleSize = 30,
  reassignThreshold = 0,
  mergeCutHeight = 0.25,
  numericLabels = TRUE,
  pamRespectsDendro = FALSE,
  saveTOMs = FALSE,
  verbose = 3
)

module_colours <- labels2colors(net$colors)

module_assignment <- data.frame(
  Gene = colnames(gene_data),
  ModuleColor = module_colours,
  ModuleLabel = net$colors
)

write.csv(module_assignment, assignment_path, row.names = FALSE)

me <- net$MEs

colnames(me) <- paste0("ME", substring(colnames(me), 3))

me_df <- data.frame(sample_id = sample_ids, me, check.names = FALSE)

write.csv(me_df, me_matrix_path, row.names = FALSE)

module_membership <- signedKME(gene_data, me)

write.csv(module_membership, membership_path, row.names = TRUE)