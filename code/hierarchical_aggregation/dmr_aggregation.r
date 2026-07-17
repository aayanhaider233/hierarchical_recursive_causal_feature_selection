suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("Usage: Rscript dmr_aggregation_02.R <cpg_matrix.csv> <dmr_rds> <output_csv>")
}

cpg_matrix <- args[1]
dmr_rds   <- args[2]
output_csv <- args[3]

beta_dt <- fread(cpg_matrix)
rownames(beta_dt) <- beta_dt[[1]]
beta_dt[[1]] <- NULL

sample_ids <- names(beta_dt)

dmr_groups <- readRDS(dmr_rds)

out <- data.table(DMR = names(dmr_groups))

for (sid in sample_ids) out[[sid]] <- NA_real_

for (i in seq_along(dmr_groups)) {
  probes <- dmr_groups[[i]]
  present <- intersect(probes, rownames(beta_dt))

  if (length(present) == 0) next

  vals <- colMeans(beta_dt[present, , drop = FALSE], na.rm = TRUE)
  out[i, (sample_ids) := as.list(vals)]
}

fwrite(out, output_csv)