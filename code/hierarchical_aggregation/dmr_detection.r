suppressPackageStartupMessages({
  library(limma)
  library(DMRcate)
  library(GenomicRanges)
  library(data.table)
  library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 4) {
  stop("Usage: Rscript dmr_detection_01.R <M_matrix.csv> <metadata.csv> <dmr_rds> <dmr_csv>")
}

M_file   <- args[1]
meta_file <- args[2]
dmr_rds  <- args[3]
dmr_csv  <- args[4]

M_dt <- fread(M_file)
rownames(M_dt) <- M_dt[[1]]
M_dt[[1]] <- NULL
M <- as.matrix(M_dt)
storage.mode(M) <- "double"

meta <- fread(meta_file)
meta <- as.data.frame(meta)
rownames(meta) <- meta[[1]]
meta[[1]] <- NULL
meta$disease <- factor(meta$disease)

if (!all(colnames(M) == rownames(meta))) {
  stop("Sample names in M and metadata do not match")
}

design <- model.matrix(~ disease, data = meta)
coef_name <- grep("^disease", colnames(design), value = TRUE)[1]

anno <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
anno$chr <- gsub("^chr", "", as.character(anno$chr))
anno <- anno[anno$chr %in% as.character(1:22), ]

fit <- lmFit(M, design)
fit <- eBayes(fit)

tt <- topTable(fit, coef = coef_name, number = Inf, adjust.method = "BH")
sig_cpgs <- rownames(tt)[tt$adj.P.Val < 0.05]

rm(fit, tt)
gc()

dmr_groups <- list()

for (chr in as.character(1:22)) {

  chr_cpgs <- intersect(sig_cpgs, rownames(anno)[anno$chr == chr])

  if (length(chr_cpgs) < 50) next

  M_chr <- M[chr_cpgs, , drop = FALSE]

  my_anno <- cpg.annotate(
    object = M_chr,
    datatype = "array",
    what = "M",
    analysis.type = "differential",
    design = design,
    coef = coef_name,
    arraytype = "450K"
  )

  dmrs <- dmrcate(my_anno)
  gr <- extractRanges(dmrs)

  if (length(gr) == 0) next

  chr_anno <- anno[anno$chr == chr, ]
  chr_gr <- GRanges(
    seqnames = paste0("chr", chr),
    ranges = IRanges(start = chr_anno$pos, width = 1)
  )
  names(chr_gr) <- rownames(chr_anno)

  for (i in seq_along(gr)) {
    hits <- findOverlaps(gr[i], chr_gr)
    cpgs <- names(chr_gr)[subjectHits(hits)]

    if (length(cpgs) > 0) {
      dmr_groups[[paste0("chr", chr, "_DMR_", i)]] <- cpgs
    }
  }

  rm(M_chr, my_anno, dmrs, gr)
  gc()
}

saveRDS(dmr_groups, dmr_rds)

if (length(dmr_groups) == 0) {
  stop("No DMRs detected")
}

dmr_df <- data.frame(
  DMR = names(dmr_groups),
  CpGs = vapply(dmr_groups, paste, collapse = ";", FUN.VALUE = character(1)),
  stringsAsFactors = FALSE
)

write.csv(dmr_df, dmr_csv, row.names = FALSE, quote = FALSE)