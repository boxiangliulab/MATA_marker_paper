#!/usr/bin/env Rscript
#
# step 3: pull genome-wide-suggestive variants (p < --pval-threshold)
# out of the alt-aligned GWAS file (step 02's output), and optionally
# clump them down to independent lead loci (most significant per
# --clump-kb window). this is what gets handed to the coloc step.
# clumping was described in the old brain_sQTL coloc.R comments but
# never actually implemented there - it is here.
#
# writes <output> (all suggestive variants) and <output>.lead_loci.tsv
# (clumped leads, skipped with --clump-kb 0).
#
# Rscript 03_select_coloc_ready_loci.R --input processed_GWAS/t2d_hg38.aligned.txt \
#   --output coloc_ready_GWAS/t2d_hg38_suggestive.txt --pval-threshold 1e-7 --clump-kb 500

suppressMessages(library(data.table))
source(file.path(
  { f <- grep("^--file=", commandArgs(FALSE), value = TRUE); if (length(f)) dirname(normalizePath(sub("^--file=", "", f))) else getwd() },
  "utils_gwas_prep.R"
))

opt <- parse_cli_args(list(
  input            = list(type = "character", required = TRUE),
  output           = list(type = "character", required = TRUE),
  `pval-threshold` = list(type = "double", default = 1e-7),
  `clump-kb`       = list(type = "double", default = 500)
))

dir.create(dirname(opt$output), showWarnings = FALSE, recursive = TRUE)

dat <- fread(opt$input)
for (col in c("variant_id", "pvalue")) {
  if (!col %in% names(dat)) stop("missing required column '", col, "'")
}
dat <- add_canonical_variant_id(dat, "variant_id")

suggestive <- dat[pvalue < opt$`pval-threshold`]
fwrite(suggestive[, !c("chr", "pos", "ref", "alt")], opt$output, sep = "\t")
message(basename(opt$input), ": ", nrow(suggestive), " / ", nrow(dat), " variants pass p < ",
        opt$`pval-threshold`, " -> ", opt$output)

if (opt$`clump-kb` > 0 && nrow(suggestive) > 0) {
  leads <- clump_loci(suggestive, window_kb = opt$`clump-kb`)
  lead_path <- paste0(opt$output, ".lead_loci.tsv")
  fwrite(leads[, !c("chr", "pos", "ref", "alt")], lead_path, sep = "\t")
  message("clumped down to ", nrow(leads), " independent lead loci (>", 2 * opt$`clump-kb`,
          "kb apart) -> ", lead_path)
}
