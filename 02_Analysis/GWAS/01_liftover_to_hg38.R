#!/usr/bin/env Rscript
#
# step 1: output of 00_harmonize_raw_gwas.R -> GRCh38 (the build the
# adipose e/sQTL mapping used). hg38 input: just re-standardizes
# variant_id. hg19 input: lifts with UCSC liftOver, drops anything that
# fails to map (mapping rate gets printed - should be >99% for SNPs).
#
# Rscript 01_liftover_to_hg38.R --input organized/t2d_hg19.txt --output processed_GWAS/t2d_hg38.txt \
#   --build hg19 --chain /path/to/hg19ToHg38.over.chain

suppressMessages(library(data.table))
source(file.path(
  { f <- grep("^--file=", commandArgs(FALSE), value = TRUE); if (length(f)) dirname(normalizePath(sub("^--file=", "", f))) else getwd() },
  "utils_gwas_prep.R"
))

opt <- parse_cli_args(list(
  input          = list(type = "character", required = TRUE),
  output         = list(type = "character", required = TRUE),
  build          = list(type = "character", default = "hg19"),
  chain          = list(type = "character", default = NA_character_),
  `liftover-bin` = list(type = "character", default = Sys.which("liftOver"))
))

if (!opt$build %in% c("hg19", "hg38")) stop("--build must be hg19 or hg38")
if (opt$build == "hg19" && is.na(opt$chain)) stop("--chain is required for --build hg19")

message("reading ", opt$input)
dat <- fread(opt$input)
for (col in c("variant_id", "beta", "se", "pvalue")) {
  if (!col %in% names(dat)) stop("missing required column '", col, "'")
}
n_in <- nrow(dat)

dat <- add_canonical_variant_id(dat, "variant_id")

if (opt$build == "hg19") {
  message("lifting hg19 -> hg38...")
  dat <- liftover_positions(dat, chain_file = opt$chain, liftover_bin = opt$`liftover-bin`)
  dat[, variant_id := paste0(chr, ":", pos, ":", ref, ":", alt)]
} else {
  message("already hg38, skipping liftover")
}

out <- dat[, !c("chr", "pos", "ref", "alt")]

dir.create(dirname(opt$output), showWarnings = FALSE, recursive = TRUE)
fwrite(out, opt$output, sep = "\t")

message("done: ", n_in, " variants in -> ", nrow(out), " retained (",
        sprintf("%.1f%%", 100 * nrow(out) / n_in), "), written to ", opt$output)
