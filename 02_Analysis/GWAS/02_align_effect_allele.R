#!/usr/bin/env Rscript
#
# step 2: make sure beta ends up "per copy of alt", same convention the
# QTL side uses. a GWAS's effect allele has no reason to land on the
# genomic alt just by chance, so this looks up the actual reference base
# at each position (via samtools faidx on the hg38 fasta) and flips
# beta's sign wherever the effect allele turns out to be the reference
# base instead - the other allele then becomes alt by elimination.
#
# anything that doesn't match the reference base at all (even trying the
# opposite strand), or is an A/T or C/G pair, goes to <output>.needs_review.txt
# instead of being silently dropped - those two can't be told apart from
# the letters alone, needs an actual look (allele frequency, etc).
#
# Rscript 02_align_effect_allele.R --input processed_GWAS/t2d_hg38.txt \
#   --fasta /data/projects/11003054/share/data/reference/hg38/GRCh38.primary_assembly.genome.fa \
#   --output processed_GWAS/t2d_hg38.aligned.txt

suppressMessages(library(data.table))
source(file.path(
  { f <- grep("^--file=", commandArgs(FALSE), value = TRUE); if (length(f)) dirname(normalizePath(sub("^--file=", "", f))) else getwd() },
  "utils_gwas_prep.R"
))

opt <- parse_cli_args(list(
  input         = list(type = "character", required = TRUE),
  fasta         = list(type = "character", required = TRUE),
  output        = list(type = "character", required = TRUE),
  `samtools-bin` = list(type = "character", default = Sys.which("samtools")),
  review        = list(type = "character", default = NA_character_)
))
if (is.na(opt$review)) opt$review <- paste0(opt$output, ".needs_review.txt")

dat <- fread(opt$input)
for (col in c("variant_id", "beta")) {
  if (!col %in% names(dat)) stop("missing required column '", col, "'")
}
dat <- add_canonical_variant_id(dat, "variant_id")
n_in <- nrow(dat)

result <- align_to_reference_fasta(dat, opt$fasta, samtools_bin = opt$`samtools-bin`)
n_not_found <- n_in - nrow(result$aligned) - nrow(result$flagged)

out_cols <- intersect(c("variant_id", "beta", "se", "pvalue", "n", "eaf"), names(result$aligned))
dir.create(dirname(opt$output), showWarnings = FALSE, recursive = TRUE)
fwrite(result$aligned[, ..out_cols], opt$output, sep = "\t")

if (nrow(result$flagged) > 0) fwrite(result$flagged, opt$review, sep = "\t")

message(n_in, " in -> ", nrow(result$aligned), " aligned, ", nrow(result$flagged),
        " flagged for review, ", n_not_found, " had no reference base returned (off-contig / bad position?)")
if (nrow(result$flagged) > 0) {
  message(nrow(result$flagged), " variants written to ", opt$review, " - go look at those by hand before trusting them")
}
