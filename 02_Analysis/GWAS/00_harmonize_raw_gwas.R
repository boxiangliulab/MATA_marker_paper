#!/usr/bin/env Rscript
#
# Step 0: raw sumstats (any columns from various GWAS resources) -> variant_id/beta/se/pvalue.
#
# Check before running: which column is the EFFECT allele. Not the same
# as ref/alt, and not guessable from the name - A1 is effect in
# PLINK/GCTA, but Allele2 is effect in SAIGE (reversed, it's in their
# docs), ALLELE1 in BOLT-LMM/regenie. Verify against the source's own
# docs and set --col-effect-allele/--col-other-allele accordingly - get
# it backwards and beta's sign flips with no error anywhere.
#
# regenie also reports LOG10P (-log10 p), not p -> --pvalue-is-neg-log10.
#
# variant_id = chr:pos:other_allele:effect_allele (effect allele last).
#
#   Rscript 00_harmonize_raw_gwas.R --input raw/bmi.txt.gz --output organized/bmi_hg19.txt \
#     --col-chr CHR --col-pos BP --col-other-allele A2 --col-effect-allele A1 \
#     --col-beta OR --effect-is-or --col-se SE --col-pvalue P --col-n N

suppressMessages(library(data.table))
source(file.path(
  { f <- grep("^--file=", commandArgs(FALSE), value = TRUE); if (length(f)) dirname(normalizePath(sub("^--file=", "", f))) else getwd() },
  "utils_gwas_prep.R"
))

opt <- parse_cli_args(list(
  input               = list(type = "character", required = TRUE),
  output              = list(type = "character", required = TRUE),
  `col-variant-id`    = list(type = "character", default = NA_character_),
  `col-chr`           = list(type = "character", default = NA_character_),
  `col-pos`           = list(type = "character", default = NA_character_),
  `col-other-allele`  = list(type = "character", default = NA_character_),
  `col-effect-allele` = list(type = "character", default = NA_character_),
  `col-beta`          = list(type = "character", required = TRUE),
  `col-se`            = list(type = "character", required = TRUE),
  `col-pvalue`        = list(type = "character", required = TRUE),
  `col-n`             = list(type = "character", default = NA_character_),
  `col-eaf`           = list(type = "character", default = NA_character_),
  `effect-is-or`      = list(type = "logical", default = FALSE),
  `pvalue-is-neg-log10` = list(type = "logical", default = FALSE)
))

has_combined_id <- !is.na(opt$`col-variant-id`)
has_split_id <- all(!is.na(c(opt$`col-chr`, opt$`col-pos`, opt$`col-other-allele`, opt$`col-effect-allele`)))
if (!has_combined_id && !has_split_id) {
  stop("need either --col-variant-id, or all of --col-chr/--col-pos/--col-other-allele/--col-effect-allele")
}

message("reading ", opt$input)
dat <- fread(opt$input)

required_cols <- c(
  if (has_combined_id) opt$`col-variant-id` else c(opt$`col-chr`, opt$`col-pos`, opt$`col-other-allele`, opt$`col-effect-allele`),
  opt$`col-beta`, opt$`col-se`, opt$`col-pvalue`,
  na.omit(c(opt$`col-n`, opt$`col-eaf`))
)
missing_cols <- setdiff(required_cols, names(dat))
if (length(missing_cols) > 0) {
  stop("missing column(s): ", paste(missing_cols, collapse = ", "),
       "\ncolumns present: ", paste(names(dat), collapse = ", "))
}

n_in <- nrow(dat)

if (has_combined_id) {
  dat[, variant_id := get(opt$`col-variant-id`)]
} else {
  # other:effect order -> beta is always "per copy of the last allele"
  dat[, variant_id := paste(get(opt$`col-chr`), get(opt$`col-pos`),
                             get(opt$`col-other-allele`), get(opt$`col-effect-allele`), sep = ":")]
}

# just fixes delimiter/chr-prefix/casing - has no idea which allele is
# the effect one, that had to already be right coming in
dat <- add_canonical_variant_id(dat, "variant_id")

beta_raw <- as.numeric(dat[[opt$`col-beta`]])
if (isTRUE(opt$`effect-is-or`)) {
  if (any(beta_raw <= 0, na.rm = TRUE)) stop("--effect-is-or set but --col-beta has values <= 0, can't log() those")
  dat[, beta := log(beta_raw)]
} else {
  dat[, beta := beta_raw]
}
dat[, se := as.numeric(get(opt$`col-se`))]

if (isTRUE(opt$`pvalue-is-neg-log10`)) {
  dat[, pvalue := 10^(-as.numeric(get(opt$`col-pvalue`)))]
} else {
  dat[, pvalue := as.numeric(get(opt$`col-pvalue`))]
}

out_cols <- c("variant_id", "beta", "se", "pvalue")
if (!is.na(opt$`col-n`)) {
  dat[, n := get(opt$`col-n`)]
  out_cols <- c(out_cols, "n")
}
if (!is.na(opt$`col-eaf`)) {
  dat[, eaf := get(opt$`col-eaf`)]
  out_cols <- c(out_cols, "eaf")
}

out <- dat[, ..out_cols]

n_na <- out[is.na(beta) | is.na(se) | is.na(pvalue), .N]
if (n_na > 0) {
  message("dropping ", n_na, " rows with missing beta/se/pvalue")
  out <- out[!(is.na(beta) | is.na(se) | is.na(pvalue))]
}

dup_id <- out[duplicated(variant_id), variant_id]
if (length(dup_id) > 0) {
  dup_path <- paste0(opt$output, ".duplicate_variant_ids.txt")
  fwrite(out[variant_id %in% dup_id], dup_path, sep = "\t")
  message(length(unique(dup_id)), " variant_id values are duplicated (", out[variant_id %in% dup_id, .N],
          " rows) -- written to ", dup_path, " for a look, not auto-resolved")
}

dir.create(dirname(opt$output), showWarnings = FALSE, recursive = TRUE)
fwrite(out, opt$output, sep = "\t")

message("done: ", n_in, " rows in -> ", nrow(out), " written to ", opt$output,
        " (coordinates untouched, run 01_liftover_to_hg38.R next)")
message("double check the effect allele call above against the source's own docs - this script just trusts what you gave it")
