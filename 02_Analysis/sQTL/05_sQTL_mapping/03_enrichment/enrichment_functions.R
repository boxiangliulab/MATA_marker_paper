#!/usr/bin/env Rscript
# Shared enrichment functions.
#

suppressPackageStartupMessages({
  library(dplyr)
})

# Read a SnpEff-annotated VCF and return one row per variant and annotation.
# The annotations live in the ANN= field of INFO, whose pipe-separated fields
# are Allele | Annotation | Impact | Gene_Name | Gene_ID | Feature_Type | ...
# A variant with several annotations yields several rows.
parse_snpeff_vcf <- function(vcf_file) {
  vcf_lines <- readLines(vcf_file)
  header_idx <- max(grep("^##", vcf_lines))
  col_names <- strsplit(gsub("^#", "", vcf_lines[header_idx + 1]), "\t")[[1]]
  data_lines <- vcf_lines[(header_idx + 2):length(vcf_lines)]

  vcf <- data.frame(do.call(rbind, strsplit(data_lines, "\t")),
                    stringsAsFactors = FALSE)
  colnames(vcf) <- col_names
  vcf$variant_id <- paste(vcf$CHROM, vcf$POS, vcf$REF, vcf$ALT, sep = ":")

  annotations <- lapply(vcf$INFO, function(info) {
    if (grepl("ANN=", info)) {
      ann_string <- sub(".*ANN=([^;]+).*", "\\1", info)
      ann_list <- strsplit(ann_string, ",")[[1]]
      ann_parsed <- lapply(ann_list, function(ann) {
        f <- strsplit(ann, "\\|")[[1]]
        field <- function(i) if (length(f) >= i) f[i] else NA
        data.frame(Allele       = field(1),
                   Annotation   = field(2),
                   Impact       = field(3),
                   Gene_Name    = field(4),
                   Gene_ID      = field(5),
                   Feature_Type = field(6),
                   stringsAsFactors = FALSE)
      })
      do.call(rbind, ann_parsed)
    } else {
      data.frame(Allele = NA, Annotation = NA, Impact = NA,
                 Gene_Name = NA, Gene_ID = NA, Feature_Type = NA,
                 stringsAsFactors = FALSE)
    }
  })

  annotated <- do.call(rbind, lapply(seq_along(annotations), function(i) {
    cbind(vcf[i, c("CHROM", "POS", "REF", "ALT", "variant_id")],
          annotations[[i]])
  }))

  rownames(annotated) <- NULL
  annotated %>% filter(!is.na(Annotation) & Annotation != "")
}

# Fisher's exact test of one annotation in the foreground vs the background set,
# returned as a one-row data frame with the log odds ratio and its 95% CI.
enrichment_fisher <- function(label, n_fg_with, n_fg_total,
                              n_bg_with, n_bg_total) {
  cont_table <- matrix(c(n_fg_with, n_fg_total - n_fg_with,
                         n_bg_with, n_bg_total - n_bg_with),
                       nrow = 2, byrow = TRUE)

  fisher_test <- tryCatch(
    fisher.test(cont_table, alternative = "two.sided"),
    error = function(e) list(estimate = NA, p.value = NA, conf.int = c(NA, NA)))

  data.frame(
    label         = label,
    n_sig         = n_fg_with,
    n_sig_total   = n_fg_total,
    prop_sig      = n_fg_with / n_fg_total,
    n_bg          = n_bg_with,
    n_bg_total    = n_bg_total,
    prop_bg       = n_bg_with / n_bg_total,
    odds_ratio    = fisher_test$estimate,
    log_odds      = log(fisher_test$estimate),
    log_ci_lower  = log(fisher_test$conf.int[1]),
    log_ci_upper  = log(fisher_test$conf.int[2]),
    pvalue        = fisher_test$p.value,
    stringsAsFactors = FALSE
  )
}

# Infinite log odds (zero counts in one cell) are not plottable.
drop_infinite <- function(df) {
  for (col in c("log_odds", "log_ci_lower", "log_ci_upper")) {
    df[[col]][is.infinite(df[[col]])] <- NA
  }
  df
}
