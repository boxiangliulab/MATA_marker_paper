library(data.table)
library(ggplot2)

# Project root is 8 levels up from this script
args0        <- commandArgs(trailingOnly = FALSE)
script_path  <- sub("--file=", "", args0[grep("--file=", args0)])
project_root <- normalizePath(file.path(dirname(normalizePath(script_path)),
                                        paste(rep("..", 8), collapse = "/")))

fp <- file.path(project_root, "data/09.eQTL/03.Saige_QTL/02.output/05.snpEff/02.output/ALL.leadSNP.variant_effects.qtlBHM_annotations.txt")
dt <- fread(fp)

# 95% CI
dt[, lo := Weight - 1.96 * Stderr]
dt[, hi := Weight + 1.96 * Stderr]

# sort ascending
dt <- dt[order(Weight)]
dt[, Annotation := factor(Annotation, levels = Annotation)]

p <- ggplot(dt, aes(x = Weight, y = Annotation)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.8, color = "grey40") +
  geom_errorbarh(aes(xmin = lo, xmax = hi),
                 height = 0.18, linewidth = 0.75, color = "firebrick") +  # end caps
  geom_point(size = 2, color = "firebrick") +
  scale_y_discrete(expand = expansion(mult = c(0.08, 0.08))) +          # padding above/below
  labs(x = "ALL lead eSNP weight", y = NULL) +
  theme(
    panel.grid.minor = element_blank(),
    axis.title.x = element_text(face = "bold"),
    axis.text.y = element_text(size = 12)
  ) + theme_classic()
ggsave(file.path(project_root, "data/09.eQTL/03.Saige_QTL/02.output/05.snpEff/02.output/ALL_lead_eSNP_weight_qtlBHM.pdf"), p, width = 5, height = 4)
