## Figure 1E. cell proportion comparison + LMM within cell group
library(dplyr)
library(glmmTMB)
library(broom.mixed)
library(dplyr)
library(emmeans)
library(ggplot2)
library(ggrepel)
library(ggpattern)
library(cowplot)
library(patchwork)
library(lmerTest)
library(tidyr)

# Set working directory to the project root before running
# setwd("path/to/your/project/root")

meta <- read.table("./08.integration_snRNA_object/83_libraries/00.annotation/02.annotation_level2/00.coordinate_per_cell_type/04.freeze_version/metadata_after_libraryQC_after_ctQC_freeze.tsv",sep = "\t",header = T)

colors_cg <- c("Immunes" = "#6baed6",
               "Vasculars" = "#df65b0",
               "Adipocytes" = "#addd8e",
               "FAPs" = "#DAE220")

immune <- c("Mast","B","Dendritic","Macrophage","Monocyte","Neutrophil","NK","T","Plasma_cell")
Vasculars <- c("Endothelial","Pericyte","SMC")
Adipocytes <- "Adipocyte"
FAPs <- c("ASPC","Mesothelial","IGFBP2+_cell")


meta$Cell_group <- "a"
meta$Cell_group[which(meta$Annotation_level1 %in% immune)] <- "Immunes"
meta$Cell_group[which(meta$Annotation_level1 %in% Vasculars)] <- "Vasculars"
meta$Cell_group[which(meta$Annotation_level1 %in% Adipocytes)] <- "Adipocytes"
meta$Cell_group[which(meta$Annotation_level1 %in% FAPs)] <- "FAPs"


########## Section 1: Cell proportion bar plot and donut plot ##########

# Bar plot
depot_col <- c("SAT" = "#BBCFC3", "VAT" = "#FBCE9E")
sat_col <- "#BBCFC3"
vat_col <- "#FBCE9E"
total_counts <- meta %>%
  group_by(Depot) %>%
  summarise(total_count = n())

cell_prop <- meta %>% group_by(Depot,Annotation_level2) %>% summarise(n_cells = n())
#cell_prop <- meta %>% group_by(Annotation_level2) %>% summarise(n_cells = n())
cell_prop <- merge(cell_prop,total_counts,by = "Depot")
#cell_prop$proportion <- round(cell_prop$n_cells / nrow(meta),digits = 3)
cell_prop$proportion <- round(cell_prop$n_cells / cell_prop$total_count,digits = 4)
#cell_prop$proportion <- cell_prop$proportion * 100
#order_levels <- cell_prop %>%
#  filter(Depot == "VAT") %>%
#  arrange(desc(proportion)) %>%
#  pull(Annotation_level2)

#cell_prop$Annotation_level2 <- factor(cell_prop$Annotation_level2, levels = order_levels)
df_plot <- cell_prop %>%
  select(Annotation_level2, Depot, proportion) %>%
  tidyr::pivot_wider(names_from = Depot, values_from = proportion)
write.table(df_plot,"./12.celltype_proportion/02.cellproportion_allcelltype_in_each_depot.tsv",sep = "\t",row.names = F,quote = F)
df <- df_plot %>%
  filter(!is.na(VAT)) %>%
  mutate(Annotation_level2 = factor(Annotation_level2,
                                    levels = Annotation_level2[order(-VAT)]))

df_long <- df %>%
  select(Annotation_level2, SAT, VAT) %>%
  pivot_longer(cols = c(SAT, VAT),
               names_to = "Depot",
               values_to = "Proportion")

df_long$Percentage <- round(df_long$Proportion * 100,4)
pdf("./12.celltype_proportion/02.cellproportion_allcelltype_in_each_depot.pdf",width = 10,height = 6)
ggplot(df_long, aes(x = Percentage,
                    y = Annotation_level2,
                    fill = Depot)) +
  geom_col(position = "stack", width = 0.7) +
  scale_fill_manual(values = depot_col) +
  labs(x = "Percentage (%)", y = "", fill = "Depot") +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.major.y = element_blank(),
    axis.text.y = element_text(size = 12),
    legend.position = "top"
  )
dev.off()


# Summarize data for the donut plot with SAT and VAT as rings
panel_a_data <- meta %>%
  group_by(Depot, Cell_group) %>%
  summarise(count = n()) %>%
  group_by(Depot) %>%
  mutate(percentage = count / sum(count) * 100) %>%
  ungroup()

# Assign a unique ring position for each Fat_type (e.g., SAT = 1, VAT = 2)
panel_a_data <- panel_a_data %>%
  mutate(ring_position = as.numeric(factor(Depot)))  # Numeric positions for each Fat_type

# Plot Panel A: Donut Plot with a smaller center and wider rings
p1<- ggplot(panel_a_data, aes(x = ring_position, y = percentage, fill = Cell_group)) +
  geom_bar(stat = "identity", color = "black", size = 0.2, width = 0.92) +
  geom_bar_pattern(
    data = subset(panel_a_data, ring_position == 2),
    stat = "identity",
    aes(pattern = Depot),
    color = NA,
    fill  = NA,
    width = 0.92,
    pattern = "circle",
    pattern_fill   = "black",
    pattern_colour = NA,
    pattern_density = 0.2,
    pattern_spacing = 0.02,
    pattern_scale   = 0.4,
    inherit.aes = TRUE
  ) +

  coord_polar(theta = "y") +
  scale_fill_manual(values = colors_cg) +
  theme_void() +
  labs(fill = "Cell Group") +
  theme(legend.position = "bottom") +
  xlim(0.4, 2.6) +

  # Only label sectors >2% to avoid overlapping
  geom_text(
    aes(label = ifelse(percentage > 2, paste0(round(percentage, 1), "%"), "")),
    position = position_stack(vjust = 0.5),
    size = 3, fontface = "bold"
  ) +

  annotate("text", x = 1, y = max(panel_a_data$percentage) + 8,
           label = "SAT", size = 5, fontface = "bold") +
  annotate("text", x = 2, y = max(panel_a_data$percentage) + 8,
           label = "VAT", size = 5, fontface = "bold")




# Summarize data for Panel B


# Calculate immune cell type proportions relative to the total cell count in each Fat_type
# Step 1: Calculate total cell count per Fat_type
total_counts <- meta %>%
  filter(Cell_group == "Immunes") %>%
  group_by(Depot) %>%
  summarise(total_count = n())

# Step 2: Filter for immune cells, calculate their count, and merge with total_counts to get proportions
panel_b_data <- meta %>%
  filter(Cell_group == "Immunes") %>%  # Only include Immune cell types
  group_by(Depot, Annotation_level2) %>%
  summarise(count = n()) %>%
  left_join(total_counts, by = "Depot") %>%  # Add total cell counts for each Fat_type
  mutate(percentage = (count / total_count) * 100) %>%  # Calculate percentage of total depot cells
  ungroup()

# Plot Panel B: Stacked Bar Plot with Immune cell types and their percentages relative to total depot cells
p2 <- ggplot(panel_b_data, aes(x = Depot, y = percentage, fill = Annotation_level2)) +
  geom_bar(stat = "identity", color = "white", width = 0.8, position = "stack") +

  #geom_bar_pattern(
  #  data = subset(panel_b_data, Depot == "VAT"),
  #  stat = "identity",
  #  aes(pattern = Depot),
  #  width = 0.8,
  #  position = "stack",
  #  color = NA,
  #  fill  = NA,
  #  pattern = "circle",
  #  pattern_fill    = "black",
  #  pattern_color   = NA,
  #  pattern_density = 0.2,
  #  pattern_spacing = 0.02,
  #  pattern_scale   = 0.20
  #) +

  scale_fill_manual(values = colors) +
  theme_minimal() +
  labs(x = "Depot", y = "Percentage", fill = "Immune Cell Type") +
  theme(
    legend.position = "bottom",
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 11, face = "bold")
  ) +
  geom_text(
    data = subset(panel_b_data, percentage > 2),
    aes(label = paste0(round(percentage, 1), "%")),
    size = 3, fontface = "bold", position = position_stack(vjust = 0.5)
  )


combined_plot <- plot_grid(p1, p2, labels = c("A", "B"), ncol = 2)

# Save combined plot to a file
ggsave("./12.celltype_proportion/01.combined_panels_AB_cell_type_proportion.pdf", plot = p1/p2, width = 3, height = 8, dpi = 300)


########## Section 2: LMM - cell proportion within each cell group ##########
# This is to fit cell proportion with depot within each cell group
# Caucasian samples are excluded for these analyses
meta <- meta[-which(meta$Race == "Caucasian"),]

# Doing test by cell groups

percentage_comparison <- function(cellgroup, meta){

  suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(tibble)
    library(ggplot2)
    library(purrr)
    library(lmerTest)
    library(broom.mixed)
    library(emmeans)
  })

  out_dir <- paste0("./08.integration_snRNA_object/83_libraries/00.annotation/Data_distribution/v3/", cellgroup)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  save_path <- paste0(out_dir, "/")

  # -----------------------------
  # Config
  # -----------------------------
  min_cells_per_donor_depot_ct <- 10   # NEW: each donorxdepotxcelltype must have >=10 nuclei
  min_total_cells_per_donor    <- 100  # keep donors with at least this total nuclei (across depots, within this cellgroup)
  min_paired_donors_per_ct     <- 10   # keep celltypes with >=10 paired donors
  depot_levels <- c("SAT","VAT")
  depot_col <- c("SAT" = "#BBCFC3", "VAT" = "#FBCE9E")

  squeeze01 <- function(p, eps = 1e-5) pmin(pmax(p, eps), 1 - eps)

  # -----------------------------
  # Subset meta
  # -----------------------------
  meta_sub <- meta %>% filter(Cell_group == cellgroup)

  # -----------------------------
  # 1) counts per donorxdepotxcelltype
  # -----------------------------
  counts <- meta_sub %>%
    count(donor_id, Depot, Annotation_level2, name = "n_cells") %>%
    mutate(
      Depot = factor(Depot, levels = depot_levels),
      donor_id = factor(donor_id),
      Annotation_level2 = factor(Annotation_level2)
    )

  # totals per donorxdepot (within this cellgroup)
  totals <- counts %>%
    group_by(donor_id, Depot) %>%
    summarise(total_cells = sum(n_cells), .groups = "drop")

  # donors with enough total nuclei in each depot(within cellgroup)
  valid_donors <- totals %>%
    group_by(donor_id,Depot) %>%
    summarise(total_cells_all = sum(total_cells), .groups = "drop") %>%
    filter(total_cells_all >= min_total_cells_per_donor) %>%
    pull(donor_id)

  # proportion table
  df_prop <- counts %>%
    left_join(totals, by = c("donor_id","Depot")) %>%
    mutate(proportion = ifelse(total_cells > 0, n_cells / total_cells, NA_real_)) %>%
    filter(donor_id %in% valid_donors)

  # -----------------------------
  # 2) join donor-level covariates
  # -----------------------------
  covars_donor <- meta_sub %>%
    distinct(donor_id, Race, Gender, BMI, Age, Self_defined_DM)

  df_prop <- df_prop %>%
    left_join(covars_donor, by = "donor_id") %>%
    mutate(
      BMI = as.numeric(BMI),
      Age = as.numeric(Age),
      Gender = factor(Gender),
      Race = factor(Race),
      Self_defined_DM = factor(Self_defined_DM),
      Depot = factor(Depot, levels = depot_levels)
    )

  # -----------------------------
  # 3) NEW filter: each donorxdepotxcelltype must have >=10 nuclei
  # -----------------------------
  df_prop <- df_prop %>% filter(n_cells >= min_cells_per_donor_depot_ct)

  # -----------------------------
  # 4) keep only paired donors (SAT & VAT) within each celltype
  # -----------------------------
  df_prop <- df_prop %>%
    group_by(donor_id, Annotation_level2) %>%
    filter(all(depot_levels %in% Depot)) %>%
    ungroup()

  # -----------------------------
  # 5) keep celltypes with >=10 paired donors
  # -----------------------------
  df_prop <- df_prop %>%
    group_by(Annotation_level2) %>%
    filter(n_distinct(donor_id) >= min_paired_donors_per_ct) %>%
    ungroup()

  # If nothing left, exit gracefully
  if (nrow(df_prop) == 0) {
    message("No data left after filtering for: ", cellgroup)
    return(invisible(NULL))
  }

  # =============================
  # A) SAT vs VAT (remove Caucasian): Paired t-test + plots with paired lines
  # =============================
  df_prop_sub <- df_prop %>% filter(Race != "Caucasian")

  # ---- paired t-test per celltype ----
  ttest_results <- list()
  for (ct in unique(df_prop_sub$Annotation_level2)) {
    df_ctw <- df_prop_sub %>%
      filter(Annotation_level2 == ct) %>%
      select(donor_id, Depot, proportion) %>%
      distinct() %>%
      pivot_wider(names_from = Depot, values_from = proportion) %>%
      drop_na(SAT, VAT)

    if (nrow(df_ctw) < min_paired_donors_per_ct) next

    test_res <- t.test(df_ctw$VAT, df_ctw$SAT, paired = TRUE)

    ttest_results[[as.character(ct)]] <- tibble(
      Annotation_level2 = as.character(ct),
      n_paired = nrow(df_ctw),
      mean_SAT = mean(df_ctw$SAT),
      mean_VAT = mean(df_ctw$VAT),
      mean_diff = mean(df_ctw$VAT - df_ctw$SAT),
      conf.low = test_res$conf.int[1],
      conf.high = test_res$conf.int[2],
      t = unname(test_res$statistic),
      p_value = test_res$p.value
    )
  }

  df_ttest <- bind_rows(ttest_results) %>%
    mutate(p_adj = p.adjust(p_value, method = "fdr"),
           label = paste0("p=", signif(p_value, 2), "\nFDR=", signif(p_adj, 2), "\nN=", n_paired))

  write.table(df_ttest, paste0(save_path, "01.paired_ttest_cellproportion_by_depot.noCaucasian.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)

  # ---- errorbar (paired t-test) ----
  pdf(paste0(save_path, "Depot_effect_CI_paired_ttest.noCaucasian.pdf"), width = 4.8, height = 4.8)
  p_ci <- ggplot(df_ttest, aes(x = mean_diff, y = reorder(Annotation_level2, mean_diff))) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2, color = "grey60") +
    geom_point(aes(color = p_adj < 0.05), size = 1.7) +
    scale_color_manual(values = c("TRUE" = "#E64B35", "FALSE" = "grey40")) +
    theme_bw() +
    labs(x = "Mean difference (VAT - SAT) with 95% CI",
         y = NULL, color = "FDR < 0.05",
         title = "Paired t-test: VAT vs SAT (proportion)")
  print(p_ci)
  dev.off()

  # ---- facet plot with paired lines + stats ----
  df_line <- df_prop_sub %>%
    select(donor_id, Annotation_level2, Depot, proportion) %>%
    distinct() %>%
    pivot_wider(names_from = Depot, values_from = proportion) %>%
    drop_na(SAT, VAT)

  df_anno <- df_prop_sub %>%
    group_by(Annotation_level2) %>%
    summarise(y = max(proportion, na.rm = TRUE) * 1.05, .groups = "drop") %>%
    mutate(Annotation_level2 = as.character(Annotation_level2)) %>%
    left_join(df_ttest %>% mutate(Annotation_level2 = as.character(Annotation_level2)),
              by = "Annotation_level2")

  pdf(paste0(save_path, "Depot_Celltype_Proportion_separate_paired.noCaucasian.pdf"), width = 7, height = 5)
  p_sep <- ggplot(df_prop_sub, aes(x = Depot, y = proportion)) +
    geom_segment(
      data = df_line,
      aes(x = "SAT", xend = "VAT", y = SAT, yend = VAT, group = donor_id),
      inherit.aes = FALSE,
      alpha = 0.25,
      linewidth = 0.25
    ) +
    geom_boxplot(aes(fill = Depot), outlier.shape = NA, width = 0.6, alpha = 0.8) +
    geom_jitter(aes(color = Depot), width = 0.12, alpha = 0.6, size = 1.0) +
    facet_wrap(~Annotation_level2, scales = "free_y") +
    scale_fill_manual(values = depot_col) +
    scale_color_manual(values = depot_col) +
    geom_text(
      data = df_anno,
      aes(x = "SAT", y = y, label = label),
      inherit.aes = FALSE,
      size = 2.8,
      lineheight = 0.95
    ) +
    theme_bw() +
    labs(x = "Depot", y = "Proportion",
         title = paste0(cellgroup, ": VAT vs SAT (paired donors; no Caucasian)")) +
    theme(strip.text = element_text(size = 10, face = "bold"),
          panel.grid = element_blank())
  print(p_sep)
  dev.off()

  # ---- merged plot with paired lines + stats ----
  ct_levels <- levels(factor(df_prop_sub$Annotation_level2))
  df_plot2 <- df_prop_sub %>%
    mutate(ct = factor(Annotation_level2, levels = ct_levels),
           idx = as.numeric(ct),
           x = ifelse(Depot == "SAT", idx - 0.2, idx + 0.2))

  df_line2 <- df_plot2 %>%
    select(donor_id, ct, idx, Depot, proportion) %>%
    distinct() %>%
    pivot_wider(names_from = Depot, values_from = proportion) %>%
    drop_na(SAT, VAT) %>%
    mutate(x_sat = idx - 0.2, x_vat = idx + 0.2)

  df_anno2 <- df_plot2 %>%
    group_by(ct, idx) %>%
    summarise(y = max(proportion, na.rm = TRUE) * 1.08, .groups="drop") %>%
    left_join(
      df_ttest %>% mutate(ct = factor(Annotation_level2, levels = ct_levels),
                          idx = as.numeric(ct)),
      by = c("ct","idx")
    )

  pdf(paste0(save_path, "Depot_Celltype_Proportion_merged_paired.noCaucasian.pdf"), width = 10, height = 4)
  p_mer <- ggplot(df_plot2, aes(x = x, y = proportion)) +
    geom_segment(
      data = df_line2,
      aes(x = x_sat, xend = x_vat, y = SAT, yend = VAT, group = interaction(ct, donor_id)),
      inherit.aes = FALSE,
      alpha = 0.20,
      linewidth = 0.20
    ) +
    geom_boxplot(
      aes(group = interaction(ct, Depot), fill = Depot),
      outlier.shape = NA,
      width = 0.32,
      alpha = 0.85
    ) +
    geom_point(aes(color = Depot), position = position_jitter(width = 0.06), alpha = 0.55, size = 0.35) +
    scale_fill_manual(values = depot_col) +
    scale_color_manual(values = depot_col) +
    scale_x_continuous(breaks = seq_along(ct_levels), labels = ct_levels) +
    geom_text(
      data = df_anno2,
      aes(x = idx, y = y, label = paste0("p=", signif(p_value,2), "\nFDR=", signif(p_adj,2))),
      inherit.aes = FALSE,
      size = 2.4,
      lineheight = 0.9
    ) +
    theme_bw() +
    labs(x = "Cell type", y = "Proportion",
         title = paste0(cellgroup, ": VAT vs SAT across cell types (paired; no Caucasian)")) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
      panel.grid = element_blank(),
      legend.position = "top"
    )
  print(p_mer)
  dev.off()

  # =============================
  # B) SAT vs VAT (no Caucasian): LMM adjusted effect (optional output)
  # =============================
  fit_results <- list()
  for (ct in unique(df_prop_sub$Annotation_level2)) {
    df_ct <- df_prop_sub %>%
      filter(Annotation_level2 == ct) %>%
      mutate(
        donor_id = factor(donor_id),
        Depot = factor(Depot, levels = depot_levels),
        Race = factor(Race),
        Gender = factor(Gender),
        Self_defined_DM = factor(Self_defined_DM)
      )

    if (n_distinct(df_ct$donor_id) < min_paired_donors_per_ct) next

    fit <- tryCatch(
      lmerTest::lmer(
        proportion ~ Depot + Age + Gender + BMI + Self_defined_DM + Race + (1|donor_id),
        data = df_ct
      ),
      error = function(e) NULL
    )
    if (is.null(fit)) next

    coef_df <- broom.mixed::tidy(fit, effects = "fixed", conf.int = TRUE) %>%
      filter(term == "DepotVAT") %>%
      mutate(celltype = as.character(ct),
             n_paired = n_distinct(df_ct$donor_id))
    fit_results[[as.character(ct)]] <- coef_df
  }

  res_lmm <- bind_rows(fit_results) %>%
    mutate(p_adj = p.adjust(p.value, method = "fdr")) %>%
    select(celltype, n_paired, term, estimate, std.error, conf.low, conf.high, statistic, p.value, p_adj)

  write.table(res_lmm, paste0(save_path,"01.LMM_cellproportion_by_depot.noCaucasian.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)

  pdf(paste0(save_path,"Depot_effect_CI_LMM.noCaucasian.pdf"), width = 4.8, height = 4.8)
  p_lmm <- ggplot(res_lmm, aes(x = estimate, y = reorder(celltype, estimate))) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2, color = "grey60") +
    geom_point(aes(color = p_adj < 0.05), size = 1.7) +
    scale_color_manual(values = c("TRUE" = "#E64B35", "FALSE" = "grey40")) +
    theme_bw() +
    labs(x = "Adjusted Depot effect (VAT vs SAT)",
         y = NULL, color = "FDR < 0.05",
         title = "LMM: adjusted VAT vs SAT (95% CI; no Caucasian)")
  print(p_lmm)
  dev.off()

  # =============================
  # C) Race x Depot (paired, logit(proportion) LMM instead of beta/logit GLMM)
  # =============================
  results_rxD <- list()
  for (cell in unique(df_prop$Annotation_level2)) {

    df_sub <- df_prop %>%
      filter(Annotation_level2 == cell) %>%
      mutate(
        donor_id = factor(donor_id),
        Depot = factor(Depot, levels = depot_levels),
        Race = factor(Race),
        Gender = factor(Gender),
        Self_defined_DM = factor(Self_defined_DM),
        logit_p = qlogis(squeeze01(proportion))
      )

    if (n_distinct(df_sub$Race) < 2 || n_distinct(df_sub$Depot) < 2) next
    if (n_distinct(df_sub$donor_id) < min_paired_donors_per_ct) next

    covars <- c("Gender", "BMI", "Age", "Self_defined_DM")
    valid_covars <- covars[sapply(covars, function(v){
      x <- df_sub[[v]]
      if (is.numeric(x)) return(sd(x, na.rm = TRUE) > 0)
      length(unique(x[!is.na(x)])) > 1
    })]

    rhs <- paste(c("Race * Depot", valid_covars, "(1|donor_id)"), collapse = " + ")
    fml <- as.formula(paste0("logit_p ~ ", rhs))

    fit <- tryCatch(lmerTest::lmer(fml, data = df_sub), error = function(e) NULL)
    if (is.null(fit)) next

    emm <- emmeans(fit, pairwise ~ Race | Depot)
    res <- as.data.frame(emm$contrasts)
    res$celltype <- as.character(cell)
    results_rxD[[as.character(cell)]] <- res
  }

  df_res <- bind_rows(results_rxD)
  write.table(df_res, paste0(save_path,"01.LMM_logitproportion_Race_by_Depot.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)

  # ---- Race plot (paired filter already in df_prop) ----
  race_col <- c(
    "Indian"   = "#55A6A2",
    "Malay"    = "#9CD17B",
    "Chinese"  = "#666D9F",
    "Caucasian"= "#1b7837"
  )

  # Only plot celltypes that survived filters
  celltypes <- levels(factor(df_prop$Annotation_level2))

  plot_celltype <- function(ct){
    df_ct <- df_prop %>% filter(Annotation_level2 == ct)
    p <- ggplot(df_ct, aes(x = Race, y = proportion, fill = Race)) +
      geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.8) +
      geom_jitter(aes(color = Race), width = 0.15, alpha = 0.7, size = 1.0) +
      scale_fill_manual(values = race_col) +
      scale_color_manual(values = race_col) +
      facet_wrap(~Depot, ncol = 2, scales = "free_y") +
      theme_bw() +
      labs(x = "Race", y = "Proportion", title = ct) +
      theme(
        strip.background = element_blank(),
        strip.text = element_text(size = 10, face = "bold"),
        axis.text.x = element_text(angle = 30, hjust = 1, size = 8),
        panel.grid = element_blank(),
        legend.position = "none"
      )
    p
  }

  # if you have patchwork loaded elsewhere, it will work; otherwise just output each page
  if (requireNamespace("patchwork", quietly = TRUE)) {
    library(patchwork)
    plots <- lapply(celltypes, plot_celltype)
    combined_plot <- wrap_plots(plots, ncol = 3) +
      plot_annotation(
        title = paste0(cellgroup, ": Cell-type proportions across Race and Depot (paired donors)"),
        theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
      )
    pdf(paste0(save_path,"Race_Depot_Celltype_Proportion_all_fixed.paired.pdf"), width = 10, height = 9)
    print(combined_plot)
    dev.off()
  } else {
    pdf(paste0(save_path,"Race_Depot_Celltype_Proportion_all_fixed.paired.pdf"), width = 8, height = 8)
    for (ct in celltypes) print(plot_celltype(ct))
    dev.off()
  }

  invisible(list(
    df_prop = df_prop,
    df_prop_noCaucasian = df_prop_sub,
    ttest = df_ttest,
    lmm_depot = res_lmm,
    lmm_race_depot = df_res
  ))
}

cell_group <- c("Vasculars","FAPs","Immunes","Adipocytes")
for (cg in cell_group){
  percentage_comparison(cg, meta)
}

##### By paired depot, plot each group proportion #####
outdir <- "./08.integration_snRNA_object/83_libraries/00.annotation/Data_distribution/v3/all/"


# 1. count: donor x depot x Cell_group
counts <- meta %>%
  count(donor_id, Depot, Cell_group, name = "n_cells")

# 2. total per donor x depot
totals <- counts %>%
  group_by(donor_id, Depot) %>%
  summarise(total_cells = sum(n_cells), .groups = "drop")

# 3. proportion
df_cg <- counts %>%
  left_join(totals, by = c("donor_id","Depot")) %>%
  mutate(prop = n_cells / total_cells) %>%
  mutate(Depot = factor(Depot, levels = c("SAT","VAT")))

df_cg <- df_cg %>%
  group_by(donor_id, Cell_group) %>%
  filter(all(c("SAT","VAT") %in% Depot)) %>%
  ungroup()

ttest_cg <- df_cg %>%
  select(donor_id, Depot, Cell_group, prop) %>%
  pivot_wider(names_from = Depot, values_from = prop) %>%
  drop_na(SAT, VAT) %>%
  group_by(Cell_group) %>%
  summarise(
    n_paired = n(),
    mean_SAT = mean(SAT),
    mean_VAT = mean(VAT),
    mean_diff = mean(VAT - SAT),
    p_value = t.test(VAT, SAT, paired = TRUE)$p.value,
    .groups = "drop"
  ) %>%
  mutate(p_adj = p.adjust(p_value, method = "fdr"))

depot_col <- c("SAT"="#BBCFC3", "VAT"="#FBCE9E")

df_line <- df_cg %>%
  select(donor_id, Depot, Cell_group, prop) %>%
  pivot_wider(names_from = Depot, values_from = prop) %>%
  drop_na(SAT, VAT)

df_anno <- df_line %>%
  group_by(Cell_group) %>%
  summarise(y = max(c(SAT, VAT), na.rm = TRUE) * 1.1, .groups="drop") %>%
  left_join(ttest_cg, by="Cell_group") %>%
  mutate(label = paste0("p=", signif(p_value,2), "\nFDR=", signif(p_adj,2), "\nN=", n_paired))

pdf(paste0(outdir,"/Cell_group_all.pdf"),width = 12,height = 10)
ggplot(df_cg, aes(x = Depot, y = prop)) +
  geom_segment(
    data = df_line,
    aes(x="SAT", xend="VAT", y=SAT, yend=VAT, group=donor_id),
    inherit.aes = FALSE,
    alpha=0.3
  ) +
  geom_boxplot(aes(fill=Depot), outlier.shape=NA, width=0.6, alpha=0.8) +
  geom_jitter(aes(color=Depot), width=0.12, alpha=0.6) +
  facet_wrap(~Cell_group, scales="free_y") +
  scale_fill_manual(values=depot_col) +
  scale_color_manual(values=depot_col) +
  geom_text(data=df_anno,
            aes(x="SAT", y=y, label=label),
            inherit.aes=FALSE, size=3) +
  theme_bw() +
  labs(title="Cell_group proportions: VAT vs SAT (paired donors)",
       x="Depot", y="Proportion")

dev.off()

######### Plot celltype enriched in cell group ######
outdir <- "./08.integration_snRNA_object/83_libraries/00.annotation/Data_distribution/v3/T2D_vs_non_T2D/"
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

# -------------------------
# Settings
# -------------------------
dm_keep <- c("Non_T2DM", "T2DM")
min_cells_per_sample_ct <- 10
min_total_cells_per_sample <- 100
min_samples_per_ct_depot <- 10

m <- meta %>%
  filter(
    !is.na(Cell_group), !is.na(Annotation_level2),
    !is.na(Self_defined_DM), Self_defined_DM %in% dm_keep,
    !is.na(donor_id), !is.na(Depot),
    !is.na(Gender), !is.na(Race), !is.na(Age), !is.na(BMI)
  ) %>%
  mutate(
    Self_defined_DM = factor(Self_defined_DM, levels = dm_keep),
    Depot = factor(Depot, levels = c("SAT","VAT")),
    Cell_group = factor(Cell_group),
    Annotation_level2 = factor(Annotation_level2),
    Gender = factor(Gender),
    Race = factor(Race),
    Age = as.numeric(Age),
    BMI = as.numeric(BMI)
  )

# -------------------------
# 1) Build df: proportion within Cell_group
#    sample = donor_id x Depot (NOT paired)
# -------------------------
counts <- m %>%
  count(donor_id, Depot, Self_defined_DM, Cell_group, Annotation_level2, name = "n_cells")

denom <- counts %>%
  group_by(donor_id, Depot, Self_defined_DM, Cell_group) %>%
  summarise(den_cells = sum(n_cells), .groups = "drop") %>%
  filter(den_cells >= min_total_cells_per_sample)

df <- counts %>%
  inner_join(denom, by = c("donor_id","Depot","Self_defined_DM","Cell_group")) %>%
  filter(n_cells >= min_cells_per_sample_ct) %>%
  mutate(proportion = n_cells / den_cells)

# join covariates back (donor-level; duplicated across depots is OK)
covars <- m %>%
  distinct(donor_id, Depot, Self_defined_DM, Gender, Race, Age, BMI)

df <- df %>%
  left_join(covars, by = c("donor_id","Depot","Self_defined_DM"))

# keep enough samples PER DEPOT, and must have both DM levels within that depot
df <- df %>%
  group_by(Depot, Cell_group, Annotation_level2) %>%
  filter(
    n() >= min_samples_per_ct_depot,
    n_distinct(Self_defined_DM) == 2
  ) %>%
  ungroup()

write.table(df,
            file = file.path(outdir, "00.proportion_table.donorXdepot_within_Cellgroup.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)



# -------------------------
# 2) Stats: within EACH Depot, test T2DM vs Non_T2DM
#    model: proportion ~ Self_defined_DM + Gender + Age + BMI + Race
#
# -------------------------
res <- df %>%
  group_by(Depot, Cell_group, Annotation_level2) %>%
  group_modify(~{
    d <- .x

    # ensure both groups
    if (n_distinct(d$Self_defined_DM) < 2) return(tibble())
    if (nrow(d) < min_samples_per_ct_depot) return(tibble())

    # keep cov with multiple levels / non-zero var
    cov_list <- c("Gender","Age","BMI","Race")
    valid_cov <- cov_list[sapply(cov_list, function(v){
      x <- d[[v]]
      if (is.numeric(x)) return(sd(x, na.rm = TRUE) > 0)
      length(unique(x[!is.na(x)])) > 1
    })]
    rhs <- paste(c("Self_defined_DM", valid_cov), collapse = " + ")
    fml <- as.formula(paste0("proportion ~ ", rhs))

    fit <- tryCatch(lm(fml, data = d), error = function(e) NULL)
    if (is.null(fit)) return(tibble())

    td <- broom::tidy(fit, conf.int = TRUE)
    term <- paste0("Self_defined_DM", dm_keep[2])  # T2DM effect vs baseline

    if (!term %in% td$term) return(tibble())

    out <- td %>%
      filter(term == !!term) %>%
      transmute(
        estimate = estimate,
        conf.low = conf.low,
        conf.high = conf.high,
        statistic = statistic,
        p_value = p.value
      )

    tibble(
      n_total = nrow(d),
      n_NonT2DM = sum(d$Self_defined_DM == dm_keep[1]),
      n_T2DM = sum(d$Self_defined_DM == dm_keep[2]),
      covars_used = ifelse(length(valid_cov) == 0, "(DM only)", paste(valid_cov, collapse = "+"))
    ) %>% bind_cols(out)
  }) %>%
  ungroup() %>%
  group_by(Depot, Cell_group) %>%
  mutate(p_adj = p.adjust(p_value, method = "fdr")) %>%
  ungroup()

write.table(res,
            file = file.path(outdir, "01.LM_T2DM_vs_NonT2DM_byDepot_Cellgroup_AnnotationL2_cov.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# -------------------------
# 3) Quick diagnostic plot (optional): per depot, per cellgroup page
# -------------------------
anno <- df %>%
  group_by(Depot, Cell_group, Annotation_level2) %>%
  summarise(y = max(proportion, na.rm = TRUE) * 1.08, .groups="drop") %>%
  left_join(res, by = c("Depot","Cell_group","Annotation_level2")) %>%
  mutate(label = paste0("p=", signif(p_value,2), "\nFDR=", signif(p_adj,2),
                        "\nN=", n_total))

pdf(file.path(outdir, "02.T2DM_vs_NonT2DM_byDepot_Cellgroup_AnnotationL2.pdf"),
    width = 11, height = 8)

for (dp in levels(df$Depot)) {
  for (cg in sort(unique(as.character(df$Cell_group)))) {
    dfg <- df %>% filter(Depot == dp, Cell_group == cg)
    if (nrow(dfg) == 0) next
    annog <- anno %>% filter(Depot == dp, Cell_group == cg)
    p <- ggplot(dfg, aes(x = Self_defined_DM, y = proportion, fill = Self_defined_DM)) +
      geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.85) +
      geom_jitter(width = 0.15, alpha = 0.6, size = 0.6) +
      facet_wrap(~Annotation_level2, scales = "free_y") +
      geom_text(
        data = annog,
        aes(x = dm_keep[2], y = y, label = label),
        inherit.aes = FALSE,
        size = 2.5,
        hjust = 1
      ) +
      theme_bw() +
      labs(
        title = paste0(dp, " | ", cg, " | T2DM vs Non_T2DM (LM; cov: Gender+Age+BMI+Race)"),
        x = "Self_defined_DM",
        y = "Proportion within Cell_group"
      ) +
      theme(
        strip.text = element_text(size = 9, face = "bold"),
        axis.text.x = element_text(angle = 20, hjust = 1),
        panel.grid = element_blank(),
        legend.position = "none"
      )
    print(p)
  }
}
dev.off()
