# ============================================================
# FIGURE 5 — SNP Dosage Class Distributions (G3 Revision)
# Three changes from original Figure6_regenerated_G3.R:
#   1. Dosage class labels: zeros→Nulliplex, simplex→Simplex, etc.
#   2. Population axis labels: D→Diploid, T→Tetraploid
#   3. Palette keys updated to match new capitalized names
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
})

# ----------------------------
# PATHS
# ----------------------------
base_dir <- "D:/Wenqian_G3_Figures"

file_s2039 <- file.path(base_dir, "S2039_prop.csv")
file_s2172 <- file.path(base_dir, "s2172_prop.csv")

stopifnot(file.exists(file_s2039), file.exists(file_s2172))

out_dir  <- file.path(base_dir, "G3_Regenerated_Figures")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

fig5_dir <- file.path(out_dir, "Figure5")
if (!dir.exists(fig5_dir)) dir.create(fig5_dir, recursive = TRUE)

out_prefix <- file.path(fig5_dir, "Figure5_SNP_dosage_class_distributions")

# ----------------------------
# CHANGE 1 — Dosage class labels (display names)
# ----------------------------
marker_display_levels <- c("Nulliplex", "Simplex", "Duplex", "Triplex", "Tetraplex")

# Column names in CSV remain unchanged
marker_cols <- c("zeros", "simplex", "duplex", "triplex", "tetraplex")

# Mapping from CSV column → display label
marker_label_map <- setNames(marker_display_levels, marker_cols)

# ----------------------------
# TLI MANUSCRIPT COLOR PALETTE
# ----------------------------
TLI_BROWN     <- "#5C3A21"
TLI_GREEN     <- "#5E8C31"
TLI_DARKGREEN <- "#2F6B2F"
TLI_BLUE      <- "#2E6F95"
TLI_PURPLE    <- "#7B3F8C"

# CHANGE 3 — Palette keys use new capitalized display names
marker_palette <- c(
  Nulliplex = TLI_BROWN,
  Simplex   = TLI_GREEN,
  Duplex    = TLI_DARKGREEN,
  Triplex   = TLI_BLUE,
  Tetraplex = TLI_PURPLE
)

se <- function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))

# ----------------------------
# READ + RESHAPE
# ----------------------------
read_prop <- function(path, study_label) {
  df <- read_csv(path, show_col_types = FALSE)
  if ("...1" %in% names(df)) df <- df %>% select(-`...1`)
  
  pop_candidates <- c("Population", "population", "pop", "ploidy", "Group", "group")
  pop_col <- intersect(pop_candidates, names(df))
  if (length(pop_col) == 0) stop("No population column found in: ", path)
  pop_col <- pop_col[1]
  
  mcols <- intersect(marker_cols, names(df))
  if (length(mcols) != length(marker_cols))
    stop("Missing marker columns in: ", path)
  
  df %>%
    rename(Population = !!pop_col) %>%
    mutate(Study = study_label) %>%
    pivot_longer(
      cols      = all_of(mcols),
      names_to  = "marker_raw",
      values_to = "value"
    ) %>%
    # CHANGE 1 — recode column names to display labels
    mutate(
      marker_type = factor(
        marker_label_map[marker_raw],
        levels = marker_display_levels
      )
    ) %>%
    # CHANGE 2 — recode D/T to Diploid/Tetraploid
    mutate(
      Population = dplyr::recode(
        as.character(Population),
        "D" = "Diploid",
        "T" = "Tetraploid"
      ),
      Population = factor(Population, levels = c("Diploid", "Tetraploid"))
    ) %>%
    filter(!is.na(Population), !is.na(marker_type))
}

summarise_prop <- function(df_long) {
  df_long %>%
    group_by(Study, Population, marker_type) %>%
    summarise(
      mean = mean(value, na.rm = TRUE),
      se   = se(value),
      .groups = "drop"
    )
}

n_by_ploidy <- function(df_long) {
  if ("ind" %in% names(df_long)) {
    df_long %>%
      distinct(Study, Population, ind) %>%
      count(Study, Population, name = "n")
  } else {
    df_long %>%
      group_by(Study, Population) %>%
      summarise(n = n_distinct(marker_type), .groups = "drop") %>%
      # safer fallback: count rows before pivot expansion
      left_join(
        df_long %>%
          group_by(Study, Population) %>%
          summarise(n = n() / length(marker_display_levels), .groups = "drop"),
        by = c("Study", "Population")
      ) %>%
      select(Study, Population, n = n.y)
  }
}

# ----------------------------
# PANEL BUILDER
# ----------------------------
plot_panel <- function(sumdat, ndf, panel_tag, panel_title,
                       show_legend = TRUE, ylim_top = NULL) {
  
  ndf2 <- ndf %>%
    mutate(label = paste0("n=", round(n)),
           y     = -0.06)
  
  if (is.null(ylim_top))
    ylim_top <- max(sumdat$mean + sumdat$se, na.rm = TRUE) * 1.12
  
  p <- ggplot(sumdat, aes(x = Population, y = mean, fill = marker_type)) +
    geom_col(
      position  = position_dodge(width = 0.90),
      width     = 0.82,
      color     = "black",
      linewidth = 0.20
    ) +
    geom_errorbar(
      aes(ymin = mean - se, ymax = mean + se),
      position  = position_dodge(width = 0.90),
      width     = 0.22,
      linewidth = 0.55,
      color     = "black"
    ) +
    geom_text(
      data         = ndf2,
      aes(x = Population, y = y, label = label),
      inherit.aes  = FALSE,
      size         = 5.0,
      fontface     = "bold",
      vjust        = 1
    ) +
    scale_fill_manual(
      values = marker_palette,
      name   = "Dosage class",          # legend title unchanged
      breaks = marker_display_levels    # explicit order in legend
    ) +
    labs(
      title = panel_title,
      x     = "Population",
      y     = "Mean SNP dosage proportion",
      fill  = "Dosage class",
      tag   = panel_tag
    ) +
    theme_classic(base_size = 15) +
    theme(
      plot.title   = element_text(size = 20, face = "bold", hjust = 0.5),
      axis.title   = element_text(size = 17, face = "bold"),
      axis.text    = element_text(size = 15),
      legend.title = element_text(size = 16, face = "bold"),
      legend.text  = element_text(size = 15),
      plot.tag     = element_text(size = 20, face = "bold"),
      plot.margin  = margin(8, 8, 18, 8)
    ) +
    coord_cartesian(ylim = c(-0.08, ylim_top))
  
  if (!show_legend) p <- p + theme(legend.position = "none")
  p
}

# ----------------------------
# RUN
# ----------------------------
s2039_long <- read_prop(file_s2039, "S2039")
s2172_long <- read_prop(file_s2172, "S2172")

s2039_sum <- summarise_prop(s2039_long)
s2172_sum <- summarise_prop(s2172_long)

s2039_n <- n_by_ploidy(s2039_long)
s2172_n <- n_by_ploidy(s2172_long)

common_ylim <- max(
  max(s2039_sum$mean + s2039_sum$se, na.rm = TRUE),
  max(s2172_sum$mean + s2172_sum$se, na.rm = TRUE)
) * 1.12

pA <- plot_panel(
  s2039_sum, s2039_n,
  panel_tag   = "(A)",
  panel_title = "S2039 Diploid vs Tetraploid",
  show_legend = TRUE,
  ylim_top    = common_ylim
)

pB <- plot_panel(
  s2172_sum, s2172_n,
  panel_tag   = "(B)",
  panel_title = "S2172 Diploid vs Tetraploid",
  show_legend = TRUE,
  ylim_top    = common_ylim
)

fig5 <- (pA / pB) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")

# ----------------------------
# SAVE
# ----------------------------
ggsave(
  paste0(out_prefix, ".png"),
  plot = fig5, width = 11, height = 10, dpi = 600, bg = "white"
)
ggsave(
  paste0(out_prefix, ".pdf"),
  plot = fig5, width = 11, height = 10, device = "pdf", bg = "white"
)
ggsave(
  paste0(out_prefix, ".tiff"),
  plot = fig5, width = 11, height = 10, dpi = 600,
  compression = "lzw", bg = "white"
)

writeLines(
  c(
    "Figure 5 regenerated outputs saved here:",
    fig5_dir,
    "",
    "Files:",
    basename(paste0(out_prefix, ".png")),
    basename(paste0(out_prefix, ".pdf")),
    basename(paste0(out_prefix, ".tiff")),
    "",
    "Inputs:",
    basename(file_s2039),
    basename(file_s2172),
    "",
    "Changes from Figure6_regenerated_G3.R:",
    "1. Dosage class labels: zeros→Nulliplex, simplex→Simplex, duplex→Duplex,",
    "   triplex→Triplex, tetraplex→Tetraplex",
    "2. Population axis labels: D→Diploid, T→Tetraploid",
    "3. Palette keys updated to match new capitalized display names",
    "",
    "Notes:",
    "- Bars: mean SNP dosage-class proportions by population and ploidy.",
    "- Error bars: standard error across individuals.",
    "- n under bars: number of individuals per ploidy group.",
    "- Exported PNG and TIFF at 600 dpi; PDF as vector."
  ),
  con = file.path(fig5_dir, "README_Figure5.txt")
)

message("Done. Figure 5 saved to: ", fig5_dir)