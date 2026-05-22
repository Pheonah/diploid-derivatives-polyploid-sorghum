# ============================================================
# Figure 7: Genomic distance matrices and MDS analyses
# S2039D/S2039T and S2172D/S2172T
# Journal-ready version for G3 resubmission
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

# -------------------------
# Input genotype files
# -------------------------

gt_s2039 <- "E:/VINCENT Momocs/RQTL_MAPPING/01_DATA/S2039/genotype/gt_s2039_all.csv"
gt_s2172 <- "E:/VINCENT Momocs/RQTL_MAPPING/01_DATA/S2172/genotype/gt_s2172_all.csv"

# -------------------------
# Output directory
# -------------------------

out_dir <- "E:/VINCENT Momocs/Wenqian_Manuscript/G3_Regenerated_Figures/Figure7_MDS"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

stopifnot(file.exists(gt_s2039), file.exists(gt_s2172))

# -------------------------
# TLI colors — consistent with Figures 5 & 6
# -------------------------

TLI_BROWN <- "#5C3A21"
TLI_GREEN <- "#5E8C31"
TLI_DARK  <- "#222222"

ploidy_palette <- c(
  Diploid    = TLI_BROWN,
  Tetraploid = TLI_GREEN
)

distance_palette <- c(
  "#F3EAD7",
  "#C7A77B",
  "#8A5A2B",
  "#4A2A12"
)

# -------------------------
# Convert genotype calls to dosage
# -------------------------

dosage_from_gt <- function(x) {
  x <- as.character(x)
  x[x == "" | x == "NA" | x == "." | x == "./." | x == ".|."] <- NA_character_
  
  out <- rep(NA_real_, length(x))
  
  is_num <- grepl("^\\d+(\\.0)?$", x)
  out[is_num] <- as.numeric(x[is_num])
  
  is_gt <- grepl("[/|]", x)
  out[is_gt] <- stringr::str_count(x[is_gt], "1")
  
  out
}

# -------------------------
# Read marker-by-individual genotype file
# -------------------------

read_gt_marker_by_individual <- function(path, family_prefix) {
  
  df <- readr::read_csv(path, show_col_types = FALSE)
  
  if ("...1" %in% names(df)) df <- df %>% select(-`...1`)
  
  id_cols <- names(df)[grepl(paste0("^", family_prefix, "[DT]_"), names(df))]
  
  if (length(id_cols) == 0)
    stop("No individual genotype columns found for ", family_prefix)
  
  G <- df %>%
    select(all_of(id_cols)) %>%
    mutate(across(everything(), dosage_from_gt)) %>%
    as.matrix()
  
  keep_markers <- rowSums(!is.na(G)) > 0
  G <- G[keep_markers, , drop = FALSE]
  
  X <- t(G)
  rownames(X) <- colnames(G)
  colnames(X) <- paste0("M", seq_len(ncol(X)))
  
  keep_inds <- rowSums(!is.na(X)) > 0
  X <- X[keep_inds, , drop = FALSE]
  
  metadata <- data.frame(
    ind = rownames(X),
    pop = ifelse(
      grepl(paste0("^", family_prefix, "D_"), rownames(X)),
      "Diploid", "Tetraploid"
    ),
    stringsAsFactors = FALSE
  )
  metadata$pop <- factor(metadata$pop, levels = c("Diploid", "Tetraploid"))
  
  message(family_prefix, ": retained ", nrow(X),
          " individuals and ", ncol(X), " markers after filtering.")
  
  list(X = X, metadata = metadata)
}

# -------------------------
# Modified Hamming distance
# -------------------------

mod_hamming <- function(X) {
  
  n <- nrow(X)
  D <- matrix(NA_real_, n, n)
  diag(D) <- 0
  rownames(D) <- colnames(D) <- rownames(X)
  
  for (i in seq_len(n)) {
    xi <- X[i, ]
    for (j in i:n) {
      xj <- X[j, ]
      ok <- !(is.na(xi) | is.na(xj))
      denom <- sum(ok)
      dij <- if (denom == 0) NA_real_ else sum(xi[ok] != xj[ok]) / denom
      D[i, j] <- dij
      D[j, i] <- dij
    }
  }
  D
}

# -------------------------
# Order individuals by ploidy, then cluster within ploidy
# -------------------------

order_by_pop_then_cluster <- function(D, pop_vec) {
  
  idxD <- which(pop_vec == "Diploid")
  idxT <- which(pop_vec == "Tetraploid")
  
  order_block <- function(idx) {
    if (length(idx) < 3) return(idx)
    DD <- D[idx, idx, drop = FALSE]
    finite_vals <- DD[is.finite(DD)]
    if (length(finite_vals) == 0) return(idx)
    mx <- max(finite_vals)
    DD[!is.finite(DD)] <- mx
    diag(DD) <- 0
    hc <- hclust(as.dist(DD), method = "average")
    idx[hc$order]
  }
  
  c(order_block(idxD), order_block(idxT))
}

# -------------------------
# Compute shared distance range across both families
# for a unified heatmap color scale
# -------------------------

# (called after both D matrices are computed — see below)

# -------------------------
# Distance heatmap — accepts global scale limits
# -------------------------

plot_dist_heatmap <- function(D, ind_order, pop_vec_ordered,
                              title, panel_tag,
                              dist_limits = NULL) {
  
  D2  <- D[ind_order, ind_order, drop = FALSE]
  df  <- expand.grid(i = seq_len(nrow(D2)), j = seq_len(ncol(D2)))
  df$val <- as.vector(D2)
  nD  <- sum(pop_vec_ordered == "Diploid")
  
  # Use global limits when provided so both heatmaps share the same scale
  if (is.null(dist_limits)) {
    dist_limits <- range(df$val, na.rm = TRUE)
  }
  
  ggplot(df, aes(x = i, y = j, fill = val)) +
    geom_raster() +
    geom_hline(yintercept = nD + 0.5, linewidth = 0.45, color = TLI_DARK) +
    geom_vline(xintercept = nD + 0.5, linewidth = 0.45, color = TLI_DARK) +
    scale_y_reverse(expand = c(0, 0)) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_fill_gradientn(
      colours  = distance_palette,
      name     = "Distance",
      limits   = dist_limits,          # FIX: shared scale across panels
      labels   = label_number(accuracy = 0.01),
      na.value = "white"
    ) +
    labs(title = title, x = NULL, y = NULL, tag = panel_tag) +
    theme_classic(base_size = 14) +
    theme(
      plot.title        = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.tag          = element_text(size = 18, face = "bold"),
      axis.text         = element_blank(),
      axis.ticks        = element_blank(),
      axis.line         = element_blank(),
      legend.title      = element_text(size = 12, face = "bold"),
      legend.text       = element_text(size = 11),
      legend.key.height = unit(0.7, "cm"),
      plot.margin       = margin(4, 8, 4, 8)
    )
}

# -------------------------
# MDS plot
# -------------------------

plot_mds <- function(D, inds, pop_vec, title, panel_tag) {
  
  D_use <- D
  finite_vals <- D_use[is.finite(D_use)]
  if (length(finite_vals) == 0) stop("Distance matrix contains no finite values.")
  mx <- max(finite_vals)
  D_use[!is.finite(D_use)] <- mx
  diag(D_use) <- 0
  
  m   <- cmdscale(as.dist(D_use), k = 2, eig = TRUE)
  pts <- as.data.frame(m$points)
  colnames(pts) <- c("MDS1", "MDS2")
  pts$ind <- inds
  pts$pop <- factor(pop_vec, levels = c("Diploid", "Tetraploid"))
  
  ggplot(pts, aes(MDS1, MDS2, color = pop)) +
    geom_point(size = 2.4, alpha = 0.85) +
    scale_color_manual(values = ploidy_palette, name = "Ploidy") +
    labs(title = title, x = "MDS axis 1", y = "MDS axis 2", tag = panel_tag) +
    theme_classic(base_size = 14) +
    theme(
      plot.title   = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.tag     = element_text(size = 18, face = "bold"),
      axis.title   = element_text(size = 14, face = "bold"),
      axis.text    = element_text(size = 12),
      legend.title = element_text(size = 12, face = "bold"),
      legend.text  = element_text(size = 11),
      plot.margin  = margin(4, 8, 4, 8)
    )
}

# -------------------------
# Generate panels for one family
# -------------------------

make_family_panels <- function(gt_path, family_prefix, family_label,
                               tags, dist_limits = NULL) {
  
  gt       <- read_gt_marker_by_individual(gt_path, family_prefix)
  X        <- gt$X
  metadata <- gt$metadata
  
  pop_vec        <- as.character(metadata$pop)
  names(pop_vec) <- rownames(X)
  
  D   <- mod_hamming(X)
  ord <- order_by_pop_then_cluster(D, pop_vec)
  pop_ord <- pop_vec[ord]
  
  p_heat <- plot_dist_heatmap(
    D, ord, pop_ord,
    title       = paste0(family_label, " distance matrix"),
    panel_tag   = tags[1],
    dist_limits = dist_limits          # pass shared scale
  )
  
  p_mds <- plot_mds(
    D,
    inds    = names(pop_vec),
    pop_vec = pop_vec,
    title   = paste0(family_label, " MDS"),
    panel_tag = tags[2]
  )
  
  list(
    heatmap         = p_heat,
    mds             = p_mds,
    distance_matrix = D,
    metadata        = metadata
  )
}

# -------------------------
# Compute distance matrices first to get shared scale
# -------------------------

gt_data_s2039 <- read_gt_marker_by_individual(gt_s2039, "S2039")
gt_data_s2172 <- read_gt_marker_by_individual(gt_s2172, "S2172")

D_s2039 <- mod_hamming(gt_data_s2039$X)
D_s2172 <- mod_hamming(gt_data_s2172$X)

# Global distance limits — both heatmaps use identical color scale
global_dist_limits <- range(
  c(D_s2039[!is.na(D_s2039)], D_s2172[!is.na(D_s2172)])
)

# -------------------------
# Build panels with shared scale
# -------------------------

s2039 <- make_family_panels(
  gt_s2039, "S2039", "S2039",
  tags        = c("(A)", "(B)"),
  dist_limits = global_dist_limits
)

s2172 <- make_family_panels(
  gt_s2172, "S2172", "S2172",
  tags        = c("(C)", "(D)"),
  dist_limits = global_dist_limits
)

# -------------------------
# FIX: 2x2 grid with guides="collect" at TOP level
# This merges duplicate Distance and Ploidy legends into one each
# -------------------------

row1 <- (s2039$heatmap | s2039$mds) +
  plot_layout(widths = c(1.1, 0.9))

row2 <- (s2172$heatmap | s2172$mds) +
  plot_layout(widths = c(1.1, 0.9))

fig7 <- (row1 / row2) +
  plot_layout(
    heights = c(1, 1),
    guides  = "collect"               # FIX: collect at top level — one Distance + one Ploidy legend
  ) &
  theme(
    legend.position = "right",
    plot.margin     = margin(6, 6, 6, 6)
  )

# -------------------------
# Export
# -------------------------

png_file  <- file.path(out_dir, "Figure7_Distance_MDS.png")
pdf_file  <- file.path(out_dir, "Figure7_Distance_MDS.pdf")
tiff_file <- file.path(out_dir, "Figure7_Distance_MDS.tiff")

ggsave(png_file,  fig7, width = 14, height = 12, dpi = 600, bg = "white")
ggsave(pdf_file,  fig7, width = 14, height = 12, device = "pdf", bg = "white")
ggsave(tiff_file, fig7, width = 14, height = 12, dpi = 600,
       compression = "lzw", bg = "white")

# -------------------------
# Save distance matrices and metadata
# -------------------------

write.csv(s2039$distance_matrix,
          file.path(out_dir, "S2039_modified_Hamming_distance_matrix.csv"),
          row.names = TRUE)
write.csv(s2172$distance_matrix,
          file.path(out_dir, "S2172_modified_Hamming_distance_matrix.csv"),
          row.names = TRUE)

readr::write_csv(s2039$metadata, file.path(out_dir, "S2039_MDS_metadata.csv"))
readr::write_csv(s2172$metadata, file.path(out_dir, "S2172_MDS_metadata.csv"))

writeLines(
  c(
    "Figure 7 regenerated outputs saved here:", out_dir, "",
    "Files:",
    basename(png_file), basename(pdf_file), basename(tiff_file), "",
    "Inputs:", basename(gt_s2039), basename(gt_s2172), "",
    "Notes:",
    "- Distance matrices computed using modified Hamming distance.",
    "- Both heatmaps share a single unified color scale (global_dist_limits).",
    "- Individuals ordered by ploidy class then clustered within ploidy.",
    "- 2x2 grid: heatmap | MDS per family row.",
    "- Single Distance legend and single Ploidy legend via guides='collect'.",
    "- TLI palette: Diploid = #5C3A21, Tetraploid = #5E8C31.",
    "- Exported PNG and TIFF at 600 dpi; PDF as vector backup."
  ),
  con = file.path(out_dir, "README_Figure7_MDS.txt")
)

message("DONE. Journal-ready Figure 7 saved to: ", out_dir)