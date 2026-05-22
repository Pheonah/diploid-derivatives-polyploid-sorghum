# ============================================================
# FIGURE 3 — S2039 Family Composite (G3 Revision v3)
# Fixes from v2:
#   - Single unified legend (duplicate removed)
#   - Chr6 shading corrected to proper chromosome row in Panel A
#   - Panel C contrast improved (darker orange/blue, lighter gray)
#   - Panel C halepense threshold lowered to 0.50
#   - Panel B chr6 shading alpha reduced to 0.06
#   - Caption updated with strong mechanistic interpretive sentence
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(patchwork)
})

# ----------------------------
# PATHS
# ----------------------------
base_dir  <- "D:/Wenqian_G3_Figures"
geno_file <- file.path(base_dir, "gt_s2039_all.csv")
out_dir   <- file.path(base_dir, "G3_Regenerated_Figures", "Figure3_v3")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
stopifnot(file.exists(geno_file))

# ----------------------------
# TLI COLOR PALETTE
# ----------------------------
TLI_ORANGE    <- "#C96A1B"
TLI_BLUE      <- "#2E6F95"
TLI_GRAY      <- "grey65"
TLI_CREAM     <- "#E9DDC7"
TLI_STEELBLUE <- "steelblue"
MIXED_COLOR   <- "#C0C0C0"   # lighter gray for better contrast

# Panel C uses slightly darker/richer versions for contrast
PC_ORANGE <- "#B85A10"
PC_BLUE   <- "#1E5F85"
PC_MIXED  <- "#D0D0D0"

# Chr6 distortion interval (manuscript text, S2039D)
CHR6_START <- 47.76
CHR6_END   <- 59.34

# ----------------------------
# HELPERS
# ----------------------------
dosage_from_gt <- function(x) {
  x <- as.character(x)
  x[x == "" | x == "NA"] <- NA_character_
  as.integer(ifelse(is.na(x), NA_integer_, str_count(x, "1")))
}

dip_f2_class <- function(dosage) {
  out <- rep(NA_integer_, length(dosage))
  out[dosage == 0] <- 0L
  out[dosage == 2] <- 1L
  out[dosage == 4] <- 2L
  out
}

# ----------------------------
# READ GENOTYPE DATA
# ----------------------------
gt_raw <- read_csv(geno_file, show_col_types = FALSE)
if ("...1" %in% names(gt_raw)) gt_raw <- gt_raw %>% select(-`...1`)
stopifnot(all(c("CHR", "POS") %in% names(gt_raw)))

gt_raw <- gt_raw %>%
  mutate(CHR = as.integer(CHR), POS = as.numeric(POS)) %>%
  filter(!is.na(CHR), !is.na(POS)) %>%
  arrange(CHR, POS)

dip_cols <- grep("^S2039D_", names(gt_raw), value = TRUE)
tet_cols <- grep("^S2039T_", names(gt_raw), value = TRUE)
stopifnot(length(dip_cols) > 0, length(tet_cols) > 0)

# ----------------------------
# SHARED: DIPLOID DOSAGE MATRIX
# ----------------------------
dip_mat <- as.matrix(gt_raw[, dip_cols])
dip_dos <- apply(dip_mat, 2, dosage_from_gt)

# Diploid allele frequency (scale 0-2, divide by 2)
dip_p <- apply(dip_dos, 1, function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA_real_)
  sum(x) / (length(x) * 2)
})

# ----------------------------
# PANEL A — INTROGRESSION MAP (S2039D)
# ----------------------------
run_k_map <- 27

map_df <- gt_raw %>%
  select(CHR, POS) %>%
  mutate(
    p_dip  = dip_p,
    grp    = case_when(
      p_dip < 0.25 ~ "S. bicolor",
      p_dip > 0.75 ~ "S. halepense",
      TRUE         ~ "Mixed"
    ),
    pos_mb = POS / 1e6,
    CHR    = factor(CHR, levels = 1:10)
  )

# S. bicolor-fixed runs
run_segs <- map_df %>%
  filter(!is.na(grp)) %>%
  arrange(CHR, POS) %>%
  group_by(CHR) %>%
  mutate(
    is_bic = grp == "S. bicolor",
    run_id = with(rle(is_bic), rep(seq_along(lengths), lengths))
  ) %>%
  filter(is_bic) %>%
  group_by(CHR, run_id) %>%
  filter(n() >= run_k_map) %>%
  summarise(seg_start = min(pos_mb), seg_end = max(pos_mb), .groups = "drop") %>%
  mutate(CHR = factor(CHR, levels = 1:10))

# Chr6 annotation using annotate() — avoids factor y-axis mapping bug
panel_A <- ggplot(map_df, aes(x = pos_mb, y = CHR)) +
  annotate(
    "rect",
    xmin  = CHR6_START,
    xmax  = CHR6_END,
    ymin  = 5.55,
    ymax  = 6.45,
    fill  = TLI_STEELBLUE,
    alpha = 0.18
  ) +
  geom_point(
    aes(color = grp),
    size  = 0.7,
    alpha = 0.55
  ) +
  geom_segment(
    data = run_segs,
    aes(
      x    = seg_start,
      xend = seg_end,
      y    = CHR,
      yend = CHR
    ),
    color       = TLI_ORANGE,
    linewidth   = 3.0,
    lineend     = "round",
    inherit.aes = FALSE
  ) +
  scale_color_manual(
    values = c(
      "S. bicolor"   = TLI_ORANGE,
      "Mixed"        = MIXED_COLOR,
      "S. halepense" = TLI_BLUE
    ),
    name = "Ancestry"
  ) +
  scale_x_continuous(
    name   = "Physical position (Mb)",
    breaks = seq(0, 80, 10),
    expand = expansion(mult = 0.01)
  ) +
  scale_y_discrete(name = "Chromosome") +
  labs(tag = "(A)", title = "S2039D \u2014 Introgression landscape") +
  guides(color = guide_legend(
    override.aes = list(size = 3, alpha = 1)
  )) +
  theme_classic(base_size = 13) +
  theme(
    plot.tag           = element_text(size = 14, face = "bold"),
    plot.title         = element_text(size = 13, face = "bold", hjust = 0.5),
    axis.title         = element_text(size = 12, face = "bold"),
    axis.text          = element_text(size = 11),
    legend.title       = element_text(size = 11, face = "bold"),
    legend.text        = element_text(size = 10),
    panel.grid.major.x = element_line(color = "grey92", linewidth = 0.3)
  )

# ----------------------------
# PANEL B — SEGREGATION DISTORTION (S2039D)
# ----------------------------
run_k_seg <- 4
alpha_p   <- 0.05
min_n     <- 30

dip_cls <- apply(dip_dos, 2, dip_f2_class)
if (!is.matrix(dip_cls)) dip_cls <- as.matrix(dip_cls)

pvals <- numeric(nrow(gt_raw))
for (i in seq_len(nrow(gt_raw))) {
  v <- dip_cls[i, ]
  v <- v[!is.na(v)]
  if (length(v) < min_n) { pvals[i] <- NA_real_; next }
  tab <- table(factor(v, levels = c(0L, 1L, 2L)))
  exp <- c(0.25, 0.50, 0.25) * sum(tab)
  if (any(exp == 0)) { pvals[i] <- NA_real_; next }
  pvals[i] <- suppressWarnings(
    chisq.test(as.numeric(tab), p = c(0.25, 0.50, 0.25))$p.value
  )
}

seg_df <- gt_raw %>%
  select(CHR, POS) %>%
  mutate(p = pvals, pos_mb = POS / 1e6) %>%
  filter(!is.na(p)) %>%
  arrange(CHR, POS) %>%
  group_by(CHR) %>%
  mutate(
    distorted = p < alpha_p,
    run_id    = with(rle(distorted), rep(seq_along(lengths), lengths)),
    run_len   = ave(distorted, run_id, FUN = length),
    in_run    = distorted & run_len >= run_k_seg
  ) %>%
  ungroup() %>%
  mutate(y = pmin(-log10(p), 20))  # cap at 20

panel_B <- ggplot(seg_df, aes(x = pos_mb, y = y)) +
  geom_rect(
    data = data.frame(
      CHR  = 6,
      xmin = CHR6_START,
      xmax = CHR6_END,
      ymin = 0,
      ymax = 20
    ),
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill        = TLI_STEELBLUE,
    alpha       = 0.06,
    inherit.aes = FALSE
  ) +
  geom_point(
    data  = seg_df %>% filter(!in_run),
    color = TLI_GRAY,
    size  = 0.8,
    alpha = 0.55
  ) +
  geom_point(
    data  = seg_df %>% filter(in_run),
    color = TLI_ORANGE,
    size  = 1.4,
    alpha = 0.90
  ) +
  geom_hline(
    yintercept = -log10(alpha_p),
    linetype   = "dashed",
    color      = "black",
    linewidth  = 0.7
  ) +
  facet_grid(. ~ CHR, scales = "free_x", space = "free_x") +
  scale_x_continuous(
    name   = "Physical position (Mb)",
    breaks = seq(0, 80, 20),
    expand = expansion(mult = 0.02)
  ) +
  scale_y_continuous(
    name   = expression(-log[10](italic(P))),
    limits = c(0, 20),
    breaks = seq(0, 20, 5)
  ) +
  labs(tag = "(B)", title = "S2039D \u2014 Segregation distortion") +
  guides(color = "none") +   # suppress Panel B legend — unified via patchwork
  theme_classic(base_size = 13) +
  theme(
    plot.tag           = element_text(size = 14, face = "bold"),
    plot.title         = element_text(size = 13, face = "bold", hjust = 0.5),
    axis.title         = element_text(size = 12, face = "bold"),
    axis.text.x        = element_text(size = 9, angle = 45, hjust = 1),
    axis.text.y        = element_text(size = 11),
    strip.text         = element_text(size = 11, face = "bold"),
    strip.background   = element_rect(fill = TLI_CREAM, color = "grey60"),
    panel.spacing.x    = unit(0.35, "lines"),
    panel.grid.major.x = element_line(color = "grey92", linewidth = 0.3)
  )

# ----------------------------
# PANEL C — T254 TETRAPLOID ANCESTRY (1 Mb windows)
# ----------------------------
tet_mat <- as.matrix(gt_raw[, tet_cols])
tet_dos <- apply(tet_mat, 2, dosage_from_gt)

tet_p <- apply(tet_dos, 1, function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA_real_)
  sum(x) / (length(x) * 4)
})

tet_raw_df <- gt_raw %>%
  select(CHR, POS) %>%
  mutate(p_tet = tet_p, pos_mb = POS / 1e6) %>%
  filter(!is.na(p_tet))

window_mb <- 1.0

tet_smooth <- tet_raw_df %>%
  group_by(CHR) %>%
  mutate(window = floor(pos_mb / window_mb) * window_mb + window_mb / 2) %>%
  group_by(CHR, window) %>%
  summarise(
    p_mean   = mean(p_tet, na.rm = TRUE),
    n_marker = n(),
    .groups  = "drop"
  ) %>%
  filter(n_marker >= 3) %>%
  mutate(
    grp = case_when(
      p_mean < 0.25 ~ "S. bicolor",
      p_mean > 0.50 ~ "S. halepense",   # lowered from 0.75 for better contrast
      TRUE          ~ "Mixed"
    ),
    CHR = factor(CHR, levels = 1:10)
  )

panel_C <- ggplot(tet_smooth, aes(x = window, y = CHR, color = grp)) +
  geom_point(size = 1.8, alpha = 0.88) +
  scale_color_manual(
    values = c(
      "S. bicolor"   = PC_ORANGE,
      "Mixed"        = PC_MIXED,
      "S. halepense" = PC_BLUE
    ),
    name = "Ancestry"
  ) +
  scale_x_continuous(
    name   = "Physical position (Mb)",
    breaks = seq(0, 80, 10),
    expand = expansion(mult = 0.01)
  ) +
  scale_y_discrete(name = "Chromosome") +
  labs(
    tag   = "(C)",
    title = "T254 \u2014 Tetraploid ancestral background (1 Mb windows)"
  ) +
  guides(color = "none") +   # unified via patchwork collect
  theme_classic(base_size = 13) +
  theme(
    plot.tag           = element_text(size = 14, face = "bold"),
    plot.title         = element_text(size = 13, face = "bold", hjust = 0.5),
    axis.title         = element_text(size = 12, face = "bold"),
    axis.text          = element_text(size = 11),
    panel.grid.major.x = element_line(color = "grey92", linewidth = 0.3)
  )

# ----------------------------
# COMBINE — SINGLE UNIFIED LEGEND
# ----------------------------
fig3 <- (panel_A / panel_B / panel_C) +
  plot_layout(
    heights = c(1.2, 1.2, 1.0),
    guides  = "collect"
  ) +
  plot_annotation(
    title    = "Figure 3. Genomic composition and segregation distortion \u2014 S2039 family",
    subtitle = paste0(
      "Orange segments (A): S. bicolor-enriched intervals (runs \u2265 27 markers). ",
      "Blue shading: chr6 segregation distortion region (47.76\u201359.34 Mb). ",
      "Orange points (B): runs of \u2265 4 consecutive distorted markers (P\u00a0<\u00a00.05, \u03c7\u00b2 test, 1:2:1 expectation); ",
      "dashed line: P\u00a0=\u00a00.05 threshold; y-axis capped at 20 (values above 20 truncated). ",
      "Panel C: T254 tetraploid ancestry aggregated in 1 Mb windows (halepense threshold p\u00a0>\u00a00.50). ",
      "Overlap between S. bicolor-enriched diploid intervals and segregation distortion peaks, ",
      "particularly on chr6, suggests selective retention or purging of introgressed genomic ",
      "regions during diploid extraction from the tetraploid background."
    ),
    theme = theme(
      plot.title    = element_text(size = 13, face = "bold",  hjust = 0.5),
      plot.subtitle = element_text(size = 9,  hjust = 0.5,
                                   color  = "grey25",
                                   margin = margin(b = 10))
    )
  ) &
  theme(
    legend.position = "right",
    legend.title    = element_text(size = 11, face = "bold"),
    legend.text     = element_text(size = 10)
  )

# ----------------------------
# SAVE ALL FORMATS
# ----------------------------
ggsave(
  file.path(out_dir, "Figure3_S2039_composite_v3.png"),
  plot = fig3, width = 14, height = 16, dpi = 600, bg = "white"
)
ggsave(
  file.path(out_dir, "Figure3_S2039_composite_v3.pdf"),
  plot = fig3, width = 14, height = 16, device = "pdf", bg = "white"
)
ggsave(
  file.path(out_dir, "Figure3_S2039_composite_v3.tiff"),
  plot = fig3, width = 14, height = 16, dpi = 600,
  compression = "lzw", bg = "white"
)

message("Done. Figure 3 v3 saved to: ", out_dir)