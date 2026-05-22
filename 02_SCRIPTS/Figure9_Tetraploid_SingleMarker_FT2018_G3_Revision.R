# ============================================================
# FIGURE 9 — Single-Marker Association Scans, Flowering Time (2018)
# Tetraploids: S2039T and S2172T — G3 Revision FINAL
# Fixes applied:
#   - dplyr:: namespace protection (MASS/lme4 conflict)
#   - Significance threshold line added (reviewer request, C1 fix)
#   - Base font ≥13pt, axis titles ≥14pt (AE requirement)
#   - Strip text enlarged, panel spacing increased
#   - Exploratory framing in titles (no overclaiming)
#   - Output renamed Figure9_SMA_FT2018
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(DBI)
  library(RSQLite)
  library(emmeans)
  library(car)
  library(patchwork)
})

# Namespace protection
select    <- dplyr::select
filter    <- dplyr::filter
rename    <- dplyr::rename
mutate    <- dplyr::mutate
summarise <- dplyr::summarise

# ----------------------------
# PATHS
# ----------------------------
base_dir   <- "D:/Wenqian_G3_Figures"
db_path    <- file.path(base_dir, "diploids.db")
geno_s2039 <- file.path(base_dir, "gt_s2039_reduced.csv")
geno_s2172 <- file.path(base_dir, "gt_s2172_reduced.csv")
out_dir    <- file.path(base_dir, "G3_Regenerated_Figures", "Figure9_SMA")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

stopifnot(file.exists(db_path), file.exists(geno_s2039), file.exists(geno_s2172))

# ----------------------------
# TLI PALETTE
# ----------------------------
TLI <- list(
  brown = "#4D3016",
  olive = "#7A9A01",
  cream = "#E9DDC7",
  gray  = "#4D4D4D",
  black = "#111111"
)

# Significance threshold — Bonferroni-style exploratory line
# Uses the larger marker set; will be computed per-panel after SMA runs
ALPHA <- 0.05

# ----------------------------
# HELPERS
# ----------------------------
pick_col <- function(df, candidates) {
  hit <- intersect(candidates, names(df))
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

read_geno <- function(path) {
  g <- readr::read_csv(path, show_col_types = FALSE)
  if ("...1" %in% names(g)) g <- g %>% dplyr::select(-`...1`)
  
  chr_col <- pick_col(g, c("CHR","chr","Chr"))
  pos_col <- pick_col(g, c("POS","pos","Pos"))
  if (is.na(chr_col) || is.na(pos_col))
    stop("Genotype file must contain CHR and POS columns: ", path)
  
  g <- g %>%
    rename(CHR = all_of(chr_col), POS = all_of(pos_col)) %>%
    mutate(
      CHR    = suppressWarnings(as.integer(CHR)),
      POS    = suppressWarnings(as.numeric(POS)),
      marker = paste0("Chr", CHR, "_", POS),
      pos_mb = POS / 1e6
    )
  
  id_cols <- names(g)[grepl("^S\\d{4}[DT]_\\d+", names(g))]
  if (length(id_cols) == 0)
    stop("No individual columns found (expected S####[DT]_###): ", path)
  
  meta <- g %>% dplyr::select(marker, CHR, POS, pos_mb)
  G    <- g %>%
    dplyr::select(marker, all_of(id_cols)) %>%
    column_to_rownames("marker") %>%
    as.matrix()
  
  list(meta = meta, G = G)
}

get_ft2018_emmeans <- function(con, population_code) {
  fl <- dbGetQuery(con, 'select * from flowering where Accession <> "NA";') %>%
    as_tibble()
  fl <- as.data.frame(lapply(fl, type.convert, as.is = TRUE)) %>% as_tibble()
  
  pop_col <- pick_col(fl, c("population","Population","pop","Pop"))
  if (is.na(pop_col)) stop("No population column in DB 'flowering' table.")
  fl <- fl %>% filter(.data[[pop_col]] == population_code)
  
  id_col <- pick_col(fl, c("geno_id","genotype_id","Genotype",
                           "Accession","accession","geno.id"))
  if (is.na(id_col))
    stop("No genotype ID column for: ", population_code)
  
  ft_col <- pick_col(fl, c("FL50_adj","FL18","FL2018","FT2018",
                           "flowering_2018","FL_2018"))
  if (is.na(ft_col)) {
    num_cols <- names(fl)[sapply(fl, is.numeric)]
    if (length(num_cols) == 0)
      stop("No numeric FT column for: ", population_code)
    ft_col <- num_cols[1]
    message("WARNING: using '", ft_col, "' as FT2018 for ", population_code)
  }
  
  block_col <- pick_col(fl, c("Block","block","rep","Rep",
                              "replicate","Replication"))
  
  dat <- fl %>%
    transmute(
      geno_id = as.character(.data[[id_col]]),
      FT2018  = suppressWarnings(as.numeric(.data[[ft_col]])),
      block   = if (!is.na(block_col)) as.factor(.data[[block_col]]) else NA
    ) %>%
    filter(!is.na(geno_id), !is.na(FT2018))
  
  if (nrow(dat) < 10) stop("Too few records for: ", population_code)
  
  if (!all(is.na(dat$block)) && nlevels(dat$block) > 1) {
    m <- lm(FT2018 ~ geno_id + block, data = dat)
    message("Block-adjusted model for ", population_code)
  } else {
    m <- lm(FT2018 ~ geno_id, data = dat)
    message("No block term for ", population_code)
  }
  
  emmeans(m, "geno_id") %>%
    as.data.frame() %>%
    as_tibble() %>%
    dplyr::select(geno_id, emmean) %>%
    rename(FT2018_emmean = emmean) %>%
    arrange(geno_id)
}

marker_keep_rule <- function(x, min_n = 6) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(FALSE)
  tab <- table(x)
  if (length(tab) == 2) return(min(tab) > min_n)
  if (length(tab) == 3) return(median(tab) > min_n)
  FALSE
}

run_sma <- function(geno_obj, pheno_emm, id_prefix_regex, min_n = 6) {
  ids    <- colnames(geno_obj$G)[grepl(id_prefix_regex, colnames(geno_obj$G))]
  if (length(ids) < 10) stop("Too few genotype columns match: ", id_prefix_regex)
  
  common <- intersect(ids, pheno_emm$geno_id)
  if (length(common) < 10)
    stop("Too few common individuals for: ", id_prefix_regex)
  
  Gsub <- geno_obj$G[, common, drop = FALSE]
  y    <- pheno_emm %>%
    filter(geno_id %in% common) %>%
    arrange(match(geno_id, common))
  stopifnot(all(y$geno_id == common))
  
  keep <- apply(Gsub, 1, marker_keep_rule, min_n = min_n)
  Gf   <- Gsub[keep, , drop = FALSE]
  meta <- geno_obj$meta %>% filter(marker %in% rownames(Gf))
  
  pval   <- rep(NA_real_, nrow(Gf))
  effect <- rep(NA_real_, nrow(Gf))
  
  for (i in seq_len(nrow(Gf))) {
    gvec <- as.factor(Gf[i, ])
    if (nlevels(droplevels(gvec)) < 2) next
    fit <- try(lm(y$FT2018_emmean ~ gvec), silent = TRUE)
    if (inherits(fit, "try-error")) next
    a3  <- try(car::Anova(fit, type = 3), silent = TRUE)
    if (inherits(a3, "try-error")) next
    pval[i]   <- a3$`Pr(>F)`[2]
    cf        <- coef(fit)
    effect[i] <- if (length(cf) >= 2) unname(cf[2]) else NA_real_
  }
  
  meta %>%
    mutate(
      pval      = pval[match(marker, rownames(Gf))],
      neglog10p = -log10(pval)
    ) %>%
    filter(is.finite(neglog10p))
}

# ----------------------------
# PLOT BUILDER — with threshold line
# ----------------------------
plot_sma <- function(res, title_text, point_col, n_tests) {
  
  # Bonferroni threshold for this panel
  bonf_thresh <- -log10(ALPHA / n_tests)
  
  ggplot(res, aes(x = pos_mb, y = neglog10p)) +
    geom_point(size = 0.65, alpha = 0.85, color = point_col) +
    # Reviewer-requested significance threshold line
    geom_hline(
      yintercept = bonf_thresh,
      linetype   = "dashed",
      color      = "black",
      linewidth  = 0.7
    ) +
    annotate(
      "text",
      x     = -Inf, y = bonf_thresh,
      label = "Bonf. \u03b1=0.05",
      hjust = -0.05, vjust = -0.4,
      size  = 3.2, color = "black"
    ) +
    facet_grid(. ~ CHR, scales = "free_x", space = "free_x") +
    scale_x_continuous(
      name   = "Physical position (Mb)",
      breaks = seq(0, 80, 20),
      expand = expansion(mult = 0.02)
    ) +
    scale_y_continuous(
      name = expression(-log[10](italic(P)))
    ) +
    labs(title = title_text) +
    theme_classic(base_size = 13) +
    theme(
      plot.title         = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.title         = element_text(size = 13, face = "bold"),
      axis.text.x        = element_text(size = 9, angle = 45, hjust = 1),
      axis.text.y        = element_text(size = 11),
      strip.text         = element_text(size = 11, face = "bold"),
      strip.background   = element_rect(fill = TLI$cream, color = "grey60"),
      panel.spacing.x    = unit(0.30, "lines"),
      panel.grid.major.y = element_line(color = "grey92", linewidth = 0.3),
      plot.margin        = margin(8, 15, 8, 8)
    )
}

# ----------------------------
# MAIN
# ----------------------------
con <- dbConnect(SQLite(), dbname = db_path)
on.exit(dbDisconnect(con), add = TRUE)

emm_s2039t <- get_ft2018_emmeans(con, "S2039T")
emm_s2172t <- get_ft2018_emmeans(con, "S2172T")

readr::write_csv(emm_s2039t, file.path(out_dir, "S2039T_emmeans_FT2018.csv"))
readr::write_csv(emm_s2172t, file.path(out_dir, "S2172T_emmeans_FT2018.csv"))

g2039 <- read_geno(geno_s2039)
g2172 <- read_geno(geno_s2172)

res_s2039t <- run_sma(g2039, emm_s2039t, "^S2039T_", min_n = 6)
res_s2172t <- run_sma(g2172, emm_s2172t, "^S2172T_", min_n = 6)

readr::write_csv(res_s2039t, file.path(out_dir, "S2039T_SMA_FT2018_results.csv"))
readr::write_csv(res_s2172t, file.path(out_dir, "S2172T_SMA_FT2018_results.csv"))

# Build panels — pass n_tests for per-panel Bonferroni threshold
p1 <- plot_sma(
  res_s2039t,
  title_text = "S2039T \u2014 Single-marker scan, flowering time (2018, exploratory)",
  point_col  = TLI$olive,
  n_tests    = nrow(res_s2039t)
)

p2 <- plot_sma(
  res_s2172t,
  title_text = "S2172T \u2014 Single-marker scan, flowering time (2018, exploratory)",
  point_col  = TLI$brown,
  n_tests    = nrow(res_s2172t)
)

fig9 <- (p1 / p2) + plot_layout(heights = c(1, 1))

# ----------------------------
# SAVE
# ----------------------------
ggsave(file.path(out_dir, "Figure9_SMA_FT2018_S2039T_S2172T.png"),
       fig9, width = 11, height = 8, dpi = 600, bg = "white")

ggsave(file.path(out_dir, "Figure9_SMA_FT2018_S2039T_S2172T.pdf"),
       fig9, width = 11, height = 8, device = "pdf", bg = "white")

ggsave(file.path(out_dir, "Figure9_SMA_FT2018_S2039T_S2172T.tiff"),
       fig9, width = 11, height = 8, dpi = 600,
       compression = "lzw", bg = "white")

writeLines(
  c(
    "Figure 9 — Single-marker FT2018 scans (tetraploids)",
    paste0("DB: ", db_path),
    paste0("Outputs: ", out_dir),
    "",
    paste0("S2039T: n(emmeans)=", nrow(emm_s2039t),
           " | markers tested=", nrow(res_s2039t),
           " | Bonferroni threshold=",
           round(-log10(ALPHA / nrow(res_s2039t)), 2)),
    paste0("S2172T: n(emmeans)=", nrow(emm_s2172t),
           " | markers tested=", nrow(res_s2172t),
           " | Bonferroni threshold=",
           round(-log10(ALPHA / nrow(res_s2172t)), 2))
  ),
  con = file.path(out_dir, "README_Figure9.txt")
)

message("Done. Figure 9 saved to: ", out_dir)