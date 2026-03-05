## ============================================================
## PACKAGES (shared)
## ============================================================
pkgs_cran <- c("Peptides","dplyr","ggplot2","ggrepel","stringr","readr","scales")
pkgs_bioc <- c("Biostrings")

to_install_cran <- setdiff(pkgs_cran, rownames(installed.packages()))
if (length(to_install_cran)) install.packages(to_install_cran)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
to_install_bioc <- pkgs_bioc[!sapply(pkgs_bioc, requireNamespace, quietly = TRUE)]
if (length(to_install_bioc)) BiocManager::install(to_install_bioc, ask = FALSE, update = FALSE)

## ============================================================
## SETTINGS
## ============================================================

invisible(lapply(c(pkgs_cran, pkgs_bioc), library, character.only = TRUE))

setwd("/path/to/directory")
fasta_file <- "BTN1A1_reconstructed.fasta"
mf_file    <- "/path/to/MF_percent.csv"
out_dir    <- "/path/to/outdir"
gene_name  <- "BTN1A1"

# Extracellular region on the *gapped* alignment (1-based, inclusive)
#extra_cols <- c(72, 292) # for BTN1A1; change as needed

# Net charge settings
pH <- 6.7
pK <- "Lehninger"

## ============================================================
## INPUTS (shared)
## ============================================================
# FASTA / alignment
aas <- Biostrings::readAAStringSet(fasta_file)
ids <- names(aas)
seqs_raw <- as.character(aas)

# Milk-fat table (expects at least: id or species; MF_percent; optional group)
mf <- read.csv(mf_file, stringsAsFactors = FALSE) %>%
  mutate(
    id         = if ("id" %in% names(.))        as.character(id)        else NA_character_,
    species    = if ("species" %in% names(.))   as.character(species)   else NA_character_,
    group      = if ("group" %in% names(.))     as.character(group)     else NA_character_,
    MF_percent = suppressWarnings(as.numeric(MF_percent))
  )

## ============================================================
## HELPERS
## ============================================================
# Clean sequence: keep only 20 standard AAs, remove gaps/stops
std <- strsplit("ACDEFGHIKLMNPQRSTVWY","")[[1]]
clean_seq <- function(s) {
  s <- toupper(s)
  s <- gsub("[\\-\\.\\*]", "", s)
  aa <- strsplit(s, "")[[1]]
  paste0(aa[aa %in% std], collapse = "")
}

# KD scale + mean
KD <- c(A=1.8,C=2.5,D=-3.5,E=-3.5,F=2.8,G=-0.4,H=-3.2,I=4.5,
        K=-3.9,L=3.8,M=1.9,N=-3.5,P=-1.6,Q=-3.5,R=-4.5,S=-0.8,
        T=-0.7,V=4.2,W=-0.9,Y=-1.3)
mean_KD <- function(s) {
  if (is.null(s) || nchar(s) == 0) return(NA_real_)
  aa <- strsplit(s, "")[[1]]
  mean(as.numeric(KD[aa]), na.rm = TRUE)
}

# Slice gapped alignment by columns, then clean to AA-only
get_aln_segment <- function(s_raw, cols) {
  L <- nchar(s_raw)
  start <- max(1L, cols[1]); end <- min(L, cols[2])
  if (start > end) return("")
  sub_raw <- substr(s_raw, start, end)
  clean_seq(sub_raw)
}

# Species extractor from FASTA id (before pipe or first token)
extract_species <- function(x) {
  dplyr::coalesce(
    stringr::str_extract(x, "^[^|]+"),
    stringr::str_extract(x, "^[^\\s]+")
  )
}

# Simple underscore normalization for species-joins
norm_sp <- function(x) gsub(" ", "_", x)

# High-contrast, colorblind-safe palette
#  Tl-bright color palette, consistent across all proteins
TOL_BRIGHT <- c(
  "Afrotherians" = "#1477AA",  # blue
  "Artiodactyls" = "#FE6677",  # coral
  "Carnivores"   = "#428833",  # green
  "Cetaceans"    = "#CCBB44",  # mustard
  "Marsupials"   = "#66CCEE",  # light blue
  "Pinnipeds"    = "#CA3377",  # magenta
  "Primates"     = "#9071A1",  # mauve
  "Rodents"      = "#6072B2",  # steel blue
  "Monotremes"   = "#333333"  # dark grey
)
TOL_FILL <- scales::alpha(TOL_BRIGHT, 0.18)

# Decide join key: prefer 'id' if mf has it; else 'species'
join_with_mf <- function(df) {
  if ("id" %in% names(mf) && any(!is.na(mf$id))) {
    out <- df %>%
      dplyr::left_join(mf %>% dplyr::select(id, MF_percent, group, species), by = "id") %>%
      dplyr::mutate(species = dplyr::coalesce(species.x, species.y)) %>%
      dplyr::select(-species.x, -species.y)
  } else {
    out <- df %>%
      dplyr::mutate(species = norm_sp(species)) %>%
      dplyr::left_join(mf %>% dplyr::mutate(species = norm_sp(species)) %>%
                         dplyr::select(species, MF_percent, group),
                       by = "species")
  }
  out
}

collapse_to_species <- function(df, cols_mean) {
  df %>%
    dplyr::group_by(species) %>%
    dplyr::summarise(
      MF_percent = dplyr::first(MF_percent),
      dplyr::across(dplyr::all_of(cols_mean), ~ mean(.x, na.rm = TRUE), .names = "{.col}"),
      group = dplyr::first(group),
      n_seqs = dplyr::n(),
      .groups = "drop"
    ) %>%
    dplyr::filter(!is.na(MF_percent))
}

# Plotting with convex hulls per group
plot_xy_hulls <- function(
    data, yvar, ylab,
    gene_label      = NULL,
    gene_size       = 5,   # new argument
    color_legend    = "Group",
    add_stats       = TRUE,
    y_expand_frac   = 0.06,
    gene_pos        = c("topleft","topright","bottomleft","bottomright"),
    stats_pos       = c("bottomright","topleft","topright","bottomleft")
) {
  gene_pos  <- match.arg(gene_pos)
  stats_pos <- match.arg(stats_pos)
  
  stopifnot(yvar %in% names(data))
  d <- data %>%
    dplyr::filter(!is.na(MF_percent), !is.na(.data[[yvar]])) %>%
    dplyr::mutate(y = .data[[yvar]],
                  group = factor(group, levels = names(TOL_BRIGHT)))
  
  # --- Stats (R^2 and p for slope in lm) ---
  stats_lab <- NULL
  if (add_stats && nrow(d) >= 3) {
    fit <- try(lm(y ~ MF_percent, data = d), silent = TRUE)
    if (!inherits(fit, "try-error")) {
      s   <- summary(fit)
      r2  <- s$r.squared
      pvl <- coef(s)["MF_percent","Pr(>|t|)"]
      stats_lab <- sprintf("R\u00B2 = %.3f; p = %.3f", r2, pvl)
    }
  }
  
  # --- Convex hulls per group (≥3 unique points) ---
  hulls <- d %>%
    dplyr::filter(!is.na(group)) %>%
    dplyr::group_by(group) %>%
    dplyr::group_modify(~{
      u <- unique(.x[, c("MF_percent", "y")])
      if (nrow(u) < 3) return(tibble::tibble(MF_percent = numeric(0), y = numeric(0)))
      idx <- chull(u$MF_percent, u$y)
      tibble::tibble(MF_percent = u$MF_percent[idx], y = u$y[idx])
    }) %>%
    dplyr::ungroup()
  
  # --- Axis ranges + small padding on Y ---
  xr <- range(d$MF_percent, na.rm = TRUE)
  yr <- range(d$y,          na.rm = TRUE)
  dy <- diff(yr); if (!is.finite(dy) || dy == 0) dy <- abs(yr[1]) * 0.1 + 1e-6
  y_limits <- c(yr[1] - y_expand_frac * dy, yr[2] + y_expand_frac * dy)
  
  # --- Annotation anchor helpers ---
  anchor_xy <- function(which,
                        x_inset_left  = 0.12,   # <- increase if text clips on the left
                        x_inset_right = 0.12,   # <- used for right-side labels
                        top_pad_mult  = 0.18,   # must match scale_y_continuous top expand
                        bot_pad_mult  = 0.02,   # must match bottom expand
                        top_frac      = 0.95,
                        bot_frac      = 0.55) {
    dx <- diff(xr); if (!is.finite(dx) || dx == 0) dx <- abs(xr[1]) * 0.1 + 1e-6
    dy <- diff(yr); if (!is.finite(dy) || dy == 0) dy <- 1
    
    x_left  <- xr[1] + x_inset_left  * dx
    x_right <- xr[2] - x_inset_right * dx
    
    y_top_pad    <- yr[2] + top_pad_mult * dy * top_frac
    y_bottom_pad <- yr[1] - bot_pad_mult * dy * bot_frac
    
    switch(which,
           topleft     = c(x_left,  y_top_pad),
           topright    = c(x_right, y_top_pad),
           bottomleft  = c(x_left,  y_bottom_pad),
           bottomright = c(x_right, y_bottom_pad)
    )
  }
  xr <- range(d$MF_percent, na.rm = TRUE)
  yr <- range(d$y,          na.rm = TRUE)
  
  # gene label positioning
  gene_xy  <- anchor_xy(gene_pos,  x_inset_left = 0.16, top_frac = 0.95)
  stats_xy <- anchor_xy(stats_pos, x_inset_right = 0.12)
  
  
  ggplot(d, aes(x = MF_percent, y = y, color = group, label = species)) +
    { if (nrow(hulls)) geom_polygon(
      data = hulls,
      aes(x = MF_percent, y = y, fill = group, group = group),
      inherit.aes = FALSE, color = NA
    ) } +
    geom_point(size = 3, alpha = 0.95) +
    geom_smooth(method = "lm", se = TRUE, linewidth = 0.8,
                linetype = "dashed", color = "grey30") +
    ggrepel::geom_text_repel(size = 3, max.overlaps = Inf,
                             box.padding = 0.4, point.padding = 0.2) +
    scale_color_manual(values = TOL_BRIGHT, drop = FALSE, na.translate = FALSE) +
    scale_fill_manual(values  = TOL_FILL,   drop = FALSE, na.translate = FALSE, guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.12))) +  # 2% bottom, 12% top
    labs(x = "Milk fat (%)", y = ylab, color = color_legend) +
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_blank(),      # no formal title
          legend.position = "right") +
    { if (!is.null(gene_label)) annotate(
      "text",
      x = gene_xy[1], y = gene_xy[2],
      label = gene_label, size = gene_size, fontface = "bold",
      hjust = 0
    ) } +
    { if (!is.null(stats_lab)) annotate(
      "text", x = stats_xy[1], y = stats_xy[2],
      label = stats_lab, size = 4
    ) }
}


## ============================================================
## BLOCK A — WHOLE PROTEIN
## ============================================================

# Compute hydropathy & net charge for each sequence
df_whole <- lapply(seq_along(seqs_raw), function(i) {
  s0 <- toupper(seqs_raw[i])
  s  <- clean_seq(s0)
  
  ch <- if (nchar(s) > 0) Peptides::charge(s, pH = pH, pKscale = pK) else NA_real_
  
  data.frame(
    id               = ids[i],
    species          = extract_species(ids[i]),
    hydropathy_KD_mean = mean_KD(s),
    net_charge          = ch,
    charge_per_res      = if (!is.na(ch)) ch / nchar(s) else NA_real_,
    stringsAsFactors = FALSE
  )
}) %>% dplyr::bind_rows()

norm_key <- function(x) {
  x <- tolower(trimws(x))
  x <- gsub("[[:space:]]+", "", x)   # remove spaces
  x <- gsub("_", "", x)              # remove underscores
  x
}

mf2 <- mf %>%
  mutate(
    join_key = norm_key(id),
    MF_percent = as.numeric(MF_percent)
  )

df_whole2 <- df_whole %>%
  mutate(
    join_key = norm_key(species)
  )

# join with milk fat table
merged_whole <- df_whole2 %>%
  left_join(mf2 %>% select(join_key, MF_percent, group), by = "join_key")


# Collapse to species-level averages
dat_sp_whole <- collapse_to_species(
  merged_whole,
  cols_mean = c("hydropathy_KD_mean", "net_charge", "charge_per_res")
)

# Make group factor consistent with fixed color palette
dat_sp_whole$group <- factor(dat_sp_whole$group, levels = names(TOL_BRIGHT))

# Plots — WHOLE
gg_whole_hydro  <- plot_xy_hulls(
  dat_sp_whole,
  yvar = "hydropathy_KD_mean",
  ylab = "Hydrophobicity (Kyte–Doolittle mean)",
  gene_label = gene_name,              # inside-plot label
  gene_pos   = "topleft",              # change if you like
  stats_pos  = "bottomright"
)

gg_whole_hydro + theme_bw()

gg_whole_charge <- plot_xy_hulls(
  dat_sp_whole,
  yvar = "net_charge",
  ylab = sprintf("Net charge (pH %.1f, %s)", pH, pK),
  gene_label = gene_name,
  gene_pos   = "topleft",
  stats_pos  = "bottomright"
)


print(gg_whole_hydro)
print(gg_whole_charge)

# Save CSV
readr::write_csv(merged_whole, file.path(out_dir, sprintf("%s_properties.csv", gene_name)))

## ============================================================
## BLOCK B - optional — EXTRACELLULAR SEGMENT (alignment cols)
## ============================================================
df_extra <- lapply(seq_along(seqs_raw), function(i) {
  s0  <- toupper(seqs_raw[i])
  sub <- get_aln_segment(s0, extra_cols)   # de-gapped extracellular segment
  
  ch_e <- if (nchar(sub) > 0) Peptides::charge(sub, pH = pH, pKscale = pK) else NA_real_
  data.frame(
    id = ids[i],
    species = extract_species(ids[i]),
    hydropathy_KD_mean_extra = mean_KD(sub),
    net_charge_extra         = ch_e,
    charge_per_res_extra     = if (!is.na(ch_e)) ch_e / nchar(sub) else NA_real_,
    length_region            = nchar(sub),
    stringsAsFactors = FALSE
  )
}) %>% dplyr::bind_rows()

merged_extra <- join_with_mf(df_extra)

dat_sp_extra <- collapse_to_species(
  merged_extra,
  cols_mean = c("hydropathy_KD_mean_extra", "net_charge_extra", "charge_per_res_extra")
)

# Plots — EXTRACELLULAR
gg_extra_hydro  <- plot_xy_hulls(
  dat_sp_extra,
  yvar = "hydropathy_KD_mean_extra",
  ylab = "Hydrophobicity (Kyte–Doolittle mean)",
  gene_label = paste0(gene_name, " extracellular"),
  gene_size  = 5,                     # smaller title text
  gene_pos   = "topleft",
  stats_pos  = "bottomright"
)

gg_extra_charge <- plot_xy_hulls(
  dat_sp_extra,
  yvar = "net_charge_extra",
  ylab = sprintf("Net charge (pH %.1f, %s)", pH, pK),
  gene_label = paste0(gene_name, " extracellular"),
  gene_size  = 5,                     # smaller title text
  gene_pos   = "topleft",
  stats_pos  = "bottomright"
)

print(gg_extra_hydro)
print(gg_extra_charge)

# Save CSV
readr::write_csv(merged_extra, file.path(out_dir, sprintf("%s_properties_EXTRA.csv", gene_name)))
