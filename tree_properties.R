# =============================
# 0) Libraries
# =============================
suppressPackageStartupMessages({
  library(ape)
  library(ggtree)
  library(Biostrings)
  library(dplyr)
  library(stringr)
  library(ggplot2)
})

# =============================
# 1) Load inputs
# =============================
setwd("/path/to/wd")

tree <- read.tree("BTN1A1_dated_tree.txt")

anc_seqs <- readAAStringSet(
  "/path/to/BTN1A1_reconstructed.fasta"
)

traits_file <- "/path/to/BTN1A1_properties.csv"
has_traits <- file.exists(traits_file)
if (has_traits) traits <- read.csv(traits_file, stringsAsFactors = FALSE)

# =============================
# 2) Ancestral sequences table
# =============================
anc_df <- tibble(label = names(anc_seqs), seq = as.character(anc_seqs)) %>%
  mutate(node = as.integer(str_remove(label, "node_"))) %>%
  filter(!is.na(node))

# =============================
# 3) Node label lookup (tips + internal nodes)
# =============================
n_tips <- length(tree$tip.label)

tip_labels  <- setNames(tree$tip.label, 1:n_tips)
node_labels <- setNames(paste0("node_", (n_tips + 1):(n_tips + tree$Nnode)),
                        (n_tips + 1):(n_tips + tree$Nnode))

all_labels <- c(tip_labels, node_labels)

# =============================
# 4) Build edge table + mutations
# =============================
edges <- as.data.frame(tree$edge)
colnames(edges) <- c("parent", "child")

edges <- edges %>%
  mutate(
    parent_label = all_labels[as.character(parent)],
    child_label  = all_labels[as.character(child)]
  ) %>%
  left_join(anc_df, by = c("parent" = "node")) %>% dplyr::rename(parent_seq = seq) %>%
  left_join(anc_df, by = c("child"  = "node")) %>% dplyr::rename(child_seq  = seq)

get_mutations <- function(pseq, cseq) {
  if (is.na(pseq) || is.na(cseq)) return(character(0))
  n <- min(nchar(pseq), nchar(cseq))
  out <- character(0)
  for (i in seq_len(n)) {
    p <- substr(pseq, i, i)
    c <- substr(cseq, i, i)
    if (p != c) out <- c(out, paste0(p, i, c))
  }
  out
}

edges$mutations <- mapply(get_mutations, edges$parent_seq, edges$child_seq, SIMPLIFY = FALSE)

# =============================
# 5) Clean traits ONCE (contains BOTH net charge + hydro + milk fat)
# =============================
if (has_traits) {
  traits_clean <- traits %>%
    transmute(
      label          = as.character(id),              # must match tips + node_##
      net_charge     = as.numeric(net_charge),
      hydrophobicity = as.numeric(hydropathy_KD_mean),
      milk_fat       = as.numeric(MF_percent)
    )
}

# =============================
# 6) Compute deltas on edges (KEEP SEPARATE COLUMNS)
# =============================
if (has_traits) {
  edges <- edges %>%
    left_join(traits_clean %>% select(label, net_charge),
              by = c("parent_label" = "label")) %>%
    dplyr::rename(parent_net_charge = net_charge) %>%
    left_join(traits_clean %>% select(label, net_charge),
              by = c("child_label" = "label")) %>%
    dplyr::rename(child_net_charge = net_charge) %>%
    mutate(
      delta_charge = child_net_charge - parent_net_charge,
      delta_charge_label = if_else(!is.na(delta_charge) & abs(delta_charge) >= 1,
                                   paste0("ΔNC=", ifelse(delta_charge > 0, "+", "−"), abs(round(delta_charge, 1))),
                                   NA_character_)
    ) %>%
    left_join(traits_clean %>% select(label, hydrophobicity),
              by = c("parent_label" = "label")) %>%
    dplyr::rename(parent_hydro = hydrophobicity) %>%
    left_join(traits_clean %>% select(label, hydrophobicity),
              by = c("child_label" = "label")) %>%
    dplyr::rename(child_hydro = hydrophobicity) %>%
    mutate(
      delta_hydro = child_hydro - parent_hydro,
      delta_hydro_label = if_else(!is.na(delta_hydro) & abs(delta_hydro) >= 0.2,   # <-- set your threshold here
                                  paste0("ΔH=", ifelse(delta_hydro > 0, "+", "−"), abs(round(delta_hydro, 2))),
                                  NA_character_)
    )
}

hydro_thr <- 0.05  # try 0.08 or 0.05 if you want more labels

edges <- edges %>%
  mutate(
    delta_hydro_label = dplyr::if_else(
      !is.na(delta_hydro) & abs(delta_hydro) >= hydro_thr,
      paste0("ΔH=", ifelse(delta_hydro > 0, "+", "−"), abs(round(delta_hydro, 2))),
      NA_character_
    )
  )

# =============================
# 7) One mutation summary table (optional)
# =============================
mutation_summary <- edges %>%
  mutate(
    branch_node = child,
    mut_label = vapply(mutations, function(m) {
      if (length(m) == 0) NA_character_
      else if (length(m) > 1) paste0(m[1], ", ...")
      else m[1]
    }, character(1))
  ) %>%
  distinct(branch_node, .keep_all = TRUE)

# =============================
# 8) Build ggtree data once
# =============================
p <- ggtree(tree)

# internal-node ages only; store in vector indexed by node id (as character)
node_ages <- branching.times(tree)
n_total <- length(tree$tip.label) + tree$Nnode
all_node_times <- rep(0, n_total)
names(all_node_times) <- as.character(seq_len(n_total))
all_node_times[names(node_ages)] <- node_ages

tree_data <- p$data %>%
  mutate(label = ifelse(isTip, label, paste0("node_", node))) %>%
  left_join(if (has_traits) traits_clean else tibble(label = character()),
            by = "label") %>%
  left_join(mutation_summary %>% select(branch_node, mut_label),
            by = c("node" = "branch_node")) %>%
  mutate(age = all_node_times[as.character(node)])

# Join delta labels for plotting at nodes
if (has_traits) {
  node_delta <- edges %>%
    transmute(
      node = child,
      delta_charge_label,
      delta_hydro_label
    ) %>%
    distinct(node, .keep_all = TRUE)
  
  tree_data <- tree_data %>%
    left_join(node_delta, by = "node")
}



# =============================
# 9) Binning + milk fat categories (do ONCE)
# =============================
if (has_traits) {
  tree_data <- tree_data %>%
    mutate(
      milk_fat_cat = case_when(
        is.na(milk_fat)   ~ "Missing",
        milk_fat < 15     ~ "≤15%",
        milk_fat <= 30    ~ "30%",
        milk_fat <= 45    ~ "45%",
        milk_fat > 45     ~ ">45%"
      ),
      milk_fat_cat = factor(milk_fat_cat, levels = c("≤15%", "30%", "45%", ">45%", "Missing")),
      
      net_charge_bin = case_when(
        net_charge <= 2.5 ~ "2.5",
        net_charge <= 5 ~ "5",
        net_charge <= 7.5  ~ "7.5",
        net_charge <= 10.5  ~ "10.5",
        TRUE              ~ "Unbinned"
      ),
      net_charge_bin = factor(net_charge_bin, levels = rev(c("2.5", "5", "7.5", "10.5", "Unbinned"))),
      
      hydro_bin = case_when(
        hydrophobicity <= -0.3  ~ "-0.3",
        hydrophobicity <= -0.2 ~ "-0.2",
        hydrophobicity <= -0.1  ~ "-0.1",
        hydrophobicity <= -0  ~ "0",
        TRUE                    ~ "Unbinned"
      ),
      hydro_bin = factor(hydro_bin, levels = rev(c("-0.3", "-0.2", "-0.1", "0", "Unbinned")))
    )
}

# =============================
# 10) Plot functions (reuse)
# =============================
size_values <- c("≤15%" = 0.3, "30%" = 1, "45%" = 1.5, ">45%" = 2.2, "Missing" = 0.5)

plot_base <- function(color_var, color_scale, title) {
  p %<+% tree_data +
    aes(color = .data[[color_var]], size = milk_fat_cat, linetype = milk_fat_cat) +
    geom_tree() +
    geom_tiplab(size = 3, align = TRUE, linetype = "dotted") +
    theme_tree2() +
    color_scale +
    scale_size_manual(values = size_values, name = "Milk Fat %") +
    scale_linetype_manual(
      values = c("≤15%"="solid","30%"="solid","45%"="solid",">45%"="solid","Missing"="dotted"),
      name = "Milk Fat %"
    ) +
    labs(title = title) +
    scale_x_continuous(
      name = "Time (Ma)",
      breaks = seq(0, ceiling(max(tree_data$x, na.rm = TRUE)), by = 20),
      labels = function(x) round(rev(x), 0)
    )
}

# --- Net charge plot ---
base_plot_charge <- plot_base(
  color_var = "net_charge_bin",
  color_scale = scale_color_manual(
    name = "Net Charge",
    values = c(
      "2.5" = "#191BFA",
      "5" = "#487BCD",
      "7.5"  = "#AA1166EE",
      "10.5"  = "#FF8C69",
      "Unbinned" = "grey70"
    ),
    drop = TRUE
  ),
  title = "CD36 extracellular — Net charge"
)

if (has_traits) {
  base_plot_charge <- base_plot_charge +
    geom_label(
      data = tree_data %>% filter(!isTip, !is.na(delta_charge_label)),
      aes(x = x - 12, y = y, label = delta_charge_label),
      fill = "gold", size = 2.5, hjust = 0
    )
}

print(base_plot_charge)

# --- Hydrophobicity plot ---
base_plot_hydro <- plot_base(
  color_var = "hydro_bin",
  color_scale = scale_color_manual(
    name = "Hydrophobicity (KD mean, binned)",
    values = c(
      "-0.3" = "#08306B",
      "-0.2" = "#487BCD",
      "-0.1" = "#AA1166EE",
      "0" = "#FF8C69",
      "Unbinned" = "grey70"
    ),
    drop = TRUE
  ),
  title = "BTN1A1 extracellular — Hydrophobicity"
)

if (has_traits) {
  base_plot_hydro <- base_plot_hydro +
    geom_label(
      data = tree_data %>% filter(!isTip, !is.na(delta_hydro_label)),
      aes(x = x - 12, y = y, label = delta_hydro_label),
      fill = "gold", size = 2.5, hjust = 0
    )
}

print(base_plot_hydro)

# =============================
# 11) Export mutation log (keep both deltas, no confusion)
# =============================
mutation_log <- edges %>%
  mutate(
    parent = as.character(parent),
    child  = as.character(child),
    mutation_string = vapply(mutations, function(m) paste(m, collapse = ", "), character(1))
  ) %>%
  select(parent, parent_label, child, child_label,
         delta_charge, delta_charge_label,
         delta_hydro,  delta_hydro_label,
         mutation_string)

write.csv(mutation_log, "asr_mutation_log.csv", row.names = FALSE)
