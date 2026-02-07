## ============================================================
## Embedding comparison: Full vs PC-only vs ncRNA-only
## - Graph metrics from PCAsub (for SNN construction)
## - cLISI & silhouettes for PLOTTING in densMAP space
## - PCA silhouettes (Euclidean & cosine) retained for TABLE
## ============================================================

## ---- Libraries ----
library(SingleCellExperiment)
library(scran)
library(scater)
library(dplyr)
library(tidyr)
library(tibble)
library(igraph)
library(ggplot2)
library(viridis)
library(cowplot)
library(ggtext)
library(FNN)
library(rtracklayer)
library(cluster)
library(forcats)
library(patchwork)  # if you prefer patchwork for plots

## ---- Load data ----
salmon_quants_hg38       <- readRDS("~/Documents/manuscripts/storm_seq/fte_analysis/cryo_fte_novasseq_grch38_salmon_norm_se_with_embeddings.rds")

## ---- Collapse clusters (labels for the FULL object) ----
to_collapsed <- function(x) {
  as_tibble(x) |>
    mutate(cluster_cat = case_when(
      value == 1 ~ "Ciliated Intermediate",
      value == 2 ~ "Secretory Intermediate",
      value == 3 ~ "UCFP",
      value == 4 ~ "Ciliated",
      value == 5 ~ "Branch",
      value == 6 ~ "Immune",
      value == 7 ~ "Secretory",
      value == 8 ~ "Secretory"
    )) |>
    pull(cluster_cat)
}
salmon_quants_hg38$cluster.collapse <- to_collapsed(salmon_quants_hg38$cluster)

## ---- Gene sets: protein-coding and non-coding ----
ens102_gtf       <- rtracklayer::import("~/Documents/manuscripts/storm_seq/fte_analysis/Homo_sapiens.GRCh38.102.gtf.gz")
ens102_gtf.genes <- ens102_gtf[ens102_gtf$type == "gene",]
ens102_gtf.genes.pc <- ens102_gtf.genes[ens102_gtf.genes$gene_biotype == "protein_coding",]
ens102_gtf.genes.nc <- ens102_gtf.genes[ens102_gtf.genes$gene_biotype != "protein_coding",]

## ---- Derive PC-only and ncRNA-only SCEs (patient 1) ----
sce_full <- salmon_quants_hg38
sce_pc   <- sce_full[rownames(sce_full) %in% ens102_gtf.genes.pc$gene_id,]
sce_nc   <- sce_full[rownames(sce_full) %in% ens102_gtf.genes.nc$gene_id,]

## ---- Utility: PCA + PCAsub + SNN clustering + densMAP ----
run_pipeline <- function(sce, pca_name, pcasub_name, dens_name, cluster_col, hvgs_prop = 0.1, ncomp = 20, k = 10, seed = 1988) {
  # HVGs
  dec      <- scran::modelGeneVar(sce)
  top.hvgs <- scran::getTopHVGs(dec, prop = hvgs_prop)
  
  # PCA -> PCAsub
  sce <- scater::runPCA(sce, subset_row = top.hvgs, ncomponents = ncomp, name = pca_name)
  npcs <- metadata(scran::getClusteredPCs(reducedDim(sce, pca_name)))$chosen
  reducedDim(sce, pcasub_name) <- reducedDim(sce, pca_name)[, seq_len(npcs), drop = FALSE]
  
  # SNN graph & clusters on PCAsub
  g <- scran::buildSNNGraph(sce, use.dimred = pcasub_name, k = k)
  cl <- igraph::cluster_walktrap(g)$membership
  colData(sce)[[cluster_col]] <- factor(cl)
  
  # densMAP from PCA (standard practice)
  set.seed(seed)
  dm <- densvis::densmap(reducedDim(sce, pca_name), n_components = 3L, n_neighbors = 15L, metric = "euclidean")
  reducedDim(sce, dens_name) <- dm
  
  sce
}

## ---- Run pipelines ----
sce_full <- run_pipeline(sce_full, pca_name="PCA",     pcasub_name="PCAsub",    dens_name="densMAP",    cluster_col="cluster")
sce_pc   <- run_pipeline(sce_pc,   pca_name="PCA_pc",  pcasub_name="PCAsub_pc", dens_name="densMAP_pc", cluster_col="cluster_pc")
sce_nc   <- run_pipeline(sce_nc,   pca_name="PCA_nc",  pcasub_name="PCAsub_nc", dens_name="densMAP_nc", cluster_col="cluster_nc")

## ---- Graph metrics from SNNs on PCAsub ----
build_graph <- function(sce, use.dimred, k = 10, type = c("snn","knn")) {
  type <- match.arg(type)
  g <- if (type == "snn") scran::buildSNNGraph(sce, use.dimred = use.dimred, k = k)
  else scran::buildKNNGraph(sce, use.dimred = use.dimred, k = k)
  V(g)$name <- colnames(sce)
  as.undirected(g, mode = "collapse", edge.attr.comb = list(weight = "sum", "ignore"))
}

compute_metrics_from_graph <- function(g, labels, use_weights = TRUE) {
  stopifnot(igraph::is_igraph(g))
  vnames <- V(g)$name; if (is.null(vnames)) { vnames <- as.character(seq_len(vcount(g))); V(g)$name <- vnames }
  lab <- as.character(labels); if (is.null(names(lab))) names(lab) <- vnames
  lab <- lab[vnames]; stopifnot(!anyNA(lab))
  
  edf <- as_data_frame(g, what = "edges")
  if (is.numeric(edf$from)) { edf$from <- vnames[edf$from]; edf$to <- vnames[edf$to] }
  edf$from_cl <- lab[edf$from]; edf$to_cl <- lab[edf$to]
  edf$w <- if ("weight" %in% edge_attr_names(g) && use_weights) edf$weight else 1
  
  per_cell_edge <- bind_rows(
    transmute(edf, cell = from, nbr = to,   same = as.integer(from_cl == to_cl), w = w),
    transmute(edf, cell = to,   nbr = from, same = as.integer(from_cl == to_cl), w = w)
  )
  pc_sum <- per_cell_edge |>
    group_by(cell) |>
    summarize(deg_w = sum(w), same_w = sum(w * same), purity = ifelse(deg_w > 0, same_w/deg_w, NA_real_), .groups="drop")
  
  per_cell <- tibble(cell=vnames, cluster=lab, degree=igraph::degree(g),
                     purity = pc_sum$purity[match(vnames, pc_sum$cell)],
                     bridging = 1 - pc_sum$purity[match(vnames, pc_sum$cell)])
  
  per_cluster_cells <- per_cell |>
    group_by(cluster) |>
    summarize(n_cells=n(), med_purity=median(purity, na.rm=TRUE),
              mean_purity=mean(purity, na.rm=TRUE), med_bridge=median(bridging, na.rm=TRUE), .groups="drop")
  
  per_cluster_edges <- bind_rows(
    transmute(edf, cluster = from_cl, w, inter = as.integer(from_cl != to_cl)),
    transmute(edf, cluster = to_cl,   w, inter = as.integer(from_cl != to_cl))
  ) |>
    group_by(cluster) |>
    summarize(incident_w=sum(w), inter_w=sum(w*inter),
              inter_frac_w = ifelse(incident_w>0, inter_w/incident_w, NA_real_), .groups="drop")
  
  list(
    per_cell    = per_cell,
    per_cluster = left_join(per_cluster_cells, per_cluster_edges, by = "cluster"),
    global_inter_frac_w = with(edf, sum(w * (from_cl != to_cl)) / sum(w))
  )
}

## ---- Labels (reference = FULL collapsed labels) ----
ref_labels <- setNames(as.character(sce_full$cluster.collapse), colnames(sce_full))
interests  <- c("Ciliated Intermediate","Secretory Intermediate","Branch","UCFP")
emb_lvls   <- c("Full","PC-only","ncRNA-only")

## ---- Build graphs on PCAsub and compute metrics ----
g_full <- build_graph(sce_full, "PCAsub", k = 10)
g_pc   <- build_graph(sce_pc,   "PCAsub_pc", k = 10)
g_nc   <- build_graph(sce_nc,   "PCAsub_nc", k = 10)

lab_full <- ref_labels[V(g_full)$name]
lab_pc   <- ref_labels[V(g_pc)$name]
lab_nc   <- ref_labels[V(g_nc)$name]

m_full <- compute_metrics_from_graph(g_full, lab_full, use_weights = TRUE)
m_pc   <- compute_metrics_from_graph(g_pc,   lab_pc,   use_weights = TRUE)
m_nc   <- compute_metrics_from_graph(g_nc,   lab_nc,   use_weights = TRUE)

graph_tbl <- bind_rows(
  m_full$per_cluster %>% mutate(embedding = "Full"),
  m_pc$per_cluster   %>% mutate(embedding = "PC-only"),
  m_nc$per_cluster   %>% mutate(embedding = "ncRNA-only")
) %>%
  filter(cluster %in% interests) %>%
  select(cluster, embedding, n_cells, med_purity, inter_frac_w) %>%
  mutate(embedding = factor(embedding, levels = emb_lvls),
         cluster   = factor(cluster,   levels = interests)) %>%
  arrange(cluster, embedding)

## ---- Silhouette helpers ----
dist_euclid <- function(mat) stats::dist(mat)
cosine_dist <- function(mat) { X <- as.matrix(mat); X <- X / sqrt(rowSums(X^2) + 1e-12); as.dist(1 - tcrossprod(X)) }

# densMAP silhouettes vs reference labels (for plotting)
sil_ref_dmap <- function(sce, dens_name, ref_labels, interests, embedding_name){
  emb <- reducedDim(sce, dens_name)
  d   <- dist_euclid(emb)
  labs <- ref_labels[colnames(sce)]
  sil  <- silhouette(as.integer(factor(labs)), d)
  bind_rows(lapply(interests, function(lbl){
    idx <- which(labs == lbl)
    if (!length(idx)) tibble(cluster = lbl, embedding = embedding_name, mean_sil_ref = NA_real_, pct_sil_pos_ref = NA_real_)
    else {
      sw <- sil[idx, "sil_width"]
      tibble(cluster = lbl, embedding = embedding_name,
             mean_sil_ref = mean(sw, na.rm = TRUE),
             pct_sil_pos_ref = mean(sw > 0, na.rm = TRUE) * 100)
    }
  }))
}

# PCA silhouettes (Euclidean & cosine) vs reference labels (for table only)
sil_ref_pca <- function(sce, pcasub_name, ref_labels, interests, embedding_name, dist_fun, metric_tag){
  emb <- reducedDim(sce, pcasub_name)
  d   <- dist_fun(emb)
  labs <- ref_labels[colnames(sce)]
  sil  <- silhouette(as.integer(factor(labs)), d)
  bind_rows(lapply(interests, function(lbl){
    idx <- which(labs == lbl)
    if (!length(idx)) tibble(cluster = lbl, embedding = embedding_name, metric = metric_tag, mean_val = NA_real_)
    else tibble(cluster = lbl, embedding = embedding_name, metric = metric_tag,
                mean_val = mean(sil[idx, "sil_width"], na.rm = TRUE))
  }))
}

## ---- Compute silhouettes ----
sil_dmap <- bind_rows(
  sil_ref_dmap(sce_full, "densMAP",    ref_labels, interests, "Full"),
  sil_ref_dmap(sce_pc,   "densMAP_pc", ref_labels, interests, "PC-only"),
  sil_ref_dmap(sce_nc,   "densMAP_nc", ref_labels, interests, "ncRNA-only")
) %>% mutate(embedding = factor(embedding, levels = emb_lvls),
             cluster   = factor(cluster,   levels = interests))

sil_pca_euc <- bind_rows(
  sil_ref_pca(sce_full, "PCAsub",    ref_labels, interests, "Full",    dist_euclid, "pca_euc"),
  sil_ref_pca(sce_pc,   "PCAsub_pc", ref_labels, interests, "PC-only", dist_euclid, "pca_euc"),
  sil_ref_pca(sce_nc,   "PCAsub_nc", ref_labels, interests, "ncRNA-only", dist_euclid, "pca_euc")
)

sil_pca_cos <- bind_rows(
  sil_ref_pca(sce_full, "PCAsub",    ref_labels, interests, "Full",    cosine_dist, "pca_cos"),
  sil_ref_pca(sce_pc,   "PCAsub_pc", ref_labels, interests, "PC-only", cosine_dist, "pca_cos"),
  sil_ref_pca(sce_nc,   "PCAsub_nc", ref_labels, interests, "ncRNA-only", cosine_dist, "pca_cos")
)

sil_pca_wide <- bind_rows(sil_pca_euc, sil_pca_cos) |>
  tidyr::pivot_wider(id_cols = c(cluster, embedding),
                     names_from = metric, values_from = mean_val,
                     values_fill = NA_real_) |>
  rename(mean_sil_ref_pca_euc = pca_euc,
         mean_sil_ref_pca_cos = pca_cos) %>%
  mutate(embedding = factor(embedding, levels = emb_lvls),
         cluster   = factor(cluster,   levels = interests))

## ---- cLISI on densMAP (lower = better) ----
cLISI_by_label_knn <- function(emb, labels, k = 10, eps = 1e-12) {
  stopifnot(nrow(emb) == length(labels))
  labs <- as.character(labels); levs <- unique(labs)
  nn <- FNN::get.knn(emb, k = k)$nn.index
  LISI <- apply(nn, 1, function(idx) {
    counts <- tabulate(match(labs[idx], levs), nbins = length(levs))
    p <- counts / (sum(counts) + eps); p <- p[p > 0]
    H <- -sum(p * log(p + eps))
    val <- exp(H); if (!is.finite(val)) val <- NA_real_
    max(1, min(val, k))
  })
  tibble(label=labs, LISI=LISI) |>
    group_by(label) |>
    summarize(cLISI_med = median(LISI, na.rm=TRUE),
              cLISI_q25 = quantile(LISI, 0.25, na.rm=TRUE),
              cLISI_q75 = quantile(LISI, 0.75, na.rm=TRUE), .groups="drop")
}

tbl_full <- cLISI_by_label_knn(reducedDim(sce_full, "densMAP"),    sce_full$cluster.collapse, k = 10) %>% mutate(embedding="Full",    cluster=label)
tbl_pc   <- cLISI_by_label_knn(reducedDim(sce_pc,   "densMAP_pc"), sce_pc$cluster.collapse,   k = 10) %>% mutate(embedding="PC-only", cluster=label)
tbl_nc   <- cLISI_by_label_knn(reducedDim(sce_nc,   "densMAP_nc"), sce_nc$cluster.collapse,   k = 10) %>% mutate(embedding="ncRNA-only", cluster=label)

cLISI_tbl <- bind_rows(tbl_full, tbl_pc, tbl_nc) |>
  select(cluster, embedding, cLISI_med, cLISI_q25, cLISI_q75) |>
  filter(cluster %in% interests) |>
  mutate(embedding=factor(embedding, levels=emb_lvls),
         cluster  =factor(cluster,   levels=interests))

## ---- Join tables: graph + densMAP silhouettes + PCA silhouettes + cLISI ----
final_tbl2 <- graph_tbl %>%
  left_join(sil_dmap,    by = c("cluster","embedding")) %>%
  left_join(sil_pca_wide,by = c("cluster","embedding")) %>%
  left_join(cLISI_tbl,   by = c("cluster","embedding")) %>%
  arrange(cluster, embedding)

## ---- Bootstrap CI for med_purity (per cluster/embedding) ----
boot_ci <- function(x, R=999) {
  x <- x[is.finite(x)]
  if (length(x) < 5) return(c(NA_real_, NA_real_))
  qs <- replicate(R, median(sample(x, replace=TRUE)))
  quantile(qs, c(0.025, 0.975), na.rm=TRUE)
}

per_cell_all <- bind_rows(
  m_full$per_cell %>% mutate(embedding="Full"),
  m_pc$per_cell   %>% mutate(embedding="PC-only"),
  m_nc$per_cell   %>% mutate(embedding="ncRNA-only")
)

ci_tbl <- per_cell_all %>%
  filter(cluster %in% interests) %>%
  group_by(cluster, embedding) %>%
  summarize(med_purity_lo = boot_ci(purity)[1],
            med_purity_hi = boot_ci(purity)[2],
            .groups="drop") %>%
  mutate(embedding = factor(embedding, levels = emb_lvls),
         cluster   = factor(cluster,   levels = interests))

final_tbl_big <- final_tbl2 %>%
  left_join(ci_tbl, by = c("cluster","embedding")) %>%
  arrange(cluster, embedding)

## ============================================================
## PLOTTING (densMAP silhouettes; plus CI strip)
## ============================================================

# Heatmaps: purity, mixing, densMAP silhouette, cLISI
hm_tbl <- final_tbl_big

p_purity <- ggplot(hm_tbl, aes(embedding, cluster, fill = med_purity)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.2f", med_purity)), size = 3.8) +
  viridis::scale_fill_viridis(name = "Median purity", option = "D", limits = c(0,1)) +
  labs(title = "Local cluster purity (higher is better)", x=NULL, y=NULL) +
  theme_minimal(base_size = 11) + theme(panel.grid = element_blank(),
                                        plot.title = element_text(face="bold", hjust=0.5))

p_mix <- ggplot(hm_tbl, aes(embedding, cluster, fill = inter_frac_w)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.2f", inter_frac_w)), size = 3.8) +
  viridis::scale_fill_viridis(name = "Inter-cluster\nedge fraction", option = "C", limits = c(0,1), direction = -1) +
  labs(title = "Cross-cluster mixing (lower is better)", x=NULL, y=NULL) +
  theme_minimal(base_size = 11) + theme(panel.grid = element_blank(),
                                        plot.title = element_text(face="bold", hjust=0.5))

p_sil <- ggplot(hm_tbl, aes(embedding, cluster, fill = mean_sil_ref)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.2f", mean_sil_ref)), size = 3.8) +
  scale_fill_gradient2(name = "Mean silhouette", low = "steelblue", mid = "white", high = "firebrick",
                       limits = c(-1,1), midpoint = 0) +
  labs(title = "Cluster compactness (silhouette, more neg. suggests misassignment)", x=NULL, y=NULL) +
  theme_minimal(base_size = 11) + theme(panel.grid = element_blank(),
                                        plot.title = element_text(face="bold", hjust=0.5))

p_cLISI <- ggplot(hm_tbl, aes(embedding, cluster, fill = cLISI_med)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.2f", cLISI_med)), size = 3.8) +
  viridis::scale_fill_viridis(name = "Score", option = "C", limits = c(1,4), direction = -1) +
  labs(title = "Neighborhood label mixing (cLISI, k = 10, lower is better)", x=NULL, y=NULL) +
  theme_minimal(base_size = 11) + theme(panel.grid = element_blank(),
                                        plot.title = element_text(face="bold", hjust=0.5))

# Bootstrap CI strip for med_purity (transitional clusters)
transitional <- interests
p_purity_ci <- final_tbl_big %>%
  filter(cluster %in% transitional) %>%
  ggplot(aes(x = embedding, y = med_purity,
             ymin = med_purity_lo, ymax = med_purity_hi,
             color = embedding)) +
  geom_pointrange(position = position_dodge(width = 0.5), size = 0.7) +
  scale_color_viridis_d(option = "D", end = 0.85, name = "Embedding") +
  coord_cartesian(ylim = c(0,1)) +
  facet_wrap(~ cluster, nrow = 1) +
  labs(y = "Median neighborhood\npurity (95% CI)", x = NULL,
       title = "Purity estimate stability (bootstrap CI)") +
  theme_bw(12) +
  theme(
    legend.position = "right",                  # single shared legend
    legend.title = element_text(face = "bold"),
    strip.background = element_rect(fill = "white", color = NA),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5),
    strip.text = element_text(face = "bold", size = 12),
    axis.text.x = element_blank(),               # remove repeated x labels
    axis.ticks.x = element_blank()
  )

## ---- Assemble minimal panel ----
top_row   <- plot_grid(p_purity, p_mix, ncol = 2)
mid_row   <- plot_grid(p_sil, ncol = 1)
row_cLISI <- plot_grid(p_cLISI, ncol = 1)
row_ci    <- plot_grid(p_purity_ci, ncol = 1)

panel_minimal <- plot_grid(
  top_row,
  mid_row,
  row_cLISI,
  row_ci,
  ncol = 1,
  rel_heights = c(1, 1.05, 1, 0.9)
)
panel_minimal

## ---- Save (optional) ----
# ggsave("embedding_comparison_panel_minimal.pdf", panel_minimal, width=10, height=8.5, useDingbats=FALSE)
# ggsave("embedding_comparison_panel_minimal.png", panel_minimal, width=10, height=8.5, dpi=300)

## ---- Final table for export (includes PCA silhouettes) ----
# Columns include:
# cluster | embedding | n_cells | med_purity | inter_frac_w |
# mean_sil_ref (densMAP) | pct_sil_pos_ref (densMAP) |
# mean_sil_ref_pca_euc | mean_sil_ref_pca_cos |
# cLISI_med | cLISI_q25 | cLISI_q75 | med_purity_lo | med_purity_hi
final_tbl_big
# write.csv(final_tbl_big, "ncRNA_PC_full_transitional_metrics.csv", row.names = FALSE)