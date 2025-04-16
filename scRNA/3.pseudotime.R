library(Seurat)
library(SCP)
library(tidyverse)

setwd("~/data_HD/Project/TS/10XData/Figure1/")

merged_object = read_rds("combinedReannotation.rds")

# ----------------------------pseudotime------------------------------
sample.list = list()
for (state in unique(merged_object$sample)) {
  sample.list[[state]] = subset(merged_object, subset = state == sample)
  sample.list[[state]] <- RunSlingshot(srt = sample.list[[state]], group.by = "celltype",
                                       reduction = "UMAP", start = "naive")
  SCP5::CellDimPlot(sample.list[[state]], group.by = "celltype", reduction = "UMAP", 
                    lineages = paste0("Lineage", 1:2), lineages_span = 0.1)+
    ggtitle(paste0("Trajectory of ", state))
  ggsave(paste0("../../figs/fig1/pseudotime/trajectory_", state, ".pdf"), width = 8, height = 5)
}

# -------------------gene expression alteration-------------------------
merged_object = RunDynamicFeatures(srt = merged_object, 
                                   lineages = paste0("Lineage", 1:lineage_num), 
                                   n_candidates = 200)
# naive2IgD-
ht_m <- DynamicHeatmap(
  srt = merged_object, lineages = c("Lineage1"),
  use_fitted = TRUE, n_split = 2, r.sq = 0.2,
  species = "Homo_sapiens", db = c("GO_BP"), anno_terms = TRUE,
  heatmap_palette = "viridis", cell_annotation = "celltype",
  separate_annotation = list("celltype"), separate_annotation_palette = c("Paired"),
  pseudotime_label = 25, pseudotime_label_color = "red",
  height = 3, width = 2,
)
ht_m$plot

target = names(ht_m$feature_split)

sample.list = list()
ht_IgD = list()
for (state in unique(merged_object$sample)) {
  sample.list[[state]] = subset(merged_object, subset = state == sample)
  sample.list[[state]] <- RunSlingshot(srt = sample.list[[state]], group.by = "celltype",
                                       reduction = "UMAP", start = "naive")
  sample.list[[state]] <- RunDynamicFeatures(
    srt = sample.list[[state]], lineages = c("Lineage1", "Lineage2"),
    n_candidates = length(target), minfreq = 0, features = target)
  lineage = ifelse(
    sample.list[[state]]@tools[["Slingshot_celltype_umap"]]@metadata[["lineages"]][["Lineage1"]][3] == "IgD-",
    c("Lineage1"), c("Lineage2")
  )
  ht_IgD[[state]] <- DynamicHeatmap(
    srt = sample.list[[state]], lineages = lineage, features = target,
    use_fitted = TRUE, n_split = 2, #feature_split = ht_m$feature_split,
    min_expcells = 0, r.sq = 0, dev.expl = 0, padjust = 1,
    species = "Homo_sapiens", db = c("GO_BP"), anno_terms = TRUE,
    heatmap_palette = "viridis", cell_annotation = "celltype",
    separate_annotation = list("celltype"), separate_annotation_palette = c("Paired"),
    pseudotime_label = 25, pseudotime_label_color = "red",
    height = 3, width = 2,
  )
  ht_IgD[[state]]$plot + ggtitle(paste0(state, " naive2IgD-"))
  ggsave(paste0("../../figs/fig1/pseudotime/heatmap_l1_", state, ".pdf"), width = 10, height = 5)
}

# naive2plasma
ht_m <- DynamicHeatmap(
  srt = merged_object, lineages = c("Lineage2"),
  use_fitted = TRUE, n_split = 2, r.sq = 0.2,
  species = "Homo_sapiens", db = c("GO_BP"), anno_terms = TRUE,
  heatmap_palette = "viridis", cell_annotation = "celltype",
  separate_annotation = list("celltype"), separate_annotation_palette = c("Paired"),
  pseudotime_label = 25, pseudotime_label_color = "red",
  height = 3, width = 2,
)
ht_m$plot

target = names(ht_m$feature_split)

ht_plasma = list()
for (state in unique(merged_object$sample)) {
  sample.list[[state]] = subset(merged_object, subset = state == sample)
  sample.list[[state]] <- RunSlingshot(srt = sample.list[[state]], group.by = "celltype",
                                       reduction = "UMAP", start = "naive")
  sample.list[[state]] <- RunDynamicFeatures(
    srt = sample.list[[state]], lineages = c("Lineage1", "Lineage2"),
    n_candidates = length(target), minfreq = 0, features = target)
  lineage = ifelse(
    sample.list[[state]]@tools[["Slingshot_celltype_umap"]]@metadata[["lineages"]][["Lineage1"]][3] == "IgD-",
    c("Lineage2"), c("Lineage1")
  )
  ht_plasma[[state]] <- DynamicHeatmap(
    srt = sample.list[[state]], lineages = lineage, features = target,
    use_fitted = TRUE, n_split = 2, feature_split = ht_m$feature_split,
    min_expcells = 0, r.sq = 0, dev.expl = 0, padjust = 1,
    species = "Homo_sapiens", db = c("GO_BP"), anno_terms = TRUE,
    heatmap_palette = "viridis", cell_annotation = "celltype",
    separate_annotation = list("celltype"), separate_annotation_palette = c("Paired"),
    pseudotime_label = 25, pseudotime_label_color = "red",
    height = 3, width = 2,
  )
  ht_plasma[[state]]$plot + ggtitle(paste0(state, " naive2plasma"))
  ggsave(paste0("../../figs/fig1/pseudotime/heatmap_l2_", state, ".pdf"), width = 10, height = 5)
}

