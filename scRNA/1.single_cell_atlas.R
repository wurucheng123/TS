library(Seurat)
library(SCP)
library(tidyverse)

setwd("~/data_HD/Project/TS/10XData/Figure1/")

merged_object = read_rds("combinedReannotation.rds")

CellDimPlot(
  srt = merged_object, group.by = c("clusters", "celltype", "sample"),
  reduction = "UMAP", theme_use = "theme_blank"
)
ggsave("../../figs/fig1/combine_umap.pdf", width = 12, height = 3)

p1 = CellStatPlot(srt = merged_object, stat.by = "celltype", group.by = "sample")
p2 = CellStatPlot(srt = merged_object, stat.by = "clusters", 
                  group.by = "sample")
p2+p1
ggsave("../../figs/fig1/combine_bar.pdf", width = 12, height = 4)

ht <- GroupHeatmap(
  srt = merged_object,
  features = c(
    "CD27", # Naive
    "MZB1", # plasma
    "IGHD", # IgD+
    "IGHG1", # IgD-
    "IGHM", # 
    "TCL1A", "AC079767.4", "AIM2", "TNFRSF13B", "NEIL1", "FGR",
    "PRDM1"
  ),
  group.by = c("clusters", "celltype"),
  heatmap_palette = "YlOrRd",
  show_row_names = TRUE, row_names_side = "left",
  add_dot = TRUE, add_reticle = TRUE
)
print(ht$plot)
ggsave("../../figs/fig1/bcell_dotplot.pdf", width = 12, height = 5)

