library(Seurat)
library(SCP)
library(tidyverse)
library(ggrepel)
library(scatterpie)

setwd("~/data_HD/Project/TS/10XData/Figure1/")
X_genes = read_tsv("~/ref/human/bed/hg19_chrX_genes.txt", col_names = F)
X_genes = X_genes$X1

Y_genes = read_tsv("~/ref/human/bed/hg19_chrY_genes.txt", col_names = F)
Y_genes = Y_genes$X1

merged_object = read_rds("combinedReannotation.rds")

# --------------------------cor heatmap-------------------------------
cell.list = list()
p.final = NULL
for (state in levels(merged_object$celltype)) {
  cell.list[[state]] = subset(merged_object, subset = state == celltype)
  df = AverageExpression(cell.list[[state]], group.by = "sample", 
                             assays = "RNA", slot = "data")
  df = as.data.frame(df[[1]])
  cor = cor(df, method = "spearman")
  dat = melt(cor)
  colnames(dat) = c("x", "y", "spearman's r")
  
  dat = dat[as.numeric(dat$x) < as.numeric(dat$y),]
  
  samples = colnames(cor)
  
  ggplot(dat, aes(x=x, y=y))+
    geom_tile(aes(fill=`spearman's r`))+
    coord_equal(clip = "off")+
    geom_text(data = data.frame(x=1:length(samples), y=0:(length(samples)-1),
                                label=samples),
              aes(x=x, y=y, label=label), angle=90)+
    theme_void()+
    scale_fill_gradient(low = "white",
                         high = "red",
                        limits = c(0.9, 1)) -> p
  p+theme(legend.position = "none") -> p1
  p+theme(legend.position = "bottom")+
    guides(fill=guide_colorbar(title.position = "top",
                               title.hjust = 0.5, barwidth = 20)) -> p2
  ggdraw()+
    draw_plot(ggplotify::as.ggplot(p1, angle = -45),
              width = 0.7, height = 0.7,
              hjust = -0.2, vjust = -0.1)+
    annotate(geom = "text", x=0.5, y=0.15, label=state)+
    annotation_custom(grob = ggpubr::get_legend(p2),
                      xmin = 0.1, xmax = 0.9,
                      ymin = 0.2, ymax = 0.3) -> p3
  if (is.null(p.final)) {
    p.final = p3
  } else {p.final = p.final|p3}
}
p.final
ggsave("../../figs/fig2/cor_heatmap.pdf", width = 20, height = 5)


# --------------------------scatter plot-------------------------
scatter_plot = function(celltype, children, target = NULL) {
  X_genes = read_tsv("~/ref/human/bed/hg19_chrX_genes.txt", col_names = F)
  X_genes = X_genes$X1
  
  Y_genes = read_tsv("~/ref/human/bed/hg19_chrY_genes.txt", col_names = F)
  Y_genes = Y_genes$X1
  pure.Y = setdiff(Y_genes, X_genes)
  
  for (cell in celltype) {
    obj = subset(merged_object, subset = celltype == cell)
    Idents(obj) = "sample"
    datasets = AverageExpression(obj, group.by = "sample", 
                                 assays = "RNA", slot = "data")
    datasets = as.data.frame(datasets[[1]])
    
    datasets = datasets |> filter(!rownames(datasets) %in% pure.Y)
    
    cor_mat = cor(datasets, method = "spearman")
    p_mat = cor.mtest(datasets, method = "spearman", conf.level = 0.95)$p
    
    # scatter plot
    datasets = log2(datasets + 1)
    datasets$chr = "autosome"
    datasets$chr[rownames(datasets) %in% X_genes] = "chrX"
    datasets = na.omit(datasets)
    
    for (child in children) {
      for (parent in c("F", "M")) {
        prt = paste0(child, parent)
        if (is.null(target)) {
          keys = datasets |> mutate(diff = !!sym(prt) - !!sym(child)) |> top_n(n = 10, wt = abs(diff))
          datasets[[paste0("diff_", prt)]] = abs(datasets[[prt]] - datasets[[child]])
          
          markers = FindMarkers(obj, ident.1 = child, ident.2 = prt, logfc.threshold = 0)
          datasets[[paste0("log2FC_", prt)]] = markers[rownames(datasets), "avg_log2FC"]
          datasets[[paste0("p_", prt)]] = markers[rownames(datasets), "p_val_adj"]
        } else {
          keys = datasets[target, ]
        }
        
        ggscatter(datasets, x = child, y = prt, add = "reg.line",
                  color = "chr", conf.int = T, 
                  add.params = list(color = "#eb2e2a", fill = "lightgray"), 
                  cor.coef = F)+
          theme(
            panel.background = element_rect(colour = "black", size = 1),
            axis.ticks = element_line(size = 1),
            axis.text = element_text(size = 16, face = "bold"),
            axis.title = element_text(size = 20, face = "bold"),
            plot.title = element_text(hjust = 0, size = 20, face = "bold")
          )+
          geom_abline(slope = 1, intercept = 0, color = "blue")+
          scale_color_manual(values = c("chrX" = "#ef7156", "autosome" = "#d3d3d3")) + # 根据染色体自定义颜色
          ggtitle(paste0(cell, "\nRs = ", 
                         round(cor_mat[child, prt], digits = 4),
                         "; p = ", sprintf("%.1e", p_mat[child, prt])))+
          geom_text_repel(data=keys,aes(x=!!sym(child),y=!!sym(prt),label=rownames(keys)), size = 3)+
          geom_point(data=keys, aes(x=!!sym(child), y=!!sym(prt)), shape=21, color="green", fill=NA, size=5)
        if (is.null(target)) {
          ggsave(paste0("../../figs/fig2/correlation/cor_exp_", cell, "_", child, "_", parent, "_scatter.pdf"),
                 width = 5, height = 5.5)
        } else {
          ggsave(paste0("../../figs/fig2/target/cor_exp_", cell, "_", child, "_", parent, "_scatter.pdf"),
                 width = 5, height = 5.5)
        }
      }
    }
    if (is.null(target)) {
      datasets$genes = rownames(datasets)
      write_tsv(datasets, paste0("../../figs/fig2/correlation/", cell, ".tsv"))
    }
  }
}

cells = c("naive", "IgD-")
children = c("TS2", "TS3")
scatter_plot(cells, children)