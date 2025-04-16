library(Seurat)
library(tidyverse)
library(RColorBrewer)
library(data.table)
library(ggpubr)
setwd("~/data_HD/Project/TS/10XData/Figure1")

IgD = subset(merged_object, subset = celltype == "IgD-")
naive = subset(merged_object, subset = celltype == "naive")

plot_samples = merged_object$sample
color_ct=c(brewer.pal(12, "Set3")[-c(2,3,9,12)],"#b3b3b3",
           brewer.pal(5, "Set1")[2],
           brewer.pal(3, "Dark2")[1],
           "#fc4e2a","#fb9a99","#f781bf","#e7298a")
names(color_ct)=levels(plot_samples)

# --------------------------XIST and TSIX---------------------------
plot_genes = c("XIST", "TSIX")

### IgD- B cells
vln.df = as.data.frame(GetAssayData(IgD, assay = "RNA", slot = "data")[plot_genes,])
vln.df$gene = rownames(vln.df)
vln.df = melt(vln.df, id="gene")
colnames(vln.df)[c(2,3)] = c("CB", "exp")

anno = IgD[[]][,c("CB", "sample")]
vln.df = inner_join(vln.df, anno, by="CB")
vln.df$gene = factor(vln.df$gene, levels = plot_genes)

my_comparisons <- list(c("TS2", "TS2M"), c("TS2F", "TS2M"),
                      c("TS3", "TS3M"), c("TS3F", "TS3M"))

vln.df |> ggplot(aes(sample, exp))+
  geom_violin(aes(fill=sample), scale = "width")+
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.signif", size = 3)+
  facet_grid(gene~.,scales = "free_y", switch = "y")+
  scale_fill_manual(values = color_ct)+
  scale_y_continuous(expand = c(0,0))+
  theme_bw()+
  theme(
    panel.grid = element_blank(),
    
    axis.title.x.bottom = element_blank(),
    axis.ticks.x.bottom = element_blank(),
    axis.text.x.bottom = element_text(angle = 45,hjust = 1,vjust = NULL,color = "black",size = 14),
    axis.title.y.left = element_blank(),
    axis.ticks.y.left = element_blank(),
    axis.text.y.left = element_blank(),
    
    legend.position = "none",
    
    panel.spacing.y = unit(0, "cm"),
    strip.text.y.left = element_text(angle=0,size = 14,hjust = 0),
    strip.background.y = element_blank(),
    strip.placement = "outside"
  )+
  ggtitle("Expression of key genes in X chromosome in IgD- B cells")
ggsave("../../figs/fig2/IgD-/X_xist_tsix_exp.pdf", width = 5, height = 4)

### naive B cells
vln.df = as.data.frame(GetAssayData(naive, assay = "RNA", slot = "data")[plot_genes,])
vln.df$gene = rownames(vln.df)
vln.df = melt(vln.df, id="gene")
colnames(vln.df)[c(2,3)] = c("CB", "exp")

anno = naive[[]][,c("CB", "sample")]
vln.df = inner_join(vln.df, anno, by="CB")
vln.df$gene = factor(vln.df$gene, levels = plot_genes)

my_comparisons <- list(c("TS2", "TS2M"), c("TS2F", "TS2M"),
                      c("TS3", "TS3M"), c("TS3F", "TS3M"))

vln.df |> ggplot(aes(sample, exp))+
  geom_violin(aes(fill=sample), scale = "width")+
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.signif", size = 3)+
  facet_grid(gene~.,scales = "free_y", switch = "y")+
  scale_fill_manual(values = color_ct)+
  scale_y_continuous(expand = c(0,0))+
  theme_bw()+
  theme(
    panel.grid = element_blank(),
    
    axis.title.x.bottom = element_blank(),
    axis.ticks.x.bottom = element_blank(),
    axis.text.x.bottom = element_text(angle = 45,hjust = 1,vjust = NULL,color = "black",size = 14),
    axis.title.y.left = element_blank(),
    axis.ticks.y.left = element_blank(),
    axis.text.y.left = element_blank(),
    
    legend.position = "none",
    
    panel.spacing.y = unit(0, "cm"),
    strip.text.y.left = element_text(angle=0,size = 14,hjust = 0),
    strip.background.y = element_blank(),
    strip.placement = "outside"
  )+
  ggtitle("Expression of key genes in X chromosome in naive B cells")
ggsave("../../figs/fig2/naive/X_xist_tsix_exp.pdf", width = 5, height = 4)
