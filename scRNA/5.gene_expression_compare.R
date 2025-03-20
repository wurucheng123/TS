library(Seurat)
library(tidyverse)
library(RColorBrewer)

merged_object = read_rds("combinedReannotation.rds")

plot_genes = c("SLC25A6", "CD99", "ZFX", "KDM5C")
plot_samples = gsub("BCELLS", "", merged_object$orig.ident)
plot_samples = as.factor(plot_samples)
merged_object$sample = plot_samples

color_ct=c(brewer.pal(12, "Set3")[-c(2,3,9,12)],"#b3b3b3",
           brewer.pal(5, "Set1")[2],
           brewer.pal(3, "Dark2")[1],
           "#fc4e2a","#fb9a99","#f781bf","#e7298a")
names(color_ct)=levels(plot_samples)

### IgD- B cells
vln.df = as.data.frame(GetAssayData(IgD, assay = "RNA", slot = "data")[plot_genes,])
vln.df$gene = rownames(vln.df)
vln.df = melt(vln.df, id="gene")
colnames(vln.df)[c(2,3)] = c("CB", "exp")

anno = IgD[[]][,c("CB", "sample")]
vln.df = inner_join(vln.df, anno, by="CB")
vln.df$gene = factor(vln.df$gene, levels = plot_genes)

vln.df |> ggplot(aes(sample, exp))+
  geom_violin(aes(fill=sample), scale = "width")+
  facet_grid(gene~.,scales = "free_y")+
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
    strip.text.y = element_text(angle=0,size = 14,hjust = 0),
    strip.background.y = element_blank()
  )+
  ggtitle("Expression of key genes in X chromosome in IgD- B cells")
ggsave("../../figs/fig2/IgD-/X_key_gene_exp.pdf", width = 5, height = 4)

### naive B cells
vln.df = as.data.frame(GetAssayData(naive, assay = "RNA", slot = "data")[plot_genes,])
vln.df$gene = rownames(vln.df)
vln.df = melt(vln.df, id="gene")
colnames(vln.df)[c(2,3)] = c("CB", "exp")

anno = naive[[]][,c("CB", "sample")]
vln.df = inner_join(vln.df, anno, by="CB")
vln.df$gene = factor(vln.df$gene, levels = plot_genes)

vln.df |> ggplot(aes(sample, exp))+
  geom_violin(aes(fill=sample), scale = "width")+
  facet_grid(gene~.,scales = "free_y")+
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
    strip.text.y = element_text(angle=0,size = 14,hjust = 0),
    strip.background.y = element_blank()
  )+
  ggtitle("Expression of key genes in X chromosome in naive B cells")
ggsave("../../figs/fig2/naive/X_key_gene_exp.pdf", width = 5, height = 4)