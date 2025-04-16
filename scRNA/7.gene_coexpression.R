library(Seurat)
library(tidyverse)
library(qqman)
library(data.table)
library(ggsci)
library(ggrepel)
library(ggpubr)
library(RcolorBrewer)

setwd("~/data_HD/Project/TS/10XData/Figure1/")

my_theme =   theme_classic(base_line_size = 1)+
  theme(plot.title = element_text(size = 20,
                                  colour = "black",
                                  hjust = 0.5),
        axis.title.y = element_text(size = 15, 
                                    color = "black",
                                    face = "bold", 
                                    vjust = 1.9, 
                                    hjust = 0.5, 
                                    angle = 90),
        axis.title.x = element_text(size = 15, 
                                    color = "black",
                                    face = "bold", 
                                    vjust = 1.9, 
                                    hjust = 0.5, 
                                    angle = 90),
        legend.title = element_text(color="black", # 修改图例的标题
                                    size=15, 
                                    face="bold"),
        legend.text = element_text(color="black", # 设置图例标签文字
                                   size = 10, 
                                   face = "bold"),
        axis.text.x = element_text(size = 13, # 修改X轴上字体大小，
                                   color = "black",
                                   face = "bold",
                                   vjust = 0.5,
                                   hjust = 0.5,
                                   angle = 0),
        axis.text.y = element_text(size = 13, # 修改y轴上字体大小，
                                   color = "black",
                                   face = "bold", 
                                   vjust = 0.5,
                                   hjust = 0.5,
                                   angle = 0)
        #  face取值：plain普通，bold加粗，italic斜体bold.italic斜体加粗
  )

color_ct=c(brewer.pal(12, "Set3")[-c(2,3,9,12)],"#b3b3b3",
           brewer.pal(5, "Set1")[2],
           brewer.pal(3, "Dark2")[1],
           "#fc4e2a","#fb9a99","#f781bf","#e7298a")

merged_object = read_rds("combinedReannotation.rds")

IgD = subset(merged_object, subset = celltype == "IgD-")
naive = subset(merged_object, subset = celltype == "naive")

# ----------------------------boxplot-----------------------------------
naive.datasets = read_tsv("../../figs/fig2/correlation/naive_info.tsv")
IgD.datasets = read_tsv("../../figs/fig2/correlation/IgD-_info.tsv")

celltypes = list(
  naive = naive.datasets,
  `IgD-` = IgD.datasets
)

lapply(names(celltypes), function(cell){
  box.data = celltypes[[cell]] |> 
    mutate(special = ifelse(is.na(special), chr, special)) |>
    select(log2FC_TS2F, log2FC_TS2M, log2FC_TS3F, log2FC_TS3M, special) |>
    melt(id.vars = "special") |>
    mutate(variable = gsub("log2FC_", "", variable))
  colnames(box.data) = c("position", "sample", "exp_diff")
  box.data$position[box.data$position == "chrX"] = "other chrX"
  
  ggplot(box.data, aes(x=sample, y=exp_diff, color=position))+
    geom_boxplot(aes(fill=position),
                 alpha=0.1)+
    geom_jitter(position = position_jitterdodge(jitter.height = 0.75,
                                                jitter.width = 0.1,
                                                dodge.width = 0.75),
                alpha = 0)+
    scale_color_manual(values = pal_npg('nrc')(9)[c(1:length(unique(box.data$position)))])+
    scale_fill_manual(values = pal_npg('nrc')(9)[c(1:length(unique(box.data$position)))])+
    stat_compare_means(aes(group = position),
                       label = "p.signif",
                       show.legend = F,
                       hide.ns = T)+
    my_theme+
    labs(x="", y="log2FC")+
    scale_y_continuous(limits = c(-2.5, 2.5))+
    ggtitle(paste0("Differential expression between children and parents in ", cell))
  ggsave(paste0("../../figs/fig4/boxplot/boxplot_diff_", cell, ".pdf"), width = 9, height = 3)
})

# --------------------------get gene position------------------------------
gene_pos = read_tsv("~/ref/human/bed/UCSC_refseq_hg19.tsv.gz", col_names = T)
gene_pos.merge = gene_pos |> dplyr::filter(grepl("^chr[0-9XY]+$", gene_pos$chrom)) |> 
  group_by(chrom, name2) |> 
  summarise(chrom, start = min(txStart), end = max(txEnd), name2, score, strand) |>
  distinct() |> dplyr::select(chrom, start, end, name2, score, strand) |>
  arrange(chrom, start, end) |> ungroup()

# -----------------------------manhattan plot-------------------------------
gene_pos.man = gene_pos.merge
gene_pos.man$chrom = case_when(
  gene_pos.man$chrom == "chrX" ~ 23,
  gene_pos.man$chrom == "chrY" ~ 24,
  .default = as.numeric(gsub("chr", "", gene_pos.man$chrom))
)
gene_pos.man = gene_pos.man |> dplyr::select("chrom", "name2") |> arrange(chrom) |> 
  unique()
auto = intersect(rownames(merged_object), autosome_genes)

source("~/r_scrips/spearman_critical_value.R")

object.list = list(merged_object, naive, IgD)
names(object.list) = c("all", "naive", "IgD")
my_genes = c("SLC25A6", "CD99", "ZFX", "KDM5C")
samples = "all"

lapply(my_genes, function(gene) {
  print(paste0("Processing gene ", gene, " ..."))
  dir.create(paste0("../../figs/fig2/manhattan_volcano/", gene))
  lapply(samples, function(sp) {
    print(paste0("Processing sample ", sp, " ..."))
    dir.create(paste0("../../figs/fig2/manhattan_volcano/", gene, "/", sp))
    lapply(names(object.list), function(cl) {
      print(paste0("Processing cell ", cl, " ..."))
      dir.create(paste0("../../figs/fig2/manhattan_volcano/", gene, "/", sp, "/", cl))
      
      if (sp == "all") {
        obj = object.list[[cl]]
      } else {
        obj = subset(object.list[[cl]], sample == sp)
      }
      exp.df = GetAssayData(obj, assay = "RNA", slot = "data") |>
        as.data.frame()
      exp.target.df = exp.df[c(gene, autosome_genes), ] |> na.omit() |> t() |> 
        as.data.frame()
      
      exp.cor.df = exp.target.df |> select(-!!sym(gene)) |>
        summarise(across(everything(), ~{
          test = cor.test(exp.target.df[[gene]], ., method = "spearman")
          tibble(Correlation = test$estimate, P_value = test$p.value)
        }, .names = "{.col}")) |>
        pivot_longer(cols = everything(), names_to = "Gene", values_to = "Results") |>
        tidyr::unnest(Results)
      exp.cor.df = na.omit(exp.cor.df)
      exp.cor.df = exp.cor.df |> mutate(adj_p = p.adjust(P_value, method = "bonferroni"))
      exp.cor.df$`-log10p` = -log10(exp.cor.df$adj_p)
      
      exp.cor.df$state = case_when(
        exp.cor.df$adj_p < 0.05 & exp.cor.df$Correlation > 0 ~ "pos_cor",
        exp.cor.df$adj_p < 0.05 & exp.cor.df$Correlation < 0 ~ "neg_cor",
        .default = "no_cor"
      )
      exp.cor.df$`-log10p`[is.infinite(exp.cor.df$`-log10p`)] = 
        max(exp.cor.df$`-log10p`[!is.infinite(exp.cor.df$`-log10p`)])
      # ---------------------Volcano Plot-----------------------
      # 过滤上下调top5的基因：
      top_m <- exp.cor.df %>% top_n(n = 10, wt = abs(Correlation)) %>%  # 如果你想改成top10，把5替换成对应的数字就好
        filter(state !='no_cor')
      net = top_m |> mutate(Source = gene, Target = Gene) |> 
        select(Source, Target, Correlation)
      net = inner_join(net, gene_pos.merge, by = join_by(Target == name2))
      write_tsv(net, paste0("../../figs/fig2/manhattan_volcano/", gene, "/", sp, "/", cl, "/table4net_top10_", gene, "_", sp, "_", cl, ".tsv"))
      
      top_50 <- exp.cor.df %>% top_n(n = 50, wt = abs(Correlation)) %>%  # 如果你想改成top10，把5替换成对应的数字就好
        filter(state !='no_cor')
      net = top_50 |> mutate(Source = gene, Target = Gene) |> 
        select(Source, Target, Correlation)
      net = inner_join(net, gene_pos.merge, by = join_by(Target == name2))
      write_tsv(net, paste0("../../figs/fig2/manhattan_volcano/", gene, "/", sp, "/", cl, "/table4net_top50_", gene, "_", sp, "_", cl, ".tsv"))
      
      ylim = max(exp.cor.df$`-log10p`) * 1.1
      xlim = max(abs(exp.cor.df$Correlation)) * 1.1
      ggplot(exp.cor.df, aes(x=Correlation, y=`-log10p` )) +
        geom_point(aes(color=state)) +
        geom_abline(intercept = 0,slope = 0,linetype="dashed") +
        theme_bw() +
        ylim(0,ylim) +
        xlim(-xlim,xlim) +
        scale_color_manual(values=c("blue" ,"gray", "red"))+
        geom_text_repel(data=top_m,aes(x=Correlation,y=`-log10p`,label=Gene), size = 3)
      ggsave(paste0("../../figs/fig2/manhattan_volcano/", gene, "/", sp, "/", cl, "/correlation_vlocano_", gene, "_", sp, "_", cl, ".pdf"), width = 5, height = 4)
      
      man.df = exp.cor.df |> select(Gene, adj_p)
      gene_pos.man = gene_pos.man |> subset(name2 %in% man.df$Gene)
      man.df = inner_join(man.df, gene_pos.man, by = c("Gene" = "name2"))
      colnames(man.df) = c("SNP", "P", "CHR")
      man.df = man.df |> group_by(CHR) |> mutate(BP = row_number()) |> ungroup()
      man.df$P[man.df$P == 0] = min(man.df$P[man.df$P > 0])
      
      cutoff = man.df |> top_n(n = -10, wt = P) |> summarise(max(P))
      
      pdf(paste0("../../figs/fig2/manhattan_volcano/", gene, "/", sp, "/", cl, "/manhatan_p_", gene, "_", sp, "_", cl, ".pdf"), width = 8, height = 4)
      manhattan(man.df,
                main = paste0("Manhattan Plot-correlation-p value-", sp, "-", gene, "-", cl),
                ylim = c(0,ylim),
                cex = 0.8,
                cex.axis = 0.9,
                col = color_ct,
                suggestiveline = F, genomewideline = -log10(0.05),
                annotatePval = cutoff$`max(P)`, annotateTop = F)
      dev.off()
      
      man.df = exp.cor.df |> select(Gene, Correlation)
      man.df = inner_join(man.df, gene_pos.man, by = c("Gene" = "name2"))
      colnames(man.df) = c("SNP", "P", "CHR")
      man.df = man.df |> group_by(CHR) |> mutate(BP = row_number())
      
      r_critical = spearman_critical_value(df = ncol(obj) - 2, alpha = 0.05/nrow(exp.cor.df))
      
      pdf(paste0("../../figs/fig2/manhattan_volcano/", gene, "/", sp, "/", cl, "/manhatan_r_", gene, "_", sp, "_", cl, ".pdf"), width = 8, height = 4)
      manhattan(man.df,
                main = paste0("Manhattan Plot-correlation-spearman r-", sp, "-", gene, "-", cl),
                ylim = c(-xlim,xlim), logp = F,
                cex = 0.8,
                cex.axis = 0.9,
                col = color_ct,
                suggestiveline = -r_critical, genomewideline = r_critical)
      dev.off()
    })
  })
})