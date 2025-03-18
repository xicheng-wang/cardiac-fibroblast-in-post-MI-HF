# setwd('e:/work/hd/scRNA/')

### 补充做个拟时序，可以和前面的内容放在一起
scRNAsub=readRDS('./scRNA_fibroblast.RDS')


## 拟时序分析 先关掉再开    suggest 4.1.3
library(Seurat)
#没有monocle要先安装 BiocManager::install
library(monocle)
library(tidyverse)
library(patchwork)
data=as.matrix(scRNAsub@assays$RNA@counts)
# data=as.matrix(scRNAsub@assays$RNA@scale.data)
data <- as(data, 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = scRNAsub@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)


dim(data)

## 以下代码一律不得修改
mycds <- newCellDataSet(data,
                        phenoData = pd,
                        featureData = fd,
                        expressionFamily = negbinomial.size())

mycds <- estimateSizeFactors(mycds)
mycds <- estimateDispersions(mycds, cores=4, relative_expr = TRUE)


##使用monocle选择的高变基因，不修改
disp_table <- dispersionTable(mycds)
disp.genes <- subset(disp_table, mean_expression >= 0.5 &
                       dispersion_empirical >= 1 * dispersion_fit)$gene_id
mycds <- setOrderingFilter(mycds, disp.genes)
plot_ordering_genes(mycds)


#降维
mycds <- reduceDimension(mycds, max_components = 2, method = 'DDRTree')
#排序
mycds <- orderCells(mycds)


table(mycds$tissue_type)

#State轨迹分布图
plot1 <- plot_cell_trajectory(mycds, color_by = "tissue_type")+
  facet_wrap(~tissue_type, nrow = 1)+
  scale_color_manual(values = c("#324c63", "#a363ae"))+theme(legend.position = "none")
plot1

##Pseudotime轨迹图
plot2 <- plot_cell_trajectory(mycds, color_by = "Pseudotime")+theme(legend.position = "right")+
  scale_color_gradientn(colours = c('white','red'))
plot2

col_vector <- readRDS('./col_vector.RDS')
plot4 <- plot_cell_trajectory(mycds, color_by = "seurat_clusters",cell_size = 0.5)+
  scale_color_manual(values = col_vector)+theme(legend.position = "right")
dev.off()
plot4


dir.create("./figures")
library(patchwork)
pdf(file="./figures/拼图.pdf", width=14, height=4)
# print(plot1|plot2|plot4)

((plot1|plot2|plot4)+plot_layout(widths = c(2,1,1)))
dev.off()



# ?plot_cell_trajectory

seurat_obj=readRDS('./hd_day2_hdWGCNA/hdWGCNA_object.rds')



rm(list=ls())
options(stringsAsFactors = F)
library(hdWGCNA)
library(tidyverse)

# # 可以根据自己的需要扩大个数
# hub_df <- GetHubGenes(seurat_obj, n_hubs = 50)
# ### !!!!!
# gene=hub_df[hub_df$module %in% c('blue','pink'),]
# gene=gene$gene_name

# gene
# gene=readRDS('./hd_day2_hdWGCNA/figures/module_turquoise_150genes.RDS')
# gene

markers <- readRDS('./hd_day2_hdWGCNA/figures/markers.RDS')

table(markers$cluster)
markers <- markers[markers$cluster %in% c("Cardioprotection-related_fibroblasts"),]
markers$gene





# marker13 <- readRDS('../82 OptimalML-101种组合算法建模 二分类/Figures/venn.RDS')
# marker13
hub_df <- readRDS('./hd_day2_hdWGCNA/figures/day2_hdWGCNA_hub_df.RDS')
hub_df <- hub_df[hub_df$module %in% c("turquoise"),]
hub_df1 <- hub_df[hub_df$gene_name %in% markers$gene, ]


library(openxlsx)
write.xlsx(hub_df1, "./figures/hub_df_turquoise_markers_genes.xlsx",  rowNames =T)
saveRDS(hub_df1, file = "./figures/hub_df_turquoise_markers_genes.RDS")




gene=readRDS('./hd_day2_hdWGCNA/figures/markers.RDS')

gene <- readRDS('figures/hub_df_turquoise_markers_genes.RDS')
gene <- gene$gene_name
gene

library(RColorBrewer)
my_pseudotime_cluster <- plot_pseudotime_heatmap(mycds[gene,],
                                                 # num_clusters = 2, # add_annotation_col = ac,
                                                 show_rownames = TRUE,
                                                 hmcols = colorRampPalette(rev(brewer.pal(9, "PiYG")))(62),
                                                 return_heatmap = TRUE)

my_pseudotime_cluster

# plot4 <- plot_pseudotime_heatmap(HSMM_1[gene_to_cluster,], num_clusters=5,
                                 # show_rownames=T, return_heatmap=T)

ggsave(my_pseudotime_cluster,filename = "figures/monocle2_heatmap.pdf",units = "cm",width = 10,height = 10)



##### 准备富集和PPI 需要的gene list
library(hdWGCNA)
library(tidyverse)
# seurat_obj=readRDS('./hd_day2_hdWGCNA/hdWGCNA_object.rds')
#
# # 可以根据自己的需要扩大个数
# hub_df <- GetHubGenes(seurat_obj, n_hubs = 150)
# table(hub_df$module)
# gene=hub_df[hub_df$module %in% c('pink','blue'),]
# gene=gene$gene_name
# gene

### 人不用运行！！！
library(Hmisc)
gene=toupper(gene)
gene=capitalize(gene)

gene

write.csv(gene,file ='./figures/hdWGCNA_gene_for_PPI.csv',quote = F)

gc()

##############
FeaturePlot(scRNAsub,features = 'CCL2',split.by = 'tissue_type', order = T)
FeaturePlot(scRNAsub,features = 'THBS1',split.by = 'tissue_type', order = T)
FeaturePlot(scRNAsub,features = 'FGF7',split.by = 'tissue_type', order = T)


FeaturePlot(scRNAsub,features = 'CCL2',split.by = 'tissue_type')
FeaturePlot(scRNAsub,features = 'THBS1',split.by = 'tissue_type')
FeaturePlot(scRNAsub,features = 'FGF7',split.by = 'tissue_type')

library("Nebulosa")
plot_density(scRNAsub, "FGF7",reduction = "umap")
plot_density(scRNAsub, "THBS1",reduction = "umap")
plot_density(scRNAsub, "CCL2",reduction = "umap")



#https://metascape.org/gp/index.html
#http://genemania.org/


##### 药物预测，分子对接准备，提上来讲###########

#https://ctdbase.org/
# 药物-蛋白-疾病
dev.off()

##在卒中中，ARG1的抑制剂

# setwd('e:/work/hd')
chemi=read.csv(file ='CTD_383_diseases_20230101092921.csv',check.names = F)
chemi=chemi[chemi$`Disease Name`=='Stroke',]
chemi=chemi$`Inference Network`
chemi=stringr::str_split(chemi,pattern = '\\|')
chemi=chemi[[1]]

intera=read.csv(file ='CTD_383_ixns_20230101093249.csv',check.names = F)

intera=intera[intera$`Chemical Name` %in% chemi,]


### https://www.rcsb.org/
### https://pubchem.ncbi.nlm.nih.gov/
## https://www.dockeasy.cn/DockCompound


# parameter_file AD4_parameters.dat
