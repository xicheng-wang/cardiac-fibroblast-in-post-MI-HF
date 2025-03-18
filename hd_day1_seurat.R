
# setwd('e:/work/hd/scRNA/')
#setwd('./scRNA/')
gc()
#BiocManager::install('Seurat',update = F,ask = F)
library(Seurat)

table(PHx_hepatocytes$orig.ident)

# ###1. 普通表达矩阵——标准表达矩阵（行为基因，列为样本（细胞））
# ###2.三联文件
sham1=Read10X(data.dir = './phx48h/Adult/')
sham2=Read10X(data.dir = './phx48h/48h/')
MCAO1=Read10X(data.dir = './phx48h/PHx_0h_onlyhepatocytes/')  #6
MCAO2=Read10X(data.dir = './phx48h/PHx_48h_onlyhepatocytes/')
#
# #sham3=read.table('...',row.names = 1,check.names = F,sep = '\t')
#
sham1=CreateSeuratObject(sham1,project='Adult_0h_1', min.cells = 5,
                         min.features = 1000 )
sham2=CreateSeuratObject(sham2,project='PHx_48h_1',min.cells = 5,
                         min.features = 1000)
MCAO1=CreateSeuratObject(MCAO1,project='Adult_0h_2',min.cells = 5,
                         min.features = 1000)
MCAO1
MCAO2=CreateSeuratObject(MCAO2,project='PHx_48h_2',min.cells = 5,
                         min.features = 1000)
#
#
# # 合并
scRNA=merge(sham1,c(sham2,MCAO1,MCAO2))



rm(MCAO1,MCAO2,sham1,sham2)
gc()
scRNA
table(scRNA$orig.ident)
##### 取一列命名为组织类型 #####
scRNA$tissue_type=stringr::str_remove(scRNA$orig.ident,pattern = '_[1-2]')
table(scRNA$tissue_type)
scRNA$tissue_type=factor(scRNA$tissue_type,levels = c('Adult_0h','PHx_48h'))
length(scRNA@active.ident)

saveRDS(scRNA, file = "PHX_0h_48h.RDS")

# pbmc <- readRDS("./")


pbmc <- readRDS('scRNA_anno_正确origident.RDS')
pbmc <- readRDS('../scRNA_harmony.RDS')

pbmc <- ifnb.data
Idents(pbmc) <- "orig.ident"
# Idents(pbmc) <- "tissue_type"
table(Idents(pbmc))

rownames(pbmc)

# pbmc <- OB20220701_result

# mito_genes=rownames(pbmc)[grep("^MT-", rownames(pbmc))]
# mito_genes=rownames(pbmc)[grep("^Mt-", rownames(pbmc))]
mito_genes=rownames(pbmc)[grep("^mt-", rownames(pbmc))]
mito_genes #13个线粒体基因
pbmc=PercentageFeatureSet(pbmc, "^mt-", col.name = "percent_mito")
# pbmc=PercentageFeatureSet(pbmc, "^MT-", col.name = "percent_mito")
fivenum(pbmc@meta.data$percent_mito)
#计算核糖体基因比例
# ribo_genes=rownames(pbmc)[grep("^RP[sl]", rownames(pbmc),ignore.case = T)]
ribo_genes=rownames(pbmc)[grep("^Rps|^Rpl", rownames(pbmc),ignore.case = T)]
ribo_genes

### 计算线粒体等基因的比例
# pbmc <- PercentageFeatureSet(pbmc, pattern = "^MT-", col.name = "pMT")
pbmc <- PercentageFeatureSet(pbmc, pattern = "^mt-", col.name = "pMT")
# pbmc <- PercentageFeatureSet(pbmc, pattern = "^HBA|^HBB", col.name = "pHB")
pbmc <- PercentageFeatureSet(pbmc, pattern = "^Hba|^Hbb", col.name = "pHB")
# pbmc <- PercentageFeatureSet(pbmc, pattern = "^RPS|^RPL", col.name = "pRP")
pbmc <- PercentageFeatureSet(pbmc, pattern = "^Rps|^Rpl", col.name = "pRP")

qcparams <- c("nFeature_RNA", "nCount_RNA", "pMT", "pHB", "pRP")
for (i in seq_along(qcparams)){
  print(VlnPlot(object = pbmc, features = qcparams[i], group.by = "orig.ident", pt.size = 0))
}
for (i in seq_along(qcparams)){
  print(RidgePlot(object = pbmc, features = qcparams[i], group.by = "orig.ident"))
}
### 批量看计算结果，很不错
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "pMT"), pt.size = 0, group.by = "orig.ident", ncol = 1, log = T)
VlnPlot(pbmc, features = c("pMT"), pt.size = 0, group.by = "orig.ident", ncol = 1, log = T)
library(ggplot2)
# ggsave("QC.pdf", path = "./", width = 20, height = 20, units = "cm")




#####################预处理  质控 批量出图

# 加载R包
library(Seurat)
library(dplyr)
library(reticulate)
library(sctransform)
library(cowplot)
library(ggplot2)
library(viridis)
library(tidyr)
library(magrittr)
library(reshape2)
library(readxl)
library(progeny)
library(readr)
library(stringr)

#设定阈值
nFeature_lower <- 500
nFeature_upper <- 10000
nCount_lower <- 1000
nCount_upper <- 100000
pMT_lower <- 0
pMT_upper <- 30
pHB_lower <- 0
pHB_upper <- 5

theme_set(theme_cowplot())

#设定配色方案
use_colors <- c(
  Tumor = "brown2",
  Normal = "deepskyblue2",
  G1 = "#46ACC8",
  G2M = "#E58601",
  S = "#B40F20",
  Epithelial = "seagreen",
  Immune = "darkgoldenrod2",
  Stromal = "steelblue",
  p018 = "#E2D200",
  p019 = "#46ACC8",
  p023 = "#E58601",
  p024 = "#B40F20",
  p027 = "#0B775E",
  p028 = "#E1BD6D",
  p029 = "#35274A",
  p030 = "#F2300F",
  p031 = "#7294D4",
  p032 = "#5B1A18",
  p033 = "#9C964A",
  p034 = "#FD6467")



# 过滤
# 定义两个过滤作图函数
qc_std_plot_helper <- function(x) x +
  # scale_color_viridis() +
  # scale_colour_gradient(low="#D54900",high="#4A4487") + #蓝橙色
  # scale_colour_gradient(low="#2e8cfa",high="#ff6074") +#prism8的颜色，挺好的
  # scale_color_gradientn(colours = rev(RColorBrewer::brewer.pal(8,"RdYlBu")))+
  scale_color_gradientn(colours = rev(RColorBrewer::brewer.pal(9,"PRGn")))+
  geom_point(size = 0.01, alpha = 0.3)

qc_std_plot <- function(pbmc) {
  qc_data <- as_tibble(FetchData(pbmc, c("nCount_RNA", "nFeature_RNA", "pMT", "pHB", "pRP")))
  plot_grid(

    qc_std_plot_helper(ggplot(qc_data, aes(log2(nCount_RNA), log2(nFeature_RNA), color = pMT))) +
      geom_hline(yintercept = log2(nFeature_lower), color = "red", linetype = 2) +
      geom_hline(yintercept = log2(nFeature_upper), color = "red", linetype = 2) +
      geom_vline(xintercept = log2(nCount_lower), color = "red", linetype = 2) +
      geom_vline(xintercept = log2(nCount_upper), color = "red", linetype = 2),
    qc_std_plot_helper(ggplot(qc_data, aes(log2(nCount_RNA), log2(nFeature_RNA), color = pHB))) +
      geom_hline(yintercept = log2(nFeature_lower), color = "red", linetype = 2) +
      geom_hline(yintercept = log2(nFeature_upper), color = "red", linetype = 2) +
      geom_vline(xintercept = log2(nCount_lower), color = "red", linetype = 2) +
      geom_vline(xintercept = log2(nCount_upper), color = "red", linetype = 2),
    qc_std_plot_helper(ggplot(qc_data, aes(log2(nCount_RNA), log2(nFeature_RNA), color = pRP))) +
      geom_hline(yintercept = log2(nFeature_lower), color = "red", linetype = 2) +
      geom_hline(yintercept = log2(nFeature_upper), color = "red", linetype = 2) +
      geom_vline(xintercept = log2(nCount_lower), color = "red", linetype = 2) +
      geom_vline(xintercept = log2(nCount_upper), color = "red", linetype = 2),

    qc_std_plot_helper(ggplot(qc_data, aes(log2(nCount_RNA), pMT, color = nFeature_RNA))) +
      geom_hline(yintercept = pMT_lower, color = "red", linetype = 2) +
      geom_hline(yintercept = pMT_upper, color = "red", linetype = 2) +
      geom_vline(xintercept = log2(nCount_lower), color = "red", linetype = 2) +
      geom_vline(xintercept = log2(nCount_upper), color = "red", linetype = 2),
    qc_std_plot_helper(ggplot(qc_data, aes(log2(nCount_RNA), pHB, color = nFeature_RNA))) +
      geom_hline(yintercept = pHB_lower, color = "red", linetype = 2) +
      geom_hline(yintercept = pHB_upper, color = "red", linetype = 2) +
      geom_vline(xintercept = log2(nCount_lower), color = "red", linetype = 2) +
      geom_vline(xintercept = log2(nCount_upper), color = "red", linetype = 2),
    qc_std_plot_helper(ggplot(qc_data, aes(log2(nCount_RNA), pRP, color = nFeature_RNA))) +
      geom_vline(xintercept = log2(nCount_lower), color = "red", linetype = 2) +
      geom_vline(xintercept = log2(nCount_upper), color = "red", linetype = 2),


    qc_std_plot_helper(ggplot(qc_data, aes(log2(nFeature_RNA), pMT, color = nCount_RNA))) +
      geom_hline(yintercept = pMT_lower, color = "red", linetype = 2) +
      geom_hline(yintercept = pMT_upper, color = "red", linetype = 2) +
      geom_vline(xintercept = log2(nFeature_lower), color = "red", linetype = 2) +
      geom_vline(xintercept = log2(nFeature_upper), color = "red", linetype = 2),
    qc_std_plot_helper(ggplot(qc_data, aes(log2(nFeature_RNA), pHB, color = nCount_RNA))) +
      geom_hline(yintercept = pHB_lower, color = "red", linetype = 2) +
      geom_hline(yintercept = pHB_upper, color = "red", linetype = 2) +
      geom_vline(xintercept = log2(nFeature_lower), color = "red", linetype = 2) +
      geom_vline(xintercept = log2(nFeature_upper), color = "red", linetype = 2),
    qc_std_plot_helper(ggplot(qc_data, aes(log2(nFeature_RNA), pRP, color = nCount_RNA))) +
      geom_vline(xintercept = log2(nFeature_lower), color = "red", linetype = 2) +
      geom_vline(xintercept = log2(nFeature_upper), color = "red", linetype = 2),

    qc_std_plot_helper(ggplot(qc_data, aes(pRP, pMT, color = nCount_RNA))) +
      geom_hline(yintercept = pMT_lower, color = "red", linetype = 2) +
      geom_hline(yintercept = pMT_upper, color = "red", linetype = 2),
    qc_std_plot_helper(ggplot(qc_data, aes(pRP, pMT, color = nFeature_RNA))) +
      geom_hline(yintercept = pMT_lower, color = "red", linetype = 2) +
      geom_hline(yintercept = pMT_upper, color = "red", linetype = 2),


    ggplot(gather(qc_data, key, value), aes(key, value)) +
      geom_violin() +
      facet_wrap(~key, scales = "free", ncol = 3),

    ncol = 3, align = "hv"
    # labels = "auto"
    # labels = "AUTO" #大写的ABCD
  )
}

# ?plot_grid
## 过滤前
pbmc
pbmc_unfiltered <- pbmc

pbmc_unfiltered
# 作图，有点慢，但是一般可以出来
qc_std_plot(pbmc_unfiltered)

dir.create("./filter_qc")
ggsave2("./filter_qc/before_filtering.pdf", path = "./", width = 40, height = 33, units = "cm")
ggsave2("./filter_qc/before_filtering.png", path = "./", width = 40, height = 33, units = "cm")


pbmc

## 过滤
pbmc <- subset(pbmc_unfiltered, subset = nFeature_RNA > nFeature_lower &
                 nFeature_RNA < nFeature_upper & nCount_RNA > nCount_lower &
                 nCount_RNA < nCount_upper & pMT < pMT_upper & pHB < pHB_upper)


pbmc
#有点慢，同样可以出
qc_std_plot(pbmc)
ggsave2("./filter_qc/after_filtering2.pdf", path = "./", width = 40, height = 33, units = "cm")
ggsave2("./filter_qc/after_filtering2.png", path = "./", width = 40, height = 33, units = "cm")

# 看过滤了多少
pbmc_unfiltered
pbmc

table(pbmc@meta.data$orig.ident)


# save(pbmc,pbmc_unfiltered, file = "./filter_qc/iPSCs.Rdata")
i <- length(rownames(pbmc@meta.data))
# dir.create("./iPSCs_所有数据")
# save(pbmc,pbmc_unfiltered, file = paste0("./filter_qc/","GSE193852DCM_filter_",i,".Rdata"))
# save(pbmc,pbmc_unfiltered, file = paste0("./filter_qc/","GSE193852DCM_filter_",i,".Rdata"))
pbmc_unfiltered
saveRDS(pbmc_unfiltered,file = paste0("./filter_qc/","pbmc_unfilter_",i,".RDS"))
saveRDS(pbmc,file = paste0("./filter_qc/","pbmc_filter_",i,".RDS"))
# save(pbmc, file = paste0("./filter_qc/","pbmc_filter_",i,".Rdata"))

pbmc


scRNA <- pbmc
scRNA
##### 降维聚类 #####
scRNA <- NormalizeData(scRNA, normalization.method = "LogNormalize", scale.factor = 10000)
scRNA<- FindVariableFeatures(scRNA, selection.method = "vst", nfeatures = 2000)
scale.genes <-  VariableFeatures(scRNA)
scRNA <- ScaleData(scRNA, features = scale.genes)
scRNA<- RunPCA(scRNA, features = VariableFeatures(scRNA))
DimPlot(scRNA, reduction = "pca", group.by = "orig.ident")

### harmony去批次
#BiocManager::
library(harmony)
scRNA_harmony <- RunHarmony(scRNA, group.by.vars = "orig.ident")
DimPlot(scRNA_harmony, reduction = "harmony", group.by = "orig.ident")

ElbowPlot(scRNA_harmony,reduction = 'harmony')

# 一定要指定“harmony”！
scRNA <- FindNeighbors(scRNA_harmony, dims = 1:5, reduction = "harmony")
scRNA <- FindClusters(scRNA)
scRNA <- RunUMAP(scRNA, dims = 1:5,reduction = 'harmony')








scRNA <- ifnb.data
scRNA

# 去批次成功
DimPlot(scRNA,split.by = 'tissue_type')

rm(scRNA_harmony)

#BiocManager
library(SingleR)
# 人用下面！！
refdata <- SingleR::HumanPrimaryCellAtlasData()
# refdata <- SingleR::MouseRNAseqData()
library(Seurat)
testdata <- GetAssayData(scRNA, slot="data")
clusters <- scRNA@meta.data$seurat_clusters
cellpred <- SingleR(test = testdata, ref = refdata,
                    # labels =refdata$label.fine,
                    labels =refdata$label.main,
                    method = "cluster", clusters = clusters,
                    assay.type.test = "logcounts", assay.type.ref = "logcounts")

celltype = data.frame(ClusterID=rownames(cellpred), celltype=cellpred$labels, stringsAsFactors = F)
scRNA@meta.data$celltype = "NA"
for(i in 1:nrow(celltype)){
  scRNA@meta.data[which(scRNA@meta.data$seurat_clusters == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}


Idents(scRNA)=scRNA$celltype

library(RColorBrewer)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
#处理后有73种差异还比较明显的颜色，基本够用
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors,
                           rownames(qual_col_pals)))


table(scRNA$orig.ident)


library(stringr)
result=str_split(scRNA$orig.ident,'-',simplify = T)
result=as.data.frame(result)
table(result$V1)
result$V1

scRNA$tissue_type <- result$V1
scRNA$tissue_type2 <- paste0(result$V1,'-',result$V2)

dev.off()
p_dim <- DimPlot(scRNA,split.by = 'tissue_type',cols = col_vector[1:19],label=T,repel = T)+NoLegend()
p_dim


p_dim <- DimPlot(scRNA,split.by = 'tissue_type2',cols = col_vector[1:19],label=T,repel = T)+NoLegend()

p_dim


##########保存
saveRDS(scRNA,file ='scRNA_anno.RDS')
gc()


##########################################  singleR注释效果差，换azimuth        ############################
############################ ##########################################



##########################################  azimuth             ############################
############################ ##########################################
pbmc1

# https://app.azimuth.hubmapconsortium.org/app/human-heart
# 一般不能超过2.5w细胞数

scRNA <- readRDS('../作曲家/scRNA_anno.RDS')
table(scRNA$tissue_type)
pbmc_normal <-subset(scRNA, tissue_type %in% c("N"))

pbmc_normal
a <- pbmc_normal@assays$RNA@counts
a
saveRDS(a, file = "pmbc_count.rds")



scRNA <- readRDS('../作曲家/scRNA_anno.RDS')
table(scRNA$tissue_type, scRNA$orig.ident)
table(scRNA$tissue_type2)
pbmc_icm3 <-subset(scRNA, tissue_type2 %in% c("ICM-3"))

pbmc_icm3
a <- pbmc_icm3@assays$RNA@counts
a
saveRDS(a, file = "pmbc_count_icm3.rds")



scRNA <- readRDS('作曲家/scRNA_anno.RDS')
table(scRNA$tissue_type, scRNA$orig.ident)
table(scRNA$tissue_type2)
pbmc_icm3 <-subset(scRNA, tissue_type2 %in% c("ICM-2","ICM-1"))

pbmc_icm3
a <- pbmc_icm3@assays$RNA@counts
a
saveRDS(a, file = "pmbc_count_icm1and2.rds")

https://app.azimuth.hubmapconsortium.org/app/human-heart



pbmc
scRNA
predictions <- read.delim('../azimuth_pred_all.tsv', row.names= 1)


scRNA <- AddMetaData(scRNA,
                     metadata = predictions)
colnames(scRNA@meta.data)
DimPlot(scRNA, reduction="umap", label = T,
        group.by="predicted.celltype.l2")

table(scRNA@meta.data$predicted.celltype.l2)
table(scRNA@meta.data$predicted.celltype.l2)
table(scRNA@meta.data$celltype)



scRNA$celltype <-scRNA@meta.data$predicted.celltype.l2

scRNA$celltype <- gsub("Arterial Endothelial","Endothelial_cells",scRNA$celltype)
scRNA$celltype <- gsub("Capillary Endothelial","Endothelial_cells",scRNA$celltype)
scRNA$celltype <- gsub("Venous Endothelial","Endothelial_cells",scRNA$celltype)
scRNA$celltype <- gsub("Macrophage","Macrophages",scRNA$celltype)
scRNA$celltype <- gsub("Monocyte/cDC","Macrophages",scRNA$celltype)
scRNA$celltype <- gsub("Ventricular Cardiomycoyte","Cardiomycoytes",scRNA$celltype)
scRNA$celltype <- gsub("B","B_cells",scRNA$celltype)
scRNA$celltype <- gsub("Endocardial","Endocardial_cells",scRNA$celltype)
scRNA$celltype <- gsub("Fibroblast","Fibroblasts",scRNA$celltype)
scRNA$celltype <- gsub("Lymphatic Endothelial","Lymphatic_Endothelial_cells",scRNA$celltype)
scRNA$celltype <- gsub("NK","NK_cells",scRNA$celltype)
scRNA$celltype <- gsub("Pericyte","Pericytes",scRNA$celltype)
scRNA$celltype <- gsub("Smooth Muscle","Smooth_Muscle_cells",scRNA$celltype)
scRNA$celltype <- gsub("T","T_cells",scRNA$celltype)
scRNA$celltype <- gsub("Mast","Mast_cells",scRNA$celltype)
scRNA$celltype <- gsub("Mesothelial","Mesothelial_cells",scRNA$celltype)
scRNA$celltype <- gsub("Neuronal","Neuronal_cells",scRNA$celltype)
scRNA$celltype <- gsub("","",scRNA$celltype)
scRNA$celltype <- gsub("","",scRNA$celltype)
scRNA$celltype <- gsub("","",scRNA$celltype)
scRNA$celltype <- gsub("","",scRNA$celltype)

DimPlot(scRNA, reduction="umap", label = T, repel = T,
        group.by="celltype")



dev.off()
table(scRNA$tissue_type)
scRNA$tissue_type <- gsub("N","Normal",scRNA$tissue_type)
scRNA$tissue_type <- gsub("Normalormal","Normal",scRNA$tissue_type)
scRNA$tissue_type <- gsub("ICM","HF",scRNA$tissue_type)
scRNA$tissue_type <- factor(scRNA$tissue_type, levels = c("Normal", "HF"))


table(scRNA$tissue_type2)
scRNA$tissue_type2 <- gsub("N","Normal",scRNA$tissue_type2)
# scRNA$tissue_type2 <- gsub("Normalormal","Normal",scRNA$tissue_type2)
scRNA$tissue_type2 <- gsub("ICM","HF",scRNA$tissue_type2)
scRNA$tissue_type2 <- factor(scRNA$tissue_type2, levels = c("Normal-1", "HF-1","HF-2","HF-3"))



p_dim1 <- DimPlot(scRNA,split.by = 'tissue_type',cols = col_vector[10:24],label=T,repel = T)+NoLegend()
p_dim1

levels(scRNA)



library(RColorBrewer)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
#处理后有73种差异还比较明显的颜色，基本够用
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors,
                           rownames(qual_col_pals)))

col_vector <- col_vector[10:24]
col_vector
# col_vector <- c(col_vector[1:10], col_vector[15],col_vector[11:14])
col_vector
saveRDS(col_vector, file = "col_vector.RDS")

col_vector <- readRDS('col_vector.RDS')
p_dim2 <- DimPlot(scRNA,split.by = 'tissue_type2',cols = col_vector,label=F,repel = T)+theme_bw()+NoLegend()

p_dim2






saveRDS(scRNA, file = "../作曲家/scRNA_anno.RDS")






scRNA=readRDS('scRNA_anno.RDS')



# 读取单细胞数据集，读取自己的


#数据是从服务器拿的
# scRNA=readRDS('scRNA_anno_衰老心脏_用这个.RDS')


library(Seurat)


### 画比例图###################
library(reshape2)
library(ggplot2)
library(dplyr)
pB2_df <- table(scRNA@meta.data$celltype,scRNA@meta.data$tissue_type) %>% melt()
colnames(pB2_df) <- c("Cluster","Sample","Number")
pB2_df$Cluster <- factor(pB2_df$Cluster)


pB2_df[1:3, 1:3]
data_cellcount1 <- pB2_df[1:length(unique(scRNA@meta.data$celltype)),]
data_cellcount2 <- pB2_df[(length(unique(scRNA@meta.data$celltype))+1):nrow(pB2_df),]
data_cellcount <- cbind(data_cellcount1, data_cellcount2)
data_cellcount$differcount <- data_cellcount[,3]-data_cellcount[,6]


factor(scRNA$celltype)
a <- factor(unique(scRNA$celltype))
a
pB2_df$Cluster <- factor(pB2_df$Cluster,levels = a)
# ## 如果要颜色一样
# pB2_df$Cluster <- factor(pB2_df$Cluster,levels = c('Edo...',
#                                                   'B_cell',
#                                                   'Monocyte',
#                                                   'T_cells',
#                                                   'Epithelial_cells',
#                                                   'Smooth_muscle_cells',
#                                                   'Fibroblasts',
#                                                   'Endothelial_cells'))

library(RColorBrewer)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
#处理后有73种差异还比较明显的颜色，基本够用
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
col_vector

sample_color <- col_vector[1:9]


col_vector <- readRDS('col_vector.RDS')
pB4 <- ggplot(data = pB2_df, aes(x =Number, y = Sample, fill =  Cluster)) +
  geom_bar(stat = "identity", width=0.8,position="fill")+
  scale_fill_manual(values=col_vector[1:20]) +
  theme_bw()+
  theme(panel.grid =element_blank()) +
  labs(x="",y="Ratio")+
  ####用来将y轴移动位置
  theme(axis.text.y = element_text(size=12, colour = "black"))+
  theme(axis.text.x = element_text(size=12, colour = "black"))+
  theme(
    axis.text.x.bottom = element_text(hjust = 1, vjust = 1, angle = 45)
  )
pB4


table(scRNA$orig.ident)
table(scRNA$tissue_type2)
# pB2_df <- table(scRNA@meta.data$celltype,scRNA@meta.data$orig.ident) %>% melt()
pB2_df <- table(scRNA@meta.data$celltype,scRNA@meta.data$tissue_type2) %>% melt()
colnames(pB2_df) <- c("Cluster","Sample","Number")
pB2_df$Cluster <- factor(pB2_df$Cluster)


pB2_df[1:3, 1:3]
data_cellcount1 <- pB2_df[1:length(unique(scRNA@meta.data$celltype)),]
data_cellcount2 <- pB2_df[(length(unique(scRNA@meta.data$celltype))+1):nrow(pB2_df),]
data_cellcount <- cbind(data_cellcount1, data_cellcount2)
data_cellcount$differcount <- data_cellcount[,3]-data_cellcount[,6]


factor(scRNA$celltype)
a <- factor(unique(scRNA$celltype))
a
pB2_df$Cluster <- factor(pB2_df$Cluster,levels = a)
# ## 如果要颜色一样
# pB2_df$Cluster <- factor(pB2_df$Cluster,levels = c('Edo...',
#                                                   'B_cell',
#                                                   'Monocyte',
#                                                   'T_cells',
#                                                   'Epithelial_cells',
#                                                   'Smooth_muscle_cells',
#                                                   'Fibroblasts',
#                                                   'Endothelial_cells'))

# library(RColorBrewer)
# qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
# #处理后有73种差异还比较明显的颜色，基本够用
# col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
col_vector
#
# sample_color <- col_vector[1:9]

pB5 <- ggplot(data = pB2_df, aes(x =Number, y = Sample, fill =  Cluster)) +
  geom_bar(stat = "identity", width=0.8,position="fill")+
  scale_fill_manual(values=col_vector[1:20]) +
  theme_bw()+
  theme(panel.grid =element_blank()) +
  labs(x="",y="Ratio")+
  ####用来将y轴移动位置
  theme(axis.text.y = element_text(size=12, colour = "black"))+
  theme(axis.text.x = element_text(size=12, colour = "black"))+
  theme(
    axis.text.x.bottom = element_text(hjust = 1, vjust = 1, angle = 45)
  )
pB5

############### 继续做柱状图###################

x <- table(scRNA@meta.data$celltype,scRNA@meta.data$tissue_type)
x3= t(t(x)/rowSums(t(x)))
x4 = as.data.frame(as.table(t(x3)))
colnames(x4) = c("sample","celltype","Freq")
library(stringr)
x4$group = str_remove(x4$sample,pattern = '_[1-2]')


top<-function(x){
  return(mean(x)+sd(x)/sqrt(length(x)))
}
bottom<-function(x){
  return(mean(x)-sd(x)/sqrt(length(x)))
}
table(x4$group)

x4$group

#### 根据个人更改！！！
dose_sham <-x4[which(x4$group=="Normal"),]
dose_MCAO<-x4[which(x4$group=="HF"),]
############################################
dose_sham=na.omit(dose_sham)
p5 <- ggplot(data=dose_sham,aes(x=celltype,y=Freq,fill=celltype))+
  stat_summary(geom = "bar",fun = "mean",
               position = position_dodge(0.9))+
  stat_summary(geom = "errorbar",
               fun.min = bottom,
               fun.max = top,
               position = position_dodge(0.9),
               width=0.2)+
  scale_y_continuous(expand = expansion(mult = c(0,0.1)))+
  theme_bw()+
  theme(panel.grid = element_blank())+
  labs(x="Celltype",y="Proportion")+
  geom_point(data=dose_sham,aes(celltype,Freq),size=3,pch=19)+
  theme(
    axis.text.x.bottom = element_text(hjust = 1, vjust = 1, angle = 45)
  ) +ggtitle("Normal")+scale_fill_manual(values = col_vector[1:20])

p5 <- p5+NoLegend()
p5
p6 <- ggplot(data=dose_MCAO,aes(x=celltype,y=Freq,fill=celltype))+
  stat_summary(geom = "bar",fun = "mean",
               position = position_dodge(0.9))+
  stat_summary(geom = "errorbar",
               fun.min = bottom,
               fun.max = top,
               position = position_dodge(0.9),
               width=0.2)+
  scale_y_continuous(expand = expansion(mult = c(0,0.1)))+
  theme_bw()+
  theme(panel.grid = element_blank())+
  labs(x="Celltype",y="Proportion")+
  geom_point(data=dose_MCAO,aes(celltype,Freq),size=3,pch=19)+
  theme(
    axis.text.x.bottom = element_text(hjust = 1, vjust = 1, angle = 45)
  ) +ggtitle("HF")+scale_fill_manual(values = col_vector)


p6

library(patchwork)

p_dim1 <- DimPlot(scRNA,split.by = 'tissue_type',cols = col_vector,label=F,repel = T)+NoLegend()+theme_void()
p_dim1

dir.create("./figures")
pB4 <- pB4+NoLegend()
pdf(file="./figures/day1_比例图2.pdf",width=16,height=5)
((pB4|pB5|p5|p6))
# (p5+p6+pB5))#+plot_annotation(tag_level = "A")&theme(plot.tag = element_text(size = 12))
dev.off()


dir.create("./figures")
# pdf(file="./figures/day1_dimplot.pdf",width=12,height=8)
# ((p_dim1|p_dim2)/
#     (p5+p6+pB5))#+plot_annotation(tag_level = "A")&theme(plot.tag = element_text(size = 12))
# dev.off()




p_dim1 <- DimPlot(scRNA,
                  # split.by = 'tissue_type',
                  cols = col_vector,label=F,repel = T,
                  pt.size = 0.1)+NoLegend()+theme_void()
p_dim1

library(patchwork)
pdf(file="./figures/day1_dimplot.pdf",width=12,height=3)
(p_dim1|p_dim2)+plot_layout(widths = c(1,3))
dev.off()







# ((p3|p1)+plot_layout(widths = c(1,3)))





# 读取单细胞数据集，读取自己的
# setwd('xxxx')
# scRNA=readRDS('scRNA_anno_衰老心脏_用这个.RDS')
library(Seurat)

unique(scRNA$celltype)


# 提单核细胞亚群,重新降维聚类
scRNAsub=subset(scRNA,celltype %in% 'Fibroblasts')


scRNAsub<- FindVariableFeatures(scRNAsub, selection.method = "vst", nfeatures = 2000)
scale.genes <-  VariableFeatures(scRNAsub)
scRNAsub <- ScaleData(scRNAsub, features = scale.genes)
scRNAsub<- RunPCA(scRNAsub, features = VariableFeatures(scRNAsub))
DimPlot(scRNAsub, reduction = "pca", group.by = "orig.ident")
ElbowPlot(scRNAsub)


## 重新harmony
library(harmony)
set.seed(1000)
scRNAsub <- RunHarmony(scRNAsub, group.by.vars = "orig.ident")
DimPlot(scRNAsub, reduction = "harmony", group.by = "orig.ident")

ElbowPlot(scRNAsub,reduction = 'harmony')


scRNAsub <- FindNeighbors(scRNAsub, reduction = 'harmony',dims = 1:10)
scRNAsub <- FindClusters(scRNAsub)
scRNAsub <- RunUMAP(scRNAsub, reduction = 'harmony',dims = 1:10)




scRNAsub <- readRDS('./scRNA_fibroblast.RDS')
# 瞄一眼，小亚群,0,2,4!!
DimPlot(scRNAsub, label=TRUE,split.by = 'tissue_type')

dev.off()
library(reshape2)
library(ggplot2)
pB2_df <- table(scRNAsub@meta.data$seurat_clusters,scRNAsub@meta.data$tissue_type) %>% melt()
colnames(pB2_df) <- c("Cluster","Sample","Number")
pB2_df$Cluster <- factor(pB2_df$Cluster)



pB2_df[1:3, 1:3]
data_cellcount1 <- pB2_df[1:length(unique(scRNAsub@meta.data$seurat_clusters)),]
data_cellcount2 <- pB2_df[(length(unique(scRNAsub@meta.data$seurat_clusters))+1):nrow(pB2_df),]
data_cellcount <- cbind(data_cellcount1, data_cellcount2)
data_cellcount$differcount <- data_cellcount[,3]-data_cellcount[,6]






library(RColorBrewer)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
#处理后有73种差异还比较明显的颜色，基本够用
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

sample_color <- col_vector[1:10]

pB4 <- ggplot(data = pB2_df, aes(x =Number, y = Sample, fill =  Cluster)) +
  geom_bar(stat = "identity", width=0.8,position="fill")+
  scale_fill_manual(values=col_vector[1:20]) +
  theme_bw()+
  theme(panel.grid =element_blank()) +
  labs(x="",y="Ratio")+
  ####用来将y轴移动位置
  theme(axis.text.y = element_text(size=12, colour = "black"))+
  theme(axis.text.x = element_text(size=12, colour = "black"))+
  theme(
    axis.text.x.bottom = element_text(hjust = 1, vjust = 1, angle = 45)
  )
pB4

saveRDS(scRNAsub,file ='scRNA_fibroblast.RDS')




#######################cellchat
## 细胞通讯分析，了解亚群
# 先定义细胞
scRNA <- readRDS("scRNA_anno.RDS")
scRNAsub <- readRDS("scRNA_fibroblast.RDS")

scRNAsub$celltype=ifelse(scRNAsub$seurat_clusters %in% c(1,6),'Cardioprotective_fibroblasts',
                         'Other_fibroblasts')


Idents(scRNAsub) <- "celltype"
VlnPlot(scRNAsub, features = "F3",split.by = "tissue_type")

saveRDS(scRNAsub,file ='scRNA_fibroblast.RDS')




pbmc.markers <- FindAllMarkers(object = scRNAsub,
                               only.pos = TRUE,
                               min.pct = 0.25,
                               logfc.threshold = 0.25)
dim(pbmc.markers)
sig.markers=pbmc.markers[(abs(as.numeric(as.vector(pbmc.markers$avg_log2FC)))>0.25 &
                            as.numeric(as.vector(pbmc.markers$p_val_adj))<0.05 &
                            as.numeric(as.vector(pbmc.markers$pct.1))> 0.5  ),]
# &  as.numeric(as.vector(pbmc.markers$pct.2))< 0.6 因为是相同的细胞，不需要对方完全不表达
length(unique(sig.markers$gene))
table(sig.markers$cluster)
write.table(sig.markers,file="./figures/markersallgenes_cluster.xls",sep="\t",row.names=F,quote=F)
save(pbmc.markers, sig.markers, file = "./figures/DEGs.Rdata")







table(scRNA$celltype)
# 排除内皮的其他细胞
scRNAother=subset(scRNA, celltype != 'Fibroblasts')
table(scRNAother$celltype)
scRNAother# # 抽取2000细胞
# set.seed(123456)
# a=sample(1:ncol(scRNAother),2000)
# scRNAother=scRNAother[,a]
# scRNAother

table(scRNAsub$celltype)
scRNA_chat=merge(scRNAsub,c(scRNAother))
scRNA_chat
## 用4.1.3的R
#devtools::install_github("sqjin/CellChat")
library(CellChat)




# 选取病鼠，如果内存不够，再进一步减少细胞数，例如随机抽2000个

# scRNA_chat <- subset(scRNA_chat, tissue_type=='Normal')


# ## 抽200细胞
# sampler
set.seed(123456)
a=sample(1:ncol(scRNA_chat),5000)
a
scRNA_chat=scRNA_chat[,a]


meta =scRNA_chat@meta.data # a dataframe with rownames containing cell mata data

data_input <- as.matrix(scRNA_chat@assays$RNA@data)
#data_input=data_input[,rownames(meta)]
identical(colnames(data_input),rownames(meta))
dir.create("./cellchat")
# save(scRNA_chat, file = "./cellchat/scRNA_chat_10000samples.Rdata")

meta =scRNA_chat@meta.data # a dataframe with rownames containing cell mata data

# data_input <- as.matrix(scRNA_chat@assays$RNA@data)
#data_input=data_input[,rownames(meta)]
# identical(colnames(data_input),rownames(meta))

library(CellChat)
library(devtools)
library(usethis)

unique(meta$celltype)
cellchat <- createCellChat(object = data_input, meta = meta, group.by = "celltype")

# CellChatDB <- CellChatDB.mouse
#人
CellChatDB <- CellChatDB.human
groupSize <- as.numeric(table(cellchat@idents))
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
cellchat@DB <- CellChatDB.use

cellchat@DB

dplyr::glimpse(CellChatDB$interaction)
##配体-受体分析
# 提取数据库支持的数据子集
cellchat <- subsetData(cellchat)
# 识别过表达基因
cellchat <- identifyOverExpressedGenes(cellchat)
# 识别配体-受体对
cellchat <- identifyOverExpressedInteractions(cellchat)
# 将配体、受体投射到PPI网络
cellchat <- projectData(cellchat, PPI.human)

unique(cellchat@idents)

# 等待
cellchat <- computeCommunProb(cellchat)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)

df.net<- subsetCommunication(cellchat)

write.csv(df.net,file ='./cellchat/cellchat.csv',quote=F)
#returns a data frame consisting of all the inferred cell-cell communications at the level of ligands/receptors. Set slot.name = "netP" to access the the inferred communications at the level of signaling pathways

#df.net <- subsetCommunication(cellchat, sources.use = c(1,2), targets.use = c(4,5))
#gives the inferred cell-cell communications sending from cell groups 1 and 2 to cell groups 4 and 5.

#df.net <- subsetCommunication(cellchat, signaling = c("WNT", "TGFb"))

cellchat <- aggregateNet(cellchat)
dev.off()
groupSize <- as.numeric(table(cellchat@idents))
groupSize

pdf(file="./cellchat/1.p_bubble_count_and_weight.pdf",width=12,height=7, onefile = T)
par(mfrow = c(1,2), xpd=TRUE)
p_bubble_count <- netVisual_circle(cellchat@net$count, vertex.weight = groupSize,
                                   weight.scale = T, label.edge= T,
                                   title.name = "Number of interactions")
p_bubble_weight <- netVisual_circle(cellchat@net$weight, vertex.weight = groupSize,
                                    weight.scale = T, label.edge= T,
                                    title.name = "Interaction weights/strength")
dev.off()



pdf(file="./cellchat/1.p_bubble_count_and_weight2.pdf",width=12,height=7, onefile = T)
par(mfrow = c(1,2), xpd=TRUE)
p_bubble_count <- netVisual_circle(cellchat@net$count, vertex.weight = groupSize,
                                   weight.scale = T, label.edge= F,
                                   title.name = "Number of interactions")
p_bubble_weight <- netVisual_circle(cellchat@net$weight, vertex.weight = groupSize,
                                    weight.scale = T, label.edge= F,
                                    title.name = "Interaction weights/strength")
dev.off()


table(scRNA_chat$celltype)
celltypes <- unique(scRNA_chat$celltype)

celltypes
celltypes <- celltypes[-c(4,9)]
celltypes
table(scRNA_chat$Mitocell)

p_bubble= netVisual_bubble(cellchat,
                           sources.use = c('Other_fibroblasts','Cardioprotective_fibroblasts'),
                           # targets.use = c('Astrocytes','Epithelial cells', 'Microglia','Monocytes','Oligodendrocytes'),
                           targets.use = celltypes,
                           remove.isolate = FALSE)+coord_flip()
p_bubble

pdf(file="./cellchat/2.p_高低对比气泡图.pdf",width=25,height=7)
p_bubble
dev.off()

# p_bubble= netVisual_bubble(cellchat,
#                            targets.use = c('Other_hepatocytes','Regeneration-related_hepatocytes'),
#                            # targets.use = c('Astrocytes','Epithelial cells', 'Microglia','Monocytes','Oligodendrocytes'),
#                            sources.use = celltypes,
#                            remove.isolate = FALSE)+coord_flip()
# p_bubble
#
# pdf(file="./cellchat/2.p_高低对比气泡图.pdf",width=12,height=7)
# p_bubble
# dev.off()


#必须把上一个图关掉
dev.off()
table(scRNA_chat$celltype)

path <- CellChatDB[["interaction"]][["pathway_name"]]
path <- as.data.frame(path)
path <- unique(path$path)
path
grep("MIF",path)
celltypes
celltypes
netVisual_aggregate(cellchat,
                    signaling = 'MIF',
                    sources.use = c('Other_hepatocytes',
                                    'Regeneration-related_hepatocytes'),
                    # targets.use = c('Astrocytes','Epithelial cells', 'Microglia','Monocytes',
                    #                 'Oligodendrocytes')
                    targets.use = celltypes)

?netVisual_aggregate



#必须把上一个图关掉
dev.off()
p_path_PTN <- netVisual_aggregate(cellchat, signaling = 'PTN',
                                  sources.use = c('sepsis_highHepatocytes','sepsis_lowHepatocytes'),
                                  targets.use = c('Astrocytes','Epithelial cells', 'Microglia','Monocytes','Oligodendrocytes'))

pdf(file="./cellchat/3.p_PTN通路.pdf",width=7,height=7)
p_path_PTN
dev.off()

cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

netAnalysis_signalingRole_network(cellchat, signaling = 'PTN', width = 8, height = 2.5, font.size = 10)




h1=netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
h2=netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
h1 + h2

pdf(file="./cellchat/4.netAnalysis_signalingRole_heatmap.pdf",width=8,height=8)
h1 + h2
dev.off()


save(cellchat,file = "./cellchat/cellchat_计算完成.Rdata")
gc()

load(file = "./cellchat/cellchat_计算完成.Rdata")
# 这一步自己的电脑跑直接奔溃？？可能要用服务器, 先保存后重新打开也是可以的
# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
gg1 <- netAnalysis_signalingRole_scatter(cellchat)
gg1
# ?netAnalysis_signalingRole_scatter

pdf(file="./cellchat/5.netAnalysis_signalingRole_散点图.pdf",width=5,height=4)
gg1
dev.off()

