run_monocle2_step1_xc <- function(scRNA){
  dir.create("./figures")
  # setwd('figures/')


  ## 拟时序分析 ------------------------------
  library(Seurat)

  #没有monocle要先安装 BiocManager::install()
  #BiocManager::install('monocle',update = F,ask = F)

  #没有的包先安装
  library(BiocGenerics)
  library(monocle)
  library(tidyverse)
  library(patchwork)


  rm(list=ls())
  options(stringsAsFactors = F)
  # scRNA=readRDS('../2. scissors表型分析/figures/3DMSCs41181_包括scissors信息.RDS')
  # DimPlot(scRNA,group.by = "scissor")
  # table(scRNA$scissor)
  # table(pbmc$scissor)
  # table(pbmc$orig.ident)
  # DimPlot(pbmc,group.by = "scissor")

  # scRNAsub=subset(scRNA,celltype=='Epithelial_cells')
  #scRNAsub=subset(scRNAsub,tissue_type=='Tumor')
  scRNAsub=scRNA
  data=as.matrix(scRNAsub@assays$RNA@counts)

  data <- as(data, 'sparseMatrix')
  pd <- new('AnnotatedDataFrame', data = scRNAsub@meta.data)
  fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
  fd <- new('AnnotatedDataFrame', data = fData)
  ## 以下代码一律不得修改 ！
  mycds <- newCellDataSet(data,
                          phenoData = pd,
                          featureData = fd,
                          expressionFamily = negbinomial.size())

  mycds <- estimateSizeFactors(mycds)
  mycds <- estimateDispersions(mycds, cores=4, relative_expr = TRUE)
  saveRDS(mycds, file = "mycds.RDS")
  ##使用monocle选择的高变基因，不修改
  disp_table <- dispersionTable(mycds)
  disp.genes <- subset(disp_table, mean_expression >= 0.1 &
                         dispersion_empirical >= 1 * dispersion_fit)$gene_id
  mycds <- setOrderingFilter(mycds, disp.genes)
  plot_ordering_genes(mycds)

  #降维
  mycds <- reduceDimension(mycds, max_components = 2, method = 'DDRTree')
  saveRDS(mycds, file = "mycds.RDS")
  gc()

  # #尝试减少分支
  # mycds <- reduceDimension(mycds, max_components = 2,
  #                          num_dim = 6,
  #                          residualModelFormulaStr = "~tissue_type", #去除样本影响
  #                          reduction_method =  'DDRTree')

  #排序，报错请用4.1.3的R，并重装monocle
  mycds <- orderCells(mycds)
  saveRDS(mycds, file = "mycds.RDS")


  # return(p_top)
}




run_monocle2_step2_xc <- function(mycds){
  dir.create("./figures")
  # setwd('figures/')


  ## 拟时序分析 ------------------------------
  library(Seurat)

  #没有的包先安装
  library(BiocGenerics)
  library(monocle)
  library(tidyverse)
  library(patchwork)


  # rm(list=ls())
  # options(stringsAsFactors = F)

  # mycds <- readRDS('mycds.RDS')
  #State轨迹分布图
  # dev.off()
  plot1 <- plot_cell_trajectory(mycds, cell_size = 0.5, color_by = "State")
  plot1

  ##Pseudotime轨迹图
  plot2 <- plot_cell_trajectory(mycds, cell_size = 0.5, color_by = "Pseudotime")+
    scale_color_gradientn(values = seq(0,1,0.2),
                          # colours = c('blue','cyan','green','yellow','orange','red'))
                          colours = c('white','red'))
  plot2


  #State轨迹分布图
  plot3 <- plot_cell_trajectory(mycds, color_by = "tissue_type")+
    facet_wrap(~tissue_type, nrow = 1)+
    scale_color_manual(values = c("#324c63", "#a363ae"))+theme(legend.position = "none")
  plot3

  # # mycds$tissue_type <- mycds$scissor
  plot3 <-  plot_cell_trajectory(mycds, cell_size = 0.5,color_by = "tissue_type")+
    scale_color_manual(values = c("#324c63", "#a363ae"))
  plot3

  # cols <- readRDS('../col_ThreeD_心衰专用色.RDS')
  cols <- readRDS('../col_vector.RDS')
  plot4 <- plot_cell_trajectory(mycds, cell_size = 0.5,color_by = "seurat_clusters")+scale_color_manual(values = cols)
  plot4


  ##合并出图
  plotc <- plot1|plot2|plot3|plot4
  plotc

  dir.create("./pseudotime")
  pdf(file="./pseudotime//6_scRNA_endo_pseudotime.pdf",width = 11,height =4)
  print(plotc)
  dev.off()


  library(viridis)
  cds <- mycds
  pdf(file="./pseudotime/cellType.trajectory pseudotime.pdf",width=3,height=3.5)
  plot_cell_trajectory(cds, color_by = "Pseudotime",
                       cell_size = 0.5,
                       show_branch_points = F,
                       use_color_gradient = F)+
    scale_color_gradientn(values = seq(0,1,0.2),
                          colours = c('blue','cyan','green','yellow','orange','red'))
  # scale_fill_gradientn(colors =   viridis(10, alpha = 1, begin = 0.5, end = 1, direction = -1, option = "B"))
  # colorRampPalette(c("#3C8DAD", "#FFEF78","#FF5C58"))(15)
  # scale_color_manual(values = c("navy","white","firebrick3"))
  dev.off()

  mycol <- ggsci::pal_simpsons()(15)
  mycol <- cols

  pdf(file="./pseudotime/cellType.trajectory 分开3.pdf",width=12,height=4)
  plot_cell_trajectory(cds, color_by = "tissue_type") +
    facet_wrap(~tissue_type, nrow = 1)+
    # scale_color_manual(values = mycol)
    # scale_color_manual(values = c("grey", "blue", "red"))
    scale_color_manual(values = c("#324c63", "#a363ae"))
  dev.off()

  #分成三种疾病状态
  pdf(file="./pseudotime/cellType.trajectory 分开3 state.pdf",width=15,height=6)
  plot_cell_trajectory(cds, color_by = "tissue_type") +
    facet_wrap(~tissue_type, nrow = 1)
  dev.off()





  ####################################################################添加“树形图”
  # 添加“树形图”
  # 基本展示完了，那如何得到文章中的结果呢？
  # 此处需要plot_complex_cell_trajectory函数添加“树形图”即可
  p1 <- plot_cell_trajectory(cds, x = 1, y = 2, color_by = "tissue_type") +
    theme(legend.position='none',panel.border = element_blank()) + #去掉第一个的legend
    # scale_color_manual(values = mycol)
    # scale_color_manual(values = c("grey","skyblue","red"))
    scale_color_manual(values = c("grey","blue","red"))
  p1
  p_tree <- plot_complex_cell_trajectory(cds, x = 1, y = 2,
                                     color_by = "tissue_type")+
    # scale_color_manual(values = mycol) +
    # scale_color_manual(values = c("grey","blue","red")) +
    scale_color_manual(values = c("#324c63", "#a363ae")) +
    theme(legend.title = element_blank())
  p_tree
  # p_tree <- p_tree+coord_flip()

  library(patchwork)
  pdf(file="./pseudotime/树形图.pdf",width=8,height=14)
  p1 / p_tree
  dev.off()

  pdf(file="./pseudotime/树形图.pdf",width=5,height=5)
  p_tree#+coord_flip()
  dev.off()


  ####################################################################山峦图
  library(monocle)
  library(tidyverse)
  library(ggridges)
  library(RColorBrewer)
  library(scales)

  cds$Cluster <- cds$tissue_type
  plotdf=pData(cds)

  p_ridge <- ggplot(plotdf, aes(x=Pseudotime,y=Cluster,fill=Cluster))+
    geom_density_ridges(scale=1) +
    geom_vline(xintercept = c(5,10),linetype=2)+
    scale_y_discrete("")+
    theme_minimal()+
    theme(
      panel.grid = element_blank()
    )
  print(p_ridge)
  ggsave("./pseudotime/山峦图.pdf",width = 13,height = 7,units = "cm")



  p_ridge2 <- ggplot(plotdf, aes(x=Pseudotime,y=Cluster,fill = stat(x))) +
    geom_density_ridges_gradient(scale=1) +
    geom_vline(xintercept = c(5,30),linetype=5)+
    scale_fill_gradientn(name="Pseudotime",colors = colorRampPalette(rev(brewer.pal(9, "PRGn")))(62))+
    # scale_fill_gradientn(name="Pseudotime",colors = colorRampPalette(rev(brewer.pal(10, "Spectral")))(99))+
    scale_y_discrete("")+
    theme_minimal()+
    theme(
      panel.grid = element_blank()
    )
  print(p_ridge2)
  ggsave("./pseudotime/tmp2山峦图.pdf",width = 13,height = 7,units = "cm")




  ##########################################               ############################
  ############################ ##########################################



  gene <- readRDS('../figures/hub_df_turquoise_markers_genes.RDS')
  gene <- gene$gene_name
  gene


  ## 查看关键差异基因在拟时序中所处的时间
  # load('../hd_day5 scoring配体打分/GT_significant.Rdata')

  # genes <- readRDS('../82 OptimalML-101种组合算法建模 二分类/Figures/venn.RDS')
  # genes
  # genes <- gsub("HLAC", "HLA-C", genes)
  genes <- gene
  # dev.off()
  # my_pseudotime_cluster <- plot_pseudotime_heatmap(mycds[rownames(df2),],
  my_pseudotime_cluster <- plot_pseudotime_heatmap(mycds[genes,],
                                                   # num_clusters = 2, # add_annotation_col = ac,
                                                   show_rownames = TRUE,
                                                   return_heatmap = TRUE)

  my_pseudotime_cluster



  # gene <- rownames(df2)
  # gene <- genes
  # gene
  table(mycds$tissue_type)
  ac <- data.frame( group = mycds$tissue_type)
  rownames(ac) <- colnames(mycds)
  nrow(ac)  #只能100行！

  ##随机挑选2000个细胞做
  a=sample(1:nrow(ac),100,replace = F)
  a
  ac=ac[a,]

  ac <- as.data.frame(ac)

  my_pseudotime_cluster <- plot_pseudotime_heatmap(mycds[gene,],
                                                   # num_clusters = 2,
                                                   add_annotation_col = ac,
                                                   hmcols = colorRampPalette(rev(brewer.pal(9, "PRGn")))(62),
                                                   show_rownames = TRUE,
                                                   return_heatmap = TRUE)


  # ?plot_pseudotime_heatmap



  library(RColorBrewer)
  dir.create("./pseudotime")
  #绘制与分支相关的基因热图
  pdf(file="./pseudotime/cellType.trajectory 分支热图.pdf",width=5,height=4)
  my_branched_heatmap <- plot_genes_branched_heatmap (mycds[gene,],
                                                      branch_point = 1,
                                                      num_clusters = 2,
                                                      cores = 4,
                                                      hmcols = colorRampPalette(rev(brewer.pal(9, "PRGn")))(62),
                                                      use_gene_short_name = TRUE,
                                                      # add_annotation_col = ac,
                                                      branch_colors = c("#979797","#F05662","#7990C8"), #pre-branch. Cell fate 1, Cell fate2
                                                      show_rownames = T,return_heatmap = TRUE)
  # my_branched_heatmap$ph
  dev.off()



  # library(export)
  # graph2pdf("分支差异分析热图.pdf")


  ##########################################  获得monocle2热图 pheatmap             ############################
  ############################ ##########################################


  #我们按照这个函数的实现，获得需要的原始表达值矩阵。很多参数都是按照默认的参数来，其它参数可以自行调整。
  cds <- mycds
  newdata <- data.frame(Pseudotime = seq(min(pData(cds)$Pseudotime),
                                         max(pData(cds)$Pseudotime),
                                         length.out = 100))

  Time_genes <- gene
  m <- genSmoothCurves(cds[Time_genes,], trend_formula = '~sm.ns(Pseudotime, df=3)',
                       relative_expr = T, new_data = newdata)

  #remove genes with no expression in any condition

  m=m[!apply(m,1,sum)==0,]

  m=log10(m+1)

  # Row-center the data.

  m=m[!apply(m,1,sd)==0,]

  m=Matrix::t(scale(Matrix::t(m),center=TRUE))

  m=m[is.na(row.names(m)) == FALSE,]

  m[is.nan(m)] = 0

  m[m>3] = 3

  m[m<-3] = -3

  heatmap_matrix <- m
  #这样我们就可以画个基本的heatmap了

  row_dist <- as.dist((1 - cor(Matrix::t(heatmap_matrix)))/2)

  row_dist[is.na(row_dist)] <- 1
  library(pheatmap)
  p1 <- pheatmap(m,

                 useRaster = F,
                 cluster_cols=FALSE,
                 cluster_rows=T,

                 show_rownames=T,

                 show_colnames=F,

                 clustering_method = "ward.D2",

                 clustering_distance_rows=row_dist,

                 cutree_rows=4,

                 border_color = NA,
                 filename='pseudotime/monocle2_geneheatmap.pdf',
                 color=colorRampPalette(c("navy","white","firebrick3"))(100)

  )


  # ====添加行注释和列注释=====

  annotation_col = data.frame(
    pseudotime = scales::rescale(newdata$Pseudotime, to = c(-1, 1)))
  row.names(annotation_col) <- colnames(m)



  annotation_row <- data.frame(Cluster=factor(cutree(p1$tree_row, 4)))

  row.names(annotation_row) <- rownames(m)



  rowcolor <- c("#85B22E","#E29827","#922927",'#57C3F3')

  names(rowcolor) <- c("1","2","3","4") #类型颜色

  ann_colors <- list(pseudotime=viridis(100),

                     Cluster=rowcolor) #颜色设置



  pheatmap(m,

           useRaster = T,

           cluster_cols=FALSE,

           cluster_rows=T,

           show_rownames=T,

           show_colnames=F,

           clustering_method = "ward.D2",

           clustering_distance_rows=row_dist,

           cutree_rows=4,

           border_color = NA,

           # filename=NA,
           filename='pseudotime/monocle2_geneheatmap2.pdf',

           color=colorRampPalette(c("navy","white","firebrick3"))(100),

           annotation_col = annotation_col,

           annotation_colors=ann_colors,

           annotation_row = annotation_row

  )


  ##合并出图
  plotc <- plot1|plot2|plot3|plot4
  plotc

  dev.off()
  p_ridge

  dir.create("./pseudotime")
  pdf(file="./pseudotime//6_scRNA_endo_pseudotime拼图.pdf",width = 15,height =9)
  print(plotc/(p_ridge|p_tree|plot_spacer()))
  dev.off()




##########################################     hepatology里的画法          ############################
  ############################ ##########################################


  ## 4.2 特定基因轨迹图----
  s.genes <- c("AEBP1","ANXA2","CCL19","CCL20","COL1A1","COL1A2","COL3A1",
               "COL4A1","COL4A2","COL5A1","CXCL10","DPT","EFEMP1",
               "FBLN5","IGFBP7","LAMC3","LGALS3","LOXL4","LTBP2","LUM",
               "MFAP4","MGP","MMP2","MMP7","PCOLCE2","PDGFD","S100A6",
               "SPP1","THBS1","THBS2","TIMP1","VCAN","CTGF","CLEC4M")

  s.genes %in% rownames(HSMM_1)
  cds_subset <- HSMM_1[s.genes,]
  ### 4.2.1 by state ----
  # 点图（抖动）
  p1 <- plot_genes_jitter(HSMM_1[s.genes,], grouping = "celltype", color_by = "celltype")
  # 小提琴图
  p2 <- plot_genes_violin(HSMM_1[s.genes,], grouping = "celltype", color_by = "celltype")
  # 伪时间图
  p3 <- plot_genes_in_pseudotime(HSMM_1[s.genes,], color_by = "celltype")
  plotc <- p1|p2|p3
  plotc
  ggsave(plotc,filename = "monocle/new.tumor.Pseudotime2.pdf",width = 16,height = 25)

  ## 4.3 拟时相关基因聚类热图----

  # monocle包的differentialGeneTest()函数可以按条件进行差异分析，将相关参数设为
  # fullModelFormulaStr = "~sm.ns(Pseudotime)"时，可以找到与拟时先关的差异基因。
  # 为了提高效率，可直接使用cluster差异基因或者高变基因，进行筛选，再绘制热图。以高变基因为例：
  #disp_table <- dispersionTable(HSMM_1)
  #disp.genes <- subset(disp_table, mean_expression >= 0.5&dispersion_empirical >= 1*dispersion_fit)
  #disp.genes <- as.character(disp.genes$gene_id)
  diff_test <- differentialGeneTest(HSMM_1[s.genes,], cores = 4,
                                    fullModelFormulaStr = "~sm.ns(Pseudotime)")
  #找到与拟时先关的差异基因
  sig_gene_names <- row.names(subset(diff_test))


  diff_test %>% arrange(qval) %>% head(50) %>% select(gene_short_name) -> gene_to_cluster
  gene_to_cluster <- gene_to_cluster[,1]
  gene_to_cluster

  #以qval为指标，挑选显著性
  plot4 <- plot_pseudotime_heatmap(HSMM_1[gene_to_cluster,], num_clusters=5,
                                   show_rownames=T, return_heatmap=T)

  ggsave(plot4,filename = "monocle/new.tumor.Pheatmap_by_state.pdf",units = "cm",width = 20,height = 10)

  table(pData(HSMM_1)$celltype)

  cds <- HSMM_1
  phe <- pData(cds)
  counts = Biobase::exprs(cds)
  dim(counts)




  n=t(scale(t( counts[gene_to_cluster,] ))) # 'scale'可以对log-ratio数值进行归一化
  n[n>2]=2
  n[n< -2]= -2
  n[1:4,1:4]

  colnames(phe)
  # ac=phe[,c(1,20,22)]#这里是对热图进行分组的信息，从colnames(phe)中选

  ac=phe[,c("Pseudotime", "orig.ident2","tissue_type")]#这里是对热图进行分组的信息，从colnames(phe)中选
  head(ac)
  rownames(ac)=colnames(n)
  dim(n)
  n[is.na(n)]<-0
  pheatmap::pheatmap(n,show_colnames =F,
                     show_rownames = T,
                     annotation_col=ac)

  od=order(ac$Pseudotime)
  p=pheatmap::pheatmap(n[,od],show_colnames =F,
                       show_rownames = T,cluster_cols = F,
                       annotation_col=ac[od,]
  )
  library(ggplotify)
  g = as.ggplot(p)
  ggsave(g,filename = "monocle/Monocle-pheatmap.pdf",width = 10,height = 8)

  # saveRDS(mycds, file = "./pseudotime/mycds.RDS")
  return(p_top)
}
