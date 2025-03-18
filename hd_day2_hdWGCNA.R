run_hdWGCNA_step1_xc <- function(scRNA,genesforhdwgcna){
  dir.create("./figures")
  # setwd('figures/')
  # 推荐本地安装
  # devtools::install_local('hdWGNCA.zip')
  #BiocManager::install('harmony',update=F,ask=F)
  library(hdWGCNA)
  #加载单细胞分析包
  library(Seurat)
  #加载作图包
  library(tidyverse)
  library(cowplot)
  library(patchwork)
  #加载共表达网络分析包
  library(WGCNA)
  gc()

  # 工作目录
  # setwd('e:/writing/benke/hdWGCNA/')
  #设置随机种子
  set.seed(12345)

  # 读取单细胞数据集，读取自己的
  # scRNA=readRDS('scRNA_anno.RDS')
  # scRNA=readRDS('../FLC_genescoreMatrix_samples_7625.RDS')

  table(scRNA$celltype)

  # 提上皮亚群,重新降维聚类

  # #之所以要提取亚群，是因为如果是全部基因的wgcna，会导致其实都是各自的marker
  # #如果我们是自己定义基因的话，应该不需要提
  # # scRNA=subset(scRNA,celltype %in% c('tumor_cells',"cholangiocytes","hepatocytes"))
  # scRNA
  # scRNA<- FindVariableFeatures(scRNA, selection.method = "vst", nfeatures = 2000)
  # scale.genes <-  VariableFeatures(scRNA)
  # scRNA <- ScaleData(scRNA, features = scale.genes)
  # scRNA<- RunPCA(scRNA, features = VariableFeatures(scRNA))
  # DimPlot(scRNA, reduction = "pca", group.by = "orig.ident")
  # ElbowPlot(scRNA)
  # scRNA <- FindNeighbors(scRNA, dims = 1:5)
  # scRNA <- FindClusters(scRNA)
  # scRNA <- RunUMAP(scRNA, dims = 1:5)
  #
  #
  # source('rscript/run_harmony.R')
  #
  # # 瞄一眼，小亚群
  # DimPlot(scRNA, label=TRUE, group.by = "orig.ident")
  # scRNA <- run_harmony_xc(scRNA, resolution = 0.2)
  #
  # dev.off()
  #
  # scRNA
  # DimPlot(scRNA)
  # DimPlot(scRNA, label=TRUE, group.by = "orig.ident")


  #
  # #过滤出至少在5%的细胞中表达的基因
  # scRNA <- SetupForWGCNA(
  #   scRNA,
  #   gene_select = "fraction", # the gene selection approach
  #   fraction = 0.05, # fraction of cells that a gene needs to be expressed in order to be included
  #   wgcna_name = "Bio_com" # the name of the hdWGCNA experiment
  # )


  # genes <- read.table('../116 怎么提取KEGG数据库特定通路的基因列表/kinase_genes.txt')[,1]
  # genes <- intersect(genes, rownames(scRNA))
  # genesforhdwgcna <- rownames(scRNA)
  genesforhdwgcna
  genes <- genesforhdwgcna

  #过滤出至少在5%的细胞中表达的基因
  scRNA <- SetupForWGCNA(
    scRNA,
    gene_select = "custom", # the gene selection approach
    gene_list = genes, # the gene selection approach
    fraction = 0.05, # fraction of cells that a gene needs to be expressed in order to be included
    wgcna_name = "Bio_com" # the name of the hdWGCNA experiment
  )


  #构建metacells!!这一步非常重要，WGCNA对数据的稀疏性非常敏感，与普通转录组的WGCNA分析相比
  # 单细胞的稀疏矩阵的解决方法是WGCNA流程的问题核心
  table(scRNA$celltype)
  table(scRNA$orig.ident)
  # construct metacells  in each group
  scRNA<- MetacellsByGroups(
    seurat_obj = scRNA,k=20,
    max_shared = 10,
    # group.by一般关注的是组织类型和细胞类型!这边组织类型是orig.ident，CT正常，PR疾病
    group.by = c("celltype",'orig.ident'), # 也可以选择别的groupby
    ident.group = 'celltype' # set the Idents of the metacell seurat object
  )

  # table(scRNA)
  # normalize metacell expression matrix:
  scRNA <- NormalizeMetacells(scRNA)
  metacell_obj <- GetMetacellObject(scRNA)


  #转置表达矩阵
  # 安全起见，另起一个对象，以角质细胞细胞为例
  seurat_obj  <- SetDatExpr(
    scRNA,
    # group_name = "Bio_com", # 选择感兴趣的细胞类群！
    group_name = unique(scRNA$celltype), # 选择感兴趣的细胞类群！
    group.by='celltype' # the metadata column containing the cell type info. This same column should have also been used in MetacellsByGroups
  )


  #选择softpower
  seurat_obj <- TestSoftPowers(
    seurat_obj,
    setDatExpr = FALSE, # 这边不用转置了，前面已转置
  )


  # plot the results:
  plot_list <- PlotSoftPowers(seurat_obj)

  # assemble with patchwork
  wrap_plots(plot_list, ncol=2)


  # dev.off()
  dir.create("./figures")
  pdf(file="./figures/day2_hdWGCNA_softpower.pdf",width=8,height=6)
  # assemble with patchwork
  p1 <- wrap_plots(plot_list, ncol=2)
  print(p1)
  dev.off()




  #查看powerTable
  power_table <- GetPowerTable(seurat_obj)
  head(power_table)

  #记得保存上面hdWGNCA关键分析过程！！！
  saveRDS(seurat_obj, file='hdWGCNA_object.rds')
  return(seurat_obj)
}




run_hdWGCNA_step2_xc <- function(seurat_obj,softpower){
  dir.create("./figures")
  # setwd('figures/')
  # 推荐本地安装
  # devtools::install_local('hdWGNCA.zip')
  #BiocManager::install('harmony',update=F,ask=F)
  library(hdWGCNA)
  #加载单细胞分析包
  library(Seurat)
  #加载作图包
  library(tidyverse)
  library(cowplot)
  library(patchwork)
  #加载共表达网络分析包
  library(WGCNA)
  gc()

  # 工作目录
  # setwd('e:/writing/benke/hdWGCNA/')
  #设置随机种子
  set.seed(12345)

  # 读取单细胞数据集，读取自己的
  # scRNA=readRDS('scRNA_anno.RDS')
  # scRNA=readRDS('../FLC_genescoreMatrix_samples_7625.RDS')

  table(scRNA$celltype)

  unique(seurat_obj$celltype)
  #构建共表达网络
  softpower=5  # 根据自己的改，选5也没问题

  # construct co-expression network:
  seurat_obj <- ConstructNetwork(
    seurat_obj, soft_power=softpower,
    # group.by='celltype',
    # group_name= unique(seurat_obj$scissor)[2],
    setDatExpr = F)
  unique(seurat_obj$scissor)
  # ?ConstructNetwork




  dir.create("./figures")
  pdf(file="./figures/day2_hdWGCNA_Dendrogram.pdf",width=8,height=6)
  #可视化WGCNA网络
  PlotDendrogram(seurat_obj, main='hdWGCNA Dendrogram')
  dev.off()



  #(可选)获取TOM矩阵，可用来进行其他高级分析
  TOM <- GetTOM(seurat_obj)


  #计算模块协调特征
  #记得scale一下 or else harmony throws an error:
  seurat_obj <- Seurat::ScaleData(
    seurat_obj,
    features = GetWGCNAGenes(seurat_obj),

  )
  # 计算ME，根据组织类型分组
  # harmony必须biocManager安装，不可以用github安装！！！
  library(harmony)

  seurat_obj
  table(seurat_obj$orig.ident)

#   #报错了
#   seurat_obj <- ModuleEigengenes(
#     seurat_obj,
#     # group.by.vars="orig.ident" #harmony对象
#     group.by.vars="tissue_type2" #harmony对象
#   )
#
# ?ModuleEigengenes


  seurat_obj <- ModuleEigengenes(
    seurat_obj)
  seurat_obj <- ModuleConnectivity(seurat_obj)

  # plot genes ranked by kME for each module
  #可视化每个模块中，按照kME打分的基因PlotKMEs(seurat_obj, ncol=5)



  # 获取hub genes
  hub_df <- GetHubGenes(seurat_obj, n_hubs = 150)
  head(hub_df)

  # hub_df <- GetHubGenes(seurat_obj, n_hubs = 10000)
  head(hub_df)
  table(hub_df$module)

  head(hub_df)
  unique(hub_df$module)
  unique(hub_df$gene_name)
  saveRDS(hub_df, file = "./figures/day2_hdWGCNA_hub_df.RDS")

  hub_df <- readRDS('./figures/day2_hdWGCNA_hub_df.RDS')





  PlotKMEs(seurat_obj)

  pdf(file="./figures/day2_hdWGCNA_PlotKMEs.pdf",width=10,height=8)
  PlotKMEs(seurat_obj)
  dev.off()



  #记得保存上面hdWGNCA关键分析过程！！！
  saveRDS(seurat_obj, file='hdWGCNA_object.rds')


  # seurat_obj <- readRDS('hdWGCNA_object.rds')

  dev.off()


  ####------------一些可视化-----------------------
  ## 模块间的相关性
  library(igraph)
  library(qgraph)
  # install.packages("qgraph")
  # BiocManager::install("qgraph")
  # 载入保存的

  # seurat_obj=readRDS('./hdWGCNA_object.rds')

  # 画模块间相关性图
  ModuleCorrelogram(seurat_obj, sig.level = 0.001, pch.cex=2)



  pdf(file="./figures/day2_hdWGCNA_模块间相关性图.pdf",width=7,height=7)
  ModuleCorrelogram(seurat_obj, sig.level = 0.001, pch.cex=2)
  dev.off()

  # # 由于识别到了每个模块的hub基因，可以去计算hub基因打分
  # # compute gene scoring for the top 25 hub genes by kME for each module
  # # (方法一)with Seurat method
  # seurat_obj <- ModuleExprScore(
  #   seurat_obj,
  #   n_genes = 25,
  #   method='Seurat'
  # )

  # compute gene scoring for the top 25 hub genes by kME for each module
  # (方法二)with UCell method #推荐这种方法
  # 由于Ucell刚刚更新，所以4.1.x的同学请用本地安装,依赖包自行安装
  # devtools::install_local("./UCell-1.3.zip")
  library(UCell)
  seurat_obj <- ModuleExprScore(
    seurat_obj,
    n_genes = 25,
    method='UCell'
  )

  # featureplot
  # 瞄一眼
  DimPlot(scRNA, label=TRUE,split.by = 'celltype')

  plot_list <- ModuleFeaturePlot(
    seurat_obj,
    features='hMEs', # plot the hMEs
    order=TRUE ,# order so the points with highest hMEs are on top
  )

  # stitch together with patchwork
  p1 <- wrap_plots(plot_list)
  pdf(file="./figures/day2_hdWGCNA_模块评分UCell.pdf",width=12,height=12)
  p1 <- wrap_plots(plot_list)
  print(p1)
  dev.off()



  ### dotplot
  # get hMEs from seurat object
  MEs <- GetMEs(seurat_obj, harmonized=TRUE)
  mods <- colnames(MEs); mods <- mods[mods != 'grey']

  # add hMEs to Seurat meta-data:
  seurat_obj@meta.data <- cbind(seurat_obj@meta.data, MEs)


  # plot with Seurat's DotPlot function
  p <- DotPlot(seurat_obj, features=mods,group.by = "celltype")

  # flip the x/y axes, rotate the axis labels, and change color scheme:
  p <- p +
    coord_flip() +
    RotatedAxis() +
    scale_color_gradient2(high='red', mid='grey95', low='blue')
  p_dot <- p +
    coord_flip() +
    RotatedAxis() +
    scale_color_gradient2(high='#a363ae', mid='grey95', low='#638fc7')#+labs(x="")

  # plot output
  p
  p_dot


  pdf(file="./figures/day2_hdWGCNA_模块评分dotplot.pdf",width=4,height=5)
  print(p_dot)
  dev.off()






  # plot with Seurat's DotPlot function
  p <- DotPlot(seurat_obj, features=mods,group.by = "seurat_clusters")


  # flip the x/y axes, rotate the axis labels, and change color scheme:
  p <- p +
    coord_flip() +
    RotatedAxis() +
    scale_color_gradient2(high='red', mid='grey95', low='blue')
  p_dot <- p +
    coord_flip() +
    RotatedAxis() +
    scale_color_gradient2(high='#a363ae', mid='grey95', low='#638fc7')#+labs(x="")

  # plot output
  p
  p_seu <- p
  p_dot
  table(seurat_obj$seurat_clusters)

  pdf(file="./figures/day2_hdWGCNA_模块评分dotplot_seurat_clusters.pdf",width=7,height=5)
  print(p)
  dev.off()

  pdf(file="./figures/day2_hdWGCNA_模块评分dotplot_拼图.pdf",width=12,height=5.5)
  print((p_seu|p_dot)+plot_layout(widths = c(3, 1)))
  dev.off()


  # table(seurat_obj$turquoise)
  colnames(seurat_obj@meta.data)
  table(seurat_obj$metacell_grouping)
  # modules <- colnames(seurat_obj@meta.data)[13:21]
  modules <- as.character(unique(hub_df$module))
  modules
  library(ggpubr)
  dev.off()
  VlnPlot(seurat_obj, group.by = "celltype", features = modules, pt.size = 0,
          ncol = 5)&
    geom_boxplot(width=.1,col="black",fill="white", outlier.alpha = 0)&
    stat_compare_means()&labs(x="", y="module score")

  # VlnPlot(seurat_obj, group.by = "seurat_clusters", features = modules, pt.size = 0)&
  #   geom_boxplot(width=.1,col="black",fill="white", outlier.alpha = 0)&
  #   stat_compare_means()&labs(x="", y="module score")


  library(Nebulosa)
  a<-plot_density(seurat_obj, modules,
                  reduction = "umap",
                  # pal = "magma")
                  # pal = "cividis")
                  # pal = "viridis")
                  # pal = "inferno")
                  pal = "plasma")

  a


  dir.create("./figures")
  pdf(file="./figures/day2_hdWGCNA_module_expression.pdf",width=20,height=20)
  a
  dev.off()

  table(hub_df$module)
  aa <- hub_df[hub_df$module %in% c("turquoise"),]
  aa <- aa$gene_name
  aa
  length(aa)

  saveRDS(aa, file = paste0("./figures/module_","turquoise","_",length(aa),"genes.RDS"))

  aa1 <- as.data.frame(aa)

  #作曲家
  library(RColorBrewer)
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  #处理后有73种差异还比较明显的颜色，基本够用
  col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors,
                             rownames(qual_col_pals)))[30:33]

  library(ggpubr)
  selected_module <- "turquoise"
  p1 <- VlnPlot(seurat_obj, group.by = "celltype", features = selected_module, pt.size = 0,
                cols = col_vector)&
    geom_boxplot(width=.1,col="black",fill="white", outlier.alpha = 0)&
    stat_compare_means()&NoLegend()
  p1
  p_green <- p1

  # p1 <- VlnPlot(seurat_obj, group.by = "celltype", features = "blue", pt.size = 0,
  #               cols = col_vector)&
  #   geom_boxplot(width=.1,col="black",fill="white", outlier.alpha = 0)&
  #   stat_compare_means()&NoLegend()
  # p1
  # p_blue <- p1


  col_vector <- readRDS('../col_vector.RDS')
  p2 <- VlnPlot(seurat_obj, group.by = "seurat_clusters", features = selected_module, pt.size = 0,
                cols = col_vector)&
    geom_boxplot(width=.1,col="black",fill="white", outlier.alpha = 0)&
    stat_compare_means(label.x.npc = 0.3)&NoLegend()
  p2

  library(Nebulosa)
  p3 <- plot_density(seurat_obj, selected_module,
                     reduction = "umap",
                     # pal = "magma")
                     # pal = "cividis")
                     # pal = "viridis")
                     # pal = "inferno")
                     pal = "plasma")

  p3

  # p4 <- DimPlot(scRNA, label=TRUE,split.by = 'tissue_type',cols = col_vector)
  # p4

  # (p4+p_dot)/(p1|p2|p3)

  library(RColorBrewer)
  #
  #
  # DoHeatmap(subset(scRNA,downsample=50,),features = aa,
  #           group.by = 'seurat_clusters',
  #           assay='RNA',slot = 'data',
  #           group.colors =c('#313c63','#b42e20','#ebc03e','#377b4c',
  #                           '#7bc7cd','#5d84a4'),lines.width = 10)+
  #   # scale_color_gradientn(values = seq(0,1,0.2),
  #   #                       colours = c('blue','cyan','green','yellow','orange','red'))
  #   # scale_fill_gradientn(colors=c('white','firebrick3'),na.value = 'white')
  #   scale_fill_gradientn(colors=c('blue','cyan','green','yellow','orange','red'),na.value = 'white')


  # DoHeatmap(scRNA, features = aa)

  # dir.create("./figures")
  # pdf(file="./figures/拼图1.pdf",width=12,height=4)
  # (p4+p_dot)
  # dev.off()
  #
  # pdf(file="./figures/拼图2.pdf",width=13,height=4)
  # # (p4+p_dot)/
  # ((p3|p1|p2)+plot_layout(widths = c(1,1,3)))
  # dev.off()






  ##########################################    富集分析           ############################
  ############################ ##########################################
  #加载seurat数据和包
  # single-cell analysis package
  # setwd('e:/writing/benke/hdWGCNA/')
  library(Seurat)

  # plotting and data science packages
  library(tidyverse)
  library(cowplot)
  library(patchwork)

  # co-expression network analysis packages:
  library(WGCNA)
  library(hdWGCNA)

  # gene enrichment packages
  # install.packages('enrichR')
  library(enrichR)

  # BiocManager::install('GeneOverlap',update = F,ask = F)
  library(GeneOverlap)

  # using the cowplot theme for ggplot
  theme_set(theme_cowplot())

  # set random seed for reproducibility
  set.seed(12345)


  # load the Zhou et al snRNA-seq dataset
  # seurat_obj <- readRDS('hdWGCNA_object.rds')

  #GO富集分析
  # enrichr databases to test
  dbs <- c('GO_Biological_Process_2021','GO_Cellular_Component_2021','GO_Molecular_Function_2021')

  # 富集分析，会逐个模块分析
  seurat_obj <- RunEnrichr(
    seurat_obj,
    dbs=dbs, # character vector of enrichr databases to test
    max_genes = 150 # number of genes per module to test
  )


  # retrieve the output table

  dir.create("./figures")
  # make GO term plots作图，在文件夹下生成！
  EnrichrBarPlot(
    seurat_obj,
    outdir = "./figures/enrichr_plots", # name of output directory
    n_terms = 10, # number of enriched terms to show (sometimes more show if there are ties!!!)
    plot_size = c(5,7), # width, height of the output .pdfs
    logscale=TRUE # do you want to show the enrichment as a log scale?
  )

  enrich_df <- GetEnrichrTable(seurat_obj)
  saveRDS(enrich_df, file = './figures/enrichr_plots/enrich_df.RDS')
  saveRDS(enrich_df, file = './figures/enrichr_plots/enrich_df非常重要！.RDS')

  # load('figures/enrichr_plots/enrich_df.RDS')
  enrich_df <- readRDS('figures/enrichr_plots/enrich_df.RDS')

  library(openxlsx)
  write.xlsx(enrich_df, "./figures/enrichr_plots/enrich_df非常重要.xlsx" , rowNames =T)


  # table(enrich_df$module)
  # rt_bar = enrich_df[enrich_df$module %in% "turquoise",]
  # colnames(rt_bar)
  # labels=rt_bar[order(rt_bar$pvalue,decreasing =T),"Description"]
  # rt_bar$Description = factor(rt_bar$Description,levels=labels)
  # rt_bar <- rt_bar[c(1:20),]
  # #绘制
  # p_bar_kegg=ggplot(data=rt_bar)+geom_bar(aes(x=Description, y=Count, fill=pvalue),
  #                                         stat='identity')+
  #   coord_flip() + scale_fill_gradient(low="#DE5E57", high = "#4291A8")+
  #   xlab("Description") + ylab("Gene count") +
  #   # theme(axis.text.x=element_text(color="black", size=10),
  #   # axis.text.y=element_text(color="black", size=10)) +
  #   scale_y_continuous(expand=c(0, 0)) +
  #   scale_x_discrete(expand=c(0,0))+
  #   theme_bw()+ ggtitle("KEGG enrichment analysis")
  # # p_bar_kegg
  #
  #
  # #气泡图
  # pdf(file="bubble_kegg.pdf",width = 8,height = 7)
  # p_dot_kegg <- dotplot(kk, showCategory = 20, color="pvalue")+
  #   scale_colour_gradient(low="#c53832",high="#1a8192")+ ggtitle("KEGG enrichment analysis")
  # print(p_dot_kegg)
  # dev.off()
  #



  dir.create("./figures")
  pdf(file="./figures/hd_WGCNA_EnrichrDotPlot拼图.pdf",width=12,height=25)
  #气泡图
  # GO_Biological_Process_2021
  p <- EnrichrDotPlot(
    seurat_obj,
    # mods = c("turquoise","black"), # use all modules (this is the default behavior)
    # mods = c("turquoise","black"), # use all modules (this is the default behavior)
    database = "GO_Biological_Process_2021", # this has to be one of the lists we used above!!!
    n_terms=6
    # number of terms for each module
  )
  print(p)
  dev.off()


  #气泡图
  # GO_Biological_Process_2021
  p1 <- EnrichrDotPlot(
    seurat_obj,
    mods = c("turquoise"), # use all modules (this is the default behavior)
    # mods = c("turquoise"), # use all modules (this is the default behavior)
    database = "GO_Biological_Process_2021", # this has to be one of the lists we used above!!!
    n_terms=15
    # number of terms for each module
  )

  dir.create("./figures")
  # pdf(file="./figures/hd_WGCNA_EnrichrDotPlot_tuiquoise.pdf",width=16,height=4)
  pdf(file="./figures/hd_WGCNA_EnrichrDotPlot_blue_green.pdf",width=7,height=7)
  p1#+coord_flip()
  dev.off()


  EnrichrDotPlot(
    seurat_obj,
    mods = c("turquoise"), # use all modules (this is the default behavior)
    database = "GO_Biological_Process_2021", # this has to be one of the lists we used above!!!
    n_terms=6
    # number of terms for each module
  )


  #气泡图
  # GO_Cellular_Component_2021
  EnrichrDotPlot(
    seurat_obj,
    mods = c("turquoise"), # use all modules (this is the default behavior)
    database = "GO_Cellular_Component_2021", # this has to be one of the lists we used above!!!
    n_terms=6
    # number of terms for each module
  )

  #气泡图
  # GO_Biological_Process_2021
  EnrichrDotPlot(
    seurat_obj,
    mods = c("turquoise"), # use all modules (this is the default behavior)
    database = "GO_Molecular_Function_2021", # this has to be one of the lists we used above!!!
    n_terms=6
    # number of terms for each module
  )

  #差异基因重叠分析
  ## 这个分析帮助我们看到，哪些模块可能是相似的
  # compute cell-type marker genes with Seurat:
  # 常规方法计算差异基因/特征基因
  Idents(seurat_obj) <- seurat_obj$celltype
  table(Idents(seurat_obj))
  markers <- Seurat::FindAllMarkers(
    seurat_obj,
    only.pos = TRUE,
    logfc.threshold=log2(1.5)
  )
  saveRDS(markers, file = "./figures/markers.RDS")
  markers <- readRDS('figures/markers.RDS')

  table(markers$cluster)
  markers <- markers[markers$cluster %in% c("Cardioprotection-related_fibroblasts"),]
  markers$gene


  # marker13 <- readRDS('../82 OptimalML-101种组合算法建模 二分类/Figures/venn.RDS')
  # marker13
  hub_df <- hub_df[hub_df$module %in% c("turquoise"),]
  hub_df1 <- hub_df[hub_df$gene_name %in% markers$gene, ]

  library(openxlsx)
  write.xlsx(hub_df1, "./figures/hub_df_turquoise_markers_genes.xlsx",  rowNames =T)


  head(markers)
  marker13
  # marker <- data.frame(avg_log2FC=2,
  #                      gene=marker13,
  #                      p_val_adj=0,
  #                      cluster="ML_13_genes")
  #
  #
  # markers1 <- markers[,c("avg_log2FC","gene","p_val_adj","cluster")]
  #
  # markers2 <- rbind(markers1, marker)
  markers

  markers <- readRDS('figures/markers.RDS')

  table(markers$cluster)

  # BiocManager::install('GeneOverlap',update = F,ask = F)
  library(GeneOverlap)
  # compute marker gene overlaps
  overlap_df <- OverlapModulesDEGs(
    seurat_obj,
    # deg_df = markers,
    deg_df = markers,
    fc_cutoff = 1 # log fold change cutoff for overlap analysis
  )

  #条形图
  # overlap barplot, produces a plot for each cell type
  plot_list <- OverlapBarPlot(overlap_df)


  dir.create("./figures")
  pdf(file="./figures/OverlapBarPlot_拼图.pdf",width=6,height=3.5)
  # stitch plots with patchwork
  p2 <- wrap_plots(plot_list, ncol=2)
  print(p2)
  dev.off()


  #气泡图
  # plot odds ratio of the overlap as a dot plot
  p1 <- OverlapDotPlot(
    overlap_df,
    plot_var = 'odds_ratio') +
    ggtitle('Overlap of modules & cell-type markers')

  dir.create("./figures")
  pdf(file="./figures/OverlapDotPlot.pdf",width=9.5,height=2.5)
  print(p1)
  dev.off()


  library(patchwork)
  pdf(file="./figures/OverlapDotPlot拼图.pdf",width=9.5,height=7)
  print(p1/p2)
  dev.off()




  #网络可视化
  # network analysis & visualization package:
  # network analysis & visualization package:
  library(igraph)

  #可视化每个模块的网络图
  ModuleNetworkPlot(seurat_obj)

  #组合网络图，在文件夹下生成
  # hubgene network
  pdf(file="./figures/HubGeneNetworkPlot.pdf",width=8,height=8)
  HubGeneNetworkPlot(
    seurat_obj,
    n_hubs = 3, n_other=5,
    edge_prop = 0.75,
    mods = "all"
  )
  dev.off()


  #UMAP可视化
  ## 利用hub基因，重新UMAP，如此可以获得分群明显的基因分布图
  g <- HubGeneNetworkPlot(seurat_obj,  return_graph=TRUE)
  g
  seurat_obj <- RunModuleUMAP(
    seurat_obj,
    n_hubs = 10, # number of hub genes to include for the UMAP embedding
    n_neighbors=15, # neighbors parameter for UMAP
    min_dist=0.1 # min distance between points in UMAP space
  )


  # get the hub gene UMAP table from the seurat object
  umap_df <- GetModuleUMAP(seurat_obj)
  saveRDS(umap_df,file = './figures/umap_df.RDS')

  # plot with ggplot
  p1 <- ggplot(umap_df, aes(x=UMAP1, y=UMAP2)) +
    geom_point(
      color=umap_df$color, # color each point by WGCNA module
      size=umap_df$kME*2 # size of each point based on intramodular connectivity
    )  +labs(title = "gene distribution")
  # umap_theme()
  pdf(file="./figures/GetModuleUMAP.pdf",width=6,height=6)
  print(p1)
  dev.off()


  ModuleUMAPPlot(
    seurat_obj,
    edge.alpha=0.25,
    sample_edges=TRUE,
    edge_prop=0.1, # proportion of edges to sample (20% here)
    label_hubs=2 ,# how many hub genes to plot per module?
    keep_grey_edges=FALSE

  )



  #监督UMAP
  g <- ModuleUMAPPlot(seurat_obj,  return_graph=TRUE)
  g
  # run supervised UMAP:
  seurat_obj <- RunModuleUMAP(
    seurat_obj,
    n_hubs = 10,
    n_neighbors=15,
    min_dist=0.1,
    supervised=TRUE,
    target_weight=0.5
  )

  # get the hub gene UMAP table from the seurat object
  umap_df <- GetModuleUMAP(seurat_obj)
  saveRDS(umap_df,file = './figures/umap_df supervised UMAP.RDS')

  # plot with ggplot
  p2 <- ggplot(umap_df, aes(x=UMAP1, y=UMAP2)) +
    geom_point(
      color=umap_df$color, # color each point by WGCNA module
      size=umap_df$kME*2 # size of each point based on intramodular connectivity
    ) +labs(title = "gene distribution") #+theme(legend.position = "top")
  # umap_theme()

  pdf(file="./figures/GetModuleUMAP supervised UMAP.pdf",width=6,height=6)
  print(p2)
  dev.off()






  ##########################################   相关性分析            ############################
  ############################ ##########################################
  library(viridis)

  score <- readRDS('../hd_day5 scoring配体打分/5种score打分 13genes.RDS')
  seurat_obj


  marker13 <- readRDS('../82 OptimalML-101种组合算法建模 二分类/Figures/venn.RDS')
  marker13 <- gsub('HLAC','HLA-C',marker13)
  marker13 <- list(marker13)
  names(marker13) <- "ML_13_genes"



  seurat_obj <- AddModuleScore(seurat_obj, marker13)

  colnames(seurat_obj@meta.data)[ncol(seurat_obj@meta.data)]
  colnames(seurat_obj@meta.data)[ncol(seurat_obj@meta.data)] <- "ML_13_genes"
  seurat_obj$tissue_type <- seurat_obj$scissor
  table(seurat_obj$tissue_type)
  cd2 <- FeatureScatter(seurat_obj, "green", "ML_13_genes",
                        group.by = "tissue_type", pt.size = 0.5, cols = c("red", "blue"))+

    # geom_point(size = 0)+
    # ggsci::scale_color_igv()+
    theme_bw()+
    # facet_wrap(~tissue_type_celltype) +
    geom_smooth(method="lm", color="grey") +
    ggsci::scale_color_aaas()+ggpubr::stat_cor(method = "pearson", label.x.npc = 0,
                                              label.y.npc = 1)+
    labs(title = "Correlation analysis")+
    theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))

  cd2
  cd2 <- cd2+NoLegend()

  colnames(seurat_obj@meta.data)
  cd3 <- FeatureScatter(seurat_obj, "turquoise", "ML_13_genes",
                        group.by = "tissue_type", pt.size = 0.5, cols = c("red", "blue"))+

    # geom_point(size = 0)+
    # ggsci::scale_color_igv()+
    theme_bw()+
    # facet_wrap(~tissue_type_celltype) +
    geom_smooth(method="lm", color="grey") +
    ggsci::scale_color_aaas()+ggpubr::stat_cor(method = "pearson", label.x.npc = 0,
                                               label.y.npc = 1)+
    labs(title = "Correlation analysis")+
    theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))

  cd3
  cd3 <- cd3+NoLegend()

  cd4 <- FeatureScatter(seurat_obj, "blue", "ML_13_genes",
                        group.by = "tissue_type", pt.size = 0.5, cols = c("red", "blue"))+

    # geom_point(size = 0)+
    # ggsci::scale_color_igv()+
    theme_bw()+
    # facet_wrap(~tissue_type_celltype) +
    geom_smooth(method="lm", color="grey") +
    ggsci::scale_color_aaas()+ggpubr::stat_cor(method = "pearson", label.x.npc = 0,
                                               label.y.npc = 1)+
    labs(title = "Correlation analysis")+
    theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))

  cd4

  cd1 <- FeatureScatter(seurat_obj, "brown", "ML_13_genes",
                        group.by = "tissue_type", pt.size = 0.5, cols = c("red", "blue"))+

    # geom_point(size = 0)+
    # ggsci::scale_color_igv()+
    theme_bw()+
    # facet_wrap(~tissue_type_celltype) +
    geom_smooth(method="lm", color="grey") +
    ggsci::scale_color_aaas()+ggpubr::stat_cor(method = "pearson", label.x.npc = 0,
                                               label.y.npc = 1)+
    labs(title = "Correlation analysis")+
    theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))

  cd1

  library(patchwork)
  pdf(file="./figures/相关性分析.pdf", width=12, height=4)
  print(cd2|cd3|cd4)
  dev.off()


  pdf(file="./figures/相关性分析13genes.pdf", width=8, height=7.5)
  print((cd2|cd1)/(cd3|cd4))
  dev.off()




  cd22 <- FeatureScatter(seurat_obj, "green", "F3",
                         group.by = "tissue_type", pt.size = 0.5, cols = c("red", "blue"))+

    # geom_point(size = 0)+
    # ggsci::scale_color_igv()+
    theme_bw()+
    # facet_wrap(~tissue_type_celltype) +
    geom_smooth(method="lm", color="grey") +
    ggsci::scale_color_aaas()+ggpubr::stat_cor(method = "pearson", label.x.npc = 0,
                                               label.y.npc = 1)+
    labs(title = "Correlation analysis")+
    theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))

  cd22
  cd22 <- cd22+NoLegend()

  colnames(seurat_obj@meta.data)
  cd33 <- FeatureScatter(seurat_obj, "turquoise", "F3",
                         group.by = "tissue_type", pt.size = 0.5, cols = c("red", "blue"))+

    # geom_point(size = 0)+
    # ggsci::scale_color_igv()+
    theme_bw()+
    # facet_wrap(~tissue_type_celltype) +
    geom_smooth(method="lm", color="grey") +
    ggsci::scale_color_aaas()+ggpubr::stat_cor(method = "pearson", label.x.npc = 0,
                                               label.y.npc = 1)+
    labs(title = "Correlation analysis")+
    theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))

  cd33
  cd33 <- cd33+NoLegend()

  cd44 <- FeatureScatter(seurat_obj, "blue", "F3",
                         group.by = "tissue_type", pt.size = 0.5, cols = c("red", "blue"))+

    # geom_point(size = 0)+
    # ggsci::scale_color_igv()+
    theme_bw()+
    # facet_wrap(~tissue_type_celltype) +
    geom_smooth(method="lm", color="grey") +
    ggsci::scale_color_aaas()+ggpubr::stat_cor(method = "pearson", label.x.npc = 0,
                                               label.y.npc = 1)+
    labs(title = "Correlation analysis")+
    theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))

  cd44

  cd33
  library(patchwork)
  pdf(file="./figures/相关性分析2-F3.pdf", width=12, height=4)
  print(cd22|cd33|cd44)
  dev.off()




  cd5 <- FeatureScatter(seurat_obj, "ML_13_genes","F3",
                        group.by = "tissue_type", pt.size = 0.5, cols = c("red", "blue"))+

    # geom_point(size = 0)+
    # ggsci::scale_color_igv()+
    theme_bw()+
    # facet_wrap(~tissue_type_celltype) +
    geom_smooth(method="lm", color="grey") +
    ggsci::scale_color_aaas()+ggpubr::stat_cor(method = "pearson", label.x.npc = 0,
                                               label.y.npc = 1)+
    labs(title = "Correlation analysis")+
    theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))

  cd5

  library(patchwork)
  pdf(file="./figures/相关性分析3-F3.pdf", width=8, height=8)
  print((cd22|cd5)/(cd33|cd44))
  dev.off()



  return(seurat_obj)

  # setwd('../')
}





##########################################   转录因子            ############################
############################ ##########################################









MotifScan<- function (seurat_obj, species_genome, pfm, EnsDb, wgcna_name = NULL)
{
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  motif_df <- data.frame(motif_name = purrr::map(1:length(pfm),
                                                 function(i) {
                                                   pfm[[i]]@name
                                                 }) %>% unlist, motif_ID = purrr::map(1:length(pfm), function(i) {
                                                   pfm[[i]]@ID
                                                 }) %>% unlist)
  gene.promoters <- ensembldb::promoters(EnsDb, filter = ~gene_biotype ==
                                           "protein_coding") %>% subset(seqnames %in% c(1:100))
  gene.coords <- ensembldb::genes(EnsDb, filter = ~gene_biotype ==
                                    "protein_coding") %>% subset(seqnames %in% c(1:100))
  gene.promoters$symbol <- gene.coords$symbol[match(gene.promoters$gene_id,
                                                    names(gene.coords))]
  gene.promoters <- keepSeqlevels(gene.promoters, value = levels(droplevels(seqnames(gene.promoters))))
  old_levels <- levels(seqnames(gene.promoters))
  new_levels <- ifelse(old_levels %in% c("X", "Y"), old_levels,
                       paste0("chr", old_levels))
  gene.promoters <- renameSeqlevels(gene.promoters, new_levels)
  genome(seqinfo(gene.promoters)) <- species_genome
  my_promoters <- GRanges(seqnames = droplevels(seqnames(gene.promoters)),
                          IRanges(start = start(gene.promoters), end = end(gene.promoters)),
                          symbol = gene.promoters$symbol, genome = species_genome)
  print("Matching motifs...")
  motif_ix <- motifmatchr::matchMotifs(pfm, my_promoters, genome = species_genome)
  tf_match <- motifmatchr::motifMatches(motif_ix)
  rownames(tf_match) <- my_promoters$symbol
  motif_df <- motif_df[motif_df$motif_ID %in% colnames(tf_match), ] # removes motif IDs not present in tf_match debug add
  colnames(tf_match) <- motif_df$motif_name
  gene_list <- rownames(seurat_obj)
  gene_list <- gene_list[gene_list %in% rownames(tf_match)]
  tf_match <- tf_match[gene_list, ]
  print("Getting TF target genes...")
  tfs <- motif_df$motif_name
  tf_targets <- list()
  n_targets <- list()
  for (cur_tf in tfs) {
    tf_targets[[cur_tf]] <- names(tf_match[, cur_tf][tf_match[,
                                                              cur_tf]])
    n_targets[[cur_tf]] <- length(tf_targets[[cur_tf]])
  }
  n_targets <- unlist(n_targets)
  motif_df$n_targets <- n_targets
  seurat_obj <- SetMotifMatrix(seurat_obj, tf_match)
  seurat_obj <- SetMotifs(seurat_obj, motif_df)
  seurat_obj <- SetMotifTargets(seurat_obj, tf_targets)
  seurat_obj <- SetPFMList(seurat_obj, pfm)
  seurat_obj
}


# require(stats4)
# library(stats4)
library(DirichletMultinomial)
library(TFBSTools)
library(hdWGCNA)
library(JASPAR2020)
library(EnsDb.Hsapiens.v86)
BiocManager::install("TFBSTools")
# BiocManager::install("JASPAR2020")
# get the pfm from JASPAR2020 using TFBSTools
pfm_core <- TFBSTools::getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)
saveRDS(pfm_core, file = "pfm_core.RDS")

seurat_obj <- readRDS('hdWGCNA_object.rds')
pfm_core <- readRDS('pfm_core.RDS')
library(stringr)
# run the motif scan with these settings for the mouse dataset
seurat_obj <- MotifScan(
  seurat_obj,
  species_genome = 'hg38',
  pfm = pfm_core,
  EnsDb = EnsDb.Hsapiens.v86
)
dim(GetMotifMatrix(seurat_obj))


library(GeneOverlap)
# OverlapModulesMotifs
# overlap between modules & TF target genes:
seurat_obj<- OverlapModulesMotifs(seurat_obj)

# look at the overlap data
head(GetMotifOverlap(seurat_obj))


library(dplyr)
library(ggplot2)
library(Seurat)

# BiocManager::install("ggseqlogo")
library(ggseqlogo)
# plot the top TFs overlapping with

MotifOverlapBarPlot
MotifOverlapBarPlot(
  seurat_obj,
  motif_font = 'xkcd_regular',
  plot_size=c(5,6)
)

source('run_MotifOverlapBarPlot_xc.R')
run_MotifOverlapBarPlot_xc(
  seurat_obj,
  motif_font = 'xkcd_regular',
  plot_size=c(5,6)
)

library(UCell)
seurat_obj <- MotifTargetScore(
  seurat_obj,
  method='UCell'
)


df <- GetMotifOverlap(seurat_obj)
df

# cur_df <- df %>% subset(tf == 'SOX9')
cur_df <- df %>% subset(tf == 'HIF1A')

# unique(rownames(df))

plot_var <- 'odds_ratio'
p <- cur_df %>%
  ggplot(aes(y=reorder(module, odds_ratio), x=odds_ratio)) +
  geom_bar(stat='identity', fill=cur_df$color) +
  geom_vline(xintercept = 1, linetype='dashed', color='gray') +
  geom_text(aes(label=Significance), color='black', size=3.5, hjust='center') +
  ylab('') +
  xlab("Odds Ratio") +
  ggtitle("HIF1A overlap") +
  theme(
    plot.title = element_text(hjust = 0.5)
  )

library(patchwork)
pdf(file="./figures/转录因子与模块之间的关系.pdf", width=4.5, height=3.5)
print(p)
dev.off()

# png(paste0(fig_dir, 'HIF1A_motif_overlap_or.png'), width=3, height=4, units='in', res=400)
p


