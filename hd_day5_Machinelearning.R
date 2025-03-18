# ### 确保自己所在的文件夹
# # seurat_obj=readRDS('./scRNA/hdWGCNA_object.rds')
# seurat_obj=readRDS('./hd_day2_hdWGCNA/hdWGCNA_object.rds')
library(hdWGCNA)
library(tidyverse)
# # 可以根据自己的需要扩大个数
# hub_df <- GetHubGenes(seurat_obj, n_hubs = 50)
#
# gene=hub_df[hub_df$module %in% c('pink','blue'),]
# gene=gene$gene_name
# gene


gene <- readRDS('./scRNA/figures/DEGS.RDS')
gene <- readRDS('./hd_day2_hdWGCNA/figures/DEGs与blue的交集.RDS')
gene <- readRDS('./figures/hub_df_turquoise_markers_genes.RDS')
# table(gene$cluster)

table(gene$module)
gene <- gene$gene
gene

# gene <- gene[gene$module %in% c("blue"),]
# gene <- gene$gene
# gene


# # 以下是为了人基因转鼠基因，人单细胞不用跑下面三行
# # ！！！！！！人物种不要多跑！！！！！！！！！！！
# library(Hmisc)
# gene=toupper(gene)
# gene=capitalize(gene)





load(file = 'GSE57338包括pdata用这个.Rdata')


train <- exprSet1
group_list1 <- group_list


load(file = 'GSE145154.Rdata')
test <- exprSet1
group_list2 <- group_list
dim(test)
group_list2

# ## 读取day1中整理好的文件
# rt=read.table('merge.normalzie.txt', header=T, sep="\t", check.names=F, row.names=1)
#
# load('GSE16561.Rdata')
# pdata1=pdata
# group_list1=c(rep('IS',39),rep('CT',24))
# colnames(pdata1)
#
# load('GSE58294.Rdata')
# pdata2=pdata
# group_list2=c(rep('CT',23),rep('IS',69))
#
#
# colnames(rt)
# ### 确定好训练和验证集，需要尝试
# train=rt[,grep('GSE58294',colnames(rt))]
# test=rt[,grep('GSE16561',colnames(rt))]

train=train[rownames(train) %in% gene,]
test=test[rownames(test) %in% gene,]
dim(test)

group1 <- data.frame(group=group_list1)
rownames(group1) <- colnames(train)
group2 <- data.frame(group=group_list2)
rownames(group2) <- colnames(test)

pheatmap::pheatmap(train, scale = "row",cluster_cols = F, annotation_col = group1,
                   color = colorRampPalette(c("#324c63","white","#ff40ff"))(50))
pheatmap::pheatmap(test, scale = "row",cluster_cols = F, annotation_col = group2,
                   color = colorRampPalette(c("#324c63","white","#ff40ff"))(50))


train=as.data.frame(t(train))

###！！！！！！！！！！！！
train$group=group_list1
table(train$group)
train$group=ifelse(train$group =='treat',1,0)

test=as.data.frame(t(test))
test$group=group_list2
table(test$group)
test$group=ifelse(test$group=='HF',1,0)

####### 单因素逻辑回归
rt=train[,colSums(train)!=0]
outTab=data.frame()
pFilter=0.05
for(gene in colnames(rt[,1:(ncol(rt)-1)])){
  set.seed(123)
  glm=glm(group ~ rt[,gene], data = rt,family= binomial(link='logit'))
  glmSummary = summary(glm)
  OR=exp(glm$coefficients)[2]
  OR_CI=exp(confint(glm,level = 0.95))
  OR.95L=OR_CI[2,1]
  OR.95H=OR_CI[2,2]
  glmP=glmSummary$coefficients[,"Pr(>|z|)"][2]
  if(glmP<pFilter){
    outTab=rbind(outTab,
                 cbind(gene=gene,
                       OR=OR,
                       OR.95L=OR.95L,
                       OR.95H=OR.95H,
                       pvalue=glmP))
  }
}

length(outTab$gene)

outTab$gene

library(openxlsx)


file.remove('./figures/单因素筛选基因列表.xlsx')
write.xlsx(outTab, "./figures/单因素筛选基因列表.xlsx",  rowNames =T)


rt=rt[,c('group',outTab$gene)]


# rt=rt[,c('group',gene)]


library(caret)
folds <- createMultiFolds(y=rt$group,k=10,times=10)###10折200次划分,这里为节省时间做10次


#
# library(dplyr)
# library(glmnet)
# #设定循环次数
# gene = c()
# for(i in 1:100){
#   #设定训练组和验证组
#   train<- rt[ folds[[i]],]
#   test <- rt[-folds[[i]],]
#   x = as.matrix(train[,c(-1)])
#   y = train[,c(1)]
#   #训练组建模
#   fitCV <- cv.glmnet(x, y, family = "binomial",
#                      type.measure = "class",
#                      nfolds = 5)
#   fit <- glmnet(x, y, family = "binomial", alpha = 1) # make sure alpha = 1
#   # get the coef
#   coef.min = coef(fitCV, s = "lambda.min")
#   # coef.min
#   cv_coef = coef.min %>% as.matrix() %>% as.data.frame %>% filter(abs(.)>0)
#   gene = c(gene,row.names(cv_coef))
# }
#
# sort_gene = data.frame(table(gene)) %>%
#   arrange(desc(Freq)) %>%
#   filter(Freq>10)
#
#
# sort_gene


###### LASSO回归降维
library(dplyr)
library(glmnet)
x=as.matrix(rt[,2:ncol(rt)])
y=unlist(rt$group)
set.seed(12345)
fit <- glmnet(x,y,alpha=1,family="binomial")
plot(fit,xvar="lambda",label=F)
set.seed(123456)
cvfit <- cv.glmnet(x,y,alpha=1,family="binomial",nfolds = 10)

plot(cvfit)
coef =coef(fit,s = cvfit$lambda.min)
index = which(coef !=0)
actCoef = coef[index]
lassoGene = row.names(coef)[index]
geneCoef = cbind(Gene=lassoGene,Coef=actCoef)
geneCoef
genes=geneCoef[,1][2:nrow(geneCoef)]

genes

save(genes,file ='./figures/lasso_genes.Rdata')

genes




load('figures/lasso_genes.Rdata')
gene <- genes
train=train[rownames(train) %in% gene,]
test=test[rownames(test) %in% gene,]
dim(test)

group1 <- data.frame(group=group_list1)
rownames(group1) <- colnames(train)
group2 <- data.frame(group=group_list2)
rownames(group2) <- colnames(test)

pheatmap::pheatmap(train, scale = "row",cluster_cols = F, annotation_col = group1,
                   show_colnames = F, height = 4,

                   filename = "./figures/hd_day5_lassogenes_heatmap.pdf",

                   color = colorRampPalette(c("#324c63","white","#ff40ff"))(50))

dev.off()
pheatmap::pheatmap(test, scale = "row",cluster_cols = F, annotation_col = group2,
                   color = colorRampPalette(c("#324c63","white","#ff40ff"))(50))







source('./57.lasso and ROC  一键lasso分析/一键lasso分析/run_lasso_xc.R')


run_lasso_binomial_xc(rt)






## 训练集
train=train[,c('group',genes)]
## 标化
train2=train[,2:ncol(train)]
train2=apply(train2,1, scale)
rownames(train2)=colnames(train)[2:ncol(train)]
train2=as.data.frame(t(train2))
train2$group=train$group
train=train2

# 验证集
test=test[,c('group',genes)]
##标化
test2=test[,2:ncol(test)]
test2=apply(test2,1, scale)
rownames(test2)=colnames(test)[2:ncol(test)]
test2=as.data.frame(t(test2))
test2$group=test$group
test=test2

#install.packages('mlr3verse')
## 主要这个包
#devtools::install_github("mlr-org/mlr3extralearners",)


## 这个可以不装
#devtools::install_local('d:/R/mlr3extralearners-main.zip')


library(mlr3verse)


train$group=factor(train$group)
test$group=factor(test$group)

dim(test)


# 创建一个任务

task_train = as_task_classif(
  train,
  target = "group"
)

task_test= as_task_classif(
  test,
  target = "group"
)

learners = list(
  learner_logreg = lrn("classif.log_reg", predict_type = "prob",
                       predict_sets = c("train", "test")),
  learner_lda = lrn("classif.lda", predict_type = "prob",
                    predict_sets = c("train", "test")),
  learner_svm = lrn("classif.svm", predict_type = "prob",
                    predict_sets = c("train", "test")),
  learner_nb = lrn("classif.naive_bayes", predict_type = "prob",
                   predict_sets = c("train", "test")),
  learner_knn = lrn("classif.kknn", scale = FALSE,
                    predict_type = "prob",predict_sets = c("train", "test")),
  learner_rpart = lrn("classif.rpart",
                      predict_type = "prob",predict_sets = c("train", "test")),
  learner_rf = lrn("classif.ranger", num.trees = 1000,
                   predict_type = "prob",predict_sets = c("train", "test"))
)


resampling_outer = rsmp("repeated_cv",#K次N折
                        folds = 5,#K折
                        repeats =10)

task_train$col_roles$stratum = task_train$target_names

design = benchmark_grid(
  tasks = task_train,
  learners =learners,
  resamplings = resampling_outer
)

bmr = benchmark(design, store_models = FALSE) ## 耗时较长

measures <- msrs(c("classif.auc","classif.acc","classif.bbrier"))

bmr_res <- bmr$aggregate(measures)

measures = list(
  msr("classif.auc", predict_sets = "train", id = "auc_train"),
  msr("classif.auc", id = "auc_test")
)

tab = bmr$aggregate(measures)
tab_1 = tab[,c('learner_id','auc_train','auc_test')]
print(tab_1)

tab2 = bmr$aggregate(msrs(c('classif.auc', 'classif.sensitivity','classif.specificity',
                             'classif.fnr', 'classif.fpr')))
tab2 = tab2[,c('learner_id','classif.auc','classif.sensitivity','classif.specificity',
               'classif.fnr', 'classif.fpr')]
print(tab2)

res1 <- bmr$score(measures = msr("classif.auc"))

data <- res1 %>%
  dplyr::select(learner_id,classif.auc)
data

library(ggplot2)
library(ggpubr)
###生成图形

p1<-ggplot(data, aes(x=`learner_id`, y=`classif.auc`, fill=factor(`learner_id`))) +
  geom_boxplot(
    alpha=0.8, #图形透明度
    color="#3D48DB", #边框颜色
    notch=F, #是否显示凹槽
    notchwidth=0.8, #凹槽宽度
    outlier.colour = "#FF6600", #离群点边界颜色
    outlier.fill = "black", #离群点填充颜色
    outlier.size = 1.5 #离群点大小
  ) +
  geom_jitter(color="#00000064", size=0.1)+
  stat_boxplot(geom = "errorbar",width=0.15,aes(color=learner_id))

library(ggpubr)

p2<-p1+
  scale_fill_brewer(palette="Dark2")+ #选择调色板，可选调色板包括："Set1", "Set2", "Set3", "Accent", "Blues", "Paired", "Pastel1", "Pastel2", "BuGn", "BuPu", "Dark2", "GnBu", "Greens", "Greys", "Oranges", "OrRd", "PuBu", "PuBuGn", "PuRd", "Purples", "RdPu", "Reds","YlGn", "YlGnBu", "YlOrBr", "YlOrRd", "BrBG", "PiYG", "PRGn", "PuOr", "RdBu", "RdGy", "RdYlBu", "RdYlGn", "Spectral"
  ggtitle("Machine Learn model comparison")+ #图片标题
  xlab("Learn Model")+ #X轴标签
  ylab("Mean AUC")+ #Y轴标签
  guides(fill=guide_legend(title = "group"))+ #图例标题
  theme(
    axis.title=element_text(size=10,face="bold",color="black"), #坐标轴标签大小、颜色
    axis.text = element_text(size=10,face="bold",color="black"), #坐标刻度大小、颜色
    # axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = .5),
    plot.title = element_text(size=20,hjust = 0.5,face="bold",color = "black"), #标题大小、颜色
    panel.grid.major = element_line(colour = "#E8DADA",size=0.5,linetype=1), #网格线颜色、粗细、类型
    axis.line = element_line(colour = "black"), #坐标轴颜色
    panel.background = element_rect(fill="white"), #背景填充颜色
    legend.title = element_text(size=15,face="bold",color="black"), #图例标题大小、颜色
    legend.text = element_text(size=15,color="black"), #图例大小、颜色
    legend.position = "right" #图例位置，可选项为："none", "right","left","top","bottom"
  )+rotate_x_text(45)+
  scale_y_continuous(limits = c(0.4, 1.0),breaks = seq(0.4,1.0,by=0.05))

p2

pdf(file="./figures/hd_day4_mlr3_所有模型对比 p2.pdf", width=7.5, height=5)
# autoplot(prediction,type = "roc",auc=T)

p2
dev.off()



pdf(file="./figures/hd_day4_mlr3_所有模型.pdf", width=5, height=5)
# autoplot(prediction,type = "roc",auc=T)

autoplot(bmr,type = "roc")+
  scale_color_discrete() +
  theme_bw()
dev.off()

#### 选择*支持向量机*模型

# final_model=learners$learner_svm
final_model=learners$learner_nb
# final_model=learners$learner_logreg
# final_model=learners$learner_rf
#SVM训练
final_model$train(task_train)


#测试！第二个数据集独立外部测试
prediction <- final_model$predict(task_test)

pdf(file="./figures/hd_day4_mlr3_所有模型.pdf", width=5, height=5)
autoplot(prediction,type = "roc",auc=T)
dev.off()



prediction_tab=as.data.table(prediction)

pdf(file="./figures/hd_day4_mlr3.pdf", width=4.5, height=4.5)
pROC::plot.roc(test$group,prediction_tab$prob.1,print.auc=T,main="learner_nb")
dev.off()

pROC::plot.roc(test$group,prediction_tab$prob.1,print.auc=T)



names(learners)


plist <- list()




# par(mfrow=c(1,7))


dev.off()
for (i in names(learners)) {

  pdf(file=paste0("./figures/拼图mlr3_",i,'.pdf'), width=5, height=5)
  #### 选择*支持向量机*模型

  # final_model=learners$learner_svm
  # final_model=learners$learner_logreg
  final_model=learners[[i]]
  #SVM训练
  final_model$train(task_train)


  #测试！第二个数据集独立外部测试
  prediction <- final_model$predict(task_test)


  # autoplot(prediction,type = "roc",auc=T)+theme_bw()

  prediction_tab=as.data.table(prediction)
  pROC::plot.roc(test$group,prediction_tab$prob.1,print.auc=T,
                 xlim=c(0,1), ylim=c(0,1),
                 main=i)
  dev.off()
}
dev.off()


# ?plot.roc
#
# library(cowplot)
# pdf(file="./figures/拼图.pdf", width=15, height=4)
# p_top <- plot_grid(plotlist=plist, ncol = 5)
# print(p_top)
# dev.off()





