library(IOBR)
load('lasso_genes.Rdata')
genes
gene <- genes


# rt=read.table('merge.normalzie.txt', header=T, sep="\t", check.names=F, row.names=1)
#
# ## 测试集
# load('GSE16561.Rdata')
# pdata1=pdata
# group_list1=c(rep('IS',39),rep('CT',24))
#
#
# ## 训练集
# load('GSE58294.Rdata')
# pdata2=pdata
# group_list2=c(rep('CT',23),rep('IS',69))
#
#
# train=rt[,grep('GSE58294',colnames(rt))]
# test=rt[,grep('GSE16561',colnames(rt))]



load(file = 'GSE57338包括pdata用这个.Rdata')

mergenormalize <- readRDS('merge.normalzie.RDS')
colnames(mergenormalize)

group_list
train <- mergenormalize[,1:length(group_list)]
dim(train)
group_list1 <- group_list


load(file = 'GSE145154.Rdata')
colnames(mergenormalize)
(ncol(mergenormalize)-length(group_list))
test <- mergenormalize[,(ncol(mergenormalize)-length(group_list)+1):ncol(mergenormalize)]
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

# train=train[rownames(train) %in% gene,]
# test=test[rownames(test) %in% gene,]
# dim(test)



# train=as.data.frame(t(train))
#
# ###！！！！！！！！！！！！
# train$group=group_list1
# table(train$group)
# train$group=ifelse(train$group =='treat',1,0)
#
# test=as.data.frame(t(test))
# test$group=group_list2
# table(test$group)
# test$group=ifelse(test$group=='HF',1,0)



############# 训练集######################################
mcp_train=IOBR::deconvo_mcpcounter(eset = train)
# mcp_train=IOBR::deconvo_cibersort(eset = train, arrays = T)
# ?IOBR::deconvo_cibersort
rownames(mcp_train)=mcp_train$ID

mcp_train$ID=NULL


# progeny的14条通路
# library(progeny)
# # ?progeny
# train1 <- as.matrix(train)
# pbmc <- progeny(train1, scale=FALSE, organism="Human", top=500, perm=1,
#                 return_assay = TRUE)
# mcp_train <- pbmc



exp_train=train[rownames(train) %in% genes,]
mcp_train=as.data.frame(t(mcp_train))

# mcp_train=mcp_train[! rownames(mcp_train) %in% c("Fibroblasts_MCPcounter","Endothelial_cells_MCPcounter" ),]

### 训练集制作第一个#######
mcp_a=mcp_train[,1]
exp_a=exp_train[,1]

#我们这里只有8个细胞类型

termCounts = 10
#8个细胞类型和13个基因的矩阵
table_a=matrix(rep(1, length(genes)*termCounts),ncol =termCounts,nrow =length(genes))


rownames(table_a)=rownames(exp_train)
colnames(table_a)=rownames(mcp_train)


for (i in 1:ncol(table_a)) {
  for (j in 1:nrow(table_a)) {
    table_a[j,i]=mcp_a[i] / exp_a[j] #除法
  }
}

table_a
pheatmap::pheatmap(table_a,cluster_rows = F,cluster_cols = F,display_numbers = T,
                   color = c("white", "red"))

colnames(table_a) <- gsub("_MCPcounter","", colnames(table_a))
p1 <- pheatmap::pheatmap(table_a,cluster_rows = F,cluster_cols = F,display_numbers = T,
                   main = "MCPcounter\ntraining corhort",
                   # color = c("white", "red"))
                   color = colorRampPalette(c("white", "purple"))(50),
                   filename = './figures/hd_day6_train.pdf')

p1 <- ggplotify::as.ggplot(p1)
p1


### 训练集制作第二到最后一个#########
ncol(mcp_train) #所有病人
for (k in 2:ncol(mcp_train)) {
  mcp_k=mcp_train[,k]
  exp_k=exp_train[,k]
  table_k=matrix(rep(1,length(genes)*termCounts),ncol =termCounts,nrow = nrow(exp_train))
  rownames(table_k)=rownames(exp_train)
  colnames(table_k)=rownames(mcp_train)
  for (i in 1:ncol(table_k)) {
    for (j in 1:nrow(table_k)) {
      table_k[j,i]=mcp_k[i] / exp_k[j]
    }
  }
  table_a=rbind(table_a,table_k)
}


library(keras)


x_train <- array_reshape(table_a, dim = c(ncol(mcp_train), length(genes), termCounts))

pheatmap::pheatmap(x_train[60,,],cluster_rows = F,cluster_cols = F)




############# 测试集######################################
mcp_test=IOBR::deconvo_mcpcounter(eset = test)
rownames(mcp_test)=mcp_test$ID


mcp_test$ID=NULL





# library(progeny)
# # ?progeny
# train1 <- as.matrix(test)
# pbmc <- progeny(train1, scale=FALSE, organism="Human", top=500, perm=1,
#                 return_assay = TRUE)
# mcp_test <- pbmc








exp_test=test[rownames(test) %in% genes,]
mcp_test=as.data.frame(t(mcp_test))
# mcp_test=mcp_test[! rownames(mcp_test) %in% c("Fibroblasts_MCPcounter","Endothelial_cells_MCPcounter" ),]


### 训练集制作第一个#######
mcp_a=mcp_test[,1]
exp_a=exp_test[,1]

table_a=matrix(rep(1, length(genes)*termCounts),ncol =termCounts,nrow =length(genes))


rownames(table_a)=rownames(exp_test)
colnames(table_a)=rownames(mcp_test)


for (i in 1:ncol(table_a)) {
  for (j in 1:nrow(table_a)) {
    table_a[j,i]=mcp_a[i] / exp_a[j]
  }
}
pheatmap::pheatmap(table_a,cluster_rows = F,cluster_cols = F)


pheatmap::pheatmap(table_a,cluster_rows = F,cluster_cols = F,display_numbers = T,
                   color = c("white", "red"))

colnames(table_a) <- gsub("_MCPcounter","", colnames(table_a))
p2 <- pheatmap::pheatmap(table_a,cluster_rows = F,cluster_cols = F,display_numbers = T,
                         main = "MCPcounter\ntesting corhort",
                         # color = c("white", "red"))
                         color = colorRampPalette(c("white", "purple"))(50),
                         filename = './figures/hd_day6_test.pdf')

p2 <- ggplotify::as.ggplot(p2)
p2

dev.off()

library(patchwork)
p1+p2

pdf(file="./figures/hd_day6_拼图 深度学习.pdf", width=10, height=5)
print(p1+p2)
dev.off()



### 训练集制作第二到最后一个#########
for (k in 2:ncol(mcp_test)) {
  mcp_k=mcp_test[,k]
  exp_k=exp_test[,k]
  table_k=matrix(rep(1,length(genes)*termCounts),ncol =termCounts,nrow = nrow(exp_test))
  rownames(table_k)=rownames(exp_test)
  colnames(table_k)=rownames(mcp_test)
  for (i in 1:ncol(table_k)) {
    for (j in 1:nrow(table_k)) {
      table_k[j,i]=mcp_k[i] / exp_k[j]
    }
  }
  table_a=rbind(table_a,table_k)
}


library(keras)
citation('keras')


x_test <- array_reshape(table_a, dim = c(ncol(mcp_test), length(genes), termCounts))

pheatmap::pheatmap(x_test[40,,],cluster_rows = F,cluster_cols = F)




###############卷积，注意修改数目 ######################！！！！！！！！！！！

### 基因个数你不一定是11个
length(genes)

dim(x_train)
x_train <- array_reshape(x_train, dim = c(231,length(genes), termCounts, 1))

dim(x_test)
x_test <- array_reshape(x_test, dim = c(67, length(genes), termCounts, 1))

table(group_list1)
table(group_list2)
y_train <- ifelse(group_list1=='treat',1,0)
y_test <- ifelse(group_list2=='HF',1,0)

set.seed(123456)
model <- keras_model_sequential()
model %>%
  layer_conv_2d(filters = 32, kernel_size = c(3, 3), padding = 'same',  input_shape = c(length(genes),
                                                                                        termCounts, 1)) %>%
  layer_activation('relu') %>%
  # layer_max_pooling_2d(pool_size=c(2, 2), strides=c(2, 2)) %>%
  layer_conv_2d(filters = 16, kernel_size = c(2, 2),
                dilation_rate = c(1,1), activation = 'softplus', padding = 'same') %>%
  layer_max_pooling_2d(pool_size=c(2, 2)) %>%
  layer_flatten() %>%
  layer_dense(64, activation = 'relu') %>%
  layer_dropout(0.5) %>%
  layer_dense(1, activation = 'sigmoid')

summary(model)

model %>% compile(
  loss = 'binary_crossentropy',
  optimizer ='adam',
  metrics = c('accuracy')
)

#这里就是在训练模型
history <- model %>% fit(
  x_train, y_train,
  epochs = 200
)

p1 <- plot(history)
p1

dir.create("./figures")
pdf(file="./figures/hd_day6_深度学习拼图_模型训练图.pdf", width=4, height=4)
print(plot(history))
dev.off()

saveRDS(model, file = "./hd_day6_deep_learning_model.RDS")

model <- readRDS('./hd_day6_deep_learning_model.RDS')

prob=model %>% predict(x_train)
pROC::plot.roc(factor(y_train),as.numeric(prob),print.auc=T,print.thres=T)

prob_test=model %>% predict(x_test)
pROC::plot.roc(factor(y_test),as.numeric(prob_test),
               print.auc=T,print.thres=T)


#
###在最佳cutoff下 ，看测试集表现
# prob_bi=ifelse(prob_test>0.750,1,0)
prob_bi=ifelse(prob_test>0.5,1,0)
prob_test
prob_bi=ifelse(prob_test>0.385,1,0)


prob_bi=ifelse(prob_test>0.053,1,0)  #ROC曲线上的那个点的值

table(prob_bi,y_test)


pdf(file="./figures/hd_day6_深度学习拼图.pdf", width=7, height=4)
par(mfrow=c(1,2))
# plot(history)
pROC::plot.roc(factor(y_train),as.numeric(prob),print.auc=T,print.thres=T,
               main="training set")
pROC::plot.roc(factor(y_test),as.numeric(prob_test),
               print.auc=T,print.thres=T,
               main="testing set")
dev.off()
# ?plot.roc




# mn_pred2 <- predict(mn_res, test_expr)
prob_test=model %>% predict(x_test)
mn_pred2 <- prob_bi

factor(y_test)
y_test
as.numeric(prob_test)
##########################################    多分类模型的混淆矩阵           ############################               #
########################### ##########################################
mn_pred2
test_y <- y_test

dir.create("./figures")
library(randomForest)
library(caret)
library(pROC)
library(caret)
library(multiROC)
# 计算测试集的混淆矩阵
mn_pred2 <- as.numeric(mn_pred2)
mn_pred2 <- as.factor(mn_pred2)
as.factor(test_y)
confusion_mn <- confusionMatrix(mn_pred2, as.factor(test_y))
print(confusion_mn)


# 绘制混淆矩阵表格
#################################################
confusion_table <- confusion_mn$table
# 将混淆矩阵转换为适合ggplot2的数据格式
confusion_matrix_df <- as.data.frame.matrix(confusion_table)
confusion_matrix_df$actual <- rownames(confusion_matrix_df)
library(tidyr)
confusion_matrix_df <- gather(confusion_matrix_df, key = "predicted", value = "count", -actual)
# 添加每组人数的文本标签
confusion_matrix_df$count_label <- paste("n =", confusion_matrix_df$count)

# 绘制混淆矩阵热力图
p_confusion <- ggplot(data = confusion_matrix_df, aes(x = predicted, y = actual, fill = count)) +
  geom_tile() +
  # scale_fill_gradient(low = "white", high = "steelblue") +
  scale_fill_gradient(low = "white", high = "#a363ae") +
  geom_text(aes(label = count_label), size = 3, color = "black") +
  labs(title = "Confusion Matrix", x = "Actual", y = "Predicted") +
  theme_minimal()
p_confusion

# 提取区分结局为1的模型性能评价指标：精确率、召回率、灵敏度、特异性和F1score
# precision <- confusion_mn$byClass[1, "Precision"]
# recall <- confusion_mn$byClass[1, "Recall"]
# f1_score <- confusion_mn$byClass[1, "F1"]
# specificity <- confusion_mn$byClass[1, "Specificity"]
# sensitivity <- confusion_mn$byClass[1, "Sensitivity"]



aa <- confusion_mn$byClass
aa <- as.data.frame(aa)
# 输出评估结果
# cat("神经网络模型的评估结果：\n")
# cat("精确率:", precision, "\n")
# cat("召回率:", recall, "\n")
# cat("特异度:", specificity, "\n")
# cat("灵敏度:", sensitivity, "\n")
# cat("F1score:", f1_score, "\n")
# # 将指标值存储到数据框中
# metric_df <- data.frame(Metric = c("precision", "recall", "F1_score", "specificity", "sensitivity"),
#                         Value = c(precision, recall, f1_score, specificity, sensitivity))



metric_df <- as.data.frame(aa[c(1,2,5,6,7),])

metric_df <- data.frame(Metric = c("sensitivity", "specificity", "precision", "recall", "F1_score" ),
                        Value = as.numeric(aa[c(1,2,5,6,7),]))


# 绘制模型性能评价的柱状图
p_perform <- ggplot(data = metric_df, aes(x = Metric, y = Value)) +
  geom_bar(stat = "identity", fill = "#b072c3") +
  geom_text(aes(label = round(Value, 2)), vjust = 0.5,
            hjust = -0.05,
            # size = 5,
            color = "black") +
  coord_flip() +
  labs(x = NULL, y = "value", title = "Model performance evaluation",
       subtitle = NULL, caption = NULL) +
  # theme_classic() +
  theme_dark() +
  theme(plot.title = element_text(hjust = 0.5), #size = 18, face = "bold"
        # plot.subtitle = element_text(hjust = 0.5, size = 14),
        # plot.caption = element_text(hjust = 1, size = 12,
        #                             margin = margin(t = 10, r = 0, b = 0, l = 0)),
        # axis.text.x = element_text(size = 14, face = "bold", color = "black"),
        # axis.text.y = element_text(size = 14, face = "bold", color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

p_perform



library(patchwork)
dir.create("./figures")
pdf(file="./figures/hd_day6_深度学习拼图_矩阵.pdf", width=10, height=3.5)
print(p1+p_confusion+p_perform)
dev.off()
# p_confusion+p_perform+p_roc+p_pr

#
# tiff("./figures/拼图_multinomial LASSO models.tiff", width = 13, height = 12, units = 'in', res = 300)
# print(p_confusion+p_perform+p_roc+p_pr)
# dev.off()
#
#
# tiff("./figures/拼图_multinomial LASSO models_roc.tiff", width = 6, height = 5, units = 'in', res = 300)
# print(p_roc)
# dev.off()








####################### 神经网络NN，注意修改数目##########################
x_train <- array_reshape(x_train, dim = c(231, length(genes), termCounts))
x_test <- array_reshape(x_test, dim = c(67, length(genes), termCounts))

#loading keras library
library(keras)
#loading the keras inbuilt mnist dataset
# converting a 2D array into a 1D array for feeding into the MLP and normalising the matrix
train_x <- array(x_train, dim = c(dim(x_train)[1], prod(dim(x_train)[-1])))
test_x <- array(x_test, dim = c(dim(x_test)[1], prod(dim(x_test)[-1])))
#converting the target variable to once hot encoded vectors using keras inbuilt function
train_y=y_train
test_y=y_test
#defining a keras sequential model
model <- keras_model_sequential()
#defining the model with 1 input layer[784 neurons], 1 hidden layer[784 neurons] with dropout rate 0.4 and 1 output layer[10 neurons]
#i.e number of digits from 0 to 9
model %>%
  layer_dense(units = length(genes)*termCounts, input_shape = length(genes)*termCounts) %>%
  layer_dropout(rate=0.4)%>%
  layer_activation(activation = 'relu') %>%
  layer_dense(units = 1) %>%
  layer_activation(activation = 'sigmoid')
#compiling the defined model with metric = accuracy and optimiser as adam.
model %>% compile(
  loss = 'binary_crossentropy',
  optimizer = 'adam',
  metrics = c('accuracy')
)
#fitting the model on the training dataset


history <- model %>% fit(
  train_x, train_y,
  epochs = 200,batch_size = 10
)


# train_x
# train_y
# x_train
# y_train
# history <- model %>% fit(
#   x_train, y_train,
#   epochs = 200
# )

plot(history)

prob=model %>% predict(train_x)
dev.off()
pROC::plot.roc(factor(train_y),as.numeric(prob),print.auc=T,print.thres=T)

prob_bi=ifelse(prob>0.5,1,0)

table(prob_bi,y_train)


prob_test=model %>% predict(test_x)
pROC::plot.roc(factor(test_y),as.numeric(prob_test),
               print.auc=T,print.thres=T,smooth=T)

prob_bi=ifelse(prob_test>0.5,1,0)

table(prob_bi,y_test)





p1 <- plot(history)
p1

dir.create("./figures/NN")

pdf(file="./figures/NN/hd_day6_深度学习拼图_模型训练图.pdf", width=4, height=4)
print(plot(history))
dev.off()

saveRDS(model, file = "./figures/NN/deep_learning_model.RDS")

prob=model %>% predict(train_x)
pROC::plot.roc(factor(y_train),as.numeric(prob),print.auc=T,print.thres=T)

prob_test=model %>% predict(test_x)
pROC::plot.roc(factor(y_test),as.numeric(prob_test),
               print.auc=T,print.thres=T)


#
###在最佳cutoff下 ，看测试集表现
prob_bi=ifelse(prob_test>0.874,1,0)
# prob_bi=ifelse(prob_test>0.5,1,0)

table(prob_bi,y_test)


pdf(file="./figures/NN/hd_day6_深度学习拼图NN.pdf", width=7, height=4)
par(mfrow=c(1,2))
# plot(history)
pROC::plot.roc(factor(y_train),as.numeric(prob),print.auc=T,print.thres=T,
               main="train set")
pROC::plot.roc(factor(y_test),as.numeric(prob_test),
               print.auc=T,print.thres=T,
               main="test set")
dev.off()
?plot.roc




# mn_pred2 <- predict(mn_res, test_expr)
prob_test=model %>% predict(test_x)
mn_pred2 <- prob_bi

factor(y_test)
y_test
as.numeric(prob_test)
##########################################    多分类模型的混淆矩阵           ############################               #
########################### ##########################################
mn_pred2
test_y <- y_test

dir.create("./figures")
library(randomForest)
library(caret)
library(pROC)
library(caret)
library(multiROC)
# 计算测试集的混淆矩阵
mn_pred2 <- as.numeric(mn_pred2)
mn_pred2 <- as.factor(mn_pred2)
as.factor(test_y)
confusion_mn <- confusionMatrix(mn_pred2, as.factor(test_y))
print(confusion_mn)


# 绘制混淆矩阵表格
#################################################
confusion_table <- confusion_mn$table
# 将混淆矩阵转换为适合ggplot2的数据格式
confusion_matrix_df <- as.data.frame.matrix(confusion_table)
confusion_matrix_df$actual <- rownames(confusion_matrix_df)
library(tidyr)
confusion_matrix_df <- gather(confusion_matrix_df, key = "predicted", value = "count", -actual)
# 添加每组人数的文本标签
confusion_matrix_df$count_label <- paste("n =", confusion_matrix_df$count)

# 绘制混淆矩阵热力图
p_confusion <- ggplot(data = confusion_matrix_df, aes(x = predicted, y = actual, fill = count)) +
  geom_tile() +
  # scale_fill_gradient(low = "white", high = "steelblue") +
  scale_fill_gradient(low = "white", high = "#a363ae") +
  geom_text(aes(label = count_label), size = 3, color = "black") +
  labs(title = "Confusion Matrix", x = "Actual", y = "Predicted") +
  theme_minimal()
p_confusion

# 提取区分结局为1的模型性能评价指标：精确率、召回率、灵敏度、特异性和F1score
# precision <- confusion_mn$byClass[1, "Precision"]
# recall <- confusion_mn$byClass[1, "Recall"]
# f1_score <- confusion_mn$byClass[1, "F1"]
# specificity <- confusion_mn$byClass[1, "Specificity"]
# sensitivity <- confusion_mn$byClass[1, "Sensitivity"]



aa <- confusion_mn$byClass
aa <- as.data.frame(aa)
# 输出评估结果
# cat("神经网络模型的评估结果：\n")
# cat("精确率:", precision, "\n")
# cat("召回率:", recall, "\n")
# cat("特异度:", specificity, "\n")
# cat("灵敏度:", sensitivity, "\n")
# cat("F1score:", f1_score, "\n")
# # 将指标值存储到数据框中
# metric_df <- data.frame(Metric = c("precision", "recall", "F1_score", "specificity", "sensitivity"),
#                         Value = c(precision, recall, f1_score, specificity, sensitivity))



metric_df <- as.data.frame(aa[c(1,2,5,6,7),])

metric_df <- data.frame(Metric = c("sensitivity", "specificity", "precision", "recall", "F1_score" ),
                        Value = as.numeric(aa[c(1,2,5,6,7),]))


# 绘制模型性能评价的柱状图
p_perform <- ggplot(data = metric_df, aes(x = Metric, y = Value)) +
  geom_bar(stat = "identity", fill = "#b072c3") +
  geom_text(aes(label = round(Value, 2)), vjust = 0.5,
            hjust = -0.05,
            # size = 5,
            color = "black") +
  coord_flip() +
  labs(x = NULL, y = "value", title = "Model performance evaluation",
       subtitle = NULL, caption = NULL) +
  # theme_classic() +
  theme_dark() +
  theme(plot.title = element_text(hjust = 0.5), #size = 18, face = "bold"
        # plot.subtitle = element_text(hjust = 0.5, size = 14),
        # plot.caption = element_text(hjust = 1, size = 12,
        #                             margin = margin(t = 10, r = 0, b = 0, l = 0)),
        # axis.text.x = element_text(size = 14, face = "bold", color = "black"),
        # axis.text.y = element_text(size = 14, face = "bold", color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

p_perform



library(patchwork)
dir.create("./figures")
pdf(file="./figures/NN/hd_day6_深度学习拼图_矩阵NN.pdf", width=12, height=4)
print(p1+p_confusion+p_perform)
dev.off()



