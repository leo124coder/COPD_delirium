#加载包
library(tidyverse)
library(forcats)
library(mice)
library(glmnet)
#数据预处理
copd_v <- read.csv("copd_ven.csv")
copd_v$sven_hours <- as.numeric(copd_v$sven_hours)
copd_v <- na.omit(copd_v)
copd_v$ventilation_status <- factor(copd_v$ventilation_status)
copd_v <- copd_v %>%
  pivot_wider(names_from = ventilation_status, values_from = sven_hours, values_fill = list(sven_hours = 0))
copd <- copd3
copd4 <- read.csv("copd.csv")
copd <- copd[, -c(1:4, 11:12, 27:28, 31:32, 37:40, 45:55, 82:91)]
copd$insurance <- factor(copd$insurance,levels = c("Other","Medicare","Medicaid"),labels = c("1","2","3"))
copd$marital_status <- factor(copd$marital_status)
copd_v$stay_id <- as.numeric(copd_v$stay_id)
copdall <- copd %>% full_join(copd_v,by="stay_id")
copdall$gender <- factor(copdall$gender,labels = c("1","2"),levels = c("F","M"))
copdall$delirium <- factor(copdall$delirium,labels = c("0","1"),levels = c("0","1"))

#将种族变量重编码
copd_all <- copdall %>% mutate(race=fct_collapse(race,
                                                 other=c("AMERICAN INDIAN/ALASKA NATIVE","ASIAN","ASIAN - CHINESE","ASIAN - SOUTH EAST ASIAN","OTHER",
                                                         "NATIVE HAWAIIAN OR OTHER PACIFIC ISLANDER",
                                                         "PATIENT DECLINED TO ANSWER","PORTUGUESE","SOUTH AMERICAN","UNABLE TO OBTAIN","UNKNOWN"),
                                                 white=c("WHITE","WHITE - EASTERN EUROPEAN","WHITE - OTHER EUROPEAN","WHITE - RUSSIAN"),
                                                 black=c("BLACK/AFRICAN","BLACK/AFRICAN AMERICAN","BLACK/CAPE VERDEAN","BLACK/CARIBBEAN ISLAND"),
                                                 hispanic=c("HISPANIC/LATINO - COLUMBIAN","HISPANIC/LATINO - CUBAN","HISPANIC/LATINO - DOMINICAN",
                                                            "HISPANIC/LATINO - GUATEMALAN","HISPANIC/LATINO - PUERTO RICAN")))
copd_all$stay_id <- NULL

#缺失值插补
set.seed(1)
imp_data1 <- mice(copd_all, #数据集
                  method = "rf", #采用pmm插补
                  m=5, # 5次插补
                  printFlag = FALSE #不显示历史记录
)
x1 <- complete(imp_data1,action=01)#提取第一次插补的完整数据
x1 <- x1 %>% mutate(marital_status=fct_collapse(marital_status,
                                                married=c("","MARRIED"),
                                                divorced="DIVORCED",
                                                single="SINGLE",
                                                widowed="WIDOWED"))
#变量筛选方法一：LASSO回归+交叉验证
aa <- x1
x <- data.matrix(aa[, c(1:60,62:67)])
y <- data.matrix(aa[, 61])

#2.lasso回归
fit<- glmnet(x, y,family="binomial",alpha=1)

#3.可视化结果
plot(fit, xvar="lambda", label=T)

#4.使用交叉验证来确定最佳的lambda。
lassoCV <- cv.glmnet(x, y,
                     family="binomial",
                     alpha = 1)
plot(cv.fit)

#5.显示两个惩罚值（调优系数）下的变量及其回归系数
ridge.coef1 <- predict(fit, s=cv.fit$lambda.1se, type = "coefficients");ridge.coef1
ridge.coef2 <- predict(fit, s=cv.fit$lambda.min, type = "coefficients");ridge.coef2

# 提取cv.glmnet数据
cv.df <- data.frame(lambda =lassoCV$lambda,#交叉验证中的lambda 
                    mse =lassoCV$cvm,#mse
                    sd=lassoCV$cvsd)#sd
p5<-ggplot(cv.df, aes(log(lambda), mse)) +
  geom_point() +#点图
  geom_errorbar(aes(ymin = mse - sd, ymax = mse + sd), width = 0.1) +#添加误差棒
  scale_x_continuous(name = "Log lambda") +
  scale_y_continuous(name = "Mean Squared Error") +
  ggtitle("10-fold Cross-validation using Lasso Regression")+
  geom_vline(xintercept = log(cv.fit$lambda.min), linetype = "dashed", color = "#E41A1C",linewidth=1) +
  geom_vline(xintercept = log(cv.fit$lambda.1se), linetype = "dashed", color = "#377EB8",linewidth=1) +
  theme_bw()
p5 
#变量筛选法二之最优子集
library(leaps)
library(scales)
library(gt)
leaps <- regsubsets(delirium ~ ., data = x1, really.big = TRUE)
broom::tidy(leaps) %>%     ##将最优子集结果可视化
  select(-`(Intercept)`) %>%
  rownames_to_column(var = "n_vars") %>%
  gt(rowname_col = "n_vars") %>%
  gt::data_color(
    columns = marital_statusdivorced:Tracheostomy,
    fn = col_numeric(
      palette = c("#fdae61", "#abdda4"),
      domain = c(0, 1)) 
  ) %>%
  gt::fmt_number(r.squared:mallows_cp, n_sigfig = 4)

##变量筛选法三之单因素分析
#批量单因素回归分析
varsU<-names(x1[,c(1:60,62:66)])#自变量

model = tibble(varsU) %>% ##建立循环
  mutate(model = map(x1[varsU], 
                     ~ glm(delirium ~ .x, data = x1, family = binomial()))) %>% 
  mutate(result = map(model, tidy),
         OR = map(model, ~ exp(coef(.x))),
         OR_ci = map(model, ~ exp(confint(.x)))) %>% 
  select(-model) %>% 
  unnest(c(result, OR, OR_ci))

model = model %>% 
  mutate(OR_ci %>% as_tibble()) %>% 
  select(-OR_ci) %>% 
  rename(LL = V1, UL = V2) %>% 
  mutate(across(term, ~ str_remove(.x, '.x'))) %>% 
  filter(if_all(term, ~ !.x=='(Intercept)')) %>% 
  mutate(`OR(95%CI)` = str_c(round(OR,2), ' (', round(LL,2), '-', round(UL,2), ')')) %>% 
  select(varsU, term, `OR(95%CI)`, p.value, OR, LL, UL, ) %>% 
  mutate(p.value = pvalue(p.value))
write.csv(model,file = "copd_model.csv")

##数据集划分与smote采样
library(caret)
library(DMwR)
###将插补后的数据的分类变量转化为哑元变量
dummy_data1=model.matrix(delirium~.,data = x1)[,-1]
data=data.frame(delirium=x1$delirium,dummy_data1)
###数据集划分训练集和测试集
set.seed(1)
train_id <- sample(1:nrow(data),0.7*nrow(data))###训练集id
train=data[train_id,]##原始训练集
test=data[-train_id,]##测试集

summary(train$delirium)##查看训练集的分类情况
summary(test$delirium)##查看验证集的分类情况

##对训练集进行数据过采样
smoteCopdRf <- SMOTE(delirium~.,train,perc.over = 250,perc.under = 330)##smoteCopdRf为过采样后的训练集
summary(smoteCopdRf$delirium)#查看采样后的训练集
#no yes 
#699 318 
###将分类变量转化为哑元变量
#dummy_data2=model.matrix(delirium~.,data = smoteCopdRf)[,-1]
#Smotedata=data.frame(delirium=smoteCopdRf$delirium,dummy_data2)##Smotedata为哑变量转换后的采样数据
###数据集划分，对过采样后的训练集进行划分，划分为训练集和验证集
set.seed(1)
train_id <- sample(1:nrow(smoteCopdRf),0.7*nrow(smoteCopdRf))###训练集id
Strain=smoteCopdRf[train_id,]##smote后的训练集
Stest=smoteCopdRf[-train_id,]##smote后的验证集
##全数据集
dataLasso <-data [,c(1:6,10,15:17,21,23,24,26,28,29,38,39,42,46,52:53,58,60,64:66,68,69,71)]
#lasso回归数据集
lassoSmoteTrain <- Strain[,c(1:6,10,15:17,21,23,24,26,28,29,38,39,42,46,52:53,58,60,64:66,68,69,71)]##lasso筛选后的smote训练集
lassoSmoteTest <- Stest[,c(1:6,10,15:17,21,23,24,26,28,29,38,39,42,46,52:53,58,60,64:66,68,69,71)]##lasso筛选后的验证集
lassoTest <- test[,c(1:6,10,15:17,21,23,24,26,28,29,38,39,42,46,52:53,58,60,64:66,68,69,71)]##lasso筛选后的测试集
#最优子集数据集
leapsSmoteTrain <- Strain[,c(1,10,15:17,20,29,39,43)]##lasso筛选后的smote训练集
leapsSmoteTest <- Stest[,c(1,10,15:17,20,29,39,43)]##lasso筛选后的验证集
leapsTest <- test[,c(1,10,15:17,20,29,39,43)]##lasso筛选后的测试集
#单因素分析数据集
logSmoteTrain <- Strain[,c(1:6,10:13,16:17,23:25,29,38:39,42,43,46,52,53,58,66,68:70)]
logSmoteTest <- Stest[,c(1:6,10:13,16:17,23:25,29,38:39,42,43,46,52,53,58,66,68:70)]
logTest <- test[,c(1:6,10:13,16:17,23:25,29,38:39,42,43,46,52,53,58,66,68:70)]
#xgboost
#最优子集数据集
leapsTrain <- train[,c(1,10,15:16,20,23,39,43,68)]##leaps筛选后的smote训练集
leapsTest <- test[,c(1,10,15:16,20,23,39,43,68)]##leaps筛选后的验证集
##xgboost建模
grid = expand.grid(
  nrounds = c(75, 100),
  colsample_bytree = 1,
  min_child_weight = 1,
  eta = c(0.01, 0.1, 0.3), #0.3 is default,
  gamma = c(0.5, 0.25),
  subsample = 0.5,
  max_depth = c(2, 3,4,5,6) 
)
#训练调优参数111111111111111111111111111111
cntrl = trainControl(
  method = "cv",
  number = 5,
  verboseIter = TRUE,
  returnData = FALSE,
  returnResamp = "final"
)
#设定参数
set.seed(1)
train.xgb.leaps = train(
  x = leapsTrain[, 2:9],
  y = ,leapsTrain[, 1],
  trControl = cntrl,
  tuneGrid = grid,
  method = "xgbTree"
)
train.xgb.leaps#显示最佳参数
paramLeaps <- list( objective = "binary:logistic",
                    booster = "gbtree",
                    eval_metric = "error",
                    eta = 0.3,
                    max_depth = 6,
                    subsample = 0.5,
                    colsample_bytree = 1,
                    gamma = 0.5
)

#数据处理
x <- as.matrix(leapsTrain[, 2:9])
y <- ifelse(leapsTrain$delirium == "1", 1, 0)
train.mat.leaps <- xgb.DMatrix(data = x, label = y)
#建模
set.seed(1)
xgb.fit.leaps <- xgb.train(params = paramLeaps, data = train.mat.leaps, nrounds = 100)
#变量重要性
impMatrix <- xgb.importance(feature_names = dimnames(x)[[2]],
                            model = xgb.fit.leaps)
impMatrix
xgb.plot.importance(impMatrix, main = "Gain by Feature")


#测试集表现
library(InformationValue)
pred <- predict(xgb.fit.leaps, x)
optimalCutoff(y, pred)
testMat <- as.matrix(leapsTest[, 2:9])
xgb.test <- predict(xgb.fit.leaps, testMat)
y.test <- ifelse(leapsTest$delirium == "1", 1, 0)
confusionMatrix(y.test, xgb.test, threshold = 0.429)
1 - misClassError(y.test, xgb.test, threshold = 0.429)
plotROC(y.test, xgb.test)##绘制ROC
#0.9284
roc1 <- roc(y.test, xgb.test)


set.seed(1)
folds <- createFolds(y=leapsTrain$delirium,k=10)###分成10份
fold_test <- leapsTrain[folds[[9]],]#取fold 1数据，建立测试集和验证集
fold_train <- leapsTrain[-folds[[9]],]#
#建立模型
x <- as.matrix(fold_train[, 2:9])
y <- ifelse(fold_train$delirium == "1", 1, 0)
train.mat.leaps <- xgb.DMatrix(data = x, label = y)
xgb.fit <- xgb.train(params = paramLeaps, data =train.mat.leaps, nrounds =100) #xgboost
pred <- predict(xgb.fit, x)
testMat <- as.matrix(fold_test[, 2:9])
xgb.test <- predict(xgb.fit, testMat)
y.test <- ifelse(fold_test$delirium == "1", 1, 0)
confusionMatrix(y.test, xgb.test, threshold = 0.429)
plotROC(y.test, xgb.test)##绘制ROC

# 0.9349 0.9475 0.9252 0.8651 0.9315 0.9484 0.9026 0.9004 0.9237 0.9401
auc <- (0.9349+0.9475+0.9252+0.8651+0.9315+0.9484+0.9026+0.9004+0.9237+0.9401)/10
auc
##0.92194


#SHAP可解释性
library(SHAPforxgboost)
object_model <- xgb.fit.leaps
traindatax <- leapsTrain[,-1]
xx1 <- as.matrix(traindatax)
shap.values <- shap.values(xgb_model = object_model,X_train =xx1)
shap_values_F <- shap.values$shap_score
shap_F <- shap.prep(xgb_model=object_model,X_train=xx1)
shap_F <- shap.prep(shap_contrib = shap_values_F,X_train = xx1)
shap.plot.summary(shap_F,scientific = TRUE)


##Randomforest建模
library(randomForest)
library(pROC)
n<-length(names(leapsTrain))     #计算数据集中自变量个数，等同n=ncol(train)
rate=1     #设置模型误判率向量初始值
for(i in 1:(n-1)){
  set.seed(1)
  rf_train<-randomForest(delirium~.,data=leapsTrain,mtry=i,ntree=1000)
  rate[i]<-mean(rf_train$err.rate)   #计算基于OOB数据的模型误判率均值
  print(rf_train)
}

rate     #展示所有模型误判率的均值
#[1] 0.1830141 0.1415931 0.1432756 0.1457645 0.1460565 0.1414199 0.1472690
plot(rate)

#寻找最优参数ntree，即指定随机森林所包含最佳决策树数目
set.seed(1)
rf_train<-randomForest(delirium~.,data=leapsTrain,mtry=6,ntree=1000)
plot(rf_train)    #绘制模型误差与决策树数量关系图

which.min(rf_train$err.rate[,1])
#[1] 572
#随机森林模型搭建
set.seed(1)
rf_train<-randomForest(delirium~.,
                       data=leapsTrain,
                       mtry=5,
                       ntree=572,
                       importance=TRUE,#是否输出因变量在模型中的重要性
                       proximity=TRUE)#是否计算模型的邻近矩阵

print(rf_train)    #展示随机森林模型简要信息

hist(treesize(rf_train))   #展示随机森林模型中每棵决策树的节点数
pred <- predict(rf_train,newdata = leapsTest)
pred_out_1<-predict(object=rf_train,newdata=leapsTest,type="prob")  #输出概率
table <- table(pred,leapsTest$delirium)
sum(diag(table))/sum(table)  #预测准确率
#0.9056277
roc2 <- roc(pred,y.test)

###交叉验证
#将变量分为10份
set.seed(1)
folds <- createFolds(y=leapsTrain$delirium,k=10)###分成10份
auc_value<-as.numeric()###建立空值
for(i in 1:10){
  fold_test <- leapsTrain[folds[[i]],] #取folds[[i]]作为测试集
  fold_train <- leapsTrain[-folds[[i]],] # 剩下的数据作为训练集
  fold_pre <- randomForest(delirium ~ .,data = fold_train,
                           mtry=6,ntree=572,importance=T)
  fold_predict <- predict(fold_pre,type='response',newdata=fold_test)
  auc_value<- append(auc_value,as.numeric(auc(as.numeric(fold_test[,1]),as.numeric(fold_predict))))
}
auc_value
#[1] 0.8545455 0.8584416 0.8389610 0.7805195 0.8220779 0.8688312 0.8090909 0.7974026 0.8207792 0.8675325
mean(auc_value)
#0.8318







##logistic回归
library(margins)
fit<-glm(delirium~.,data = leapsTrain,family = binomial) #参数“famliy”表示连接函数为二项分布
summary(fit)
(fit$null.deviance-fit$deviance)/fit$null.deviance  #根据公式计算，下面为结果
coef(fit)
confint(fit)
exp(coef(fit))  #显示几率比
effects<-margins(fit)
summary(effects)
plot(effects,main="AME and Confidence Intervals")
#模型预测
prob_train<-predict(fit,type="response")   #参数表示预测事件发生的条件概率,对数几率
pred_train<-prob_train>0.5
table<-table(predicted=pred_train,Actual=leapsTrain$delirium)
table
Accuracy<-(table[1,1]+table[2,2])/sum(table)  #准确率
Error_rate<-(table[2,1]+table[1,2])/sum(table)   #错分率
Accuracy
#0.6910714
Error_rate
#0.3089286
#测试集查看测试误差
prob_test<-predict(fit,type="response",newdata=leapsTest)
pred_test<-prob_test>0.5
table<-table(predicted=pred_test,Actual=leapsTest$delirium)
table #混淆矩阵
Accuracy<-(table[1,1]+table[2,2])/sum(table)
Error_rate<-(table[2,1]+table[1,2])/sum(table)
Accuracy
#[1] 0.6761566
Error_rate
#[1]  0.3238434
roc3 <- roc(pred_test,y.test)

#ROC曲线
library(ROCR)
pred_object<-prediction(prob_test,leapsTest$delirium)  #建立预测对象
perf<-performance(pred_object,measure = "tpr",x.measure = "fpr")
plot(perf,main="ROC curve(Test Set)",lwd=2,col="blue",
     xlab="1-Specificity",ylab="Sensitivity")
abline(0,1)
performance(pred_object, "auc")@y.values    #用预测对象计算
#[[1]]
#[1] 0.8108771

regplot(fit)
#将变量分为10份
set.seed(1)
folds <- createFolds(y=leapsTrain$delirium,k=10)###分成10份
###第二种方法
auc_value<-as.numeric()###建立空值
#建立循环
for(i in 1:10){
  fold_test <- leapsTrain[folds[[i]],] #取folds[[i]]作为测试集
  fold_train <- leapsTrain[-folds[[i]],] # 剩下的数据作为训练集
  fold_pre <- glm(delirium ~ .,
                  family = binomial(link = logit), data =fold_train )#logistic
  fold_predict <- predict(fold_pre,type='response',newdata=fold_test)
  auc_value<- append(auc_value,as.numeric(auc(as.numeric(fold_test[,1]),as.numeric(fold_predict))))
}
auc_value
#0.6359926 0.6842301 0.6730983 0.6081633 0.7439703 0.7725417 0.7499072 0.6894249 0.6768089 0.6727273
mean(auc_value)
#0.6906865





#knn建模

##library(caret)
# 设置10折交叉训练
control <- trainControl(method = 'cv',number = 10)
# knn模型训练
knn.model <- train(delirium~.,leapsTrain,
                   method = 'knn',
                   preProcess = c('center','scale'),
                   trControl = control,
                   tuneLength = 5)
knn.model
# 测试数据真实值
truth <- leapsTest$delirium
# 测试数据预测值
pred <- predict(knn.model,newdata = leapsTest)
# 计算混淆矩阵
caret::confusionMatrix(table(pred,truth))

roc4 <- roc(pred,y.test)
#k折交叉验证
library(caret)
library(pROC)
#将变量分为10份
library(kknn)
set.seed(1)
folds <- createFolds(y=leapsTrain$delirium,k=10)###将smote后的train数据集分成10份
###第二种方法
auc_value<-as.numeric()###建立空值
#建立循环
for(i in 1:10){
  fold_test <- leapsTrain[folds[[i]],] #取folds[[i]]作为测试集
  fold_train <- leapsTrain[-folds[[i]],] # 剩下的数据作为训练集
  fold_pre <- kknn(delirium~.,fold_train,fold_test,k=5)
  fold_predict <- predict(fold_pre,newdata=fold_test)
  auc_value<- append(auc_value,as.numeric(auc(as.numeric(fold_test[,1]),as.numeric(fold_predict))))
}
auc_value
##[1] 0.7857143 0.7558442 0.7506494 0.8038961 0.7298701 0.7779221 0.8038961 0.8298701 0.8194805 0.8181818
mean(auc_value)
##[1] 0.7875325

##原始数据集knn
##control <- trainControl(method = 'cv',number = 10)
# knn模型训练
##knn.model <- train(delirium~.,dataleaps,
method = 'knn',
preProcess = c('center','scale'),
trControl = control,
tuneLength = 5)
##knn.model
set.seed(1)
folds <- createFolds(y=dataleaps$delirium,k=10)###将原始数据集分成10份
###第二种方法
auc_value<-as.numeric()###建立空值
#建立循环
for(i in 1:10){
  fold_test <- dataleaps[folds[[i]],] #取folds[[i]]作为测试集
  fold_train <- dataleaps[-folds[[i]],] # 剩下的数据作为训练集
  fold_pre <- knn.model
  fold_predict <- predict(fold_pre,newdata=fold_test)
  auc_value<- append(auc_value,as.numeric(auc(as.numeric(fold_test[,1]),as.numeric(fold_predict))))
}
##auc_value
## [1] 0.5357143 0.4900990 0.4801980 0.4950495 0.4950495 0.4702970 0.4752475 0.5000000 0.5135314 0.4851485
#比较失败，可能是原始数据集的不均衡导致
##mean(auc_value)
##[1]  0.4940335

#roc绘制
plot(roc1,
     print.auc=TRUE, print.auc.x=0.4, print.auc.y=0.5,
     # 图像上输出AUC值,坐标为（x，y）
     auc.polygon=TRUE, auc.polygon.col="#fff7f7", # 设置ROC曲线下填充色
     max.auc.polygon=FALSE,  # 填充整个图像
     grid=c(0.5, 0.2), grid.col=c("black", "black"),  # 设置间距为0.1，0.2，线条颜色
     print.thres=TRUE, print.thres.cex=0.9, # 图像上输出最佳截断值，字体缩放倍数
     smooth=F, # 绘制不平滑曲线
     main="Comparison of two ROC curves", # 添加标题
     col="#FF2E63",  # 曲线颜色
     legacy.axes=TRUE)   # 使横轴从0到1，表示为1-特异度
plot.roc(roc2,
         add=T,  # 增加曲线
         col="#0074b3", # 曲线颜色为红色
         print.thres=TRUE, print.thres.cex=0.9,  # 图像上输出最佳截断值，字体缩放倍数
         print.auc=TRUE, print.auc.x=0.4,print.auc.y=0.4,
         # 图像上输出AUC值,坐标为（x，y）
         smooth = F)  # 绘制不平滑曲线
plot.roc(roc3,
         add=T,  # 增加曲线
         col="#f47720", # 曲线颜色为红色
         print.thres=TRUE, print.thres.cex=0.9,  # 图像上输出最佳截断值，字体缩放倍数
         print.auc=TRUE, print.auc.x=0.4,print.auc.y=0.35,
         # 图像上输出AUC值,坐标为（x，y）
         smooth = F)  # 绘制不平滑曲线
plot.roc(roc4,
         add=T,  # 增加曲线
         col="#459943", # 曲线颜色为红色
         print.thres=TRUE, print.thres.cex=0.9,  # 图像上输出最佳截断值，字体缩放倍数
         print.auc=TRUE, print.auc.x=0.4,print.auc.y=0.45,
         # 图像上输出AUC值,坐标为（x，y）
         smooth = F)  # 绘制不平滑曲线
library(shapviz)
library(patchwork)
x <- c("gcs_verbal","hos_days","spo2_mean","mdrd_est","dbp_mean","gcs_motor","gender2","NonInvasiveVent")
shp <- shapviz(xgb.fit.leaps, X_pred = xx1, X = traindatax)
sv_waterfall(shp, row_id = 1) +
  theme(axis.text = element_text(size = 11))
sv_importance(shp, kind = "beeswarm")
sv_dependence(shp, v = x) &
  theme_gray(base_size = 9)

