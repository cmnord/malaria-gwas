library(tidyverse)
library(data.table)
library(dplyr)
library(reshape2)
library(lattice)
library(ggplot2)
library(caret)
library(pROC) #<- not installed

algns1 = read.csv('duffy/outer/duffy_all_primates_outer1.csv')#,row.names=1)
algns2 = read.csv('duffy/outer/duffy_all_primates_outer2.csv')#,row.names=1)
algns3 = read.csv('duffy/outer/duffy_all_primates_outer3.csv')#,row.names=1)
algns4 = read.csv('duffy/outer/duffy_all_primates_outer4.csv')#,row.names=1)
algns = merge(algns1, algns2, by=c("X","Species"), all=TRUE)
#dim(algns)
#algns$X
algns = merge(algns, algns3, by=c("X", "Species"), all=TRUE)
#dim(algns)
#algns$X
algns = merge(algns, algns4, by=c("X", "Species"), all=TRUE)
#dim(algns)
#algns$X
row.names(algns)<-algns$X
algns$X<-NULL
row.names(algns1)<-algns1$X
row.names(algns2)<-algns2$X
row.names(algns3)<-algns3$X
row.names(algns4)<-algns4$X
algns1$X<-NULL
algns2$X<-NULL
algns3$X<-NULL
algns4$X<-NULL
#set.seed(1234)
#levels(algns$Species) = c('infected', 'resistant')

#index <- createDataPartition(algns$Species, p=0.8, list=FALSE)
#trainDF = algns[index,]
#testDF = algns[-index,]
#trainSet <- data.matrix(algns[index,])
#testSet <- data.matrix(algns[-index,])
for(header in colnames(algns4)){
  if (all(algns4[,header]==0) || all(algns4[,header]==1)){
    algns4[,header]<-NULL
  }
}
trainSet <- algns4
outcomeName = "Species"
predictors = colnames(trainSet)[!colnames(trainSet) %in% outcomeName]

tc <- trainControl(method='cv', number=3, verboseIter=T, returnResamp='all', summaryFunction = twoClassSummary, classProbs = TRUE)

#tc <- trainControl(method="cv", 
#                   number=5, 
#                   verboseIter = T, 
                   #savePredictions = T, 
                   #summaryFunction = twoClassSummary,
#                   classProbs=TRUE,
#                   returnResamp = "all")

tg <- expand.grid(alpha=c(0.5, 1.0), 
                  lambda = seq(0.001,1.0,by = 0.001))

#levels(trainSet[,outcomeName]) <- c('infected', 'resistant')
trainY = as.matrix(trainSet$Species)
# here we gooooooo
m1 <- train(trainSet[,predictors], 
                  as.factor(trainY), 
                  method='glmnet', 
                  metric='ROC',
                  trControl=tc,
                  tuneGrid = tg,
                  preProc = c("center", "scale")) 
                  #tuneGrid=tg)

#summary(m1)


#plot(m1)
res = m1$results
best = m1$bestTune
row = filter(res, res$lambda==best$lambda[1], res$alpha==best$alpha[1])
print(row)
#print(row$Rsquared)
#print(row$RsquaredSD)

#pred <- predict(object=m1, testSet[,predictors], type='raw')
#pred
#testDF$Species
#print(cor(pdf_all$Pred, pdf_all$Obs, method='spearman'))
imps = varImp(m1)
imp = imps$importance
imp$feats=row.names(imp)
diiii = filter(imp,imp$Overall>0)
