library(tidyverse)
library(data.table)
library(dplyr)
library(reshape2)
library(lattice)
library(ggplot2)
library(caret)
#library(pROC) <- not installed

algns = read.csv('receptor_genes/tfr1/primates/tfr1_test3_alignment.csv', row.names = 1)

set.seed(1234)

infected = c(1, 1, 1, 0, 0, 0, 0, 0, 0, 0)
algns$infected = infected

trainSet = data.matrix(algns)

outcomeName = "infected"
predictors = colnames(trainSet)[!colnames(trainSet) %in% outcomeName]

tc <- trainControl(method="cv", 
                   number=5, 
                   verboseIter = T, 
                   savePredictions = T, 
                   returnResamp = "all")

tg <- expand.grid(alpha=c(0.0, 0.5, 1.0), 
                  lambda = seq(0.001,1.0,by = 0.001))

# here we gooooooo
m1 <- train(trainSet[,predictors], 
                  as.factor(trainSet[,outcomeName]), 
                  method='glmnet', 
                  trControl=tc, 
                  tuneGrid=tg)

plot(m1)
res = m1$results
best = m1$bestTune
row = filter(res, res$lambda==best$lambda[1], res$alpha==best$alpha[1])
print(row$Rsquared)
print(row$RsquaredSD)
print(cor(pdf_all$Pred, pdf_all$Obs, method='spearman'))
imps = varImp(m1)
