############ Script for Choosing loess parameters ###########
library(caret)
library(data.table)
library(ggplot2)
library(tidyverse)
library(viridisLite)
library(stats)
######### 10 fold validation #######################

#data <- urine_msms_irt[[1]]
#setnames(data, "Retention time", "rt")
#setnames(data, "NormalizedRetentionTime", "irt")

# set span values to test
#set.seed(123)

span_values <- c(0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, seq(0.2, 0.75, by = 0.05))

#cross_val_span <- data.frame(span = span_values, R2 = numeric(length(span_values)),
#                             RMSE = numeric(length(span_values)),
#                             MAE = numeric(length(span_values)))

################### Create loess model in caret###############

CVloess <- list(type = "Regression", library = "stats", loop = NULL)

parameters <- data.frame(parameter = "span",
			 class = "numeric",
			 label = "span")

CVloess$parameters <- parameters

grid <- function(x, y, len=NULL, search = "grid"){
	library(stats)
	span_values <- expand.grid(span=seq(0.01, 0.75, length.out=len))
	span_values
}

CVloess$grid <- grid

Fit <- function(x,y, wts, param, lev, weights, classProbs, ...){
	loess(y ~ x, span = param$span, degree = 2, ...)}
CVloess$fit <- Fit

Pred <- function(modelFit, newdata, preProc = NULL, submodels = NULL){
	predict(modelFit, newdata)}

CVloess$predict <- Pred

Prob <- function(modelFit, newdata, preProc=NULL, submodels=NULL){
	preduct(modelFit, newdata, type = "probabilities")}
CVloess$prob <- Prob

Sort <- function(x){x[order(x$span), ]}
CVloess$sort <- Sort

CVloess$levels <- function(x){lev(x)}
######################### test model ################

fitControl <- trainControl(method = "repeatedcv",
                           ## 10-fold CV...
                           number = 10,
                           ## repeated ten times
                           repeats = 5)

fitGrid <- data.frame(span = span_values)

#model <- train(irt ~ rt, data = data, method = CVloess, tuneGrid = fitGrid, tuneLength = NULL,
#               trControl = fitControl)
#print(model)

#modFit <- model$results
#ggplot(modFit, aes(x = span, y = RMSE)) + geom_point() + geom_line()
#trellis.par.set(caretTheme())
#densityplot(model, pch = "|")

TrainLoess <- function(irt_data){
  model <- train(irt ~ rt, irt_data, method = CVloess, 
                 tuneGrid = fitGrid, tuneLength = NULL,
                 trControl = fitControl)
  return(model)
}


