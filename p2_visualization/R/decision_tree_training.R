#!/usr/bin/env Rscript
### author: Shengdi Li, date: Nov. 8 2023 ###
### train decision tree based on given labels ###
library(optparse)
set.seed(123)
option_list = list(
    make_option(c("-i", "--input"), type="character", default=NULL, 
                help="input table (ex. downsampling.txt)", metavar="character"),
    make_option(c("-f", "--fset"), type="character", default=NULL, 
                help="names of feature columns used for analysis", metavar="character"),
    make_option(c("-b", "--batch"), type="character", default=NULL, 
                help="names of batch(s)", metavar="character"),
    make_option(c("-g", "--group"), type="character", default=NULL, 
                help="matrix for transforming protein id into group labels", metavar="character"),
    make_option(c("-o", "--output"), type="character", default=NULL, 
		          help="output file prefix", metavar="character")
    )

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
#opt=data.frame(input="/scratch/shli/icstag/data/all_tag_libs.ds.tsv",group="/g/steinmetz/shli/gitlab/ics_tagging2/inputs/groupings.txt",batch="S8_training_run1",output="/g/steinmetz/shli/gitlab/ics_tagging2/test",fset="/g/steinmetz/shli/gitlab/ics_tagging2/inputs/feature_set_S8_GFP_relevant.txt")

library(ggplot2)
library(gridExtra)
library(tidyverse)
library(dplyr)
library(caret)
library(rpart.plot)

fit_tree <- function(df, label, maxdepth){
  ## select parameters
  df$label<-label
  df <- df[,!apply(df, 2, function(x) any(is.na(x)))]
  ## fix feature names (exclude special characters)
  colnames(df) <- make.names(colnames(df))
  ## split into training and test set (70-30)
  idx <- createDataPartition(df$label, p = 0.75)[[1]]
  training <- df[idx,]
  ## learn model, 5x CV
  control <- trainControl(method = 'repeatedcv', number = 5, repeats = 3, savePredictions = T)
  ### tree depth set here ###
  rpart.grid <- expand.grid(maxdepth=maxdepth) 
  model <- train(
    label ~ ., 
    data = training,
    method = 'rpart2',
    tuneGrid = rpart.grid,
    trControl = control
  )
  return(list(model = model, index = idx, df = df))
}
make_confusion_matrix <- function(tree_obj){
  ## get test data (all labeled cells not used for training)
  test_set <- tree_obj$df[-tree_obj$index,]
  ## predict classes for test data 
  pred <- predict(tree_obj$model, test_set)
  ## generate confusion matrix
  conf_mat <- confusionMatrix(data = pred,
                              reference = factor(test_set$label,
                                                 levels = levels(pred))) 
}
### load input table
data.table<-read.table(opt$input,header=T,sep = "\t",stringsAsFactors = F)

### load feature set
fset<-read.table(opt$fset,header=T,sep="\t",stringsAsFactors = F)
groupings<-read.table(opt$group,header=T,sep="\t",stringsAsFactors = F)
rownames(groupings)<-groupings[,1]
### extract features

data.table<-data.table[,c(fset[,1],"prot_id","batch")]
batches<-unlist(strsplit(opt$batch,","))
data.table.set1<-data.table[data.table$batch %in% batches,]
### filter out POGZ
data.table.set1<-data.table.set1[data.table.set1$prot_id %in% groupings[groupings$off_target!="Yes",1],]
## remove a cell if there's NA value
for(i in 1:ncol(data.table.set1))
{
  data.table.set1<-data.table.set1[!is.na(data.table.set1[,i]),]
}

trees<-list()
for(i in 3:ncol(groupings))
{
  trees[[names(groupings)[i]]]<-fit_tree(data.table.set1[,-((ncol(data.table.set1)-1):ncol(data.table.set1))],paste0("G",groupings[data.table.set1$prot_id,i]),10)$model$finalModel
  pdf(paste0(opt$output,".dc_tree.",names(groupings)[i],".pdf"),height = 8,width = 12)
  rpart.plot::rpart.plot(trees[[names(groupings)[i]]],main="",box.palette = 0)
  dev.off()
}

saveRDS(trees,file = paste0(opt$output,".models.RDS"))