#!/usr/bin/env Rscript
### author: Shengdi Li, date: Nov. 8 2023 ###
### predict cell label based on tree model ###
library(optparse)
library(stringr)
set.seed(123)
option_list = list(
    make_option(c("-i", "--input"), type="character", default=NULL, 
                help="input table (ex. downsampling.txt)", metavar="character"),
    make_option(c("-f", "--fset"), type="character", default=NULL, 
                help="names of feature columns used for analysis", metavar="character"),
    make_option(c("-b", "--batch"), type="character", default=NULL, 
                help="names of batch(s)", metavar="character"),
    make_option(c("-t", "--tree"), type="character", default=NULL, 
                help="tree file", metavar="character"),
    make_option(c("-o", "--output"), type="character", default=NULL, 
		          help="output file", metavar="character")
    )

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

#library(ggplot2)
#library(gridExtra)

## define function
predict_tree<-function(tbl,tr){
  out<-"x"
  for(i in 1:nrow(tbl))
  {
    j<-3
    flag<-TRUE
    while((j<=ncol(tr))&(flag))
    {
      tmp.conditions<-tr[tr[,j]!="",c(1:2,j)]
      match<-TRUE
      for(n in 1:nrow(tmp.conditions))
      {
        if((tmp.conditions[n,3]==">")&(tbl[i,tmp.conditions[n,1]]<tmp.conditions[n,2]))
        {
          match<-FALSE
        }
        else if((tmp.conditions[n,3]=="<")&(tbl[i,tmp.conditions[n,1]]>=tmp.conditions[n,2]))
        {
          match<-FALSE
        }
      }
      if(match)
      {
        out<-c(out,names(tr)[j])
        flag<-FALSE
      }
      j<-j+1
    }
    if(flag)
    {
      print(paste0("Error: ",tbl[i,]))
    }
  }
  return(out[-1])
}

#opt<-data.frame(input="/g/steinmetz/shli/gitlab/ics_tagging2/test.in",tree = "/g/steinmetz/shli/gitlab/ics_tagging2/inputs/hc_k23.txt", fset="/g/steinmetz/shli/gitlab/ics_tagging2/inputs/feature_set_S8_GFP_relevant.txt",batch="S8_training_run1,S8_calibration_run1",output="/g/steinmetz/shli/gitlab/ics_tagging2/test.out")
### load input table
data.table<-read.table(opt$input,header=T,sep = "\t",stringsAsFactors = F)

### load feature set
fset<-read.table(opt$fset,header=T,sep="\t",stringsAsFactors = F)
dctree<-read.table(opt$tree,header=T,sep="\t",stringsAsFactors = F)
### extract features

data.table<-data.table[,c(fset[,1],"prot_id","batch")]
## set1 umap with all parameters
batches<-unlist(strsplit(opt$batch,","))
data.table.set1<-data.table[data.table$batch %in% batches,]
## remove a cell if there's NA value
for(i in 1:ncol(data.table.set1))
{
  data.table.set1<-data.table.set1[!is.na(data.table.set1[,i]),]
}
### remove protein ids with too low number
#data.table.set1 <- data.table.set1[data.table.set1$prot_id %in% names(table(data.table.set1$prot_id))[table(data.table.set1$prot_id)>1000],]
#data.table.set1 <- data.table.set1[sample(nrow(data.table.set1),100),]
for(mth in c("km","hc"))
{
  for(i in 2:40)
  {
    prediction<-predict_tree(data.table.set1,dctree[,c("param","threshold",names(dctree)[(grep(paste0(mth,"_k",i,"b"),names(dctree)))])])
    data.table.set1<-cbind(data.table.set1,prediction)
    names(data.table.set1)[ncol(data.table.set1)]<-paste0(mth,"_k",i)
    print(paste0(mth,"_k",i))
  }
}
for(mth in c("lc_light","lc_comp"))
{
  prediction<-predict_tree(data.table.set1,dctree[,c("param","threshold",names(dctree)[(grep(paste0(mth,"b"),names(dctree)))])])
  data.table.set1<-cbind(data.table.set1,prediction)
  names(data.table.set1)[ncol(data.table.set1)]<-mth
}
#data.table.set1$predict<-predict_tree(data.table.set1,dctree)
result<-data.frame(prot_id="x",batch="x")
result<-cbind(result,data.frame(matrix(rep(0,ncol(dctree)-2),nrow =1 )))
names(result)[3:ncol(result)]<-names(dctree)[3:ncol(dctree)]
for(prot_id in names(table(data.table.set1$prot_id)))
{
  for(batch in names(table(data.table.set1$batch)))
  {
    if(sum((data.table.set1$prot_id==prot_id)&(data.table.set1$batch==batch))>0)
    {
      tmp.result<-data.frame(prot_id=prot_id,batch=batch)
      tmp.result<-cbind(tmp.result,data.frame(matrix(rep(0,ncol(dctree)-2),nrow =1 )))
      names(tmp.result)[3:ncol(result)]<-names(dctree)[3:ncol(dctree)]
      for(n in 3:ncol(tmp.result))
      {
        regmatch<-regmatches(names(tmp.result)[n],regexpr(".+b",names(tmp.result)[n]))
        regmatch<-substr(regmatch,1,str_length(regmatch)-1)
        tmp.result[1,n]<-sum((data.table.set1$prot_id==prot_id)&(data.table.set1$batch==batch)&(data.table.set1[,regmatch]==names(tmp.result)[n]))
      }
      result<-rbind(result,tmp.result)
    }
    print(paste0(prot_id,"_",batch))
  }
}
result<-result[-1,]
write.table(result,opt$output,col.names = T,row.names = F,quote=F,sep="\t")
