#!/usr/bin/env Rscript
### author: Shengdi Li, date: Nov. 8 2023 ###
### merge and generate downsampling table ###
library(optparse)
library(stringr)
set.seed(123)
option_list = list(
  make_option(c("-m", "--meta"), type="character", default=NULL, 
              help="path to meta sample file", metavar="character"),
  make_option(c("-f", "--folder"), type="character", default=NULL, 
              help="folder path to the table files", metavar="character"),
  make_option(c("-o", "--output"), type="character", default=NULL, 
              help="output file", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
#opt<-data.frame(meta="/g/steinmetz/shli/gitlab/ics_tagging2/inputs/samples.txt",folder = "/scratch/shli/icstag/tables/",output="/g/steinmetz/shli/gitlab/ics_tagging2/test.txt")

metatable<-read.table(opt$meta,header=T,sep="\t",stringsAsFactors = F)

result<-data.frame()
for(i in 1:nrow(metatable))
{
  if(nrow(result)>0)
  {
    tmp.in<-read.table(paste0(opt$folder,"/",metatable$id[i],".ds.tsv"),header=T,sep="\t",stringsAsFactors = F)
    tmp.in$prot_id<-metatable$protein_id[i]
    tmp.in$batch<-metatable$batch[i]
    result<-rbind(result,tmp.in)
  }
  else
  {
    result<-read.table(paste0(opt$folder,"/",metatable$id[i],".ds.tsv"),header=T,sep="\t",stringsAsFactors = F)
    result$prot_id<-metatable$protein_id[i]
    result$batch<-metatable$batch[i]
  }
}
write.table(result,opt$output,col.names=T,row.names=F,sep="\t",quote=F)