#!/usr/bin/env Rscript
### author: Shengdi Li, date: Nov. 8 2023 ###
### merge and generate downsampling table ###
library(optparse)
library(stringr)
set.seed(123)
option_list = list(
    make_option(c("-i", "--input"), type="character", default=NULL, 
                help="path to input tsv file", metavar="character"),
    make_option(c("-s", "--featureset"), type="character", default=NULL, 
                help="list of features to be included in the output with header", metavar="character"),
    make_option(c("-o", "--output"), type="character", default=NULL, 
		         help="output file", metavar="character"),
	make_option(c("-p", "--prefix"), type="character", default=NULL, 
			  	help="image file prefix", metavar="character"),
	make_option(c("-f", "--folder"), type="character", default=NULL, 
				help="image folder path", metavar="character"),
    make_option(c("-n", "--number"), type="integer", default=NULL, 
                help="number of rows to keep per class (will keep all if smaller than n)", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
#opt<-data.frame(input="/scratch/shli/icstag_2/tables/ACTN4_S8_calibration_run1.tsv",featureset = "/g/steinmetz/shli/gitlab/ics_tagging2/inputs/feature_set_S8.txt", folder = "/g/steinmetz/project/vulcan_images/Melanie/20231130_Calibration_Tag_S8/43_ACTN4/43_ACTN4_images_20240104_1405_26",prefix="43_ACTN4",output="/g/steinmetz/shli/gitlab/ics_tagging2/test.txt",number = 2000)

rand_sampling<-function(input,n){
  return(input[sample(1:nrow(input),n,replace = F),])
}

fsetdata<-read.table(opt$featureset,header=T,sep="\t",stringsAsFactors = F)
tmp.in<-read.table(opt$input,header=T,sep="\t")
image_num<-sprintf("%08d",0:(nrow(tmp.in)-1))
image_folder<-sprintf("%08d", rep(seq(0, by = 10000, length.out = ceiling(nrow(tmp.in) / 10000)),each = 10000)[1:nrow(tmp.in)])

tmp.in$image_path<-paste0(opt$folder,"/",image_folder,"/",opt$prefix,"_",image_num,".tiff")
write.table(tmp.in,paste0(opt$output,".complete.tsv"),col.names=T,sep="\t",quote=F,row.names=F)

#for(i in 1:nrow(tmp.in))
tmp.in<-tmp.in[tmp.in[,"root.Cells.Singlets.SSC.Singlets.FSC.Singlets.RM_Ecc.mNG..Dapi."]=="True",]
tmp.in<-tmp.in[,names(tmp.in) %in% c(fsetdata[,1],"image_path")]
#tmp.in<-tmp.in[,names(tmp.in) %in% fsetdata[,1]]
print(paste0(opt$input,":",nrow(tmp.in),":",ncol(tmp.in)))
if(nrow(tmp.in) > opt$number)
{
  tmp.in<-rand_sampling(tmp.in,opt$number)
}
#  tmp.in$prot_id=strsplit(sample_name,"_")[[1]][1]
#  tmp.in$batch=str_c(strsplit(sample_name,"_")[[1]][-1],collapse = "_")
write.table(tmp.in,opt$output,col.names=T,sep="\t",quote=F,row.names=F)

