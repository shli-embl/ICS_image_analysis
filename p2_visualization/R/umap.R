#!/usr/bin/env Rscript
### author: Shengdi Li, date: Nov. 8 2023 ###
### generate Umap ###
library(optparse)
set.seed(123)
option_list = list(
    make_option(c("-i", "--input"), type="character", default=NULL, 
                help="input table (ex. downsampling.txt)", metavar="character"),
    make_option(c("-f", "--fset"), type="character", default=NULL, 
                help="names of feature columns used for analysis", metavar="character"),
    make_option(c("-b", "--batch"), type="character", default=NULL, 
                help="names of batch(s)", metavar="character"),
    make_option(c("-o", "--output"), type="character", default=NULL, 
		          help="output file prefix", metavar="character")
    )

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

library(ggplot2)
library(umap)
library(gridExtra)

#opt<-data.frame(input="/scratch/shli/icstag/data/all_tag_libs.ds.tsv",fset="/g/steinmetz/shli/gitlab/ics_tagging2/inputs/feature_set_S8_GFP_relevant.txt",batch="S8_training_run1,S8_calibration_run1",output="/g/steinmetz/shli/gitlab/ics_tagging2/test.out.pdf")
### load input table
data.table<-read.table(opt$input,header=T,sep = "\t",stringsAsFactors = F)

### load feature set
fset<-read.table(opt$fset,header=T,sep="\t",stringsAsFactors = F)
### extract features

data.table<-data.table[,c(fset[,1],"prot_id","batch")]

## log-transform columns with positive flags
for(col_name in fset[fset$log_scale==1,1])
{
  data.table[,col_name]<-log2(data.table[,col_name]+1)
}
## scaling
for(col_name in fset[,1])
{
  data.table[,col_name]<-scale(data.table[,col_name])
}

## set1 umap with all parameters
batches<-unlist(strsplit(opt$batch,","))
data.table.set1<-data.table[data.table$batch %in% batches,]
## remove a cell if there's NA value
for(i in 1:ncol(data.table.set1))
{
  data.table.set1<-data.table.set1[!is.na(data.table.set1[,i]),]
}
data.table.set1<-data.table.set1[sample(1:nrow(data.table.set1),10000),]
data.table.set1.umap<-umap(data.table.set1[,-((ncol(data.table.set1)-1):ncol(data.table.set1))])

p.levels<-c("EGFR","DERL1","TUBA1C","CDH2","OCLN","DSG2","ZYX",
            "RAB11A","GOLGA5","IER3IP1","DDX6","KRT19","CLTA","ULK1","EZR","SLC16A1","EDC4","LEMD2","ACTN4","ATG3","POM121","CCR-NOT","CERK","ATP1B3","WRAP53","MAPRE3","EMD","SERF1A",
            "DDX21","NOP56","ACTB","HNRNPA1",
            "LMNA","CANX","PRDX1","ARF1","LMNB1","POGZ","HP1BP3","EPHX2","CTDNEP1","ERV3-1","RAD21","TOMM70","IDH3B","COIL","ERBIN","TOMM20","GOLGA2","OPTN")

plot_list<-list()
p.tmp<-ggplot() + geom_point(aes(x = data.table.set1.umap$layout[,1],y = data.table.set1.umap$layout[,2], color = factor(data.table.set1[,"prot_id"],levels=p.levels))) +
  labs(x="layout1",y="layout2",color="cell_tag") +
  theme_bw() +
  theme(panel.border = element_rect(colour = "black", fill=NA), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        legend.text=element_text(size=7),
        legend.title =element_blank(),
        legend.key.height= unit(0.3, 'cm'),
        legend.key.width = unit(0.3, 'cm'))
plot_list[[1]]<-ggplotGrob(p.tmp)
i<-2
for(cur_id in names(table(data.table.set1$prot_id)))
{
  p.tmp<-ggplot() + geom_point(aes(x = data.table.set1.umap$layout[,1],y = data.table.set1.umap$layout[,2])) +
    geom_point(aes(x = data.table.set1.umap$layout[which(data.table.set1[,"prot_id"] %in% cur_id),1],y = data.table.set1.umap$layout[which(data.table.set1[,"prot_id"] %in% cur_id),2],col = factor(data.table.set1[which(data.table.set1[,"prot_id"] %in% cur_id),"batch"],levels = batches))) +
    labs(x="layout1",y="layout2",title = cur_id) +
    scale_color_manual(values = c("red","#00D400")) +
    theme_bw() +
    theme(panel.border = element_rect(colour = "black", fill=NA), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
          legend.position="none")
  plot_list[[i]]<-ggplotGrob(p.tmp)
  i<-i+1
}

lay <- rbind(c(1,1,1,2:6),
             c(1,1,1,7:11),
             12:19,
             20:27,
             28:35,
             36:43,
             c(44:51))
pdf(opt$output,width = 30,height = 30,useDingbats = F)
grid.arrange(grobs = plot_list, layout_matrix = lay)
dev.off()
