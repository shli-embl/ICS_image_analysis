#!/usr/bin/env Rscript
### author: Shengdi Li, date: Nov. 8 2023 ###
### make plots and prune decision forest based on correlation matrix ###
library(optparse)
set.seed(123)
option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="input decision forest (output of the predict.R)", metavar="character"),
  make_option(c("-o", "--output"), type="character", default=NULL, 
              help="output prefix", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
#opt<-data.frame(input = "/g/steinmetz/shli/gitlab/ics_tagging2/plots/decision_forest.predict.txt", output = "/g/steinmetz/shli/gitlab/ics_tagging2/plots/decision_tree")
library(ggplot2)
library(stringr)
library(gridExtra)
library(pheatmap)
library(caret)
library(RColorBrewer)
#set.seed(123)
input<-read.table(opt$input,header=T,stringsAsFactors = F,sep="\t")
rownames(input)<-paste0(input$prot_id,"_",input$batch)
input2<-t(input[,-1:-2])
annotation_row<-data.frame(batch = input$batch)
rownames(annotation_row)<-rownames(input)
annotation_col<-data.frame(method = substr(rownames(input2),1,2))
rownames(annotation_col)<-colnames(input)[-1:-2]                       

pdf(paste0(opt$output,".forest.pdf"),height = 10,width =8)
pheatmap(input2, cluster_rows = T, cluster_cols = T,angle_col = 90,
         annotation_row = annotation_col,annotation_col = annotation_row,show_rownames = F)
dev.off()


#### pruning tree based on correlation (remove redundancy)

## use only the training pool
input<-input[input$batch=="S8_training_run1",]

rownames(input)<-input$prot_id
### remove off-target
off_target<-c("EMD","GOLGA5","IER3IP1","CERK","LMNA","POGZ")
input<-input[!input$prot_id %in% off_target,]
input2<-t(input[,-1:-2])
## remove bins with too little cells
rowsum.input<-apply(input2,1,sum)
input2<-input2[rowsum.input>1000,]

annotation_col<-data.frame(method = substr(rownames(input2),1,2))
rownames(annotation_col)<-rownames(input2)                   
#pheatmap(input2[rowsum.input>2000,], cluster_rows = T, cluster_cols = T,angle_col = 90,
#         annotation_row = annotation_col,show_rownames = F)

### transform into enrichment score
enrich.mat<-input2
#enrich.mat<-enrich.mat[rowsum.input>2000,]
#annotation_col2<-annotation_col[rowsum.input>2000,]
rowsum.input2<-apply(enrich.mat,1,sum)
for(i in 1:nrow(enrich.mat))
{
  enrich.mat[i,]<-log2((enrich.mat[i,]+1)/(rowsum.input2[i]/length(table(input$prot_id))))
}
#pheatmap(enrich.mat, cluster_rows = T, cluster_cols = T,angle_col = 90,show_rownames = F)

### get correlation matrix
remove_id<-findCorrelation(cor(t(enrich.mat)), cutoff=0.55, verbose=F,exact = T)
tbl.step1<-enrich.mat[-remove_id,]
remove_id2<-findCorrelation(cor(tbl.step1), cutoff=0.75, verbose=F,exact = T)
tbl.step2<-tbl.step1[,-remove_id2]
break_list<-seq(-4,4,by=0.1)
color_list<-colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(length(break_list))
pdf(paste0(opt$output,".pruned.pdf"))
pheatmap(tbl.step2, cluster_rows = T, cluster_cols = T,angle_col = 90,show_rownames = T,
         breaks = break_list,color = color_list)
dev.off()

###draw complete matrix on the gene axis
pdf(paste0(opt$output,".pruned_withfullgene.pdf"),width = 10,height =5)
pheatmap(tbl.step1, cluster_rows = T, cluster_cols = T,angle_col = 90,show_rownames = T,
         breaks = break_list,color = color_list)
dev.off()

### plot PCA via the pruned vs full forest space
### full
input<-read.table(opt$input,header=T,stringsAsFactors = F,sep="\t")
pca<-prcomp(input[,-1:-2])
pdf(paste0(opt$output,".pca.full.pdf"),width = 9,height = 6,useDingbats = F)
ggplot() + geom_text(aes(x = pca$x[,1],y =pca$x[,2], col = input$batch, label = input$prot_id)) +
  geom_point(aes(x = pca$x[,1],y =pca$x[,2], col = input$batch)) +
  labs(x = "PC1", y = "PC2",col = "") + scale_color_manual(values = c("#BE1E2D","black","#009444")) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.border = element_rect(colour = "black", fill=NA),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
dev.off()  

pca<-prcomp(input[,-1:-2][,-remove_id])
pdf(paste0(opt$output,".pca.ensemble.pdf"),width = 9,height = 6,useDingbats = F)
ggplot() + geom_text(aes(x = pca$x[,1],y =pca$x[,2], col = input$batch, label = input$prot_id)) +
  geom_point(aes(x = pca$x[,1],y =pca$x[,2], col = input$batch)) +
  labs(x = "PC1", y = "PC2",col = "") + scale_color_manual(values = c("#BE1E2D","black","#009444")) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.border = element_rect(colour = "black", fill=NA),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
dev.off()  

