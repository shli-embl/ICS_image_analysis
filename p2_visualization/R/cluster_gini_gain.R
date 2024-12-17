#!/usr/bin/env Rscript
### author: Shengdi Li, date: Nov. 8 2023 ###
### grouping tags based on distance matrix (PCA) ###
library(optparse)
set.seed(123)
option_list = list(
    make_option(c("-i", "--input"), type="character", default=NULL, 
                help="input table (ex. downsampling.txt)", metavar="character"),
    make_option(c("-f", "--fset"), type="character", default=NULL, 
                help="names of feature columns used for analysis", metavar="character"),
    make_option(c("-b", "--batch"), type="character", default=NULL, 
                help="names of batch(s)", metavar="character"),
    make_option(c("-l", "--local"), type="character", default=NULL, 
                help="localization of the proteins", metavar="character"),
    make_option(c("-o", "--output"), type="character", default=NULL, 
		          help="output file prefix", metavar="character")
    )

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

library(ggplot2)
library(gridExtra)

## define cloud distance calculation function
downsampling <- function(df,index,neach){
  do.call(rbind,apply(data.frame(names(table(df[,index]))),1,function(x){
    df.tmp<-df[df[,index]==x,]
    df.tmp[sample(nrow(df.tmp),neach),]
  }))
}
distmat.eu <-function(indata,index,neach){
  data.ds<-downsampling(indata,index,neach)
  dist.mat<-as.matrix(dist(data.ds[,-index]))
  comp.mat<-compare.mat(data.ds[,index])
  tags<-names(table(data.ds[,index]))
  distmat.result<-list()
  for(i in 1:(length(tags)-1))
  {
    for(j in (i+1):length(tags))
    {
      intradistance<-sum(dist.mat[data.ds[,index] %in% c(tags[i],tags[j]),data.ds[,index] %in% c(tags[i],tags[j])]*comp.mat[["intra"]][data.ds[,index] %in% c(tags[i],tags[j]),data.ds[,index] %in% c(tags[i],tags[j])])
      interdistance<-sum(dist.mat[data.ds[,index] %in% c(tags[i],tags[j]),data.ds[,index] %in% c(tags[i],tags[j])]*comp.mat[["inter"]][data.ds[,index] %in% c(tags[i],tags[j]),data.ds[,index] %in% c(tags[i],tags[j])])
      distmat.result[[paste0(tags[i],"_",tags[j])]]<-data.frame(pop1=tags[i],pop2=tags[j],intra=intradistance,inter=interdistance,ratio=interdistance/(interdistance+intradistance))
    }
  }
  do.call(rbind,distmat.result)
}
compare.mat <-function(tag){
  compare.matrix<-list()
  compare.matrix[["intra"]]<-matrix(0,nrow= length(tag),ncol=length(tag))
  compare.matrix[["inter"]]<-compare.matrix[["intra"]]
  for (i in 1:(length(tag)-1)){
    for (j in ((i+1):(length(tag)))){
      if(tag[i]==tag[j]){
        compare.matrix[["intra"]][i,j]<-1
      } else {
        compare.matrix[["inter"]][i,j]<-1
      }
    }
  }
  compare.matrix
}
get_dist_mat <- function(indata){
  tags<-names(table(c(indata[,1],indata[,2])))
  dist_mat<-matrix(NA,nrow = length(tags),ncol = length(tags))
  rownames(dist_mat)<-tags
  colnames(dist_mat)<-tags
  for(i in 1:dim(indata)[1]){
    dist_mat[indata[i,1],indata[i,2]]<-indata[i,3]
  }
  dist_mat
}
gini_index<-function(intable){
  ### always treat first column as the target class
  prop.tbl<-intable[,1]/sum(intable[,1])
  row.sum<-apply(intable,1,sum)
  prop.tbl.sum<-row.sum/sum(row.sum)
  gini1<-1-sum(prop.tbl^2)
  gini2<-1-sum(prop.tbl.sum^2)
  return(gini2-gini1)
}

perm_test<-function(intable,n=1000){
  output<-list()
  ## use fisher test
 # output[["pvalue"]]<-fisher.test(intable)$p.value
  ## use gini index
  output[["gini"]]<-gini_index(intable)
  sum_per_cluster<-apply(intable,1,sum)
  pool<-0
  for(i in 1:length(sum_per_cluster))
  {
    pool<-c(pool,rep(i,sum_per_cluster[i]))
  }
  pool<-pool[-1]
  
  perm<-0
  for(j in 1:n)
  {
    counts.new<-data.frame(class1 = rep(0,length(sum_per_cluster)),class2 = sum_per_cluster)
    distribution<-table(sample(pool,sum(intable[,1])))
    for(m in 1:length(distribution))
    {
      counts.new[names(distribution)[m],1]<-distribution[m]
    }
    counts.new[,2]<-counts.new[,2]-counts.new[,1]
#    perm<-c(perm,fisher.test(counts.new)$p.value)
    perm<-c(perm,gini_index(counts.new))
  }
  output[["permutation"]]<-perm[-1]
 # output[["permutation_p"]]<-rank(c(-output[["pvalue"]],-output[["permutation"]]))[1]/(n+1)
  output[["permutation_p"]]<-rank(c(output[["gini"]],output[["permutation"]]))[1]/(n+1)
  return(output)
}

#opt<-data.frame(input="/scratch/shli/icstag/data/all_tag_libs.ds.tsv",local = "/g/steinmetz/shli/gitlab/ics_tagging2/inputs/localization.txt", fset="/g/steinmetz/shli/gitlab/ics_tagging2/inputs/feature_set_S8_GFP_relevant.txt",batch="S8_training_run1",output="/g/steinmetz/shli/gitlab/ics_tagging2/test.out")
### load input table
data.table<-read.table(opt$input,header=T,sep = "\t",stringsAsFactors = F)

### load feature set
fset<-read.table(opt$fset,header=T,sep="\t",stringsAsFactors = F)
localization<-read.table(opt$local,header=T,sep="\t",stringsAsFactors = F)
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
### remove protein ids with too low number
data.table.set1 <- data.table.set1[data.table.set1$prot_id %in% names(table(data.table.set1$prot_id))[table(data.table.set1$prot_id)>100],]
data.pca.pr <- prcomp(data.table.set1[,-((ncol(data.table.set1)-1):ncol(data.table.set1))])
indata<-data.frame(cbind(data.pca.pr$x[,1:2],data.table.set1$prot_id))

dist_mat<-get_dist_mat(distmat.eu(indata,3,min(table(indata$V3)))[,c("pop1","pop2","ratio")])

dist_mat_sym <- dist_mat
dist_mat_sym[is.na(dist_mat_sym)] <- 0
dist_mat_sym <- dist_mat_sym + t(dist_mat_sym) + diag(0.5, nrow = nrow(dist_mat_sym), ncol = ncol(dist_mat_sym))

#set.seed(123)
hirclust <- hclust(as.dist(t(dist_mat_sym)))
#plot1 <- as.ggplot(pheatmap(dist_mat_sym[hirclust$order,hirclust$order], cluster_rows = FALSE, cluster_cols = FALSE))
cl_list<-list()
#for(k in c(5:length(table(indata$V3))))
for(k in c(2:40))
{
  cl_list[[paste0("k_",k)]]<-factor(cutree(hirclust,k = k))
}
anno<-data.frame(do.call(cbind,cl_list))
rownames(localization)<-localization$Protein_id
lc<-localization[rownames(anno),-1]
## replace NA with 0
for(i in 1:ncol(lc))
{
  lc[is.na(lc[,i]),i]<-0
}
anno<-cbind(anno,lc)
write.table(dist_mat_sym[hirclust$order,hirclust$order],paste0(opt$output,".distmat.txt"),row.names = T,col.names = T,sep="\t",quote=F)
write.table(anno[hirclust$order,],paste0(opt$output,".local.txt"),row.names = T,col.names = T,sep="\t",quote=F)

######
anno$Nuclear<-apply(anno[,c("Nuclear_membrane","Nucleoplasm","Nuclear_bodies","Nucleoli","Nucleoli_fibrillar_center","Nuclear_speckles")],1,function(x){
  if(sum(as.numeric(x))>0)
  {
    o<-1
  }
  else
  {
    o<-0
  }
  o
})
anno$Cyto<-apply(anno[,c("Cytoskeleton","Cytostol")],1,function(x){
  if(sum(as.numeric(x))>0)
  {
    o<-1
  }
  else
  {
    o<-0
  }
  o
})
anno$Membrane<-apply(anno[,c("Cell_junctions","Plasma_membrane")],1,function(x){
  if(sum(as.numeric(x))>0)
  {
    o<-1
  }
  else
  {
    o<-0
  }
  o
})
anno$Particles<-apply(anno[,c("Vesicles","Endosomes","Lysosomes","Lipid_droplets")],1,function(x){
  if(sum(as.numeric(x))>0)
  {
    o<-1
  }
  else
  {
    o<-0
  }
  o
})

anno$ER_Golgi<-apply(anno[,c("ER","Golgi")],1,function(x){
  if(sum(as.numeric(x))>0)
  {
    o<-1
  }
  else
  {
    o<-0
  }
  o
})


#write.table(paste0(opt$output,".cluster.txt"),col.names=T,row.names=T,quote=F,sep="\t")
#exclude off-target for L-index assessment
off_target<-c("EMD","GOLGA5","IER3IP1","CERK","LMNA","POGZ")
anno<-anno[!rownames(anno) %in% off_target,]
result<-data.frame(nclass=0,localization="x",L_index=0,perm = 0)
for(g in 2:40)
{
  for(l in c("Nuclear","Cyto","Membrane","Particles","Mitochondria","ER_Golgi"))
  {
    counts<-data.frame(class1=rep(0,g),class2=rep(0,g))
    for(i in 1:g)
    {
# FIX BUG 20241217, use row count instead of sum
      counts[i,1]<-sum(as.numeric(anno[anno[,paste0("k_",g)]==i,l]))
      
 #     counts[i,2]<-sum(anno[anno[,paste0("k_",g)]==i,(ncol(anno)-3):ncol(anno)])-counts[i,1]
      counts[i,2]<-nrow(anno[anno[,paste0("k_",g)]==i,])-counts[i,1]
    }
    test<-perm_test(counts)
#    tmp<-data.frame(nclass=g,localization=l,L_index=-log2(test[["pvalue"]]),perm=test[["permutation_p"]])
    tmp<-data.frame(nclass=g,localization=l,L_index=test[["gini"]],perm=test[["permutation_p"]])
    result<-rbind(result,tmp)
  }
}
result<-result[-1,]

p1<-ggplot(result[(result$nclass<50),])+geom_line(aes(x = nclass,y=L_index,col = localization)) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

p2<-ggplot(result[(result$nclass<50),])+geom_line(aes(x = nclass,y=perm,col = localization)) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

pdf(paste0(opt$output,".lindex.pdf"),width = 8,height = 8)
grid.arrange(p1,p2)
dev.off()
