#######################################Run SCENIC in R
library("AUCell")
library("RcisTarget")
library("GENIE3")
library("Seurat")
library("SCENIC")
######################################control###############################
setwd("/boot3/xugang/SCENIC/control")
load("/boot3/xugang/Rdata/Seurat_control.RData")
control_exp<-control_object@assays$RNA@counts
exprMat<-as.matrix(control_exp)
org="hgnc"
dbDir="cisTarget_databases"
myDatasetTitle="SCENIC example on human dengue"
data(defaultDbNames)
dbs <- defaultDbNames[[org]]
scenicOptions <- initializeScenic(org=org, dbDir=dbDir, dbs=dbs,datasetTitle=myDatasetTitle, nCores=10)
saveRDS(scenicOptions, file="int/scenicOptions.Rds")
genesKept <- geneFiltering(exprMat, scenicOptions)
exprMat_filtered <- exprMat[genesKept, ]
runCorrelation(exprMat_filtered, scenicOptions)
exprMat_filtered_log <- log2(exprMat_filtered+1)
runGenie3(exprMat_filtered_log, scenicOptions,nParts=100)
exprMat_log <- log2(exprMat+1)
scenicOptions@settings$verbose <- TRUE
scenicOptions@settings$nCores <- 10
scenicOptions@settings$seed <- 123
scenicOptions@settings$dbs <- scenicOptions@settings$dbs["10kb"] # Toy run settings
runSCENIC_1_coexNetwork2modules(scenicOptions)
runSCENIC_2_createRegulons(scenicOptions, coexMethod=c("top5perTarget")) # Toy run settings
runSCENIC_3_scoreCells(scenicOptions, exprMat_log)
runSCENIC_4_aucell_binarize(scenicOptions,exprMat_log)
save.image("SCENIC_control.RData")
######################################fever###############################
setwd("/boot3/xugang/SCENIC/fever")
load("/boot3/xugang/Rdata/Seurat_fever.RData")
fever_exp<-fever_object@assays$RNA@counts
exprMat<-as.matrix(fever_exp)
org="hgnc"
dbDir="cisTarget_databases"
myDatasetTitle="SCENIC example on human dengue"
data(defaultDbNames)
dbs <- defaultDbNames[[org]]
scenicOptions <- initializeScenic(org=org, dbDir=dbDir, dbs=dbs,datasetTitle=myDatasetTitle, nCores=10)
saveRDS(scenicOptions, file="int/scenicOptions.Rds")
genesKept <- geneFiltering(exprMat, scenicOptions)
exprMat_filtered <- exprMat[genesKept, ]
runCorrelation(exprMat_filtered, scenicOptions)
exprMat_filtered_log <- log2(exprMat_filtered+1)
runGenie3(exprMat_filtered_log, scenicOptions,nParts=100)
exprMat_log <- log2(exprMat+1)
scenicOptions@settings$verbose <- TRUE
scenicOptions@settings$nCores <- 10
scenicOptions@settings$seed <- 123
scenicOptions@settings$dbs <- scenicOptions@settings$dbs["10kb"] # Toy run settings
runSCENIC_1_coexNetwork2modules(scenicOptions)
runSCENIC_2_createRegulons(scenicOptions, coexMethod=c("top5perTarget")) # Toy run settings
runSCENIC_3_scoreCells(scenicOptions, exprMat_log)
runSCENIC_4_aucell_binarize(scenicOptions,exprMat_log)
save.image("SCENIC_fever.RData")
######################################severe###############################
rm(list=ls())
library("AUCell")
library("RcisTarget")
library("GENIE3")
library("Seurat")
library("SCENIC")
setwd("/boot3/xugang/SCENIC/severe")
load("/boot3/xugang/Rdata/Seurat_severe.RData")
severe_exp<-severe_object@assays$RNA@counts
exprMat<-as.matrix(severe_exp)
org="hgnc"
dbDir="cisTarget_databases"
myDatasetTitle="SCENIC example on human dengue"
data(defaultDbNames)
dbs <- defaultDbNames[[org]]
scenicOptions <- initializeScenic(org=org, dbDir=dbDir, dbs=dbs,datasetTitle=myDatasetTitle, nCores=10)
saveRDS(scenicOptions, file="int/scenicOptions.Rds")
genesKept <- geneFiltering(exprMat, scenicOptions)
exprMat_filtered <- exprMat[genesKept, ]
runCorrelation(exprMat_filtered, scenicOptions)
exprMat_filtered_log <- log2(exprMat_filtered+1)
runGenie3(exprMat_filtered_log, scenicOptions,nParts=100)
exprMat_log <- log2(exprMat+1)
scenicOptions@settings$verbose <- TRUE
scenicOptions@settings$nCores <- 10
scenicOptions@settings$seed <- 123
scenicOptions@settings$dbs <- scenicOptions@settings$dbs["10kb"] # Toy run settings
runSCENIC_1_coexNetwork2modules(scenicOptions)
runSCENIC_2_createRegulons(scenicOptions, coexMethod=c("top5perTarget")) # Toy run settings
runSCENIC_3_scoreCells(scenicOptions, exprMat_log)
runSCENIC_4_aucell_binarize(scenicOptions,exprMat_log)
save.image("SCENIC_severe.RData")

#################only for control   heatmap
rm(list = ls())
load("D:\\xg\\data\\Dengue\\Rdata\\Seurat_control.RData")
TF_all<-readRDS("D:\\xg\\data\\Dengue\\SCENIC\\control\\int\\4.1_binaryRegulonActivity.Rds")
TF_all<-TF_all[-grep("_",rownames(TF_all)),]

cli_info<-control_object@meta.data[,c("Cell_type","Diseaseseverity")]
annot_df <- cli_info[intersect(colnames(TF_all),rownames(cli_info)),]
annot_df2<-annot_df[order(annot_df$Cell_type),] 
gsmorder<-rownames(annot_df2) 
gsmorder2<-match(gsmorder,colnames(TF_all)) 
gsmorder2<-na.omit(gsmorder2) 
######################reorder cells
TF_all1<-TF_all[,gsmorder2]

library(RColorBrewer)
colsmy2 = list(Cell_type= c("B_cell"="#B7282A", "cDC"="#775948", "monocyte"="#8C6AB9", "NK_cell"="#E48328",
                            "NKT_cell"="#5D9C3D", "T_cell"="#4E78B3", "pDC"="#CE6CBD"), 
               Diseaseseverity= c('control'=brewer.pal(9, 'Reds')[2])) 

pdf("D:\\xg\\data\\Dengue\\pic\\SCENIC_control(no_extended).pdf")
plot.new()
pheatmap::pheatmap(TF_all1, border_color=NA,fontsize_row=10,
                   annotation_col=annot_df,cluster_cols=F,show_colnames=F,annotation_colors=colsmy2, 
                   color=colorRampPalette(c("white","black"))(100), breaks=seq(0, 1, length.out = 100), )
dev.off()

#########################################target of shared TF 
control_tf_target<-read.table("D:\\xg\\data\\Dengue\\result\\control_tf_target.txt",sep = "\t",header = F,as.is = T)
colnames(control_tf_target)<-c("Diseaseseverity","TF","target")
fever_tf_target<-read.table("D:\\xg\\data\\Dengue\\result\\fever_tf_target.txt",sep = "\t",header = F,as.is = T)
colnames(fever_tf_target)<-c("Diseaseseverity","TF","target")
severe_tf_target<-read.table("D:\\xg\\data\\Dengue\\result\\severe_tf_target.txt",sep = "\t",header = F,as.is = T)
colnames(severe_tf_target)<-c("Diseaseseverity","TF","target")
tf_target<-rbind(control_tf_target,rbind(fever_tf_target,severe_tf_target))
colnames(tf_target)<-c("Diseaseseverity","TF","target")
################################################
inter_target_mat<-matrix(ncol = 3)
inter_tf<-c("EOMES","IRF8","PAX5","SPI1","TBX21")
for(i in 1:5){
  control_interTF<-control_tf_target[which(control_tf_target$TF==inter_tf[i]),3]
  fever_interTF<-fever_tf_target[which(fever_tf_target$TF==inter_tf[i]),3]
  severe_interTF<-severe_tf_target[which(severe_tf_target$TF==inter_tf[i]),3]
  inter_target<-intersect(control_interTF,intersect(fever_interTF,severe_interTF))
  inter_target_res<-cbind(rep("inter_TF",length(inter_target)),cbind(rep(inter_tf[i],length(inter_target)),inter_target))
  inter_target_mat<-rbind(inter_target_mat,inter_target_res)
}

inter_target_mat<- inter_target_mat[-1,]
colnames( inter_target_mat)<-c("Diseaseseverity","TF","target")
tf_target_res<-rbind(tf_target,inter_target_mat)
write.table(tf_target_res,"D:\\xg\\data\\Dengue\\result\\interTF_target.txt",sep = "\t",quote = F,col.names = T,row.names = F)
################################
tf_target<-read.csv("D:\\xg\\data\\Dengue\\result\\tf_target.csv")
library(ggplot2)
library(dplyr)
col_cell<-c("#FCD6C8","#F58B74","#E5453B","#BE1D2C")
names(col_cell)<-c("control","fever","severe","inter")
pdf("D:\\xg\\data\\Dengue\\result\\tf_target.pdf",height = 5,width = 5)
tf_target %>% 
  ggplot(aes(x = TF, fill = Diseaseseverity)) + 
  geom_bar(position = position_fill(),width = 0.8) + 
  scale_fill_manual(values = col_cell)+ 
  theme_classic() + 
  labs(y = 'Percent') + 
  theme(legend.position="right")+scale_y_continuous(expand = c(0,0))+
  theme(axis.text.x =element_text(size=12), axis.text.y=element_text(size=12))
dev.off()
