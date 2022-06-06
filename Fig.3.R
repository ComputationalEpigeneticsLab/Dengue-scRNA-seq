###############################################only for severe
library(reshape2)
library(aplot)
library(ggplot2)
library(ggpubr)
library(Seurat)
rm(list = ls())
load("D:\\xg\\data\\Dengue\\Rdata\\Seurat_severe.RData")
exp_severe<-as.matrix(severe_object@assays[["RNA"]]@data)
marker_gene_severe<-read.table("D:\\xg\\data\\Dengue\\result\\diff_gene\\severe_markergene.txt",header = T,as.is = T,sep = "\t")
marker_gene_severe_50<-marker_gene_severe %>% group_by(Cell_type) %>% top_n(n=50,wt=avg_log2FC)
marker_gene_severe_50<-marker_gene_severe_50[order(marker_gene_severe_50$Cell_type),]
marker_gene_severe_50_dup<-marker_gene_severe_50
for(i in 1:350){
  marker_gene_severe_50_dup[i,6]<-paste0(marker_gene_severe_50_dup$Gene_symble[i],"-",i)
}
gene_name_severe<-marker_gene_severe_50$Gene_symble
marker_severe_cell_info<-severe_object@meta.data[order(severe_object@meta.data$Cell_type),]
marker_exp_severe<-exp_severe[gene_name_severe,rownames(marker_severe_cell_info)]
######
spearman_mat<-matrix(nrow = 350,ncol = 350)
colnames(spearman_mat)<-rownames(marker_exp_severe)
rownames(spearman_mat)<-rownames(marker_exp_severe)
for(i in 1:dim(marker_exp_severe)[1]){
  for(j in 1:dim(marker_exp_severe)[1]){
    x<-marker_exp_severe[i,]
    y<-marker_exp_severe[j,]
    spearman_mat[i,j]<-cor(marker_exp_severe[i,],marker_exp_severe[j,],method = "spearman")
  }
  cat(i,sep = "\n")
}
spearman_mat_res<-spearman_mat
ano_cols<-data.frame(cell_type=marker_gene_severe_50_dup$Cell_type)
rownames(ano_cols)<-marker_gene_severe_50_dup$Gene_symble
colnames(spearman_mat_res)<-rownames(ano_cols)
rownames(spearman_mat_res)<-rownames(ano_cols)

get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

upper_mat <- get_upper_tri(spearman_mat_res)
melted_cormat <- melt(upper_mat, na.rm = TRUE)


ano_cols<-data.frame(cell_type=marker_gene_severe_50_dup$Cell_type)
ano_cols$gene<-marker_gene_severe_50_dup$Gene_symble
ano_cols$gene<-factor(ano_cols$gene,levels = rev(ano_cols$gene))
ano_cols$x<-1
#install.packages("aplot")

p_cor.up <-ggplot(melted_cormat, aes(Var2, Var1)) + 
  geom_tile(aes(fill = value)) 
p1<-p_cor.up + scale_fill_gradient2(low = "#5573BF", high = "#E25259", mid = "#F4EFEF", 
                                    midpoint = 0, limit = c(-1,1), space = "Lab", 
                                    name="Pearson\nCorrelation")  +
  #theme(axis.text.x = element_text(vjust = 0.5, hjust = 0.5, angle = 90))+
  #coord_fixed(ratio=1)+
  #theme(axis.text= element_text(size = 8,family="ARL"))+
  #theme(plot.margin = unit(c(0.1,0,0,0), unit="mm"))+
  labs(x = NULL, y = NULL, title = NULL)+
  theme(plot.title = element_text(size = 13,hjust = 1,family = "ARL" ))+
  theme(legend.key.width=unit(3,'mm'),legend.key.height=unit(0.5,'cm'))+
  theme(legend.title = element_text(size = 8))+
  theme(axis.ticks.y = element_blank(),
        axis.text.y = element_blank())+
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank())+
  theme(panel.background = element_blank(),#去除背景
        panel.border = element_blank())+
  coord_flip()
#######
p2<-ggplot(ano_cols,aes(x=x,y=gene))+
  geom_tile(aes(fill=cell_type))+
  scale_x_continuous(expand = c(0,0))+
  theme(panel.background = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "right",
        legend.title = element_blank())+
  theme(legend.key.width=unit(2,'mm'),legend.key.height=unit(0.4,'cm'))+
  scale_fill_manual(values = c("#B7282A","#775948", "#8C6AB9","#E48328","#5D9C3D","#CE6CBD", "#4E78B3"))
##########
p3<-ggplot(ano_cols,aes(x=x,y=gene))+
  geom_tile(aes(fill=cell_type))+
  scale_x_continuous(expand = c(0,0))+
  theme(panel.background = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank())+
  guides(fill=F)+
  scale_fill_manual(values = c("#B7282A","#775948", "#8C6AB9","#E48328","#5D9C3D","#CE6CBD", "#4E78B3"))+
  coord_flip()

p1%>%
  insert_left(p2,width = 0.01) %>% 
  insert_top(p3,height  = 0.01)

#####
marker_gene_Severe<-read.table("D:\\xg\\data\\DEV\\Seurat\\DEG\\Severe_markergene.txt",header = T,as.is = T,sep = "\t")
marker_gene_Severe<-read.table("G:\\DEV\\Seurat\\DEG\\Severe_markergene.txt",header = T,as.is = T,sep = "\t")

marker_gene_50<-marker_gene_Severe %>% group_by(celltype) %>% top_n(n=50,wt=avg_log2FC)
#####
marker_gene<-marker_gene_50[which(marker_gene_50$celltype=="monocyte"),]
inter_row<-intersect(marker_gene$Gene_Symbol,rownames(spearman_mat))
inter_col<-intersect(marker_gene$Gene_Symbol,colnames(spearman_mat))
#####
module_gene<-intersect(inter_col,inter_row)
sample<-rownames(severe_object@meta.data[which(severe_object@meta.data$Cell_type=="monocyte"),])
##33
marker_exp<-as.matrix(severe_object@assays[["RNA"]]@data)[module_gene,sample]
####
mean_markerexp<-apply(marker_exp,2,mean)
diff_tf<-c('JUN','XBP1', 'CEBPD', 'KLF4', 'TBX21','IRF7')

for(j in 1:length(diff_tf)){
  mat<-t(rbind(mean_markerexp,as.matrix(severe_object@assays[["RNA"]]@data)[diff_tf[j],sample]))
  colnames(mat)<-c("marker_gene","TF")
  mat<-as.data.frame(mat)
  pdf(paste0("D:\\xg\\data\\DEV\\result\\pearson\\density\\severe\\pDC_",diff_tf[j],".pdf"),width = 5,height = 5)
  p<-ggplot(mat, aes(x=marker_gene, y=TF) ) +
    stat_density_2d(aes(fill = ..level..), geom = "polygon", colour=NA)+
    scale_fill_gradient(low="#FEDECF", high="#CE6CBD")+
    theme(panel.background = element_blank(),
          panel.border = element_blank(),
          axis.line=element_line(colour="black"))+
    stat_cor(data=as.data.frame(mat), method = "pearson")+ labs(title = diff_tf[j])
  
  print(p)
  dev.off()
  cat(j,"\n")
}
