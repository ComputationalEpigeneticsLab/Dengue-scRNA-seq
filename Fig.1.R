########################################figures 1####################################
library(dplyr)
library(ggplot2)
library(Seurat)
######################################control
rm(list = ls())
control_exp<-read.table("D:\\xg\\data\\Dengue\\result\\data\\control_exp.txt",header=T,as.is=T,sep="\t",stringsAsFactors = F)
cell_info_res<-read.table("D:\\xg\\data\\Dengue\\result\\data\\cell_info_res.txt",header=T,as.is=T,sep="\t",stringsAsFactors = F)
control_object<-CreateSeuratObject(counts = control_exp, project = "PBMC", min.cells = 5)

#####################################Calculate the proportion of mitochondrial genes
control_object[["percent.mt"]] <- PercentageFeatureSet(control_object, pattern = "^MT-")
VlnPlot(control_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#########################################Remove cell and gene outliers
control_object <- subset(control_object, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 10)
##########################################update meta.data
clininfo<-cell_info_res[rownames(control_object@meta.data),]
control_object@meta.data$Patient_id<-clininfo$Patient_id
control_object@meta.data$Cell_type<-clininfo$Cell_type
control_object@meta.data$Sex<-clininfo$Sex
control_object@meta.data$Diseaseseverity<-clininfo$Diseaseseverity
#########################################Normalization
control_object <- NormalizeData(object = control_object, normalization.method = "LogNormalize", scale.factor = 10000)
#########################################Find Variable Features
control_object <- FindVariableFeatures(object = control_object, selection.method = "vst", nfeatures = 2000)
##########################################scaledata and runpca
all.genes <- rownames(control_object)
control_object <- ScaleData(control_object, features = all.genes)
control_object <- RunPCA(control_object, features = VariableFeatures(object = control_object))
# VizDimLoadings(control_object,dims = 1:2,reduction = "pca")
# DimPlot(control_object,reduction = "pca")
# control_object<-RunUMAP(control_object,reduction = "pca",dims = 1:20)
# DimPlot(control_object, reduction = "umap", group.by = "Cell_type",pt.size = 1,
#         cols=c("T_cell" = "#4E78B3", "monocyte" = "#8C6AB9", "NK_cell" = "#E48328", 
#                "NKT_cell" = "#5D9C3D", "B_cell" = "#B7282A", "cDC" = "#775948","pDC"="#CE6CBD"))
set.seed(101)
control_object<-RunTSNE(control_object,reduction = "pca",dims = 1:20)
TSNEPlot(control_object,reduction = "tsne", group.by = "Cell_type",pt.size=1,
         cols=c("T_cell" = "#4E78B3", "monocyte" = "#8C6AB9", "NK_cell" = "#E48328", 
                "NKT_cell" = "#5D9C3D", "B_cell" = "#B7282A", "cDC" = "#775948","pDC"="#CE6CBD"))
############################################CD8A,CD79A
# FeaturePlot(control_object, features = "CD8A",reduction = "tsne",pt.size = 1.5)
save.image("D:\\xg\\data\\Dengue\\Rdata\\Seurat_control.RData")
############################################Fever
rm(list = ls())
fever_exp<-read.table("D:\\xg\\data\\Dengue\\result\\data\\fever_exp.txt",header=T,as.is=T,sep="\t",stringsAsFactors = F)
cell_info_res<-read.table("D:\\xg\\data\\Dengue\\result\\data\\cell_info_res.txt",header=T,as.is=T,sep="\t",stringsAsFactors = F)

fever_object<-CreateSeuratObject(counts = fever_exp, project = "PBMC", min.cells = 5)
fever_object[["percent.mt"]] <- PercentageFeatureSet(fever_object, pattern = "^MT-")
#####################################Calculate the proportion of mitochondrial genes
VlnPlot(fever_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
fever_object <- subset(fever_object, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & percent.mt < 10)
##########################################update meta.data
clininfo<-cell_info_res[rownames(fever_object@meta.data),]
fever_object@meta.data$Patient_id<-clininfo$Patient_id
fever_object@meta.data$Cell_type<-clininfo$Cell_type
fever_object@meta.data$Sex<-clininfo$Sex
fever_object@meta.data$Diseaseseverity<-clininfo$Diseaseseverity
#########################################Normalization
fever_object <- NormalizeData(object = fever_object, normalization.method = "LogNormalize", scale.factor = 10000)
#########################################Find Variable Features
fever_object <- FindVariableFeatures(object = fever_object, selection.method = "vst", nfeatures = 2000)
##########################################scaledata and runpca
all.genes <- rownames(fever_object)
fever_object <- ScaleData(fever_object, features = all.genes)
fever_object <- RunPCA(fever_object, features = VariableFeatures(object = fever_object))
# fever_object<-RunUMAP(fever_object,reduction = "pca",dims = 1:20)
# DimPlot(fever_object, reduction = "umap", group.by = "Cell_type",pt.size = 1,
#         cols=c("T_cell" = "#4E78B3", "monocyte" = "#8C6AB9", "NK_cell" = "#E48328", 
#                "NKT_cell" = "#5D9C3D", "B_cell" = "#B7282A", "cDC" = "#775948","pDC"="#CE6CBD"))
set.seed(101)
fever_object<-RunTSNE(fever_object,dims = 1:20)
TSNEPlot(fever_object,reduction = "tsne", group.by = "Cell_type",pt.size=1,
         cols=c("T_cell" = "#4E78B3", "monocyte" = "#8C6AB9", "NK_cell" = "#E48328", 
                "NKT_cell" = "#5D9C3D", "B_cell" = "#B7282A", "cDC" = "#775948","pDC"="#CE6CBD"))
# FeaturePlot(fever_object, features = "CD8A",reduction = "tsne",pt.size = 1.5)
save.image("D:\\xg\\data\\Dengue\\Rdata\\Seurat_fever.RData")

#############################################severe
rm(list = ls())
severe_exp<-read.table("D:\\xg\\data\\Dengue\\result\\data\\severe_exp.txt",header=T,as.is=T,sep="\t",stringsAsFactors = F)
cell_info_res<-read.table("D:\\xg\\data\\Dengue\\result\\data\\cell_info_res.txt",header=T,as.is=T,sep="\t",stringsAsFactors = F)
severe_object<-CreateSeuratObject(counts = severe_exp, project = "PBMC", min.cells = 5)
severe_object[["percent.mt"]] <- PercentageFeatureSet(severe_object, pattern = "^MT-")
VlnPlot(severe_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
severe_object <- subset(severe_object, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)
severe_cell_filtered<-as.matrix(severe_object@assays[["RNA"]]@counts)
clininfo<-cell_info_res[rownames(severe_object@meta.data),]
##########################################update meta.data
severe_object@meta.data$Patient_id<-clininfo$Patient_id
severe_object@meta.data$Cell_type<-clininfo$Cell_type
severe_object@meta.data$Sex<-clininfo$Sex
severe_object@meta.data$Diseaseseverity<-clininfo$Diseaseseverity
#########################################Normalization
severe_object <- NormalizeData(object = severe_object, normalization.method = "LogNormalize", scale.factor = 10000)
#########################################Find Variable Features
severe_object <- FindVariableFeatures(object = severe_object, selection.method = "vst", nfeatures = 2000)
##########################################scaledata and runpca
all.genes <- rownames(severe_object)
severe_object <- ScaleData(severe_object, features = all.genes)
##
severe_object <- RunPCA(severe_object, features = VariableFeatures(object = severe_object))
# severe_object<-RunUMAP(severe_object,reduction = "pca",dims = 1:20)

set.seed(101)
severe_object<-RunTSNE(severe_object,dims = 1:20)
TSNEPlot(severe_object,reduction = "tsne", group.by = "Cell_type",pt.size=1,
         cols=c("T_cell" = "#4E78B3", "monocyte" = "#8C6AB9", "NK_cell" = "#E48328", 
                "NKT_cell" = "#5D9C3D", "B_cell" = "#B7282A", "cDC" = "#775948","pDC"="#CE6CBD"))
# FeaturePlot(severe_object, features = "CD8A",reduction = "tsne",pt.size = 1.5)
save.image("D:\\xg\\data\\Dengue\\Rdata\\Seurat_severe.RData")

#################################################fisher test ,control
rm(list = ls())
load("G:\\Dengue\\Rdata\\Seurat_control.RData")
load("G:\\Dengue\\Rdata\\Seurat_fever.RData")
load("G:\\Dengue\\Rdata\\Seurat_severe.RData")
#####
cell_info<-rbind(control_object@meta.data,fever_object@meta.data,severe_object@meta.data)
####
length(which(cell_info$Cell_type=="T_cell"))
####
control_info<-control_object@meta.data
fever_info<-fever_object@meta.data
severe_info<-severe_object@meta.data

#######control
cell_type<-unique(control_info$Cell_type)
control_fisher_res<-data.frame(Diseaseseverity=rep("control",length(cell_type)),Cell_type=NA,or=NA,p_value=NA)
for(i in 1:length(cell_type)){
  
  a<-length(which(control_info$Cell_type==cell_type[i]))
  b<-dim(control_info)[1]-a
  c<-length(which(fever_info$Cell_type==cell_type[i]))+length(which(severe_info$Cell_type==cell_type[i]))
  d<-dim(cell_info)[1]-dim(control_info)[1]-c
  f_or<-fisher.test(matrix(c(a,b,c,d),nrow=2))[["estimate"]][["odds ratio"]]
  control_fisher_res$Cell_type[i]<-cell_type[i]
  control_fisher_res$or[i]<-f_or
  control_fisher_res$p_value[i]<-fisher.test(matrix(c(a,b,c,d),nrow=2))$p.value
}
control_fisher_res<-control_fisher_res[order(control_fisher_res$or),]
map<-data.frame(or=control_fisher_res$or)
rownames(map)<-control_fisher_res$Cell_type
#map<-map[order(map$or),]
bk<-c(seq(0,2,by=0.01),seq(2.1,4,by=0.01))
pheatmap(map,cluster_rows = F,cluster_cols = F,cellwidth=18,breaks = bk,
         color = c(colorRampPalette(colors = c("#FCD6C8","#F58B74"))(length(bk)/2),
                   colorRampPalette(colors = c("#F58B74","#E5453B"))(length(bk)/2)))
