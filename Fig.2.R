
############################### load object ####################################
library(Seurat)
library(RColorBrewer)
library(patchwork)
library(pheatmap)
load('Seurat_control.RData')# control_object
load('Seurat_fever.RData')# fever_object
load('Seurat_severe.RData')# severe_object
scedata<-c('control','fever','severe')
scedatalist=list(control_object,fever_object,severe_object)

############################### find Marker ####################################

myFindmarkergene<-function(object0,celltype0){
    library(Seurat)
    markers <- FindMarkers(object0, ident.1 = celltype0,group.by = 'Cell_type',
                       logfc.threshold = 0.25,assay='RNA',test.use = "wilcox",
                       min.cells.group = 3,
                       min.pct=0.25)#Alter
    return(markers)
}
myFindallmarkergene<-function(object0){
    celltype=unique(object0@meta.data$Cell_type)
    mymarker_all<-myFindmarkergene(object0,celltype[1])
    mymarker_all$celltype<-celltype[1]
    mymarker_all$Gene_Symbol<-rownames(mymarker_all)
    rownames(mymarker_all)<-1:dim(mymarker_all)[1]
    for(i in 2:length(celltype)){
        mymarker<-myFindmarkergene(object0,celltype[i])
        mymarker$celltype<-celltype[i]
        mymarker$Gene_Symbol<-rownames(mymarker)
        rownames(mymarker)<-1:dim(mymarker)[1]
        mymarker_all<-rbind(mymarker_all,mymarker)
        }
    mymarker_all<-mymarker_all[which(mymarker_all$p_val_adj < 0.05),]
    return(mymarker_all)
}
outmarkerfile<-function(object0,filpn){
    sce_markers<-myFindallmarkergene(object0)
    write.table(sce_markers,filpn,sep='\t',quote=F,row.names=F)
    return(sce_markers)
}
scedata<-c('control','fever','severe')
filpn=paste0('findmarker/',scedata[1],'_markergene.txt')
sce_markers_control<-outmarkerfile(control_object,filpn)
filpn=paste0('findmarker/',scedata[2],'_markergene.txt')
sce_markers_fever<-outmarkerfile(fever_object,filpn)
filpn=paste0('findmarker/',scedata[3],'_markergene.txt')
sce_markers_severe<-outmarkerfile(severe_object,filpn)


#control
library(dplyr)
top_control <- sce_markers_control %>% group_by(celltype) %>% top_n(n = 50, wt = avg_log2FC)
top_control<-top_control[order(top_control$celltype),]

#fever
library(dplyr)
top_fever <- sce_markers_fever %>% group_by(celltype) %>% top_n(n = 50, wt = avg_log2FC)
top_fever<-top_fever[order(top_fever$celltype),]

#severe
library(dplyr)
top_severe <- sce_markers_severe %>% group_by(celltype) %>% top_n(n = 50, wt = avg_log2FC)
top_severe<-top_severe[order(top_severe$celltype),]

#integration
top_all=list(top_control,top_fever,top_severe)

myneedgene<-read.csv('CellMarker_all_onlypblood.txt',head=T,sep='\t')
myneedgene<-myneedgene[order(myneedgene$CellType),]
colsmy = list(Diseaseseverity= c('control'=brewer.pal(9, 'Reds')[2],
                     'dengue_fever'=brewer.pal(9, 'Reds')[4],
                     'severe_dengue'=brewer.pal(9, 'Reds')[6]),
             Cell_type= c("B_cell"="#B7282A", "cDC"="#775948",
                           "monocyte"="#8C6AB9", "NK_cell"="#E48328", 
                           "NKT_cell"="#5D9C3D", "T_cell"="#4E78B3", 
                           "pDC"="#CE6CBD"))

for(kk in c(1,3)){
#select
scedatalist0<-scedatalist[[kk]]
scedata0<-scedata[kk]
top_all0<-top_all[[kk]]$Gene_Symbol

#annotation
library(circlize)
cellInfo<-scedatalist0@meta.data
cellInfo<-cellInfo[order(cellInfo$Cell_type),]
cellInfo2_n<-as.data.frame(cellInfo %>% group_by(Cell_type) %>% summarise(n = n()))
annot_df<-select(cellInfo,Cell_type,Diseaseseverity)
top_anno <- HeatmapAnnotation(df=annot_df,col=colsmy)

#expression
scedatalist1<-scedatalist0[['RNA']]@scale.data
scedatalist2<-scedatalist1[top_all0,]
scedatalist2<-scedatalist2[,rownames(cellInfo)]
labels_row = rownames(scedatalist2)

#gap position in heatmap
library(ComplexHeatmap)
cellInfo2_n3<-c()
for(i in 1:nrow(cellInfo2_n)){
    cellInfo2_n3<-c(cellInfo2_n3,rep(cellInfo2_n[i,1],cellInfo2_n[i,2]))
}

#Screening the known mark gene
index_at<-c()
for(i in 1:7){
labels_row2<-rep('',nrow(scedatalist2))
labels_row2[(50*(i-1)+1):(50*i)]<-labels_row[(50*(i-1)+1):(50*i)]
mymark<-strsplit(myneedgene[i,2],', ')[[1]]
index_at<-c(index_at,which(labels_row2 %in% mymark))
}


#preparations for heatmap
row_anno = rowAnnotation(foo = anno_mark(at = index_at,
                                         labels = labels_row[index_at]))
col_fun2 = colorRamp2(c(0,1,2),c('#007C9D01','#00979D01','#FFD700'))
mylist<-data.frame(Cell_type=c(rep('B_cell',50),rep('cDC',50),rep('monocyte',50),
                               rep('NK_cell',50),rep('NKT_cell',50),rep('pDC',50),rep('T_cell',50)))
mylist<-as.matrix(mylist)

colors <- structure(c('#B7282A','#775948','#8C6AB9','#E48328','#5D9C3D','#CE6CBD','#4E78B3'),
                    names=c("B_cell", "cDC",
                           "monocyte", "NK_cell", 
                           "NKT_cell", "pDC", "T_cell"
                           ))

a2<-Heatmap(mylist,col=colors,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
#             cell_fun = ,
            width = 20,#40,20,40
        show_column_names = FALSE,
        show_heatmap_legend = FALSE)

#main heatmap
a<-Heatmap(scedatalist2,
col = col_fun2,
cluster_rows = FALSE,
cluster_columns = FALSE,
show_column_names = FALSE,
show_row_names = FALSE,
column_split = cellInfo2_n3,
top_annotation = top_anno,
column_title = NULL,
row_names_side = "left",
row_names_gp = gpar(fontsize = 7),
right_annotation = row_anno,
heatmap_legend_param = list(title = "")
)
pdf(paste0('findmarker/',scedata0,'.pdf'),width = 10,height =10)
a2+a
dev.off()

}

i=2
scedatalist0<-scedatalist[[i]]
scedata0<-scedata[i]
top_all0<-top_all[[i]]$Gene_Symbol

library(circlize)
cellInfo<-scedatalist0@meta.data
cellInfo<-cellInfo[order(cellInfo$Cell_type),]
cellInfo2_n<-as.data.frame(cellInfo %>% group_by(Cell_type) %>% summarise(n = n()))
annot_df<-select(cellInfo,Cell_type,Diseaseseverity)
top_anno <- HeatmapAnnotation(df=annot_df,col=colsmy)

scedatalist1<-scedatalist0[['RNA']]@scale.data
scedatalist2<-scedatalist1[top_all0,]
scedatalist2<-scedatalist2[,rownames(cellInfo)]
labels_row = rownames(scedatalist2)

library(ComplexHeatmap)
cellInfo2_n3<-c()
for(i in 1:nrow(cellInfo2_n)){
    cellInfo2_n3<-c(cellInfo2_n3,rep(cellInfo2_n[i,1],cellInfo2_n[i,2]))
}

myneedgene2<-myneedgene[which(myneedgene$CellType != 'pDC'),]
index_at<-c()
index_save<-c()
for(i in 1:6){#fever
labels_row2<-rep('',nrow(scedatalist2))
labels_row2[(50*(i-1)+1):(50*i)]<-labels_row[(50*(i-1)+1):(50*i)]
mymark<-strsplit(myneedgene2[i,2],', ')[[1]]
index_at<-c(index_at,which(labels_row2 %in% mymark))
}

row_anno = rowAnnotation(foo = anno_mark(at = index_at,
                                         labels = labels_row[index_at]))
col_fun2 = colorRamp2(c(0,1,2),c('#007C9D01','#00979D01','#FFD700'))

mylist<-data.frame(Cell_type=c(rep('B_cell',50),rep('cDC',50),rep('monocyte',50),
                               rep('NK_cell',50),rep('NKT_cell',50),rep('T_cell',50)))

mylist<-as.matrix(mylist)

# fever
colors <- structure(c('#B7282A','#775948','#8C6AB9','#E48328','#5D9C3D','#4E78B3'),
                    names=c("B_cell", "cDC",
                           "monocyte", "NK_cell", 
                           "NKT_cell",  "T_cell"
                           ))

a2<-Heatmap(mylist,col=colors,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        width = 20,#40,20,40
        show_column_names = FALSE,
        show_heatmap_legend = FALSE)

a<-Heatmap(scedatalist2,
col = col_fun2,
cluster_rows = FALSE,
cluster_columns = FALSE,
show_column_names = FALSE,
show_row_names = FALSE,
column_split = cellInfo2_n3,
top_annotation = top_anno,
column_title = NULL,
row_names_side = "left",
row_names_gp = gpar(fontsize = 7),
right_annotation = row_anno,
heatmap_legend_param = list(title = "")
)
pdf('findmarker/fever.pdf',width = 10,height =10)
a2+a
dev.off()
