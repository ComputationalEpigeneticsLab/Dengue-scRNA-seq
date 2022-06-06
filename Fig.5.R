####https://github.com/sqjin/CellChat
rm(list=ls())
load("D:\\xg\\data\\Dengue\\Rdata\\Seurat_control.RData")

setwd("D:\\xg\\data\\Dengue\\result\\Cellchat_res\\control")
library(CellChat)
library(ggalluvial)
library(CellChat)
library(ggplot2)
library(ggalluvial)
library(svglite)
library(Seurat)
#pDC_cell<-rownames(severe_object@meta.data[which(severe_object@meta.data$Cell_type=="pDC"),])
meta<-control_object@meta.data[-which(control_object@meta.data$Cell_type=="pDC"),]
data.input  <- as.matrix(as.data.frame(control_object@assays$RNA@data)[,rownames(meta)])
#data.input<-data.input[,-which( colnames(data.input) %in% pDC_cell )]
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "Cell_type")

groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group
#####Set the ligand-receptor interaction database
CellChatDB <- CellChatDB.human# use CellChatDB.human if running on human data

showDatabaseCategory(CellChatDB)

dplyr::glimpse(CellChatDB$interaction)# Show the structure of the database

CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")# use Secreted Signaling for cell-cell communication analysis

cellchat@DB <- CellChatDB.use# set the used database in the object

cellchat <- subsetData(cellchat) # subset the expression data of signaling genes for saving computation cost

future::plan("multiprocess", workers = 4) # do parallel

cellchat <- identifyOverExpressedGenes(cellchat)

cellchat <- identifyOverExpressedInteractions(cellchat)

cellchat <- projectData(cellchat, PPI.human)

#cellchat <- computeCommunProb(cellchat)
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 10)
#######
cellchat <- computeCommunProbPathway(cellchat)

cellchat <- aggregateNet(cellchat)

cellchat@netP$pathways

setwd("D:\\xg\\data\\Dengue\\result\\Cellchat_res\\control\\part1")
df.net <- subsetCommunication(cellchat)
par(mfrow = c(1,2), xpd=TRUE)
pdf("all.cell.interactions.pdf",width = 10,height = 10)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, 
                 label.edge= F, title.name = "Number of interactions")
dev.off()
#$####3
pdf("weight_all.cell.interactions.pdf",width = 10,height = 10)
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T,
                 label.edge= F, title.name = "Interaction weights/strength")
dev.off()
######mat <- cellchat@net$weight
mat <- cellchat@net$weight

par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  pdf(paste0(rownames(mat)[i],"_netVisual_circle.pdf"),width = 10,height = 10)
  p<-netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
  print(p)
  dev.off()
}
########part2
pathways<-cellchat@netP$pathways
for(j in 1:length(pathways)){
  vertex.receiver = seq(1,4)
  result.name <- pathways[j]
  if(!dir.exists(result.name)){
    dir.create(result.name)
  }
}

###
for(i in 1:length(pathways)){
  pathways.show<-pathways[i]
  #####
  setwd(paste0("D:\\xg\\data\\Dengue\\result\\Cellchat_res\\control\\part1\\",pathways.show))
  #####
  pdf(paste0(pathways.show,"_signaling_pathway_network_hierarchy.pdf"),height = 10,width = 10)
  vertex.receiver = seq(1,4) # a numeric vector. 
  p<-netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)
  print(p)
  dev.off()
  #################
  pdf(paste0(pathways.show,"_signaling_pathway_network_circle.pdf"),height = 10,width = 10)
  par(mfrow=c(1,1))
  p<-netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
  print(p)
  dev.off()
  ######
  pdf(paste0(pathways.show,"_signaling_pathway_network_chord.pdf"),height = 10,width = 10)
  par(mfrow=c(1,1))
  p<-netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
  print(p)
  dev.off()
  ########
  pdf(paste0(pathways.show,"_signaling_network_heatmap.pdf"),height = 10,width = 10)
  par(mfrow=c(1,1))
  p<-netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
  print(p)
  dev.off()
  #######
  pdf(paste0(pathways.show,"_signaling_network_contribution.pdf"),height = 10,width = 10)
  p<-netAnalysis_contribution(cellchat, signaling = pathways.show)
  print(p)
  dev.off()
  ######
  pairLR.pathway <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)
  #####
  for(j in 1:length(pairLR.pathway)){
    LR.show <- pairLR.pathway[j,]
    vertex.receiver = seq(1,4)
    pdf(paste0(LR.show[j],"_signaling_network.pdf"),height = 10,width = 10)
    p1<-netVisual_individual(cellchat, signaling = pathways.show,  pairLR.use = LR.show, vertex.receiver = vertex.receiver)
    print(p1)
    dev.off()
    #####
    pdf(paste0(LR.show[j],"_signaling_network_circle.pdf"),height = 10,width = 10)
    p2<-netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")
    print(p2)
    dev.off()
    #####
    pdf(paste0(LR.show[j],"_signaling_network_chord.pdf"),height = 10,width = 10)
    p3<-netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "chord")
    print(p3)
    dev.off()
  }
  
}

setwd("D:\\xg\\data\\Dengue\\result\\Cellchat_res\\control\\part2")
for(k in 1:6){
  pdf(paste0("netVisual_bubble_","k=",k,".pdf"),height = 12,width = 10)
  p<-netVisual_bubble(cellchat, sources.use = k, targets.use = c(1:6), remove.isolate = FALSE)
  print(p)
  dev.off()
  ######
  pairLR.use <- extractEnrichedLR(cellchat, signaling = cellchat@netP$pathways)
  pdf(paste0("all_pathway_pairLR.use_k=",k,".pdf"),height = 12,width = 10)
  p<-netVisual_bubble(cellchat, sources.use = k, targets.use = c(1:6), pairLR.use = pairLR.use, remove.isolate = TRUE)
  print(p)
  dev.off()
  ########
  pdf(paste0("netVisual_chord_gene_k=",k,".pdf"),height = 20,width = 20)
  p<-netVisual_chord_gene(cellchat, sources.use = k, targets.use = c(1:6), lab.cex = 1,legend.pos.y = 30)
  print(p)
  dev.off()
  ######
  
  
}

###part3
setwd("D:\\xg\\data\\Dengue\\result\\Cellchat_res\\control\\part3")
result.name <- "plotGeneExpression"
if(!dir.exists(result.name)){
  dir.create(result.name)
}
for(n in (1:length(pathways))){
  
  pdf(paste0(pathways[n],"_plotGeneExpression.pdf"),height = 10,width = 10)
  p<-plotGeneExpression(cellchat, signaling = pathways[n])
  print(p)
  dev.off()
}

########
result.name <- "netAnalysis_signalingRole_network"
if(!dir.exists(result.name)){
  dir.create(result.name)
}
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
for(i in 1:length(pathways)){
  pdf(paste0(result.name,"\\",pathways[i],"_netAnalysis_signalingRole_network.pdf"),height = 10,width = 10)
  p<-netAnalysis_signalingRole_network(cellchat, signaling = pathways[i], width = 12, height = 4, font.size = 10)
  print(p)
  dev.off()
  
}

#####
pdf("netAnalysis_signalingRole_heatmap.pdf",height = 18,width = 18)
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing",width = 14, height = 18,font.size=14)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming",width = 14, height = 18,font.size = 14)
print(ht1 + ht2)

dev.off()


selectK(cellchat, pattern = "outgoing")
result.name <- "patterns"
if(!dir.exists(result.name)){
  dir.create(result.name)
}
nPatterns = 4
pdf(paste0(result.name,"\\","outgoing_CommunicationPatterns4.pdf"),width = 10,height = 10)
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns,font.size = 6,height = 8)
dev.off()
###
pdf(paste0(result.name,"\\","outgoing_netAnalysis_river4.pdf"),height = 10,width = 10)
netAnalysis_river(cellchat, pattern = "outgoing",font.size = 4)
dev.off()
###
pdf(paste0(result.name,"\\","outgoing_netAnalysis_dot4.pdf"),height = 6,width = 8)

netAnalysis_dot(cellchat, pattern = "outgoing")
dev.off()
#############

#############

selectK(cellchat, pattern = "incoming")
####
nPatterns = 5
pdf(paste0(result.name,"\\","incoming_CommunicationPatterns.pdf"),height = 10,width = 10)

cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns,font.size = 6,height = 8)
dev.off()
#####
pdf(paste0(result.name,"\\","incoming_netAnalysis_river.pdf"),height = 10,width = 10)

netAnalysis_river(cellchat, pattern = "incoming",font.size = 4)
dev.off()
#######
pdf(paste0(result.name,"\\","incoming_netAnalysis_dot.pdf"),height = 6,width = 8)

netAnalysis_dot(cellchat, pattern = "incoming")
dev.off()
#############
###Identify signaling groups based on function similarity
cellchat <- computeNetSimilarity(cellchat, type = "functional")
reticulate::py_install(packages = 'umap-learn')

cellchat <- netEmbedding(cellchat, type = "functional")
#reticulate::py_install(packages = 'umap-learn')

cellchat <- netClustering(cellchat, type = "functional")
pdf("D:\\xg\\data\\Dengue\\result\\Cellchat_res\\control\\part3\\similarity\\function.pdf",height = 10,width = 10)
netVisual_embedding(cellchat, type = "functional", label.size = 3.5,font.size = 12)
dev.off()
###########
pdf("D:\\xg\\data\\Dengue\\result\\Cellchat_res\\control\\part3\\similarity\\function_col.pdf",height = 10,width = 10)
netVisual_embeddingZoomIn(cellchat, type = "functional", nCol = 2)
dev.off()
#####Identify signaling groups based on structure similarity

cellchat <- computeNetSimilarity(cellchat, type = "structural")
cellchat <- netEmbedding(cellchat, type = "structural")
#> Manifold learning of the signaling networks for a single dataset
cellchat <- netClustering(cellchat, type = "structural")
#> Classification learning of the signaling networks for a single dataset
# Visualization in 2D-space
pdf("D:\\xg\\data\\Dengue\\result\\Cellchat_res\\control\\part3\\similarity\\structure.pdf",height = 10,width = 10)

netVisual_embedding(cellchat, type = "structural", label.size = 3.5)
dev.off()
#############
pdf("D:\\xg\\data\\Dengue\\result\\Cellchat_res\\control\\part3\\similarity\\structure_col.pdf",height = 10,width = 10)
netVisual_embeddingZoomIn(cellchat, type = "structural", nCol = 2)
dev.off()
saveRDS(cellchat, file = "D:\\xg\\data\\Dengue\\result\\Cellchat_res\\control.rds")
##########











###############Cellchat: compare control and fever
rm(list = ls())
control_cellchat<- readRDS("D:\\xg\\data\\Dengue\\result\\Cellchat_res\\control.rds")
fever_cellchat <- readRDS("D:\\xg\\data\\Dengue\\result\\Cellchat_res\\fever.rds")
#severe_cellchat <- readRDS("D:\\xg\\data\\Dengue\\result\\Cellchat_res\\severe.rds")
###control、fever
object.list <- list(control = control_cellchat, fever = fever_cellchat)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
pdf("D:\\xg\\data\\Dengue\\result\\Cellchat_res\\compare_cf\\compareInteractions.pdf")
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2
dev.off()

#####Compare the number of interactions and interaction strength among different cell populations
pdf("D:\\xg\\data\\Dengue\\result\\Cellchat_res\\compare_cf\\netVisual_diffInteraction.pdf")
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")
dev.off()
######
pdf("D:\\xg\\data\\Dengue\\result\\Cellchat_res\\compare_cf\\netVisual_heatmap.pdf")
gg1 <- netVisual_heatmap(cellchat)
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
#> Do heatmap based on a merged object
gg1 + gg2
dev.off()
#####
pdf("D:\\xg\\data\\Dengue\\result\\Cellchat_res\\compare_cf\\netVisual_circle.pdf")
weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}

dev.off()
###################################Identify signaling groups based on their functional similarity
cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
#> Compute signaling network similarity for datasets 1 2
cellchat <- netEmbedding(cellchat, type = "functional")
#> Manifold learning of the signaling networks for datasets 1 2
cellchat <- netClustering(cellchat, type = "functional")
#> Classification learning of the signaling networks for datasets 1 2
# Visualization in 2D-space
pdf("D:\\xg\\data\\Dengue\\result\\Cellchat_res\\compare_cf\\similarity_function.pdf",width = 10,height = 8)
netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 4,point.shape=c(21,0))
dev.off()
pdf("D:\\xg\\data\\Dengue\\result\\Cellchat_res\\compare_cf\\similarity_function2.pdf",width = 10,height = 8)
netVisual_embeddingPairwiseZoomIn(cellchat, type = "functional", nCol = 2)
dev.off()
######################################Identify signaling groups based on their structural similarity
cellchat <- computeNetSimilarityPairwise(cellchat, type = "structural")
#> Compute signaling network similarity for datasets 1 2
cellchat <- netEmbedding(cellchat, type = "structural")
#> Manifold learning of the signaling networks for datasets 1 2
cellchat <- netClustering(cellchat, type = "structural")
#> Classification learning of the signaling networks for datasets 1 2
# Visualization in 2D-space
pdf("D:\\xg\\data\\Dengue\\result\\Cellchat_res\\compare_cf\\similarity_structual.pdf",width = 10,height = 8)

netVisual_embeddingPairwise(cellchat, type = "structural", label.size = 3.5)
dev.off()
pdf("D:\\xg\\data\\Dengue\\result\\Cellchat_res\\compare_cf\\similarity_structual2.pdf",width = 10,height = 8)
netVisual_embeddingPairwiseZoomIn(cellchat, type = "functional", nCol = 2)
dev.off()
###############Larger distance implies larger difference of the communication networks between two datasets in terms of either functional or structure similarity.
pdf("D:\\xg\\data\\Dengue\\result\\Cellchat_res\\compare_cf\\distance.pdf")

rankSimilarity(cellchat, type = "functional")
dev.off()
##########################################Compare the overall information flow of each signaling pathway
pdf("D:\\xg\\data\\Dengue\\result\\Cellchat_res\\compare_cf\\infomation_flow.pdf")
gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE)
gg1 + gg2
dev.off()
####################################Compare outgoing (or incoming) signaling associated with each cell population
library(ComplexHeatmap)
i = 1
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 8, height = 12)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 8, height = 12)
pdf("D:\\xg\\data\\Dengue\\result\\Cellchat_res\\compare_cf\\netAnalysis_signalingRole_heatmap_outgoing.pdf",width = 12,height = 12)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()
###
pdf("D:\\xg\\data\\Dengue\\result\\Cellchat_res\\compare_cf\\netAnalysis_signalingRole_heatmap_incoming.pdf",width = 12,height = 12)

ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i], width = 8, height = 12, color.heatmap = "GnBu")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+1], width = 8,, height = 12, color.heatmap = "GnBu")

draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()
###########################################################Identify the upgulated and down-regulated signaling ligand-receptor pairs
setwd("D:\\xg\\data\\Dengue\\result\\Cellchat_res\\compare_cf\\netVisual_bubble")
for(i in c(1:6)){
  pdf(paste0("netVisual_bubble_source=",i,".pdf"),width = 10,height = 14)
  p<-netVisual_bubble(cellchat, sources.use = i, targets.use = c(1:6),  comparison = c(1, 2), angle.x = 45)
  print(p)
  dev.off()
}
#####

gg1 <- netVisual_bubble(cellchat, sources.use = 3, targets.use = c(1:6),  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in LS", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg2 <- netVisual_bubble(cellchat, sources.use = 3, targets.use = c(1:6),  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in LS", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg1 + gg2
######Compare the signaling gene expression distribution between different datasets

cellchat@meta$datasets = factor(cellchat@meta$datasets, levels = c("control", "fever")) # set factor level
pathway.show<-intersect(object.list[[1]]@netP$pathways,object.list[[2]]@netP$pathways)[-c(3,5,6,7)][-22]
setwd("D:\\xg\\data\\Dengue\\result\\Cellchat_res\\compare_cf\\plotGeneExpression")
for(i in 1:length(pathway.show)){
  pdf(paste0(pathway.show[i],".pdf"))
  p<-plotGeneExpression(cellchat, signaling = pathway.show[i], split.by = "datasets", colors.ggplot = T)
  print(p)
  dev.off()
  cat(i,"\n")
}
plotGeneExpression(cellchat, signaling = "CXCL", split.by = "datasets", colors.ggplot = T)
object.list[[1]]@netP$pathways

pathway.show01<-intersect(object.list[[1]]@netP$pathways,object.list[[2]]@netP$pathways)
#####
saveRDS(cellchat, file = "D:\\xg\\data\\Dengue\\result\\Cellchat_res\\compare_cf.rds")
