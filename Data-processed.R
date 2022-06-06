rm(list = ls())
exp_all<-read.table("D:\\xg\\data\\DEV\\expression_res.txt",header = T,as.is = T,sep = "\t",stringsAsFactors = F)
###########################cell infomation
cell_info<-read.table("D:\\xg\\data\\DEV\\single_cell\\GSE116672_clincal_info.txt",row.names=1,header = T,as.is = T,sep = "\t",stringsAsFactors = F)
##########################remove unknown cells
cell_info_res<-cell_info[-which(cell_info$Cell_type=="unknown"),]
write.table(cell_info_res,"D:\\xg\\data\\DEV\\result\\data\\cell_info_res.txt",row.names = T,col.names = T,sep = "\t",
            quote = F)
exp_known<-exp_all[,rownames(cell_info_res)]
mRNA_list<-read.table("D:\\xg\\data\\GRch38\\mRNA.txt",header = T,as.is = T,sep = "\t",stringsAsFactors = F)
exp_mRNA<-exp_known[intersect(mRNA_list$gene_name,rownames(exp_known)),]
##############################classed by Diseaseseverity
info_control<-cell_info_res[which(cell_info_res$Diseaseseverity=="control"),]
info_fever<-cell_info_res[which(cell_info_res$Diseaseseverity=="dengue_fever"),]
info_severe<-cell_info_res[which(cell_info_res$Diseaseseverity=="severe_dengue"),]

control_exp<-exp_mRNA[,intersect(rownames(info_control),colnames(exp_mRNA))]
fever_exp<-exp_mRNA[,intersect(rownames(info_fever),colnames(exp_mRNA))]
severe_exp<-exp_mRNA[,intersect(rownames(info_severe),colnames(exp_mRNA))]
#
write.table(exp_mRNA,"D:\\xg\\data\\Dengue\\result\\data\\mRNAexp_known.txt",row.names = T,col.names = T,sep = "\t",
            quote = F)

write.table(info_control,"D:\\xg\\data\\Dengue\\result\\data\\info_control.txt",row.names = T,col.names = T,sep = "\t",
            quote = F)
write.table(info_fever,"D:\\xg\\data\\Dengue\\result\\data\\info_fever.txt",row.names = T,col.names = T,sep = "\t",
            quote = F)
write.table(info_severe,"D:\\xg\\data\\Dengue\\result\\data\\info_severe.txt",row.names = T,col.names = T,sep = "\t",
            quote = F)

write.table(control_exp,"D:\\xg\\data\\Dengue\\result\\data\\control_exp.txt",row.names = T,col.names = T,sep = "\t",
            quote = F)
write.table(fever_exp,"D:\\xg\\data\\Dengue\\result\\data\\fever_exp.txt",row.names = T,col.names = T,sep = "\t",
            quote = F)
write.table(severe_exp,"D:\\xg\\data\\Dengue\\result\\data\\severe_exp.txt",row.names = T,col.names = T,sep = "\t",
            quote = F)
