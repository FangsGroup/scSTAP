#Transcriptome and proteome profile
windowsFonts(ARL = windowsFont("Arial"))

#set work path
setwd("D:/Research/R/Data")

#Bar-Figure 4A
library(ggplot2)

RNA<-read.csv("quan_RNA.csv",sep=",",header=T)
PRO<-read.csv("quan_PRO.csv",sep=",",header=T)

p<-ggplot()+
  geom_col(data=PRO,aes(x=Sample,y=Number,fill=Group),position=position_dodge(width=0.8),width=0.6,color="white")+
  geom_col(data=RNA,aes(x=Sample,y=Number,fill=Group),position=position_dodge(width=0.8),width=0.6,color="white")+
  scale_y_continuous(limits=c(0,6000),breaks=scales::pretty_breaks(n=5),sec.axis = sec_axis( ~.*5,name="Genes",breaks=scales::pretty_breaks(n=5)))+
  labs(x='Single oocytes',y='Protein groups')+
  scale_fill_manual(values = c("#C3416B","#416595"),labels= c("Protein groups","Genes"))+
  coord_fixed(ratio=0.0006)+
  theme(panel.background=element_rect(fill="white",colour="black",size=0.25),
        panel.border=element_rect(fill=NA,color="black",size=2,linetype="solid"),
        axis.line=element_line(colour="black",size=0),
        axis.ticks=element_line(colour="black",size=0.8),
        axis.ticks.length=unit(0.5,'line'),
        axis.title.x=element_text(size=22,face="plain",color="black",vjust=-0.2),
        axis.title.y=element_text(size=22,face="plain",color="black",vjust=2),
        axis.title.y.right=element_text(size=22,face="plain",color="black",vjust=2),
        axis.text.x=element_text(size=16,face="plain",color="black",vjust=0.1),
        axis.text.y=element_text(size=18,face="plain",color="black",hjust=0.2),
        axis.text.y.right=element_text(size=18,face="plain",color="black",hjust=0.2),
        legend.text=element_text(size=18,face="plain",colour="black"),
        legend.title=element_blank(),
        legend.margin=margin(0.1,0.1,0.1,0.1,'cm'),
        legend.key.size=unit(0.7,'cm'),
        legend.position=c(0.9,0.83))+
  guides(fill=guide_legend(nrow=2,byrow=T))

ggsave("Bar-T&P.png",p,width=16,height=4,dpi=600)


#For GA
#Transcriptome and proteome profile
windowsFonts(ARL = windowsFont("Arial"))

#set work path
setwd("D:/Research/R/Data")

#Bar-Figure 4A
library(ggplot2)

RNA<-read.csv("quan_RNA.csv",sep=",",header=T)
PRO<-read.csv("quan_PRO.csv",sep=",",header=T)

p<-ggplot()+
  geom_col(data=PRO,aes(x=Sample,y=Number,fill=Group),position=position_dodge(width=0.8),width=0.7,color="white")+
  geom_col(data=RNA,aes(x=Sample,y=Number,fill=Group),position=position_dodge(width=0.8),width=0.7,color="white")+
  scale_y_continuous(limits=c(0,6000),breaks=scales::pretty_breaks(n=5),sec.axis = sec_axis( ~.*5,name="Genes",breaks=scales::pretty_breaks(n=5)))+
  labs(x='Single oocytes',y='Protein groups')+
  scale_fill_manual(values = c("#C3416B","#416595"),labels= c("Protein groups","Genes"))+
  coord_fixed(ratio=0.00045)+
  theme(panel.background=element_rect(fill="white",colour="black",size=0.25),
        panel.border=element_rect(fill=NA,color="black",size=2,linetype="solid"),
        axis.line=element_line(colour="black",size=0),
        axis.ticks=element_line(colour="black",size=0.8),
        axis.ticks.length=unit(0.5,'line'),
        axis.title.x=element_text(size=26,face="plain",color="black",vjust=0),
        axis.title.y=element_text(size=26,face="plain",color="black",vjust=2),
        axis.title.y.right=element_text(size=26,face="plain",color="black",vjust=2),
        axis.text.x=element_text(size=22,face="plain",color="black",vjust=0.1),
        axis.text.y=element_text(size=22,face="plain",color="black",hjust=0.2),
        axis.text.y.right=element_text(size=22,face="plain",color="black",hjust=0.2),
        legend.text=element_text(size=22,face="plain",colour="black"),
        legend.title=element_blank(),
        legend.margin=margin(0.1,0.1,0.1,0.1,'cm'),
        legend.key.size=unit(0.7,'cm'),
        legend.position=c(0.9,0.83))+
  guides(fill=guide_legend(nrow=2,byrow=T))

ggsave("Bar-T&P.png",p,width=19,height=4,dpi=1200)




#PCA
library(dplyr)
library(Seurat)
library(SeuratData)
library(SeuratDisk)
library(cowplot)
library("Matrix")
library(AnnotationDbi)
library(clusterProfiler)
library(org.Mm.eg.db)
library(ggplot2)
library(ggrepel)

#Transcriptome-Figure 4B
RNA<-read.csv("RNA_15.csv",sep=",",header=T)
RNA_SMB<-bitr(RNA$gene_id,fromType="ENSEMBL",toType="SYMBOL",OrgDb=org.Mm.eg.db)

MATQ_RNA<-merge(RNA_SMB,RNA,by.x="ENSEMBL",by.y="gene_id")
MATQ_RNA<- MATQ_RNA%>% distinct(SYMBOL, .keep_all = T)
RNA<-MATQ_RNA[,-(1:2)]
rownames(RNA)<-MATQ_RNA$SYMBOL

mmoc<-CreateSeuratObject(counts=RNA,assay="RNA")
mmoc<-NormalizeData(mmoc, assay="RNA", normalization.method="LogNormalize")
mmoc<-ScaleData(mmoc, assay="RNA")
mmoc<-FindVariableFeatures(mmoc,selection.method="vst",nfeatures = 2000, assay="RNA")
mmoc<-RunPCA(mmoc, npcs=7, reduction.name="pcaRNA",reduction.key="pcaRNA_", assay="RNA")

mmoc@reductions$pcaRNA@stdev
#[1] 23.708072 20.015325 17.591102 13.175478 11.850388  8.827989  8.412763

mmoc<- FindNeighbors(mmoc, reduction = "pcaRNA", k.param = 7, dims=1:7, assay="RNA")
mmoc<- FindClusters(mmoc, resolution = 0.8)

pcaRNA<-data.frame(mmoc@meta.data,mmoc@reductions$pcaRNA@cell.embeddings,group=c(rep("GV", 9),rep("MII", 6)))

label<-rownames(pcaRNA)

label[2]<-""
label[3]<-""
label[6]<-""
label[10]<-""
label[11]<-""

p<-ggplot(pcaRNA)+
  geom_point(aes(x=pcaRNA_1,y=pcaRNA_2,fill=mmoc@active.ident),size=7,shape=21,stroke=0.7,color="#004F64")+
  geom_text_repel(aes(x=pcaRNA_1,y=pcaRNA_2,label=label,color=group),size = 7,box.padding = unit(0.6,"lines"),point.padding = unit(1,"lines"),
                  segment.color="transparent",show.legend = FALSE)+
  scale_y_continuous(limits=c(-60,50))+
  scale_x_continuous(limits=c(-60,50))+
  labs(x='pcaRNA_1 (23.7%)',y='pcaRNA_2 (20.0%)')+
  geom_text(x=-9,y=10,aes(label="GV2"),parse=T,size=7,color="#5D110B")+
  geom_text(x=20,y=2,aes(label="GV3"),parse=T,size=7,color="#5D110B")+
  geom_text(x=3,y=3,aes(label="GV6"),parse=T,size=7,color="#5D110B")+
  geom_text(x=10,y=22,aes(label="MII1"),parse=T,size=7,color="#004F64")+
  geom_text(x=5,y=16,aes(label="MII2"),parse=T,size=7,color="#004F64")+
  coord_fixed(ratio=1)+
  scale_fill_manual(values=c("#5179C9"),labels=c("Cluster 1"))+
  scale_color_manual(values=c("#5D110B","#004F64"))+
  ggtitle("Transcriptome")+
  theme(panel.background=element_rect(fill="white",colour="black",size=0.25),
        panel.border=element_rect(fill=NA,color="black",size=3,linetype="solid"),
        plot.title=element_text(hjust=0.5,size=40,color="black"),
        axis.line=element_line(colour="black",size=0),
        axis.ticks=element_line(colour="black",size=1.4),
        axis.ticks.length=unit(0.5,'line'),
        axis.title.x=element_text(size=28,face="plain",color="black",vjust=-0.5),
        axis.title.y=element_text(size=28,face="plain",color="black",vjust=1.7),
        axis.text.x=element_text(size=25,face="plain",color="black",vjust=-0.2),
        axis.text.y=element_text(size=25,face="plain",color="black"),
        legend.text=element_text(size=25,face="plain",colour="black"),
        legend.position=c(0.20,0.88),
        legend.title=element_blank(),
        legend.margin=margin(0.1,0.1,0.1,0.1,'cm'),
        legend.key.height=unit(1.1,'cm'),
        legend.key=element_rect(fill="white"))+
  guides(fill=guide_legend(nrow=2,byrow=T))

ggsave("pcaRNA.png",p,width=7,height=7.5,dpi=1200)

#Proteome-Figure 4C
PRO<-read.csv("PRO_15.csv",sep=",",header=T)
PRO_SMB<-bitr(PRO$PG.ProteinGroups,fromType="UNIPROT",toType="SYMBOL",OrgDb=org.Mm.eg.db)

MATQ_PRO<-merge(PRO_SMB,PRO,by.x="UNIPROT",by.y="PG.ProteinGroups")
MATQ_PRO<- MATQ_PRO%>% distinct(SYMBOL, .keep_all = T)
PRO<-MATQ_PRO[,-(1:2)]
rownames(PRO)<-MATQ_PRO$SYMBOL

mmoc[["PRO"]]<-CreateAssayObject(counts=PRO)
mmoc<-NormalizeData(mmoc, assay="PRO", normalization.method="LogNormalize")
mmoc<-ScaleData(mmoc, assay="PRO")
mmoc<-FindVariableFeatures(mmoc,selection.method="vst",nfeatures = 2000, assay="PRO")
mmoc<-RunPCA(mmoc, npcs=7, reduction.name="pcaPRO",reduction.key="pcaPRO_",assay="PRO")

mmoc@reductions$pcaPRO@stdev
#[1] 21.323544 17.187363 17.061533 11.454338 10.258102  9.847127  9.586008

mmoc<- FindNeighbors(mmoc, dims=1:7, k.param = 7, reduction = "pcaPRO", assay="PRO")
mmoc<- FindClusters(mmoc, resolution = 0.8)

pcaPRO<-data.frame(mmoc@meta.data,mmoc@reductions$pcaPRO@cell.embeddings,group=c(rep("GV", 9),rep("MII", 6)))

p<-ggplot(pcaPRO)+
  geom_point(aes(x=pcaPRO_1,y=pcaPRO_2,fill=mmoc@active.ident),size=7,shape=21,stroke=0.7,color="#5D110B")+
  geom_text_repel(aes(x=pcaPRO_1,y=pcaPRO_2,label=rownames(pcaPRO),color=group),size = 7,box.padding = unit(0.6,"lines"),point.padding = unit(1,"lines"),
                  segment.color="transparent",show.legend = FALSE)+
  scale_y_continuous(limits=c(-40,40))+
  scale_x_continuous(limits=c(-40,40))+
  labs(x='pcaPRO_1 (21.3%)',y='pcaPRO_2 (17.2%)')+
  coord_fixed(ratio=1)+
  scale_fill_manual(values=c("#E43075"),labels=c("Cluster 1"))+
  scale_color_manual(values=c("#5D110B","#004F64"))+
  ggtitle("Proteome")+
  theme(panel.background=element_rect(fill="white",colour="black",size=0.25),
        panel.border=element_rect(fill=NA,color="black",size=3,linetype="solid"),
        plot.title=element_text(hjust=0.5,size=40,color="black"),
        axis.line=element_line(colour="black",size=0),
        axis.ticks=element_line(colour="black",size=1.4),
        axis.ticks.length=unit(0.5,'line'),
        axis.title.x=element_text(size=28,face="plain",color="black",vjust=-0.5),
        axis.title.y=element_text(size=28,face="plain",color="black",vjust=1.7),
        axis.text.x=element_text(size=25,face="plain",color="black",vjust=-0.2),
        axis.text.y=element_text(size=25,face="plain",color="black"),
        legend.text=element_text(size=25,face="plain",colour="black"),
        legend.position=c(0.20,0.88),
        legend.title=element_blank(),
        legend.margin=margin(0.1,0.1,0.1,0.1,'cm'),
        legend.key.height=unit(1.1,'cm'),
        legend.key=element_rect(fill="white"))+
  guides(fill=guide_legend(nrow=2,byrow=T))

ggsave("pcaPRO.png",p,width=7,height=7.5,dpi=1200)

#DEGs-limma
library(ggbeeswarm)
library(limma)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(AnnotationDbi)
library(clusterProfiler)
library(org.Mm.eg.db)
library(Hmisc)

#Transcriptome-Figure 4E
RNA<-read.csv("RNA_15.csv",sep=",",header=T)
RNA_SMB<-bitr(RNA$gene_id,fromType="ENSEMBL",toType="SYMBOL",OrgDb=org.Mm.eg.db)

MATQ_RNA<-merge(RNA_SMB,RNA,by.x="ENSEMBL",by.y="gene_id")
MATQ_RNA<- MATQ_RNA%>% distinct(SYMBOL, .keep_all = T)
Data<-MATQ_RNA[,-(1:2)]
rownames(Data)<-MATQ_RNA$SYMBOL

Data.log<-log2(Data)

group_list <- factor(group[,2],ordered = F)
design <- model.matrix(~0+group_list)
colnames(design) <- levels(group_list)
rownames(design) <- colnames(Data.log)

fit <- lmFit(Data.log,design)

cont.matrix <- makeContrasts(contrasts = paste0(unique(group_list),collapse = "-"),levels = design)
fit2 <- contrasts.fit(fit,cont.matrix)
fit2 <- eBayes(fit2)
tmpOut <- topTable(fit2,coef = 1,n = Inf,adjust.method = "BH",sort.by = "logFC")
limma.na <- na.omit(tmpOut)

limma.na$change <- as.factor(ifelse(limma.na$adj.P.Val<0.05 & abs(limma.na$logFC)>log2(2),ifelse(limma.na$logFC>log2(2),"Down","Up"),"Non"))
limma.na$label <- ifelse(limma.na$adj.P.Val<0.05 & abs(limma.na$logFC)>log2(2),as.character(rownames(limma.na)),"")

limma.na$label[3]<-""
limma.na$label[289]<-""
limma.na$label[825]<-""

p<-ggplot(data=limma.na)+
  geom_point(aes(x = -logFC,y = -log10(adj.P.Val),fill=change),color="white",size=5,shape=21,stroke=0.2,alpha=0.65)+
  geom_text_repel(aes(x = -logFC,y = -log10(adj.P.Val),label = label),size = 7,force=2,box.padding = unit(0.9,"lines"),point.padding = unit(0.7,"lines"),segment.color = "black",show.legend = FALSE)+
  geom_vline(xintercept = c(-log2(2),log2(2)),lty = 4,col = "black",lwd = 0.8)+
  geom_hline(yintercept = -log10(0.05),lty = 4,col = "black",lwd = 0.8)+
  scale_y_continuous(limits=c(0,5))+
  scale_x_continuous(limits=c(-5,4))+
  scale_fill_manual(values = c("#0000FF","#A3A3A3","#CC0000"))+
  xlab(expression(Log[2](FoldChange)))+
  ylab(expression(-Log[10](italic(p.adjust))))+
  theme_classic(base_family="ARL")+
  coord_fixed(ratio=1.8)+
  ggtitle("Transcriptome")+
  theme(panel.background=element_rect(fill="white",colour="black",size=0.25),
        panel.border=element_rect(fill=NA,color="black",size=3,linetype="solid"),
        plot.title=element_text(hjust=0.5,size=40,color="black"),
        axis.line=element_line(colour="black",size=0),
        axis.ticks=element_line(colour="black",size=1.4),
        axis.ticks.length=unit(0.5,'line'),
        axis.title.x=element_text(size=28,face="plain",color="black",vjust=0),
        axis.title.y=element_text(size=28,face="plain",color="black",vjust=2),
        axis.text.x=element_text(size=25,face="plain",color="black",vjust=0),
        axis.text.y=element_text(size=25,face="plain",color="black"),
        legend.text=element_text(size=25,face="plain",colour="black"),
        legend.position=c(0.85,0.88),
        legend.title=element_blank(),
        legend.margin=margin(0.1,0.1,0.1,0.1,'cm'),
        legend.background=element_rect(fill="transparent"),
        legend.key.height=unit(0.5,'cm'),
        legend.key=element_rect(fill="transparent",color="transparent"))+
  guides(fill=guide_legend(nrow=3,byrow=T))

ggsave("Volcano_RNA.png",p,width=7,height=7.7,dpi=1200)

#Proteome-Figure 4F
PRO<-read.csv("PRO_15.csv",sep=",",header=T)
PRO_SMB<-bitr(PRO$PG.ProteinGroups,fromType="UNIPROT",toType="SYMBOL",OrgDb=org.Mm.eg.db)

MATQ_PRO<-merge(PRO_SMB,PRO,by.x="UNIPROT",by.y="PG.ProteinGroups")
MATQ_PRO<- MATQ_PRO%>% distinct(SYMBOL, .keep_all = T)
Data<-MATQ_PRO[,-(1:2)]
rownames(Data)<-MATQ_PRO$SYMBOL

Data.log<-log2(Data)

group<-read.csv("Group.csv",sep=",",header=T)

group_list <- factor(group[,2],ordered = F)
design <- model.matrix(~0+group_list)
colnames(design) <- levels(group_list)
rownames(design) <- colnames(Data.log)

fit <- lmFit(Data.log,design)

cont.matrix <- makeContrasts(contrasts = paste0(unique(group_list),collapse = "-"),levels = design)
fit2 <- contrasts.fit(fit,cont.matrix)
fit2 <- eBayes(fit2)
tmpOut <- topTable(fit2,coef = 1,n = Inf,adjust.method = "BH",sort.by = "logFC")
limma.na <- na.omit(tmpOut)

limma.na$change <- as.factor(ifelse(limma.na$adj.P.Val<0.05 & abs(limma.na$logFC)>log2(2),ifelse(limma.na$logFC>log2(2),"Down","Up"),"Non"))
limma.na$label <- ifelse(limma.na$adj.P.Val<0.05 & abs(limma.na$logFC)>log2(2),as.character(rownames(limma.na)),"")

limma.na$label[1]<-""
limma.na$label[6]<-""
limma.na$label[18]<-""
limma.na$label[50]<-""
limma.na$label[55]<-""
limma.na$label[88]<-""

p<-ggplot(data=limma.na)+
  geom_point(aes(x = -logFC,y = -log10(adj.P.Val),fill=change),color="white",size=5,shape=21,stroke=0.2,alpha=0.65)+
  geom_text_repel(aes(x = -logFC,y = -log10(adj.P.Val),label = label),size = 7,force=1, box.padding = unit(1,"lines"),point.padding = unit(0.7,"lines"),segment.color = "black",show.legend = FALSE)+
  geom_vline(xintercept = c(-log2(2),log2(2)),lty = 4,col = "black",lwd = 0.8)+
  geom_hline(yintercept = -log10(0.05),lty = 4,col = "black",lwd = 0.8)+
  scale_y_continuous(limits=c(0,12.6))+
  scale_x_continuous(limits=c(-5,4))+
  scale_fill_manual(values = c("#0000FF","#A3A3A3","#CC0000"))+
  xlab(expression(Log[2](FoldChange)))+
  ylab(expression(-Log[10](italic(p.adjust))))+
  coord_fixed(ratio=0.72)+
  ggtitle("Proteome")+
  theme_classic(base_family="ARL")+
  theme(panel.background=element_rect(fill="white",colour="black",size=0.25),
        panel.border=element_rect(fill=NA,color="black",size=3,linetype="solid"),
        plot.title=element_text(hjust=0.5,size=40,color="black"),
        axis.line=element_line(colour="black",size=0),
        axis.ticks=element_line(colour="black",size=1.4),
        axis.ticks.length=unit(0.5,'line'),
        axis.title.x=element_text(size=28,face="plain",color="black",vjust=0),
        axis.title.y=element_text(size=28,face="plain",color="black",vjust=1),
        axis.text.x=element_text(size=25,face="plain",color="black",vjust=0),
        axis.text.y=element_text(size=25,face="plain",color="black"),
        legend.text=element_text(size=25,face="plain",colour="black"),
        legend.position=c(0.85,0.88),
        legend.title=element_blank(),
        legend.margin=margin(0.1,0.1,0.1,0.1,'cm'),
        legend.background=element_rect(fill="transparent"),
        legend.key.height=unit(0.5,'cm'),
        legend.key=element_rect(fill="transparent",color="transparent"))+
  guides(fill=guide_legend(nrow=3,byrow=T))

ggsave("Volcano_PRO.png",p,width=7,height=7.7,dpi=1200)

#WNN
library(psych)

#UMAP-Figure 4D
mmoc<-FindMultiModalNeighbors(mmoc, reduction.list=list("pcaRNA","pcaPRO"),
                              dims.list=list(1:7,1:7), k.nn=5, knn.range=5, modality.weight.name=c("RNA.weight","PRO.weight"))

mmoc <- RunUMAP(mmoc, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
mmoc <- FindClusters(mmoc, graph.name = "wsnn", algorithm = 3, verbose = FALSE)

WNN<-data.frame(mmoc@meta.data,mmoc@reductions$wnn.umap@cell.embeddings,group=c(rep("GV", 9),rep("MII", 6)))

p<-ggplot(WNN)+
  geom_point(aes(x=wnnUMAP_1,y=wnnUMAP_2,fill=mmoc@active.ident,color=mmoc@active.ident),size=7,shape=21,stroke=0.7)+
  geom_text_repel(aes(x=wnnUMAP_1,y=wnnUMAP_2,label=rownames(WNN),color=mmoc@active.ident),size = 7,box.padding = unit(0.6,"lines"),point.padding = unit(1,"lines"),
                  segment.color="transparent",show.legend = FALSE)+
  scale_y_continuous(limits=c(-3,3))+
  scale_x_continuous(limits=c(-3,3))+
  coord_fixed(ratio=1)+
  scale_fill_manual(values=c("#EE776E","#00A9C9"),labels=c("Cluster 1","Cluster 2"))+
  scale_color_manual(values=c("#5D110B","#004F64"),labels=c("Cluster 1","Cluster 2"))+
  ggtitle("WNN")+
  theme(panel.background=element_rect(fill="white",colour="black",size=0.25),
        panel.border=element_rect(fill=NA,color="black",size=3,linetype="solid"),
        plot.title=element_text(hjust=0.5,size=40,color="black"),
        axis.line=element_line(colour="black",size=0),
        axis.ticks=element_line(colour="black",size=1.4),
        axis.ticks.length=unit(0.5,'line'),
        axis.title.x=element_text(size=28,face="plain",color="black",vjust=0),
        axis.title.y=element_text(size=28,face="plain",color="black",vjust=1.7),
        axis.text.x=element_text(size=25,face="plain",color="black",vjust=-0.2),
        axis.text.y=element_text(size=25,face="plain",color="black"),
        legend.text=element_text(size=25,face="plain",colour="black"),
        legend.position=c(0.20,0.88),
        legend.title=element_blank(),
        legend.margin=margin(0.1,0.1,0.1,0.1,'cm'),
        legend.key.height=unit(1.1,'cm'),
        legend.key=element_rect(fill="white"))+
  guides(fill=guide_legend(nrow=2,byrow=T))

ggsave("WNN.png",p,width=7,height=7.7,dpi=1200)

#DEGs-Figure 4G
PRO<-read.csv("PRO_15.csv",sep=",",header=T)
RNA<-read.csv("RNA_15.csv",sep=",",header=T)

PRO_SMB<-bitr(PRO$PG.ProteinGroups,fromType="UNIPROT",toType="SYMBOL",OrgDb=org.Mm.eg.db)
RNA_SMB<-bitr(RNA$gene_id,fromType="ENSEMBL",toType="SYMBOL",OrgDb=org.Mm.eg.db)

MATQ<-merge(RNA_SMB,PRO_SMB,by.x="SYMBOL",by.y="SYMBOL")
MATQ_RNA<-merge(MATQ,RNA,by.x="ENSEMBL",by.y="gene_id")
MATQ_RNA<-MATQ_RNA[order(MATQ_RNA[,2]),]
MATQ_PRO<-merge(MATQ,PRO,by.x="UNIPROT",by.y="PG.ProteinGroups")
MATQ_PRO<-MATQ_PRO[order(MATQ_PRO[,2]),]

Co_RNA<-MATQ_RNA[,-c(1,2,3)]
Co_PRO<-MATQ_PRO[,-c(1,2,3)]

W<-data.frame(RNA=WNN$RNA.weight,PRO=WNN$PRO.weight)
Wt<-t(W)

x=1
X=matrix(data=NA,nrow=nrow(Co_RNA),ncol=ncol(Co_RNA))
for(i in 1:nrow(Co_RNA)){
  for(j in 1:ncol(Co_RNA)){
    X[i,j]<-Co_RNA[i,j]*Wt[1,j]+Co_PRO[i,j]*Wt[2,j]
  }
}

rownames(X)<-MATQ_RNA$SYMBOL
colnames(X)<-colnames(Co_RNA)

group<-read.csv("Group.csv",sep=",",header=T)

Data.log<-log2(X)

group_list <- factor(group[,2],ordered = F)
design <- model.matrix(~0+group_list)
colnames(design) <- levels(group_list)
rownames(design) <- colnames(Data.log)

fit <- lmFit(Data.log,design)

cont.matrix <- makeContrasts(contrasts = paste0(unique(group_list),collapse = "-"),levels = design)
fit2 <- contrasts.fit(fit,cont.matrix)
fit2 <- eBayes(fit2)
tmpOut <- topTable(fit2,coef = 1,n = Inf,adjust.method = "BH",sort.by = "logFC")
limma.na <- na.omit(tmpOut)

limma.na$change <- as.factor(ifelse(limma.na$adj.P.Val<0.05 & abs(limma.na$logFC)>log2(2),ifelse(limma.na$logFC>log2(2),"Down","Up"),"Non"))
limma.na$label <- ifelse(limma.na$adj.P.Val<0.05 & abs(limma.na$logFC)>log2(2),as.character(limma.na$ID),"")

limma.na$label[157]<-""

p<-ggplot(data=limma.na)+
  geom_point(aes(x = -logFC,y = -log10(adj.P.Val),fill=change),color="white",size=5,shape=21,stroke=0.2,alpha=0.65)+
  geom_text_repel(aes(x = -logFC,y = -log10(adj.P.Val),label = label),size = 7,box.padding = unit(0.9,"lines"),point.padding = unit(0.6,"lines"),segment.color = "black",show.legend = FALSE)+
  geom_vline(xintercept = c(-log2(2),log2(2)),lty = 4,col = "black",lwd = 0.8)+
  geom_hline(yintercept = -log10(0.05),lty = 4,col = "black",lwd = 0.8)+
  geom_text(x=1.4,y=4,aes(label="Bpgm"),parse=T,size=6.7)+
  scale_y_continuous(limits=c(0,8))+
  scale_x_continuous(limits=c(-3,2.5))+
  scale_fill_manual(values = c("#0000FF","#A3A3A3","#CC0000"))+
  xlab(expression(Log[2](FoldChange)))+
  ylab(expression(-Log[10](italic(p.adjust))))+
  theme_classic(base_family="ARL")+
  coord_fixed(ratio=0.6875)+
  ggtitle("WNN")+
  theme(panel.background=element_rect(fill="white",colour="black",size=0.25),
        panel.border=element_rect(fill=NA,color="black",size=3,linetype="solid"),
        plot.title=element_text(hjust=0.5,size=40,color="black"),
        axis.line=element_line(colour="black",size=0),
        axis.ticks=element_line(colour="black",size=1.4),
        axis.ticks.length=unit(0.5,'line'),
        axis.title.x=element_text(size=28,face="plain",color="black",vjust=0),
        axis.title.y=element_text(size=28,face="plain",color="black",vjust=2),
        axis.text.x=element_text(size=25,face="plain",color="black",vjust=0),
        axis.text.y=element_text(size=25,face="plain",color="black"),
        legend.text=element_text(size=25,face="plain",colour="black"),
        legend.position=c(0.85,0.88),
        legend.title=element_blank(),
        legend.margin=margin(0.1,0.1,0.1,0.1,'cm'),
        legend.background=element_rect(fill="transparent"),
        legend.key.height=unit(0.5,'cm'),
        legend.key=element_rect(fill="transparent",color="transparent"))+
  guides(fill=guide_legend(nrow=3,byrow=T))

ggsave("Volcano-WNN.png",p,width=7,height=7.7,dpi=1200)


