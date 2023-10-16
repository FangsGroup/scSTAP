#Compared with T&T-seq
windowsFonts(ARL = windowsFont("Arial"))

#set work path
setwd("D:/Research/R/Data")

#Venn-Figure 6E
library(VennDiagram)

#Transcriptome
Data1<-read.csv("Mouse_transcript.csv",sep=",",header=T)

M<-strsplit(Data1$GENE,"\\.")

x=1
X=matrix(data=NA,nrow=dim(Data1)[1],ncol=1)
for(i in 1:dim(Data1)[1]){
  X[i,1]<-M[[i]][1]
}

Data1<-data.frame(Gene=as.character(X[,1]))

gene1<-bitr(Data1$Gene,fromType="ENSEMBL",toType="SYMBOL",OrgDb=org.Mm.eg.db)

Data2<-read.csv("RNA_total.csv",sep=",",header=T)

gene2<-bitr(Data2$gene_id,fromType="ENSEMBL",toType="SYMBOL",OrgDb=org.Mm.eg.db)

venn.diagram(list(Transcriptome_scSTAP=gene2$SYMBOL,Transcriptome_TTseq=gene1$SYMBOL),
             category.names = c("Transcriptome_scSTAP","Transcriptome_TTseq"),
             filename="VennDiagram.jpg",height=5000,width=5000,resolution=500,
             imagetype="tiff",
             main.pos=c(0.5,0.9), main.cex=4, main.fontfamily="ARL",
             fontfamily="ARL", cat.fontfamily="ARL",cat.pos = c(-20,-220),
             col="white",lwd=2,cat.dist = c(0.02,0.04),
             fill=c("#013EDF","#FF4F4F"),ext.text=F,
             alpha=0.3,cex=3,cat.cex=3.2,cat.fontface="plain",margin=0.2)

#Translatome vs Proteome
Data1<-read.csv("Mouse_translatome.csv",sep=",",header=T)

M<-strsplit(Data1$GENE,"\\.")

x=1
X=matrix(data=NA,nrow=dim(Data1)[1],ncol=1)
for(i in 1:dim(Data1)[1]){
  X[i,1]<-M[[i]][1]
}

Data1<-data.frame(Gene=as.character(X[,1]))

gene1<-bitr(Data1$Gene,fromType="ENSEMBL",toType="SYMBOL",OrgDb=org.Mm.eg.db)

Data2<-read.csv("PRO_total.csv",sep=",",header=T)
head(Data2)

gene2<-bitr(Data2$PG.ProteinGroups,fromType="UNIPROT",toType="SYMBOL",OrgDb=org.Mm.eg.db)

venn.diagram(list(Proteome_scSTAP=gene2$SYMBOL,Translatome_TTseq=gene1$SYMBOL),
             category.names = c("Proteome_scSTAP","Translatome_TTseq"),
             filename="VennDiagram-Translatome.jpg",height=5000,width=5000,resolution=500,
             imagetype="tiff",
             main.pos=c(0.5,0.9), main.cex=4, main.fontfamily="ARL",cat.pos = c(-20,-220),
             fontfamily="ARL", cat.fontfamily="ARL",
             col="white",lwd=2,
             fill=c("#013EDF","#FF4F4F"),ext.text=F,
             alpha=0.3,cex=3,cat.cex=3.2,cat.fontface="plain",margin=0.2)

#DEGs
library(AnnotationDbi)
library(clusterProfiler)
library(org.Mm.eg.db)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(psych)
library(corrplot)
library(limma)
library(PerformanceAnalytics)

#Plot-Figure S10A
Data1<-read.csv("Mouse_translatome_FC.csv",sep=",",header=T)

M<-strsplit(Data1$GENE,"\\.")

x=1
X=matrix(data=NA,nrow=dim(Data1)[1],ncol=1)
for(i in 1:dim(Data1)[1]){
  X[i,1]<-M[[i]][1]
}

Data1<-data.frame(Gene=as.character(X[,1]),logFC=Data1$log2.FoldChange.)

gene1<-bitr(Data1$Gene,fromType="ENSEMBL",toType="SYMBOL",OrgDb=org.Mm.eg.db)

TRA<-merge(Data1,gene1,by.x="Gene",by.y="ENSEMBL")

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
limma.na <- transform(limma.na,SYMBOL=rownames(limma.na))

matrix<-merge(TRA,limma.na,by.x="SYMBOL",by.y="SYMBOL")

matrix$change <- as.factor(ifelse(abs(matrix$logFC.x)>1 | abs(matrix$logFC.y)>1,
                                  ifelse(abs(matrix$logFC.x)>1,
                                         ifelse(abs(matrix$logFC.y)>1,
                                                ifelse(matrix$logFC.x*matrix$logFC.y>1,"Class I","Class II"),"Class III"),"Class IV"),"Class V"))

matrix$label <- ifelse(abs(matrix$logFC.x)>1 & abs(matrix$logFC.y)>1,as.character(matrix$SYMBOL),"")

p<-ggplot(data=matrix)+
  geom_point(aes(x = -logFC.y,y = logFC.x,fill=change),color="white",size=5,shape=21,stroke=0.2,alpha=0.65)+
  geom_text_repel(aes(x = -logFC.y,y = logFC.x,label = label),size = 6,force=2,box.padding = unit(0.7,"lines"),
                  point.padding = unit(0.3,"lines"),segment.color = "black",color="black",show.legend = FALSE)+
  geom_vline(xintercept = c(-log2(2),log2(2)),lty = 4,col = "black",lwd = 0.8)+
  geom_hline(yintercept = c(-log2(2),log2(2)),lty = 4,col = "black",lwd = 0.8)+
  scale_y_continuous(limits=c(-12.5,7.5))+
  scale_x_continuous(limits=c(-4,4))+
  scale_fill_manual(values = c("#FA7F6F","#82B0D2","#FFBE7A","#8ECFC9","#E7DAD2"))+
  xlab(expression(Proteome_log[2](FoldChange)))+
  ylab(expression(Translatome_log[2](FoldChange)))+
  theme_classic(base_family="ARL")+
  theme(panel.background=element_rect(fill="white",colour="black",size=0.25),
        panel.border=element_rect(fill=NA,color="black",size=2.4,linetype="solid"),
        axis.line=element_line(colour="black",size=0),
        axis.ticks=element_line(colour="black",size=1),
        axis.ticks.length=unit(0.5,'line'),
        axis.title.x=element_text(size=24,face="plain",color="black",vjust=0),
        axis.title.y=element_text(size=24,face="plain",color="black",vjust=1),
        axis.text.x=element_text(size=23,face="plain",color="black",vjust=0),
        axis.text.y=element_text(size=23,face="plain",color="black"),
        legend.text=element_text(size=20,face="plain",colour="black"),
        legend.position=c(0.82,0.22),
        legend.title=element_blank(),
        legend.margin=margin(0.1,0.1,0.1,0.1,'cm'),
        legend.background=element_rect(fill="transparent"),
        legend.key.height=unit(0.5,'cm'),
        legend.key=element_rect(fill="transparent",color="transparent"))+
  guides(fill=guide_legend(nrow=5,byrow=T))

ggsave("Volcano_T&T.png",p,width=6,height=6,dpi=600)

#Bar-Figure S10A
bar<-data.frame(Class_I = dim(matrix[matrix$change=="Class I",])[1],
                Class_II = dim(matrix[matrix$change=="Class II",])[1],
                Class_III = dim(matrix[matrix$change=="Class III",])[1],
                Class_IV = dim(matrix[matrix$change=="Class IV",])[1],
                Class_V = dim(matrix[matrix$change=="Class V",])[1])
bar<-t(bar)
colnames(bar)<-c("Gene_Number")
bar<-transform(bar, Group = c("Class I","Class II","Class III","Class IV","Class V"))

p<-ggplot()+
  geom_col(data=bar,aes(x=Group,y=Gene_Number,fill=Group),width=0.6,color="white",show.legend=F)+
  geom_text(data=bar,aes(x=Group,y=(Gene_Number+25),label=Gene_Number),size=5,position = position_dodge(width = 0.6))+
  labs(y='Gene number')+
  scale_y_continuous(position = "right",limits=c(0,900))+
  scale_fill_manual(values = c("#FA7F6F","#82B0D2","#FFBE7A","#8ECFC9","#E7DAD2"),labels= c("Class I","Class II","Class III","Class IV","Class V"))+
  theme(panel.background=element_rect(fill="white",colour="black",size=0.25),
        panel.border=element_rect(fill=NA,color="black",size=2.2,linetype="solid"),
        axis.line=element_line(colour="black",size=0),
        axis.ticks=element_line(colour="black",size=1),
        axis.ticks.length=unit(0.5,'line'),
        axis.title.y.right=element_text(size=24,face="plain",color="black",vjust=1.5),
        axis.title.x=element_blank(),
        axis.text.y.right=element_text(size=23,face="plain",color="black",hjust=-0.2),
        axis.text.x=element_text(size=20,face="plain",color="black",angle=90,vjust=0.4))

ggsave("Bar-T&T.png",p,width=3.8,height=6,dpi=600)

#GO-Figure S10B
GO<-matrix[matrix$change=="Class II",]

gene<-bitr(GO$SYMBOL,fromType="SYMBOL",toType="ENTREZID",OrgDb=org.Mm.eg.db)

All.params<-enrichGO(gene=gene$ENTREZID,OrgDb=org.Mm.eg.db,ont="All",pAdjustMethod="BH",pvalueCutoff=0.05,qvalueCutoff=0.05,readable=T)
All.list<-setReadable(All.params,org.Mm.eg.db,keyType="ENTREZID")
go<-as.data.frame(All.params)

mixedToFloat <- function(x){
  x <- sapply(x, as.character)
  is.integer  <- grepl("^-?\\d+$", x)
  is.fraction <- grepl("^-?\\d+\\/\\d+$", x)
  is.float <- grepl("^-?\\d+\\.\\d+$", x)
  is.mixed    <- grepl("^-?\\d+ \\d+\\/\\d+$", x)
  stopifnot(all(is.integer | is.fraction | is.float | is.mixed))
  
  numbers <- strsplit(x, "[ /]")
  
  ifelse(is.integer,  as.numeric(sapply(numbers, `[`, 1)),
         ifelse(is.float,    as.numeric(sapply(numbers, `[`, 1)),
                ifelse(is.fraction, as.numeric(sapply(numbers, `[`, 1)) /
                         as.numeric(sapply(numbers, `[`, 2)),
                       as.numeric(sapply(numbers, `[`, 1)) +
                         as.numeric(sapply(numbers, `[`, 2)) /
                         as.numeric(sapply(numbers, `[`, 3)))))
  
}


go$Description<-capitalize(go$Description)
go$Description<-factor(go$Description,levels = rev(go$Description))

go$number <- factor(rev(1:nrow(go)))

go$GeneRatio <- mixedToFloat(go$GeneRatio)

labels <- go$Description
names(labels) = rev(1:nrow(go))

p <- ggplot(data=go, aes(x=Description, y=GeneRatio)) +
  geom_point(aes(color=p.adjust,size=Count)) + coord_flip() +
  theme_test() +
  scale_color_gradientn(colours=c("#E02020","#FDE59A","#07388F"))+
  scale_x_discrete(labels=labels) +
  scale_y_continuous(limits = c(0.07,0.18))+
  xlab("GO terms") +
  ylab("Gene ratio") +
  ggtitle("GO terms entichment of class II genes")+
  theme(panel.background=element_rect(fill="#FFFFFF",colour="black",size=0.25),
        panel.border=element_rect(fill=NA,color="black",size=2,linetype="solid"),
        plot.title=element_text(colour="black",size=20,hjust=0.9),
        axis.ticks.length=unit(0.4,'line'),
        axis.ticks=element_line(colour="black",size=0.8),
        axis.line=element_line(colour="black",size=0),
        axis.title.x=element_text(size=18,face="plain",color="black",vjust=0),
        axis.title.y=element_text(size=18,face="plain",color="black",vjust=2.0),
        axis.text.x=element_text(size=14,face="plain",color="black",vjust=0),
        axis.text.y=element_text(size=16,face="plain",color="black",lineheight=0.6),
        legend.title=element_text(size=14,face="plain",colour="black"),
        legend.text=element_text(size=14,face="plain",colour="black"))

ggsave("go_point-T&T.png",p,width=9,height=3.6)

#Individual Bar-Figure S10C
#Bsg
Bar<-Matrixx[Matrixx$SYMBOL=="Bsg",]
Bar
Bar<-data.frame(TTseq_GV1_Tr=Bar$TTseq_GV1_Tr/119.7*100,TTseq_GV2_Tr=Bar$TTseq_GV2_Tr/119.7*100,TTseq_MII1_Tr=Bar$TTseq_MII1_Tr/119.7*100,TTseq_MII2_Tr=Bar$TTseq_MII2_Tr/119.7*100)
Bar<-t(Bar)
colnames(Bar)<-c("Gene")
Bar<-transform(Bar,Sample=c("TTseq_GV1","TTseq_GV2","TTseq_MII1","TTseq_MII2"),
               Group=c(rep("GV",2),rep("MII",2)))
Bar$Sample<-factor(Bar$Sample,level=c("TTseq_GV1","TTseq_GV2","TTseq_MII1","TTseq_MII2"))

p<-ggplot()+
  geom_col(data=Bar,aes(x=Sample,y=Gene,fill=Group),position=position_dodge(width=0.8),width=0.8,color="white",show.legend=F)+
  labs(y='Bsg')+
  coord_cartesian(xlim=c(0,5),ylim=c(1,110),expand=F)+
  scale_fill_manual(values = c("#F4CC66","#758837"),labels= c("GV","MII"))+
  theme(panel.background=element_rect(fill="white",colour="black",size=0.25),
        panel.border=element_rect(fill=NA,color="black",size=2,linetype="solid"),
        axis.line=element_line(colour="#595959",size=0),
        axis.ticks.y=element_line(colour="black",size=0.6),
        axis.ticks.length.y=unit(0.5,'line'),
        axis.ticks.x=element_line(colour="#595959",size=0),
        axis.ticks.length.x=unit(0,'line'),
        axis.title.y=element_text(size=28,face="plain",color="black",vjust=1.5),
        axis.title.x=element_blank(),
        axis.text.y=element_text(size=24,face="plain",color="black",hjust=1),
        axis.text.x=element_blank())

ggsave("Bar.png",p,width=2.4,height=1.8,dpi=600)

Bar<-Matrixx[Matrixx$SYMBOL=="Bsg",]
Bar
Bar<-data.frame(TTseq_GV1_Tl=Bar$TTseq_GV1_Tl/57.6*100,TTseq_GV2_Tl=Bar$TTseq_GV2_Tl/57.6*100,TTseq_MII1_Tl=Bar$TTseq_MII1_Tl/57.6*100,TTseq_MII2_Tl=Bar$TTseq_MII2_Tl/57.6*100)
Bar<-t(Bar)
colnames(Bar)<-c("Gene")
Bar<-transform(Bar,Sample=c("TTseq_GV1","TTseq_GV2","TTseq_MII1","TTseq_MII2"),
               Group=c(rep("GV",2),rep("MII",2)))
Bar$Sample<-factor(Bar$Sample,level=c("TTseq_GV1","TTseq_GV2","TTseq_MII1","TTseq_MII2"))

p<-ggplot()+
  geom_col(data=Bar,aes(x=Sample,y=Gene,fill=Group),position=position_dodge(width=0.8),width=0.8,color="white",show.legend=F)+
  coord_cartesian(xlim=c(0,5),ylim=c(1,110),expand=F)+
  scale_fill_manual(values = c("#F4CC66","#758837"),labels= c("GV","MII"))+
  theme(panel.background=element_rect(fill="white",colour="black",size=0.25),
        panel.border=element_rect(fill=NA,color="black",size=2,linetype="solid"),
        axis.line=element_line(colour="#595959",size=0),
        axis.ticks.y=element_line(colour="black",size=0),
        axis.ticks.length.y=unit(0,'line'),
        axis.ticks.x=element_line(colour="#595959",size=0),
        axis.ticks.length.x=unit(0,'line'),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.y=element_blank(),
        axis.text.x=element_blank())

ggsave("Bar.png",p,width=1.4,height=1.8,dpi=600)

Bar<-Matrixx[Matrixx$SYMBOL=="Bsg",]
Bar<-data.frame(GV1=Bar$GV1.x/245.9*100,GV2=Bar$GV2.x/245.9*100,GV3=Bar$GV3.x/245.9*100,GV4=Bar$GV4.x/245.9*100,GV5=Bar$GV5.x/245.9*100,GV6=Bar$GV6.x/245.9*100,GV7=Bar$GV7.x/245.9*100,
                GV8=Bar$GV8.x/245.9*100,GV9=Bar$GV9.x/245.9*100,MII1=Bar$MII1.x/245.9*100,MII2=Bar$MII2.x/245.9*100,MII3=Bar$MII3.x/245.9*100,MII4=Bar$MII4.x/245.9*100,MII5=Bar$MII5.x/245.9*100,MII6=Bar$MII6.x/245.9*100)
Bar<-t(Bar)
colnames(Bar)<-c("Gene")
Bar<-transform(Bar,Sample=c("scSTAP_GV1","scSTAP_GV2","scSTAP_GV3","scSTAP_GV4","scSTAP_GV5","scSTAP_GV6","scSTAP_GV7","scSTAP_GV8","scSTAP_GV9",
                            "scSTAP_MII1","scSTAP_MII2","scSTAP_MII3","scSTAP_MII4","scSTAP_MII5","scSTAP_MII6"),
               Group=c(rep("GV",9),rep("MII",6)))
Bar$Sample<-factor(Bar$Sample,level=c("scSTAP_GV1","scSTAP_GV2","scSTAP_GV3","scSTAP_GV4","scSTAP_GV5","scSTAP_GV6","scSTAP_GV7","scSTAP_GV8","scSTAP_GV9",
                                      "scSTAP_MII1","scSTAP_MII2","scSTAP_MII3","scSTAP_MII4","scSTAP_MII5","scSTAP_MII6"))

p<-ggplot()+
  geom_col(data=Bar,aes(x=Sample,y=Gene,fill=Group),position=position_dodge(width=0.8),width=0.8,color="white",show.legend=F)+
  coord_cartesian(xlim=c(0,16),ylim=c(1,110),expand=F)+
  scale_fill_manual(values = c("#CE5F56","#00789E"),labels= c("GV","MII"))+
  theme(panel.background=element_rect(fill="white",colour="black",size=0.25),
        panel.border=element_rect(fill=NA,color="black",size=2,linetype="solid"),
        axis.line=element_line(colour="#595959",size=0),
        axis.ticks.y=element_line(colour="black",size=0.0),
        axis.ticks.length.y=unit(0.0,'line'),
        axis.ticks.x=element_line(colour="#595959",size=0),
        axis.ticks.length.x=unit(0,'line'),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.y=element_blank(),
        axis.text.x=element_blank())

ggsave("Bar.png",p,width=4.2,height=1.8,dpi=600)

Bar<-Matrixx[Matrixx$SYMBOL=="Bsg",]
Bar<-data.frame(GV1=Bar$GV1.y/743.48*100,GV2=Bar$GV2.y/743.48*100,GV3=Bar$GV3.y/743.48*100,GV4=Bar$GV4.y/743.48*100,GV5=Bar$GV5.y/743.48*100,GV6=Bar$GV6.y/743.48*100,GV7=Bar$GV7.y/743.48*100,
                GV8=Bar$GV8.y/743.48*100,GV9=Bar$GV9.y/743.48*100,MII1=Bar$MII1.y/743.48*100,MII2=Bar$MII2.y/743.48*100,MII3=Bar$MII3.y/743.48*100,MII4=Bar$MII4.y/743.48*100,MII5=Bar$MII5.y/743.48*100,MII6=Bar$MII6.y/743.48*100)
Bar<-t(Bar)
colnames(Bar)<-c("Gene")
Bar<-transform(Bar,Sample=c("scSTAP_GV1","scSTAP_GV2","scSTAP_GV3","scSTAP_GV4","scSTAP_GV5","scSTAP_GV6","scSTAP_GV7","scSTAP_GV8","scSTAP_GV9",
                            "scSTAP_MII1","scSTAP_MII2","scSTAP_MII3","scSTAP_MII4","scSTAP_MII5","scSTAP_MII6"),
               Group=c(rep("GV",9),rep("MII",6)))
Bar$Sample<-factor(Bar$Sample,level=c("scSTAP_GV1","scSTAP_GV2","scSTAP_GV3","scSTAP_GV4","scSTAP_GV5","scSTAP_GV6","scSTAP_GV7","scSTAP_GV8","scSTAP_GV9",
                                      "scSTAP_MII1","scSTAP_MII2","scSTAP_MII3","scSTAP_MII4","scSTAP_MII5","scSTAP_MII6"))

p<-ggplot()+
  geom_col(data=Bar,aes(x=Sample,y=Gene,fill=Group),position=position_dodge(width=0.8),width=0.8,color="white",show.legend=F)+
  coord_cartesian(xlim=c(0,16),ylim=c(1,110),expand=F)+
  scale_fill_manual(values = c("#CE5F56","#00789E"),labels= c("GV","MII"))+
  theme(panel.background=element_rect(fill="white",colour="black",size=0.25),
        panel.border=element_rect(fill=NA,color="black",size=2,linetype="solid"),
        axis.line=element_line(colour="#595959",size=0),
        axis.ticks.y=element_line(colour="black",size=0.0),
        axis.ticks.length.y=unit(0.0,'line'),
        axis.ticks.x=element_line(colour="#595959",size=0),
        axis.ticks.length.x=unit(0,'line'),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.y=element_blank(),
        axis.text.x=element_blank())

ggsave("Bar.png",p,width=4.2,height=1.8,dpi=600)

#Cycs
Bar<-Matrixx[Matrixx$SYMBOL=="Cycs",]
Bar
Bar<-data.frame(TTseq_GV1_Tr=Bar$TTseq_GV1_Tr/181*100,TTseq_GV2_Tr=Bar$TTseq_GV2_Tr/181*100,TTseq_MII1_Tr=Bar$TTseq_MII1_Tr/181*100,TTseq_MII2_Tr=Bar$TTseq_MII2_Tr/181*100)
Bar<-t(Bar)
colnames(Bar)<-c("Gene")
Bar<-transform(Bar,Sample=c("TTseq_GV1","TTseq_GV2","TTseq_MII1","TTseq_MII2"),
               Group=c(rep("GV",2),rep("MII",2)))
Bar$Sample<-factor(Bar$Sample,level=c("TTseq_GV1","TTseq_GV2","TTseq_MII1","TTseq_MII2"))

p<-ggplot()+
  geom_col(data=Bar,aes(x=Sample,y=Gene,fill=Group),position=position_dodge(width=0.8),width=0.8,color="white",show.legend=F)+
  labs(y='Cycs')+
  coord_cartesian(xlim=c(0,5),ylim=c(1,110),expand=F)+
  scale_fill_manual(values = c("#F4CC66","#758837"),labels= c("GV","MII"))+
  theme(panel.background=element_rect(fill="white",colour="black",size=0.25),
        panel.border=element_rect(fill=NA,color="black",size=2,linetype="solid"),
        axis.line=element_line(colour="#595959",size=0),
        axis.ticks.y=element_line(colour="black",size=0.6),
        axis.ticks.length.y=unit(0.5,'line'),
        axis.ticks.x=element_line(colour="#595959",size=0),
        axis.ticks.length.x=unit(0,'line'),
        axis.title.y=element_text(size=28,face="plain",color="black",vjust=1.5),
        axis.title.x=element_blank(),
        axis.text.y=element_text(size=24,face="plain",color="black",hjust=1),
        axis.text.x=element_blank())

ggsave("Bar.png",p,width=2.4,height=1.8,dpi=600)

Bar<-Matrixx[Matrixx$SYMBOL=="Cycs",]
Bar
Bar<-data.frame(TTseq_GV1_Tl=Bar$TTseq_GV1_Tl/336.3*100,TTseq_GV2_Tl=Bar$TTseq_GV2_Tl/336.3*100,TTseq_MII1_Tl=Bar$TTseq_MII1_Tl/336.3*100,TTseq_MII2_Tl=Bar$TTseq_MII2_Tl/336.3*100)
Bar<-t(Bar)
colnames(Bar)<-c("Gene")
Bar<-transform(Bar,Sample=c("TTseq_GV1","TTseq_GV2","TTseq_MII1","TTseq_MII2"),
               Group=c(rep("GV",2),rep("MII",2)))
Bar$Sample<-factor(Bar$Sample,level=c("TTseq_GV1","TTseq_GV2","TTseq_MII1","TTseq_MII2"))

p<-ggplot()+
  geom_col(data=Bar,aes(x=Sample,y=Gene,fill=Group),position=position_dodge(width=0.8),width=0.8,color="white",show.legend=F)+
  coord_cartesian(xlim=c(0,5),ylim=c(1,110),expand=F)+
  scale_fill_manual(values = c("#F4CC66","#758837"),labels= c("GV","MII"))+
  theme(panel.background=element_rect(fill="white",colour="black",size=0.25),
        panel.border=element_rect(fill=NA,color="black",size=2,linetype="solid"),
        axis.line=element_line(colour="#595959",size=0),
        axis.ticks.y=element_line(colour="black",size=0),
        axis.ticks.length.y=unit(0,'line'),
        axis.ticks.x=element_line(colour="#595959",size=0),
        axis.ticks.length.x=unit(0,'line'),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.y=element_blank(),
        axis.text.x=element_blank())

ggsave("Bar.png",p,width=1.4,height=1.8,dpi=600)

Bar<-Matrixx[Matrixx$SYMBOL=="Cycs",]
Bar<-data.frame(GV1=Bar$GV1.x/58.9*100,GV2=Bar$GV2.x/58.9*100,GV3=Bar$GV3.x/58.9*100,GV4=Bar$GV4.x/58.9*100,GV5=Bar$GV5.x/58.9*100,GV6=Bar$GV6.x/58.9*100,GV7=Bar$GV7.x/58.9*100,
                GV8=Bar$GV8.x/58.9*100,GV9=Bar$GV9.x/58.9*100,MII1=Bar$MII1.x/58.9*100,MII2=Bar$MII2.x/58.9*100,MII3=Bar$MII3.x/58.9*100,MII4=Bar$MII4.x/58.9*100,MII5=Bar$MII5.x/58.9*100,MII6=Bar$MII6.x/58.9*100)
Bar<-t(Bar)
colnames(Bar)<-c("Gene")
Bar<-transform(Bar,Sample=c("scSTAP_GV1","scSTAP_GV2","scSTAP_GV3","scSTAP_GV4","scSTAP_GV5","scSTAP_GV6","scSTAP_GV7","scSTAP_GV8","scSTAP_GV9",
                            "scSTAP_MII1","scSTAP_MII2","scSTAP_MII3","scSTAP_MII4","scSTAP_MII5","scSTAP_MII6"),
               Group=c(rep("GV",9),rep("MII",6)))
Bar$Sample<-factor(Bar$Sample,level=c("scSTAP_GV1","scSTAP_GV2","scSTAP_GV3","scSTAP_GV4","scSTAP_GV5","scSTAP_GV6","scSTAP_GV7","scSTAP_GV8","scSTAP_GV9",
                                      "scSTAP_MII1","scSTAP_MII2","scSTAP_MII3","scSTAP_MII4","scSTAP_MII5","scSTAP_MII6"))

p<-ggplot()+
  geom_col(data=Bar,aes(x=Sample,y=Gene,fill=Group),position=position_dodge(width=0.8),width=0.8,color="white",show.legend=F)+
  coord_cartesian(xlim=c(0,16),ylim=c(1,110),expand=F)+
  scale_fill_manual(values = c("#CE5F56","#00789E"),labels= c("GV","MII"))+
  theme(panel.background=element_rect(fill="white",colour="black",size=0.25),
        panel.border=element_rect(fill=NA,color="black",size=2,linetype="solid"),
        axis.line=element_line(colour="#595959",size=0),
        axis.ticks.y=element_line(colour="black",size=0.0),
        axis.ticks.length.y=unit(0.0,'line'),
        axis.ticks.x=element_line(colour="#595959",size=0),
        axis.ticks.length.x=unit(0,'line'),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.y=element_blank(),
        axis.text.x=element_blank())

ggsave("Bar.png",p,width=4.2,height=1.8,dpi=600)

Bar<-Matrixx[Matrixx$SYMBOL=="Cycs",]
Bar<-data.frame(GV1=Bar$GV1.y/512.5*100,GV2=Bar$GV2.y/512.5*100,GV3=Bar$GV3.y/512.5*100,GV4=Bar$GV4.y/512.5*100,GV5=Bar$GV5.y/512.5*100,GV6=Bar$GV6.y/512.5*100,GV7=Bar$GV7.y/512.5*100,
                GV8=Bar$GV8.y/512.5*100,GV9=Bar$GV9.y/512.5*100,MII1=Bar$MII1.y/512.5*100,MII2=Bar$MII2.y/512.5*100,MII3=Bar$MII3.y/512.5*100,MII4=Bar$MII4.y/512.5*100,MII5=Bar$MII5.y/512.5*100,MII6=Bar$MII6.y/512.5*100)
Bar<-t(Bar)
colnames(Bar)<-c("Gene")
Bar<-transform(Bar,Sample=c("scSTAP_GV1","scSTAP_GV2","scSTAP_GV3","scSTAP_GV4","scSTAP_GV5","scSTAP_GV6","scSTAP_GV7","scSTAP_GV8","scSTAP_GV9",
                            "scSTAP_MII1","scSTAP_MII2","scSTAP_MII3","scSTAP_MII4","scSTAP_MII5","scSTAP_MII6"),
               Group=c(rep("GV",9),rep("MII",6)))
Bar$Sample<-factor(Bar$Sample,level=c("scSTAP_GV1","scSTAP_GV2","scSTAP_GV3","scSTAP_GV4","scSTAP_GV5","scSTAP_GV6","scSTAP_GV7","scSTAP_GV8","scSTAP_GV9",
                                      "scSTAP_MII1","scSTAP_MII2","scSTAP_MII3","scSTAP_MII4","scSTAP_MII5","scSTAP_MII6"))

p<-ggplot()+
  geom_col(data=Bar,aes(x=Sample,y=Gene,fill=Group),position=position_dodge(width=0.8),width=0.8,color="white",show.legend=F)+
  coord_cartesian(xlim=c(0,16),ylim=c(1,110),expand=F)+
  scale_fill_manual(values = c("#CE5F56","#00789E"),labels= c("GV","MII"))+
  theme(panel.background=element_rect(fill="white",colour="black",size=0.25),
        panel.border=element_rect(fill=NA,color="black",size=2,linetype="solid"),
        axis.line=element_line(colour="#595959",size=0),
        axis.ticks.y=element_line(colour="black",size=0.0),
        axis.ticks.length.y=unit(0.0,'line'),
        axis.ticks.x=element_line(colour="#595959",size=0),
        axis.ticks.length.x=unit(0,'line'),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.y=element_blank(),
        axis.text.x=element_blank())

ggsave("Bar.png",p,width=4.2,height=1.8,dpi=600)

#Eif1ax
Bar<-Matrixx[Matrixx$SYMBOL=="Eif1ax",]
Bar
Bar<-data.frame(TTseq_GV1_Tr=Bar$TTseq_GV1_Tr/165*100,TTseq_GV2_Tr=Bar$TTseq_GV2_Tr/165*100,TTseq_MII1_Tr=Bar$TTseq_MII1_Tr/165*100,TTseq_MII2_Tr=Bar$TTseq_MII2_Tr/165*100)
Bar<-t(Bar)
colnames(Bar)<-c("Gene")
Bar<-transform(Bar,Sample=c("TTseq_GV1","TTseq_GV2","TTseq_MII1","TTseq_MII2"),
               Group=c(rep("GV",2),rep("MII",2)))
Bar$Sample<-factor(Bar$Sample,level=c("TTseq_GV1","TTseq_GV2","TTseq_MII1","TTseq_MII2"))

p<-ggplot()+
  geom_col(data=Bar,aes(x=Sample,y=Gene,fill=Group),position=position_dodge(width=0.8),width=0.8,color="white",show.legend=F)+
  labs(y='Eif1ax')+
  coord_cartesian(xlim=c(0,5),ylim=c(1,110),expand=F)+
  scale_fill_manual(values = c("#F4CC66","#758837"),labels= c("GV","MII"))+
  theme(panel.background=element_rect(fill="white",colour="black",size=0.25),
        panel.border=element_rect(fill=NA,color="black",size=2,linetype="solid"),
        axis.line=element_line(colour="#595959",size=0),
        axis.ticks.y=element_line(colour="black",size=0.6),
        axis.ticks.length.y=unit(0.5,'line'),
        axis.ticks.x=element_line(colour="black",size=0.0),
        axis.ticks.length.x=unit(0.0,'line'),
        axis.title.y=element_text(size=28,face="plain",color="black",vjust=1.5),
        axis.title.x=element_blank(),
        axis.text.y=element_text(size=24,face="plain",color="black",hjust=1),
        axis.text.x=element_blank())

ggsave("Bar.png",p,width=2.4,height=1.8,dpi=600)

Bar<-Matrixx[Matrixx$SYMBOL=="Eif1ax",]
Bar
Bar<-data.frame(TTseq_GV1_Tl=Bar$TTseq_GV1_Tl/178.3*100,TTseq_GV2_Tl=Bar$TTseq_GV2_Tl/178.3*100,TTseq_MII1_Tl=Bar$TTseq_MII1_Tl/178.3*100,TTseq_MII2_Tl=Bar$TTseq_MII2_Tl/178.3*100)
Bar<-t(Bar)
colnames(Bar)<-c("Gene")
Bar<-transform(Bar,Sample=c("TTseq_GV1","TTseq_GV2","TTseq_MII1","TTseq_MII2"),
               Group=c(rep("GV",2),rep("MII",2)))
Bar$Sample<-factor(Bar$Sample,level=c("TTseq_GV1","TTseq_GV2","TTseq_MII1","TTseq_MII2"))

p<-ggplot()+
  geom_col(data=Bar,aes(x=Sample,y=Gene,fill=Group),position=position_dodge(width=0.8),width=0.8,color="white",show.legend=F)+
  coord_cartesian(xlim=c(0,5),ylim=c(1,110),expand=F)+
  scale_fill_manual(values = c("#F4CC66","#758837"),labels= c("GV","MII"))+
  theme(panel.background=element_rect(fill="white",colour="black",size=0.25),
        panel.border=element_rect(fill=NA,color="black",size=2,linetype="solid"),
        axis.line=element_line(colour="#595959",size=0),
        axis.ticks.y=element_line(colour="black",size=0),
        axis.ticks.length.y=unit(0,'line'),
        axis.ticks.x=element_line(colour="black",size=0.0),
        axis.ticks.length.x=unit(0.0,'line'),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.y=element_blank(),
        axis.text.x=element_blank())

ggsave("Bar.png",p,width=1.4,height=1.8,dpi=600)

Bar<-Matrixx[Matrixx$SYMBOL=="Eif1ax",]
Bar<-data.frame(GV1=Bar$GV1.x/110.26*100,GV2=Bar$GV2.x/110.26*100,GV3=Bar$GV3.x/110.26*100,GV4=Bar$GV4.x/110.26*100,GV5=Bar$GV5.x/110.26*100,GV6=Bar$GV6.x/110.26*100,GV7=Bar$GV7.x/110.26*100,
                GV8=Bar$GV8.x/110.26*100,GV9=Bar$GV9.x/110.26*100,MII1=Bar$MII1.x/110.26*100,MII2=Bar$MII2.x/110.26*100,MII3=Bar$MII3.x/110.26*100,MII4=Bar$MII4.x/110.26*100,MII5=Bar$MII5.x/110.26*100,MII6=Bar$MII6.x/110.26*100)
Bar<-t(Bar)
colnames(Bar)<-c("Gene")
Bar<-transform(Bar,Sample=c("scSTAP_GV1","scSTAP_GV2","scSTAP_GV3","scSTAP_GV4","scSTAP_GV5","scSTAP_GV6","scSTAP_GV7","scSTAP_GV8","scSTAP_GV9",
                            "scSTAP_MII1","scSTAP_MII2","scSTAP_MII3","scSTAP_MII4","scSTAP_MII5","scSTAP_MII6"),
               Group=c(rep("GV",9),rep("MII",6)))
Bar$Sample<-factor(Bar$Sample,level=c("scSTAP_GV1","scSTAP_GV2","scSTAP_GV3","scSTAP_GV4","scSTAP_GV5","scSTAP_GV6","scSTAP_GV7","scSTAP_GV8","scSTAP_GV9",
                                      "scSTAP_MII1","scSTAP_MII2","scSTAP_MII3","scSTAP_MII4","scSTAP_MII5","scSTAP_MII6"))

p<-ggplot()+
  geom_col(data=Bar,aes(x=Sample,y=Gene,fill=Group),position=position_dodge(width=0.8),width=0.8,color="white",show.legend=F)+
  coord_cartesian(xlim=c(0,16),ylim=c(1,110),expand=F)+
  scale_fill_manual(values = c("#CE5F56","#00789E"),labels= c("GV","MII"))+
  theme(panel.background=element_rect(fill="white",colour="black",size=0.25),
        panel.border=element_rect(fill=NA,color="black",size=2,linetype="solid"),
        axis.line=element_line(colour="#595959",size=0),
        axis.ticks.y=element_line(colour="black",size=0.0),
        axis.ticks.length.y=unit(0.0,'line'),
        axis.ticks.x=element_line(colour="black",size=0.0),
        axis.ticks.length.x=unit(0.0,'line'),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.y=element_blank(),
        axis.text.x=element_blank())

ggsave("Bar.png",p,width=4.2,height=1.8,dpi=600)

Bar<-Matrixx[Matrixx$SYMBOL=="Eif1ax",]
Bar<-data.frame(GV1=Bar$GV1.y/368.31*100,GV2=Bar$GV2.y/368.31*100,GV3=Bar$GV3.y/368.31*100,GV4=Bar$GV4.y/368.31*100,GV5=Bar$GV5.y/368.31*100,GV6=Bar$GV6.y/368.31*100,GV7=Bar$GV7.y/368.31*100,
                GV8=Bar$GV8.y/368.31*100,GV9=Bar$GV9.y/368.31*100,MII1=Bar$MII1.y/368.31*100,MII2=Bar$MII2.y/368.31*100,MII3=Bar$MII3.y/368.31*100,MII4=Bar$MII4.y/368.31*100,MII5=Bar$MII5.y/368.31*100,MII6=Bar$MII6.y/368.31*100)
Bar<-t(Bar)
colnames(Bar)<-c("Gene")
Bar<-transform(Bar,Sample=c("scSTAP_GV1","scSTAP_GV2","scSTAP_GV3","scSTAP_GV4","scSTAP_GV5","scSTAP_GV6","scSTAP_GV7","scSTAP_GV8","scSTAP_GV9",
                            "scSTAP_MII1","scSTAP_MII2","scSTAP_MII3","scSTAP_MII4","scSTAP_MII5","scSTAP_MII6"),
               Group=c(rep("GV",9),rep("MII",6)))
Bar$Sample<-factor(Bar$Sample,level=c("scSTAP_GV1","scSTAP_GV2","scSTAP_GV3","scSTAP_GV4","scSTAP_GV5","scSTAP_GV6","scSTAP_GV7","scSTAP_GV8","scSTAP_GV9",
                                      "scSTAP_MII1","scSTAP_MII2","scSTAP_MII3","scSTAP_MII4","scSTAP_MII5","scSTAP_MII6"))

p<-ggplot()+
  geom_col(data=Bar,aes(x=Sample,y=Gene,fill=Group),position=position_dodge(width=0.8),width=0.8,color="white",show.legend=F)+
  coord_cartesian(xlim=c(0,16),ylim=c(1,110),expand=F)+
  scale_fill_manual(values = c("#CE5F56","#00789E"),labels= c("GV","MII"))+
  theme(panel.background=element_rect(fill="white",colour="black",size=0.25),
        panel.border=element_rect(fill=NA,color="black",size=2,linetype="solid"),
        axis.line=element_line(colour="#595959",size=0),
        axis.ticks.y=element_line(colour="black",size=0.0),
        axis.ticks.length.y=unit(0.0,'line'),
        axis.ticks.x=element_line(colour="black",size=0.0),
        axis.ticks.length.x=unit(0.0,'line'),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.y=element_blank(),
        axis.text.x=element_blank())

ggsave("Bar.png",p,width=4.2,height=1.8,dpi=600)

#Psmc4
Bar<-Matrixx[Matrixx$SYMBOL=="Psmc4",]
Bar
Bar<-data.frame(TTseq_GV1_Tr=Bar$TTseq_GV1_Tr/233.8*100,TTseq_GV2_Tr=Bar$TTseq_GV2_Tr/233.8*100,TTseq_MII1_Tr=Bar$TTseq_MII1_Tr/233.8*100,TTseq_MII2_Tr=Bar$TTseq_MII2_Tr/233.8*100)
Bar<-t(Bar)
colnames(Bar)<-c("Gene")
Bar<-transform(Bar,Sample=c("TTseq_GV1","TTseq_GV2","TTseq_MII1","TTseq_MII2"),
               Group=c(rep("GV",2),rep("MII",2)))
Bar$Sample<-factor(Bar$Sample,level=c("TTseq_GV1","TTseq_GV2","TTseq_MII1","TTseq_MII2"))

p<-ggplot()+
  geom_col(data=Bar,aes(x=Sample,y=Gene,fill=Group),position=position_dodge(width=0.8),width=0.8,color="white",show.legend=F)+
  labs(y='Psmc4')+
  coord_cartesian(xlim=c(0,5),ylim=c(1,110),expand=F)+
  scale_fill_manual(values = c("#F4CC66","#758837"),labels= c("GV","MII"))+
  theme(panel.background=element_rect(fill="white",colour="black",size=0.25),
        panel.border=element_rect(fill=NA,color="black",size=2,linetype="solid"),
        axis.line=element_line(colour="#595959",size=0),
        axis.ticks.y=element_line(colour="black",size=0.6),
        axis.ticks.length.y=unit(0.5,'line'),
        axis.ticks.x=element_line(colour="#595959",size=0),
        axis.ticks.length.x=unit(0,'line'),
        axis.title.y=element_text(size=28,face="plain",color="black",vjust=1.5),
        axis.title.x=element_blank(),
        axis.text.y=element_text(size=24,face="plain",color="black",hjust=1),
        axis.text.x=element_blank())

ggsave("Bar.png",p,width=2.4,height=1.8,dpi=600)

Bar<-Matrixx[Matrixx$SYMBOL=="Psmc4",]
Bar
Bar<-data.frame(TTseq_GV1_Tl=Bar$TTseq_GV1_Tl/234.9*100,TTseq_GV2_Tl=Bar$TTseq_GV2_Tl/234.9*100,TTseq_MII1_Tl=Bar$TTseq_MII1_Tl/234.9*100,TTseq_MII2_Tl=Bar$TTseq_MII2_Tl/234.9*100)
Bar<-t(Bar)
colnames(Bar)<-c("Gene")
Bar<-transform(Bar,Sample=c("TTseq_GV1","TTseq_GV2","TTseq_MII1","TTseq_MII2"),
               Group=c(rep("GV",2),rep("MII",2)))
Bar$Sample<-factor(Bar$Sample,level=c("TTseq_GV1","TTseq_GV2","TTseq_MII1","TTseq_MII2"))

p<-ggplot()+
  geom_col(data=Bar,aes(x=Sample,y=Gene,fill=Group),position=position_dodge(width=0.8),width=0.8,color="white",show.legend=F)+
  coord_cartesian(xlim=c(0,5),ylim=c(1,110),expand=F)+
  scale_fill_manual(values = c("#F4CC66","#758837"),labels= c("GV","MII"))+
  theme(panel.background=element_rect(fill="white",colour="black",size=0.25),
        panel.border=element_rect(fill=NA,color="black",size=2,linetype="solid"),
        axis.line=element_line(colour="#595959",size=0),
        axis.ticks.y=element_line(colour="black",size=0),
        axis.ticks.length.y=unit(0,'line'),
        axis.ticks.x=element_line(colour="#595959",size=0),
        axis.ticks.length.x=unit(0,'line'),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.y=element_blank(),
        axis.text.x=element_blank())

ggsave("Bar.png",p,width=1.4,height=1.8,dpi=600)

Bar<-Matrixx[Matrixx$SYMBOL=="Psmc4",]
Bar<-data.frame(GV1=Bar$GV1.x/80.2*100,GV2=Bar$GV2.x/80.2*100,GV3=Bar$GV3.x/80.2*100,GV4=Bar$GV4.x/80.2*100,GV5=Bar$GV5.x/80.2*100,GV6=Bar$GV6.x/80.2*100,GV7=Bar$GV7.x/80.2*100,
                GV8=Bar$GV8.x/80.2*100,GV9=Bar$GV9.x/80.2*100,MII1=Bar$MII1.x/80.2*100,MII2=Bar$MII2.x/80.2*100,MII3=Bar$MII3.x/80.2*100,MII4=Bar$MII4.x/80.2*100,MII5=Bar$MII5.x/80.2*100,MII6=Bar$MII6.x/80.2*100)
Bar<-t(Bar)
colnames(Bar)<-c("Gene")
Bar<-transform(Bar,Sample=c("scSTAP_GV1","scSTAP_GV2","scSTAP_GV3","scSTAP_GV4","scSTAP_GV5","scSTAP_GV6","scSTAP_GV7","scSTAP_GV8","scSTAP_GV9",
                            "scSTAP_MII1","scSTAP_MII2","scSTAP_MII3","scSTAP_MII4","scSTAP_MII5","scSTAP_MII6"),
               Group=c(rep("GV",9),rep("MII",6)))
Bar$Sample<-factor(Bar$Sample,level=c("scSTAP_GV1","scSTAP_GV2","scSTAP_GV3","scSTAP_GV4","scSTAP_GV5","scSTAP_GV6","scSTAP_GV7","scSTAP_GV8","scSTAP_GV9",
                                      "scSTAP_MII1","scSTAP_MII2","scSTAP_MII3","scSTAP_MII4","scSTAP_MII5","scSTAP_MII6"))

p<-ggplot()+
  geom_col(data=Bar,aes(x=Sample,y=Gene,fill=Group),position=position_dodge(width=0.8),width=0.8,color="white",show.legend=F)+
  coord_cartesian(xlim=c(0,16),ylim=c(1,110),expand=F)+
  scale_fill_manual(values = c("#CE5F56","#00789E"),labels= c("GV","MII"))+
  theme(panel.background=element_rect(fill="white",colour="black",size=0.25),
        panel.border=element_rect(fill=NA,color="black",size=2,linetype="solid"),
        axis.line=element_line(colour="#595959",size=0),
        axis.ticks.y=element_line(colour="black",size=0.0),
        axis.ticks.length.y=unit(0.0,'line'),
        axis.ticks.x=element_line(colour="#595959",size=0),
        axis.ticks.length.x=unit(0,'line'),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.y=element_blank(),
        axis.text.x=element_blank())

ggsave("Bar.png",p,width=4.2,height=1.8,dpi=600)

Bar<-Matrixx[Matrixx$SYMBOL=="Psmc4",]
Bar<-data.frame(GV1=Bar$GV1.y/111.5*100,GV2=Bar$GV2.y/111.5*100,GV3=Bar$GV3.y/111.5*100,GV4=Bar$GV4.y/111.5*100,GV5=Bar$GV5.y/111.5*100,GV6=Bar$GV6.y/111.5*100,GV7=Bar$GV7.y/111.5*100,
                GV8=Bar$GV8.y/111.5*100,GV9=Bar$GV9.y/111.5*100,MII1=Bar$MII1.y/111.5*100,MII2=Bar$MII2.y/111.5*100,MII3=Bar$MII3.y/111.5*100,MII4=Bar$MII4.y/111.5*100,MII5=Bar$MII5.y/111.5*100,MII6=Bar$MII6.y/111.5*100)
Bar<-t(Bar)
colnames(Bar)<-c("Gene")
Bar<-transform(Bar,Sample=c("scSTAP_GV1","scSTAP_GV2","scSTAP_GV3","scSTAP_GV4","scSTAP_GV5","scSTAP_GV6","scSTAP_GV7","scSTAP_GV8","scSTAP_GV9",
                            "scSTAP_MII1","scSTAP_MII2","scSTAP_MII3","scSTAP_MII4","scSTAP_MII5","scSTAP_MII6"),
               Group=c(rep("GV",9),rep("MII",6)))
Bar$Sample<-factor(Bar$Sample,level=c("scSTAP_GV1","scSTAP_GV2","scSTAP_GV3","scSTAP_GV4","scSTAP_GV5","scSTAP_GV6","scSTAP_GV7","scSTAP_GV8","scSTAP_GV9",
                                      "scSTAP_MII1","scSTAP_MII2","scSTAP_MII3","scSTAP_MII4","scSTAP_MII5","scSTAP_MII6"))

p<-ggplot()+
  geom_col(data=Bar,aes(x=Sample,y=Gene,fill=Group),position=position_dodge(width=0.8),width=0.8,color="white",show.legend=F)+
  coord_cartesian(xlim=c(0,16),ylim=c(1,110),expand=F)+
  scale_fill_manual(values = c("#CE5F56","#00789E"),labels= c("GV","MII"))+
  theme(panel.background=element_rect(fill="white",colour="black",size=0.25),
        panel.border=element_rect(fill=NA,color="black",size=2,linetype="solid"),
        axis.line=element_line(colour="#595959",size=0),
        axis.ticks.y=element_line(colour="black",size=0.0),
        axis.ticks.length.y=unit(0.0,'line'),
        axis.ticks.x=element_line(colour="#595959",size=0),
        axis.ticks.length.x=unit(0,'line'),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.y=element_blank(),
        axis.text.x=element_blank())

ggsave("Bar.png",p,width=4.2,height=1.8,dpi=600)

#Rplp0
Bar<-Matrixx[Matrixx$SYMBOL=="Rplp0",]
Bar
Bar<-data.frame(TTseq_GV1_Tr=Bar$TTseq_GV1_Tr/428.5*100,TTseq_GV2_Tr=Bar$TTseq_GV2_Tr/428.5*100,TTseq_MII1_Tr=Bar$TTseq_MII1_Tr/428.5*100,TTseq_MII2_Tr=Bar$TTseq_MII2_Tr/428.5*100)
Bar<-t(Bar)
colnames(Bar)<-c("Gene")
Bar<-transform(Bar,Sample=c("TTseq_GV1","TTseq_GV2","TTseq_MII1","TTseq_MII2"),
               Group=c(rep("GV",2),rep("MII",2)))
Bar$Sample<-factor(Bar$Sample,level=c("TTseq_GV1","TTseq_GV2","TTseq_MII1","TTseq_MII2"))

p<-ggplot()+
  geom_col(data=Bar,aes(x=Sample,y=Gene,fill=Group),position=position_dodge(width=0.8),width=0.8,color="white",show.legend=F)+
  labs(y='Rplp0')+
  coord_cartesian(xlim=c(0,5),ylim=c(1,110),expand=F)+
  scale_fill_manual(values = c("#F4CC66","#758837"),labels= c("GV","MII"))+
  theme(panel.background=element_rect(fill="white",colour="black",size=0.25),
        panel.border=element_rect(fill=NA,color="black",size=2,linetype="solid"),
        axis.line=element_line(colour="#595959",size=0),
        axis.ticks.y=element_line(colour="black",size=0.6),
        axis.ticks.length.y=unit(0.5,'line'),
        axis.ticks.x=element_line(colour="#595959",size=0),
        axis.ticks.length.x=unit(0,'line'),
        axis.title.y=element_text(size=28,face="plain",color="black",vjust=1.5),
        axis.title.x=element_blank(),
        axis.text.y=element_text(size=24,face="plain",color="black",hjust=1),
        axis.text.x=element_blank())

ggsave("Bar.png",p,width=2.4,height=1.8,dpi=600)

Bar<-Matrixx[Matrixx$SYMBOL=="Rplp0",]
Bar
Bar<-data.frame(TTseq_GV1_Tl=Bar$TTseq_GV1_Tl/312.03*100,TTseq_GV2_Tl=Bar$TTseq_GV2_Tl/312.03*100,TTseq_MII1_Tl=Bar$TTseq_MII1_Tl/312.03*100,TTseq_MII2_Tl=Bar$TTseq_MII2_Tl/312.03*100)
Bar<-t(Bar)
colnames(Bar)<-c("Gene")
Bar<-transform(Bar,Sample=c("TTseq_GV1","TTseq_GV2","TTseq_MII1","TTseq_MII2"),
               Group=c(rep("GV",2),rep("MII",2)))
Bar$Sample<-factor(Bar$Sample,level=c("TTseq_GV1","TTseq_GV2","TTseq_MII1","TTseq_MII2"))

p<-ggplot()+
  geom_col(data=Bar,aes(x=Sample,y=Gene,fill=Group),position=position_dodge(width=0.8),width=0.8,color="white",show.legend=F)+
  coord_cartesian(xlim=c(0,5),ylim=c(1,110),expand=F)+
  scale_fill_manual(values = c("#F4CC66","#758837"),labels= c("GV","MII"))+
  theme(panel.background=element_rect(fill="white",colour="black",size=0.25),
        panel.border=element_rect(fill=NA,color="black",size=2,linetype="solid"),
        axis.line=element_line(colour="#595959",size=0),
        axis.ticks.y=element_line(colour="black",size=0),
        axis.ticks.length.y=unit(0,'line'),
        axis.ticks.x=element_line(colour="#595959",size=0),
        axis.ticks.length.x=unit(0,'line'),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.y=element_blank(),
        axis.text.x=element_blank())

ggsave("Bar.png",p,width=1.4,height=1.8,dpi=600)

Bar<-Matrixx[Matrixx$SYMBOL=="Rplp0",]
Bar<-data.frame(GV1=Bar$GV1.x/263.24*100,GV2=Bar$GV2.x/263.24*100,GV3=Bar$GV3.x/263.24*100,GV4=Bar$GV4.x/263.24*100,GV5=Bar$GV5.x/263.24*100,GV6=Bar$GV6.x/263.24*100,GV7=Bar$GV7.x/263.24*100,
                GV8=Bar$GV8.x/263.24*100,GV9=Bar$GV9.x/263.24*100,MII1=Bar$MII1.x/263.24*100,MII2=Bar$MII2.x/263.24*100,MII3=Bar$MII3.x/263.24*100,MII4=Bar$MII4.x/263.24*100,MII5=Bar$MII5.x/263.24*100,MII6=Bar$MII6.x/263.24*100)
Bar<-t(Bar)
colnames(Bar)<-c("Gene")
Bar<-transform(Bar,Sample=c("scSTAP_GV1","scSTAP_GV2","scSTAP_GV3","scSTAP_GV4","scSTAP_GV5","scSTAP_GV6","scSTAP_GV7","scSTAP_GV8","scSTAP_GV9",
                            "scSTAP_MII1","scSTAP_MII2","scSTAP_MII3","scSTAP_MII4","scSTAP_MII5","scSTAP_MII6"),
               Group=c(rep("GV",9),rep("MII",6)))
Bar$Sample<-factor(Bar$Sample,level=c("scSTAP_GV1","scSTAP_GV2","scSTAP_GV3","scSTAP_GV4","scSTAP_GV5","scSTAP_GV6","scSTAP_GV7","scSTAP_GV8","scSTAP_GV9",
                                      "scSTAP_MII1","scSTAP_MII2","scSTAP_MII3","scSTAP_MII4","scSTAP_MII5","scSTAP_MII6"))

p<-ggplot()+
  geom_col(data=Bar,aes(x=Sample,y=Gene,fill=Group),position=position_dodge(width=0.8),width=0.8,color="white",show.legend=F)+
  coord_cartesian(xlim=c(0,16),ylim=c(1,110),expand=F)+
  scale_fill_manual(values = c("#CE5F56","#00789E"),labels= c("GV","MII"))+
  theme(panel.background=element_rect(fill="white",colour="black",size=0.25),
        panel.border=element_rect(fill=NA,color="black",size=2,linetype="solid"),
        axis.line=element_line(colour="#595959",size=0),
        axis.ticks.y=element_line(colour="black",size=0.0),
        axis.ticks.length.y=unit(0.0,'line'),
        axis.ticks.x=element_line(colour="#595959",size=0),
        axis.ticks.length.x=unit(0,'line'),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.y=element_blank(),
        axis.text.x=element_blank())

ggsave("Bar.png",p,width=4.2,height=1.8,dpi=600)

Bar<-Matrixx[Matrixx$SYMBOL=="Rplp0",]
Bar<-data.frame(GV1=Bar$GV1.y/388.82*100,GV2=Bar$GV2.y/388.82*100,GV3=Bar$GV3.y/388.82*100,GV4=Bar$GV4.y/388.82*100,GV5=Bar$GV5.y/388.82*100,GV6=Bar$GV6.y/388.82*100,GV7=Bar$GV7.y/388.82*100,
                GV8=Bar$GV8.y/388.82*100,GV9=Bar$GV9.y/388.82*100,MII1=Bar$MII1.y/388.82*100,MII2=Bar$MII2.y/388.82*100,MII3=Bar$MII3.y/388.82*100,MII4=Bar$MII4.y/388.82*100,MII5=Bar$MII5.y/388.82*100,MII6=Bar$MII6.y/388.82*100)
Bar<-t(Bar)
colnames(Bar)<-c("Gene")
Bar<-transform(Bar,Sample=c("scSTAP_GV1","scSTAP_GV2","scSTAP_GV3","scSTAP_GV4","scSTAP_GV5","scSTAP_GV6","scSTAP_GV7","scSTAP_GV8","scSTAP_GV9",
                            "scSTAP_MII1","scSTAP_MII2","scSTAP_MII3","scSTAP_MII4","scSTAP_MII5","scSTAP_MII6"),
               Group=c(rep("GV",9),rep("MII",6)))
Bar$Sample<-factor(Bar$Sample,level=c("scSTAP_GV1","scSTAP_GV2","scSTAP_GV3","scSTAP_GV4","scSTAP_GV5","scSTAP_GV6","scSTAP_GV7","scSTAP_GV8","scSTAP_GV9",
                                      "scSTAP_MII1","scSTAP_MII2","scSTAP_MII3","scSTAP_MII4","scSTAP_MII5","scSTAP_MII6"))

p<-ggplot()+
  geom_col(data=Bar,aes(x=Sample,y=Gene,fill=Group),position=position_dodge(width=0.8),width=0.8,color="white",show.legend=F)+
  coord_cartesian(xlim=c(0,16),ylim=c(1,110),expand=F)+
  scale_fill_manual(values = c("#CE5F56","#00789E"),labels= c("GV","MII"))+
  theme(panel.background=element_rect(fill="white",colour="black",size=0.25),
        panel.border=element_rect(fill=NA,color="black",size=2,linetype="solid"),
        axis.line=element_line(colour="#595959",size=0),
        axis.ticks.y=element_line(colour="black",size=0.0),
        axis.ticks.length.y=unit(0.0,'line'),
        axis.ticks.x=element_line(colour="#595959",size=0),
        axis.ticks.length.x=unit(0,'line'),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.y=element_blank(),
        axis.text.x=element_blank())

ggsave("Bar.png",p,width=4.2,height=1.8,dpi=600)

#Rpl31
Bar<-Matrixx[Matrixx$SYMBOL=="Rpl31",]
Bar
Bar<-data.frame(TTseq_GV1_Tr=Bar$TTseq_GV1_Tr/155.2*100,TTseq_GV2_Tr=Bar$TTseq_GV2_Tr/155.2*100,TTseq_MII1_Tr=Bar$TTseq_MII1_Tr/155.2*100,TTseq_MII2_Tr=Bar$TTseq_MII2_Tr/155.2*100)
Bar<-t(Bar)
colnames(Bar)<-c("Gene")
Bar<-transform(Bar,Sample=c("GV1","GV2","MII1","MII2"),
               Group=c(rep("GV",2),rep("MII",2)))
Bar$Sample<-factor(Bar$Sample,level=c("GV1","GV2","MII1","MII2"))

p<-ggplot()+
  geom_col(data=Bar,aes(x=Sample,y=Gene,fill=Group),position=position_dodge(width=0.8),width=0.8,color="white",show.legend=F)+
  labs(y='Rpl31')+
  coord_cartesian(xlim=c(0,5),ylim=c(1,110),expand=F)+
  scale_fill_manual(values = c("#F4CC66","#758837"),labels= c("GV","MII"))+
  theme(panel.background=element_rect(fill="white",colour="black",size=0.25),
        panel.border=element_rect(fill=NA,color="black",size=2,linetype="solid"),
        axis.line=element_line(colour="#595959",size=0),
        axis.ticks.y=element_line(colour="black",size=0.6),
        axis.ticks.length.y=unit(0.5,'line'),
        axis.ticks.x=element_line(colour="black",size=0.6),
        axis.ticks.length.x=unit(0.5,'line'),
        axis.title.y=element_text(size=28,face="plain",color="black",vjust=1.5),
        axis.title.x=element_blank(),
        axis.text.y=element_text(size=24,face="plain",color="black",hjust=1),
        axis.text.x=element_text(size=14,face="plain",color="black",angle=90,vjust=0.4))

ggsave("Bar.png",p,width=2.4,height=2.29,dpi=600)

Bar<-Matrixx[Matrixx$SYMBOL=="Rpl31",]
Bar
Bar<-data.frame(TTseq_GV1_Tl=Bar$TTseq_GV1_Tl/204.56*100,TTseq_GV2_Tl=Bar$TTseq_GV2_Tl/204.56*100,TTseq_MII1_Tl=Bar$TTseq_MII1_Tl/204.56*100,TTseq_MII2_Tl=Bar$TTseq_MII2_Tl/204.56*100)
Bar<-t(Bar)
colnames(Bar)<-c("Gene")
Bar<-transform(Bar,Sample=c("GV1","GV2","MII1","MII2"),
               Group=c(rep("GV",2),rep("MII",2)))
Bar$Sample<-factor(Bar$Sample,level=c("GV1","GV2","MII1","MII2"))

p<-ggplot()+
  geom_col(data=Bar,aes(x=Sample,y=Gene,fill=Group),position=position_dodge(width=0.8),width=0.8,color="white",show.legend=F)+
  coord_cartesian(xlim=c(0,5),ylim=c(1,110),expand=F)+
  scale_fill_manual(values = c("#F4CC66","#758837"),labels= c("GV","MII"))+
  theme(panel.background=element_rect(fill="white",colour="black",size=0.25),
        panel.border=element_rect(fill=NA,color="black",size=2,linetype="solid"),
        axis.line=element_line(colour="#595959",size=0),
        axis.ticks.y=element_line(colour="black",size=0),
        axis.ticks.length.y=unit(0,'line'),
        axis.ticks.x=element_line(colour="black",size=0.6),
        axis.ticks.length.x=unit(0.5,'line'),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.y=element_blank(),
        axis.text.x=element_text(size=14,face="plain",color="black",angle=90,vjust=0.4))

ggsave("Bar.png",p,width=1.4,height=2.29,dpi=600)

Bar<-Matrixx[Matrixx$SYMBOL=="Rpl31",]
Bar<-data.frame(GV1=Bar$GV1.x/63.8*100,GV2=Bar$GV2.x/63.8*100,GV3=Bar$GV3.x/63.8*100,GV4=Bar$GV4.x/63.8*100,GV5=Bar$GV5.x/63.8*100,GV6=Bar$GV6.x/63.8*100,GV7=Bar$GV7.x/63.8*100,
                GV8=Bar$GV8.x/63.8*100,GV9=Bar$GV9.x/63.8*100,MII1=Bar$MII1.x/63.8*100,MII2=Bar$MII2.x/63.8*100,MII3=Bar$MII3.x/63.8*100,MII4=Bar$MII4.x/63.8*100,MII5=Bar$MII5.x/63.8*100,MII6=Bar$MII6.x/63.8*100)
Bar<-t(Bar)
colnames(Bar)<-c("Gene")
Bar<-transform(Bar,Sample=c("GV1","GV2","GV3","GV4","GV5","GV6","GV7","GV8","GV9",
                            "MII1","MII2","MII3","MII4","MII5","MII6"),
               Group=c(rep("GV",9),rep("MII",6)))
Bar$Sample<-factor(Bar$Sample,level=c("GV1","GV2","GV3","GV4","GV5","GV6","GV7","GV8","GV9",
                                      "MII1","MII2","MII3","MII4","MII5","MII6"))

p<-ggplot()+
  geom_col(data=Bar,aes(x=Sample,y=Gene,fill=Group),position=position_dodge(width=0.8),width=0.8,color="white",show.legend=F)+
  coord_cartesian(xlim=c(0,16),ylim=c(1,110),expand=F)+
  scale_fill_manual(values = c("#CE5F56","#00789E"),labels= c("GV","MII"))+
  theme(panel.background=element_rect(fill="white",colour="black",size=0.25),
        panel.border=element_rect(fill=NA,color="black",size=2,linetype="solid"),
        axis.line=element_line(colour="#595959",size=0),
        axis.ticks.y=element_line(colour="black",size=0.0),
        axis.ticks.length.y=unit(0.0,'line'),
        axis.ticks.x=element_line(colour="black",size=0.6),
        axis.ticks.length.x=unit(0.5,'line'),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.y=element_blank(),
        axis.text.x=element_text(size=14,face="plain",color="black",angle=90,vjust=0.4))

ggsave("Bar.png",p,width=4.2,height=2.29,dpi=600)

Bar<-Matrixx[Matrixx$SYMBOL=="Rpl31",]
Bar<-data.frame(GV1=Bar$GV1.y/417.601*100,GV2=Bar$GV2.y/417.601*100,GV3=Bar$GV3.y/417.601*100,GV4=Bar$GV4.y/417.601*100,GV5=Bar$GV5.y/417.601*100,GV6=Bar$GV6.y/417.601*100,GV7=Bar$GV7.y/417.601*100,
                GV8=Bar$GV8.y/417.601*100,GV9=Bar$GV9.y/417.601*100,MII1=Bar$MII1.y/417.601*100,MII2=Bar$MII2.y/417.601*100,MII3=Bar$MII3.y/417.601*100,MII4=Bar$MII4.y/417.601*100,MII5=Bar$MII5.y/417.601*100,MII6=Bar$MII6.y/417.601*100)
Bar<-t(Bar)
colnames(Bar)<-c("Gene")
Bar<-transform(Bar,Sample=c("GV1","GV2","GV3","GV4","GV5","GV6","GV7","GV8","GV9",
                            "MII1","MII2","MII3","MII4","MII5","MII6"),
               Group=c(rep("GV",9),rep("MII",6)))
Bar$Sample<-factor(Bar$Sample,level=c("GV1","GV2","GV3","GV4","GV5","GV6","GV7","GV8","GV9",
                                      "MII1","MII2","MII3","MII4","MII5","MII6"))

p<-ggplot()+
  geom_col(data=Bar,aes(x=Sample,y=Gene,fill=Group),position=position_dodge(width=0.8),width=0.8,color="white",show.legend=F)+
  coord_cartesian(xlim=c(0,16),ylim=c(1,110),expand=F)+
  scale_fill_manual(values = c("#CE5F56","#00789E"),labels= c("GV","MII"))+
  theme(panel.background=element_rect(fill="white",colour="black",size=0.25),
        panel.border=element_rect(fill=NA,color="black",size=2,linetype="solid"),
        axis.line=element_line(colour="#595959",size=0),
        axis.ticks.y=element_line(colour="black",size=0.0),
        axis.ticks.length.y=unit(0.0,'line'),
        axis.ticks.x=element_line(colour="black",size=0.6),
        axis.ticks.length.x=unit(0.5,'line'),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.y=element_blank(),
        axis.text.x=element_text(size=14,face="plain",color="black",angle=90,vjust=0.4))

ggsave("Bar.png",p,width=4.2,height=2.29,dpi=600)

