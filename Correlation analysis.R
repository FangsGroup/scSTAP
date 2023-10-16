#Correlation analysis
windowsFonts(ARL = windowsFont("Arial"))

#set work path
setwd("D:/Research/R/Data")

#Venn-Figure 5A
library(VennDiagram)

PRO<-read.csv("PRO_total.csv",sep=",",header=T)
PRO_SMB<-bitr(PRO$PG.ProteinGroups,fromType="UNIPROT",toType="SYMBOL",OrgDb=org.Mm.eg.db)

RNA<-read.csv("RNA_total.csv",sep=",",header=T)
RNA_SMB<-bitr(RNA$gene_id,fromType="ENSEMBL",toType="SYMBOL",OrgDb=org.Mm.eg.db)

venn.diagram(list(Proteome=PRO_SMB$SYMBOL,Transcriptome=RNA_SMB$SYMBOL),
             category.names = c("Proteome","Transcriptome"),
             filename="VennDiagram-T&P.jpg",height=5000,width=5000,resolution=500,
             imagetype="tiff",
             main.pos=c(0.5,0.9), main.cex=4, main.fontfamily="ARL",
             fontfamily="ARL", cat.fontfamily="ARL",
             col="white",lwd=2,
             fill=c("#E43075","#5179C9"),ext.text=F,
             alpha=0.3,cex=3,cat.cex=3.2,cat.fontface="plain",margin=0.2)

#Global Corr-Figure 5B
library(ggplot2)

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

GV1<-data.frame(RNA=Co_RNA[,1],PRO=Co_PRO[,1])
GV2<-data.frame(RNA=Co_RNA[,2],PRO=Co_PRO[,2])
GV3<-data.frame(RNA=Co_RNA[,3],PRO=Co_PRO[,3])
GV4<-data.frame(RNA=Co_RNA[,4],PRO=Co_PRO[,4])
GV5<-data.frame(RNA=Co_RNA[,5],PRO=Co_PRO[,5])
GV6<-data.frame(RNA=Co_RNA[,6],PRO=Co_PRO[,6])
GV7<-data.frame(RNA=Co_RNA[,7],PRO=Co_PRO[,7])
GV8<-data.frame(RNA=Co_RNA[,8],PRO=Co_PRO[,8])
GV9<-data.frame(RNA=Co_RNA[,9],PRO=Co_PRO[,9])
MII1<-data.frame(RNA=Co_RNA[,10],PRO=Co_PRO[,10])
MII2<-data.frame(RNA=Co_RNA[,11],PRO=Co_PRO[,11])
MII3<-data.frame(RNA=Co_RNA[,12],PRO=Co_PRO[,12])
MII4<-data.frame(RNA=Co_RNA[,13],PRO=Co_PRO[,13])
MII5<-data.frame(RNA=Co_RNA[,14],PRO=Co_PRO[,14])
MII6<-data.frame(RNA=Co_RNA[,15],PRO=Co_PRO[,15])

P<-rbind(GV1,GV2,GV3,GV4,GV5,GV6,GV7,GV8,GV9,MII1,MII2,MII3,MII4,MII5,MII6)

cor(P)
#RNA       PRO
#RNA 1.0000000 0.1865799
#PRO 0.1865799 1.0000000

b <- ggplot(P, aes(x = RNA, y = PRO))+
  geom_point(alpha=0.5,shape=21,stroke=0.1,fill="#A00000",color="white",size=3)+
  geom_text(x=0.2,y=0.25,aes(label="italic(R)==0.1866"),parse=T,size=6)+
  scale_x_log10(limits=c(0.1,20000))+scale_y_log10(limits=c(1,200000))+
  labs(x='Exon RPM',y='Protein intensity')+
  theme_classic(base_family="ARL")+
  coord_fixed(ratio=1)+
  theme(panel.background=element_rect(fill="white",colour="black",size=0.25),
        panel.border=element_rect(fill=NA,color="black",size=2,linetype="solid"),
        axis.line=element_line(colour="black",size=0),
        axis.ticks=element_line(colour="black",size=0.8),
        axis.ticks.length=unit(0.5,'line'),
        axis.title.x=element_text(size=20,face="plain",color="black",vjust=-0.5),
        axis.title.y=element_text(size=20,face="plain",color="black",vjust=2),
        axis.text.x=element_text(size=18,face="plain",color="black",vjust=-0.2),
        axis.text.y=element_text(size=18,face="plain",color="black"),
        legend.text=element_text(size=10,face="plain",colour="black"),
        legend.position=c(0.15,0.92),
        legend.title=element_blank(),
        legend.margin=margin(0.05,0.05,0.05,0.05,'cm'),
        legend.key.size=unit(0.5,'cm'))

ggsave("plot-global-corr.png",b,width=5.5,height=5.5,dpi=600)

#Individual Corr
library(psych)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(AnnotationDbi)
library(clusterProfiler)
library(org.Mm.eg.db)
library(Hmisc)
library(stringr)

#Bar-Figure 5C
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

x=1
X=matrix(data=NA,nrow=nrow(Co_RNA),ncol=2)
for(i in 1:nrow(Co_RNA)){
  cor<-corr.test(t(rbind(Co_RNA[i,],Co_PRO[i,])))
  X[i,1]<-cor$r[1,2]
  X[i,2]<-cor$p[1,2]
}

rownames(X)<-MATQ_RNA$SYMBOL
colnames(X)<-c("R","p.values")
X<-X[order(-X[,1]),]

P<-data.frame(R=X[,1],N=1:nrow(X))

p<-ggplot()+
  geom_col(data=P,aes(x=N,y=R,fill=R),width=0.8)+
  coord_cartesian(xlim=c(-10,1355),ylim=c(-1.1,1.1),expand=F)+
  scale_fill_gradientn(colours=c("#07388F","#07388F","#416DB5","#FDE59A","#FB4B33","#E02020","#E02020"))+
  theme_classic()+
  labs(x='RNA-Protein pairs', y='Correlation coefficient')+
  theme(panel.background=element_rect(fill="white",colour="black",size=0.25),
        panel.border=element_rect(fill=NA,color="black",size=1.6,linetype="solid"),
        axis.line=element_line(colour="black",size=0),
        axis.ticks=element_line(colour="black",size=0.8),
        axis.ticks.length=unit(0.5,'line'),
        axis.title.x=element_text(size=22,face="plain",color="black",vjust=-0.2),
        axis.title.y=element_text(size=22,face="plain",color="black",vjust=2),
        axis.text.x=element_text(size=20,face="plain",color="black",vjust=0.1),
        axis.text.y=element_text(size=20,face="plain",color="black"),
        legend.text=element_text(size=18,face="plain",colour="black"),
        legend.title=element_blank(),
        legend.position=c(0.83,0.1),
        legend.margin=margin(0.01,0.01,0.01,0.01,'cm'),
        legend.key.height=unit(0.7,'cm'),
        legend.key.width=unit(1.1,'cm'))+
  guides(fill=guide_colorbar(direction="horizontal",reverse=T))

ggsave("corr_15.png",p,width=16,height=5,dpi=1200)

#Cpeb1-Figure 5C
RNA_corr<-MATQ_RNA[MATQ_RNA$SYMBOL=="Cpeb1",]
PRO_corr<-MATQ_PRO[MATQ_PRO$SYMBOL=="Cpeb1",]

RNA_corr<-t(RNA_corr[,-c(1,2,3)])
colnames(RNA_corr)<-c("RNA")

PRO_corr<-t(PRO_corr[,-c(1,2,3)])
colnames(PRO_corr)<-c("PRO")

Corr<-cbind(RNA_corr,PRO_corr)
cor(Corr)
#RNA       PRO
#RNA 1.0000000 0.9122512
#PRO 0.9122512 1.0000000

Corr<-transform(Corr,Group=c(rep("GV",9),rep("MII",6)))

b <- ggplot(Corr, aes(x = RNA, y = PRO))+
  geom_point(aes(fill = Group, color= Group),size=4,shape=21,stroke=0.2,color="black")+
  geom_text(x=100,y=520,aes(label="Cpeb1"),parse=T,size=8)+
  geom_text(x=380,y=53,aes(label="italic(R)==0.9123"),parse=T,size=6)+
  geom_smooth(method = "lm",se=F,color="black")+
  scale_x_continuous(limits=c(0,500))+scale_y_continuous(limits=c(0,600))+
  scale_fill_manual(values = c("#EE776E","#00A9C9"),labels= c("GV","MII"))+
  scale_color_manual(values = c("#5D110B","#004F64"),labels= c("GV","MII"))+
  labs(x='Exon RPM',y='Protein int.')+
  coord_fixed(ratio=0.4167)+
  theme_classic(base_family="ARL")+
  theme(panel.background=element_rect(fill="white",colour="black",size=0.25),
        panel.border=element_rect(fill=NA,color="black",size=1.5,linetype="solid"),
        axis.line=element_line(colour="black",size=0),
        axis.ticks=element_line(colour="black",size=0.8),
        axis.ticks.length=unit(0.5,'line'),
        axis.title.x=element_text(size=17,face="plain",color="black",vjust=-0.1),
        axis.title.y=element_text(size=17,face="plain",color="black",vjust=2,hjust=0.3),
        axis.text.x=element_text(size=16,face="plain",color="black",vjust=-0.1),
        axis.text.y=element_text(size=16,face="plain",color="black"),
        legend.text=element_text(size=13,face="plain",colour="black"),
        legend.position = 'none',
        legend.title=element_blank(),
        legend.margin=margin(0.001,0.001,0.001,0.001,'cm'),
        legend.key.size=unit(0.05,'cm'))+
  guides(fill=guide_legend(nrow=2,byrow=T))

ggsave("Cpeb1.png",b,width=4,height=2,dpi=1200)

#Slc3a2-Figure 5C
RNA_corr<-MATQ_RNA[MATQ_RNA$SYMBOL=="Slc3a2",]
PRO_corr<-MATQ_PRO[MATQ_PRO$SYMBOL=="Slc3a2",]

RNA_corr<-t(RNA_corr[,-c(1,2,3)])
colnames(RNA_corr)<-c("RNA")

PRO_corr<-t(PRO_corr[,-c(1,2,3)])
colnames(PRO_corr)<-c("PRO")

Corr<-cbind(RNA_corr,PRO_corr)
cor(Corr)
#RNA       PRO
#RNA 1.0000000 0.9093293
#PRO 0.9093293 1.0000000

Corr<-transform(Corr,Group=c(rep("GV",9),rep("MII",6)))

b <- ggplot(Corr, aes(x = RNA, y = PRO))+
  geom_point(aes(fill = Group, color= Group),size=4,shape=21,stroke=0.2,color="black")+
  geom_text(x=60,y=612,aes(label="Slc3a2"),parse=T,size=8)+
  geom_text(x=230,y=61,aes(label="italic(R)==0.9093"),parse=T,size=6)+
  geom_smooth(method = "lm",se=F,color="black")+
  scale_x_continuous(limits=c(0,300))+scale_y_continuous(limits=c(0,700))+
  scale_fill_manual(values = c("#EE776E","#00A9C9"),labels= c("GV","MII"))+
  scale_color_manual(values = c("#5D110B","#004F64"),labels= c("GV","MII"))+
  labs(x='Exon RPM',y='Protein int.')+
  coord_fixed(ratio=0.2143)+
  theme_classic(base_family="ARL")+
  theme(panel.background=element_rect(fill="white",colour="black",size=0.25),
        panel.border=element_rect(fill=NA,color="black",size=1.5,linetype="solid"),
        axis.line=element_line(colour="black",size=0),
        axis.ticks=element_line(colour="black",size=0.8),
        axis.ticks.length=unit(0.5,'line'),
        axis.title.x=element_text(size=17,face="plain",color="black",vjust=-0.1),
        axis.title.y=element_text(size=17,face="plain",color="black",vjust=2,hjust=0.3),
        axis.text.x=element_text(size=16,face="plain",color="black",vjust=-0.1),
        axis.text.y=element_text(size=16,face="plain",color="black"),
        legend.text=element_text(size=13,face="plain",colour="black"),
        legend.position = 'none',
        legend.title=element_blank(),
        legend.margin=margin(0.001,0.001,0.001,0.001,'cm'),
        legend.key.size=unit(0.05,'cm'))+
  guides(fill=guide_legend(nrow=2,byrow=T))

ggsave("Slc3a2.png",b,width=4,height=2,dpi=1200)

#Myo1b-Figure 5C
RNA_corr<-MATQ_RNA[MATQ_RNA$SYMBOL=="Myo1b",]
PRO_corr<-MATQ_PRO[MATQ_PRO$SYMBOL=="Myo1b",]

RNA_corr<-t(RNA_corr[,-c(1,2,3)])
colnames(RNA_corr)<-c("RNA")

PRO_corr<-t(PRO_corr[,-c(1,2,3)])
colnames(PRO_corr)<-c("PRO")

Corr<-cbind(RNA_corr,PRO_corr)
cor(Corr)
#RNA       PRO
#RNA 1.0000000 0.6020951
#PRO 0.6020951 1.0000000

Corr<-transform(Corr,Group=c(rep("GV",9),rep("MII",6)))

b <- ggplot(Corr, aes(x = RNA, y = PRO))+
  geom_point(aes(fill = Group, color= Group),size=4,shape=21,stroke=0.2,color="black")+
  geom_text(x=10,y=525,aes(label="Myo1b"),parse=T,size=8)+
  geom_text(x=38.5,y=52.5,aes(label="italic(R)==0.6021"),parse=T,size=6)+
  geom_smooth(method = "lm",se=F,color="black")+
  scale_x_continuous(limits=c(0,50))+scale_y_continuous(limits=c(0,600))+
  scale_fill_manual(values = c("#EE776E","#00A9C9"),labels= c("GV","MII"))+
  scale_color_manual(values = c("#5D110B","#004F64"),labels= c("GV","MII"))+
  labs(x='Exon RPM',y='Protein int.')+
  coord_fixed(ratio=0.0417)+
  theme_classic(base_family="ARL")+
  theme(panel.background=element_rect(fill="white",colour="black",size=0.25),
        panel.border=element_rect(fill=NA,color="black",size=1.5,linetype="solid"),
        axis.line=element_line(colour="black",size=0),
        axis.ticks=element_line(colour="black",size=0.8),
        axis.ticks.length=unit(0.5,'line'),
        axis.title.x=element_text(size=17,face="plain",color="black",vjust=-0.1),
        axis.title.y=element_text(size=17,face="plain",color="black",vjust=2,hjust=0.3),
        axis.text.x=element_text(size=16,face="plain",color="black",vjust=-0.1),
        axis.text.y=element_text(size=16,face="plain",color="black"),
        legend.text=element_text(size=13,face="plain",colour="black"),
        legend.position = 'none',
        legend.title=element_blank(),
        legend.margin=margin(0.001,0.001,0.001,0.001,'cm'),
        legend.key.size=unit(0.05,'cm'))+
  guides(fill=guide_legend(nrow=2,byrow=T))

ggsave("Myo1b.png",b,width=4,height=2,dpi=1200)

#Bub1b-Figure 5C
RNA_corr<-MATQ_RNA[MATQ_RNA$SYMBOL=="Bub1b",]
PRO_corr<-MATQ_PRO[MATQ_PRO$SYMBOL=="Bub1b",]

RNA_corr<-t(RNA_corr[,-c(1,2,3)])
colnames(RNA_corr)<-c("RNA")

PRO_corr<-t(PRO_corr[,-c(1,2,3)])
colnames(PRO_corr)<-c("PRO")

Corr<-cbind(RNA_corr,PRO_corr)
cor(Corr)
#RNA       PRO
#RNA  1.000000 -0.751593
#PRO -0.751593  1.000000

Corr<-transform(Corr,Group=c(rep("GV",9),rep("MII",6)))

b <- ggplot(Corr, aes(x = RNA, y = PRO))+
  geom_point(aes(fill = Group, color= Group),size=4,shape=21,stroke=0.2,color="black")+
  geom_text(x=320,y=263,aes(label="Bub1b"),parse=T,size=8)+
  geom_text(x=110,y=26.3,aes(label="italic(R)==-0.7516"),parse=T,size=6)+
  geom_smooth(method = "lm",se=F,color="black")+
  scale_x_continuous(limits=c(0,400))+scale_y_continuous(limits=c(0,300))+
  scale_fill_manual(values = c("#EE776E","#00A9C9"),labels= c("GV","MII"))+
  scale_color_manual(values = c("#5D110B","#004F64"),labels= c("GV","MII"))+
  labs(x='Exon RPM',y='Protein int.')+
  coord_fixed(ratio=0.6667)+
  theme_classic(base_family="ARL")+
  theme(panel.background=element_rect(fill="white",colour="black",size=0.25),
        panel.border=element_rect(fill=NA,color="black",size=1.5,linetype="solid"),
        axis.line=element_line(colour="black",size=0),
        axis.ticks=element_line(colour="black",size=0.8),
        axis.ticks.length=unit(0.5,'line'),
        axis.title.x=element_text(size=17,face="plain",color="black",vjust=-0.1),
        axis.title.y=element_text(size=17,face="plain",color="black",vjust=2,hjust=0.3),
        axis.text.x=element_text(size=16,face="plain",color="black",vjust=-0.1),
        axis.text.y=element_text(size=16,face="plain",color="black"),
        legend.text=element_text(size=13,face="plain",colour="black"),
        legend.position = 'none',
        legend.title=element_blank(),
        legend.margin=margin(0.001,0.001,0.001,0.001,'cm'),
        legend.key.size=unit(0.05,'cm'))+
  guides(fill=guide_legend(nrow=2,byrow=T))

ggsave("Bub1b.png",b,width=4,height=2,dpi=1200)

#Oas1e-Figure 5C
RNA_corr<-MATQ_RNA[MATQ_RNA$SYMBOL=="Oas1e",]
PRO_corr<-MATQ_PRO[MATQ_PRO$SYMBOL=="Oas1e",]

RNA_corr<-t(RNA_corr[,-c(1,2,3)])
colnames(RNA_corr)<-c("RNA")

PRO_corr<-t(PRO_corr[,-c(1,2,3)])
colnames(PRO_corr)<-c("PRO")

Corr<-cbind(RNA_corr,PRO_corr)
cor(Corr)
#RNA        PRO
#RNA  1.0000000 -0.8266141
#PRO -0.8266141  1.0000000

Corr<-transform(Corr,Group=c(rep("GV",9),rep("MII",6)))

b <- ggplot(Corr, aes(x = RNA, y = PRO))+
  geom_point(aes(fill = Group, color= Group),size=4,shape=21,stroke=0.2,color="black")+
  geom_text(x=470,y=4854,aes(label="Oas1e"),parse=T,size=8)+
  geom_text(x=170,y=466,aes(label="italic(R)==-0.8266"),parse=T,size=6)+
  geom_smooth(method = "lm",se=F,color="black")+
  scale_x_continuous(limits=c(0,600))+scale_y_continuous(limits=c(0,5600))+
  scale_fill_manual(values = c("#EE776E","#00A9C9"),labels= c("GV","MII"))+
  scale_color_manual(values = c("#5D110B","#004F64"),labels= c("GV","MII"))+
  labs(x='Exon RPM',y='Protein int.')+
  coord_fixed(ratio=0.0538)+
  theme_classic(base_family="ARL")+
  theme(panel.background=element_rect(fill="white",colour="black",size=0.25),
        panel.border=element_rect(fill=NA,color="black",size=1.5,linetype="solid"),
        axis.line=element_line(colour="black",size=0),
        axis.ticks=element_line(colour="black",size=0.8),
        axis.ticks.length=unit(0.5,'line'),
        axis.title.x=element_text(size=17,face="plain",color="black",vjust=-0.1),
        axis.title.y=element_text(size=17,face="plain",color="black",vjust=2,hjust=0.3),
        axis.text.x=element_text(size=16,face="plain",color="black",vjust=-0.1),
        axis.text.y=element_text(size=16,face="plain",color="black"),
        legend.text=element_text(size=13,face="plain",colour="black"),
        legend.position = 'none',
        legend.title=element_blank(),
        legend.margin=margin(0.001,0.001,0.001,0.001,'cm'),
        legend.key.size=unit(0.05,'cm'))+
  guides(fill=guide_legend(nrow=2,byrow=T))

ggsave("Oas1e.png",b,width=4,height=2,dpi=1200)

#GO-Figure S9
X.H<-X[abs(X[,1]) > 0.666,]

gene<-row.names(X.H)
gene<-bitr(gene,fromType="SYMBOL",toType="ENTREZID",OrgDb=org.Mm.eg.db)

BP.params<-enrichGO(gene=gene$ENTREZID,OrgDb=org.Mm.eg.db,ont="BP",pAdjustMethod="BH",pvalueCutoff=0.05,qvalueCutoff=0.05,readable=T)
BP.list<-setReadable(BP.params,org.Mm.eg.db,keyType="ENTREZID")
go<-as.data.frame(BP.params)

go<-go[1:20,]

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

go$GeneRatio <- mixedToFloat(go$GeneRatio)

go$Description<-capitalize(go$Description)
go$Description<-factor(go$Description,levels = rev(go$Description))

go$number <- factor(rev(1:nrow(go)))

p <- ggplot(data=go, aes(x=Description, y=GeneRatio)) +
  geom_point(aes(color=p.adjust,size=Count)) + coord_flip() +
  theme_test() +
  scale_color_gradientn(colours=c("#E02020","#FDE59A","#07388F"))+
  scale_x_discrete(labels = function(x) str_wrap(x,width = 55))+
  scale_y_continuous(limits=c(0.02,0.13))+
  xlab("GO term") +
  ylab("Gene ratio") +
  theme(panel.background=element_rect(fill="white",colour="black",size=0.25),
        panel.border=element_rect(fill=NA,color="black",size=1.6,linetype="solid"),
        plot.title=element_text(colour="black",size=22,hjust=0.5),
        axis.ticks=element_line(colour="black",size=0.6),
        axis.ticks.length=unit(0.4,'line'),
        axis.line=element_line(colour="black",size=0.0),
        axis.title.x=element_text(size=20,face="plain",color="black"),
        axis.title.y=element_text(size=20,face="plain",color="black",vjust=2.0),
        axis.text.x=element_text(size=16,face="plain",color="black"),
        axis.text.y=element_text(size=14,face="plain",color="black",lineheight=0.6),
        legend.text=element_text(size=14,face="plain",colour="black"),
        legend.title=element_text(size=14,face="plain",colour="black")) +
  labs(title = "Biological process")+
  guides(color=guide_colorbar(order = 1),size=guide_legend(order = 2))

ggsave("go_point-BP.png",p,width=8,height=6)

CC.params<-enrichGO(gene=gene$ENTREZID,OrgDb=org.Mm.eg.db,ont="CC",pAdjustMethod="BH",pvalueCutoff=0.05,qvalueCutoff=0.05,readable=T)
CC.list<-setReadable(CC.params,org.Mm.eg.db,keyType="ENTREZID")
go<-as.data.frame(CC.params)

go<-go[1:20,]

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

go$GeneRatio <- mixedToFloat(go$GeneRatio)

go$Description<-capitalize(go$Description)
go$Description<-factor(go$Description,levels = rev(go$Description))

go$number <- factor(rev(1:nrow(go)))


p <- ggplot(data=go, aes(x=Description, y=GeneRatio)) +
  geom_point(aes(color=p.adjust,size=Count)) + coord_flip() +
  theme_test() +
  scale_color_gradientn(colours=c("#E02020","#FDE59A","#07388F"))+
  scale_x_discrete(labels = function(x) str_wrap(x,width = 50))+
  scale_y_continuous(limits=c(0.015,0.14))+
  xlab("GO term") +
  ylab("Gene ratio") +
  theme(panel.background=element_rect(fill="white",colour="black",size=0.25),
        panel.border=element_rect(fill=NA,color="black",size=1.6,linetype="solid"),
        plot.title=element_text(colour="black",size=22,hjust=0.5),
        axis.ticks=element_line(colour="black",size=0.6),
        axis.ticks.length=unit(0.4,'line'),
        axis.line=element_line(colour="black",size=0.0),
        axis.title.x=element_text(size=20,face="plain",color="black"),
        axis.title.y=element_text(size=20,face="plain",color="black",vjust=2.0),
        axis.text.x=element_text(size=16,face="plain",color="black"),
        axis.text.y=element_text(size=16,face="plain",color="black",lineheight=0.6),
        legend.text=element_text(size=14,face="plain",colour="black"),
        legend.title=element_text(size=14,face="plain",colour="black")) +
  labs(title = "Cellular component")+
  guides(color=guide_colorbar(order = 1),size=guide_legend(order = 2))

ggsave("go_point-CC.png",p,width=8,height=6)


#Fold change Corr
library(ggbeeswarm)
library(limma)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(AnnotationDbi)
library(clusterProfiler)
library(org.Mm.eg.db)
library(Hmisc)

#Plot-Figure 5D
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
limma.na.rna <- na.omit(tmpOut)

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
limma.na.pro <- na.omit(tmpOut)

limma.na.pro <- transform(limma.na.pro,SYMBOL=rownames(limma.na.pro))
limma.na.rna <- transform(limma.na.rna,SYMBOL=rownames(limma.na.rna))
matrix<-merge(limma.na.pro,limma.na.rna,by.x="SYMBOL",by.y="SYMBOL")

matrix$change <- as.factor(ifelse(abs(matrix$logFC.x)>1 | abs(matrix$logFC.y)>1,
                                  ifelse(abs(matrix$logFC.x)>1,
                                         ifelse(abs(matrix$logFC.y)>1,
                                                ifelse(matrix$logFC.x*matrix$logFC.y>1,"Class I","Class II"),"Class III"),"Class IV"),"Class V"))

matrix$label <- ifelse(abs(matrix$logFC.x)>1 & abs(matrix$logFC.y)>1, as.character(matrix$SYMBOL),"")

p<-ggplot(data=matrix)+
  geom_point(aes(x = -logFC.x,y = -logFC.y,fill=change),color="white",size=5,shape=21,stroke=0.2,alpha=0.65,show.legend = T)+
  geom_text_repel(aes(x = -logFC.x,y = -logFC.y,label = label),size = 5,force=2,box.padding = unit(0.5,"lines"),point.padding = unit(0.5,"lines"),segment.color = "black",color="black",show.legend = FALSE)+
  geom_vline(xintercept = c(-log2(2),log2(2)),lty = 4,col = "black",lwd = 0.8)+
  geom_hline(yintercept = c(-log2(2),log2(2)),lty = 4,col = "black",lwd = 0.8)+
  scale_y_continuous(limits=c(-5,4))+
  scale_x_continuous(limits=c(-5,4))+
  scale_fill_manual(values = c("#FA7F6F","#82B0D2","#FFBE7A","#8ECFC9","#E7DAD2"))+
  xlab(expression(Proteome_log[2](FoldChange)))+
  ylab(expression(Transcriptome_log[2](FoldChange)))+
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
        legend.position=c(0.18,0.83),
        legend.title=element_blank(),
        legend.margin=margin(0.1,0.1,0.1,0.1,'cm'),
        legend.background=element_rect(fill="transparent"),
        legend.key.height=unit(0.5,'cm'),
        legend.key=element_rect(fill="transparent",color="transparent"))+
  guides(fill=guide_legend(nrow=5,byrow=T))

ggsave("Volcano_Omics.png",p,width=6,height=6,dpi=1200)

#Bar-Figure 5D
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
  coord_fixed(ratio=0.0105)+
  scale_fill_manual(values = c("#FA7F6F","#82B0D2","#FFBE7A","#8ECFC9","#E7DAD2"),labels= c("Class I","Class II","Class III","Class IV","Class V"))+
  theme(panel.background=element_rect(fill="white",colour="black",size=0.25),
        panel.border=element_rect(fill=NA,color="black",size=2.2,linetype="solid"),
        axis.line=element_line(colour="black",size=0),
        axis.ticks=element_line(colour="black",size=1),
        axis.ticks.length=unit(0.5,'line'),
        axis.title.y.right=element_text(size=24,face="plain",color="black",vjust=2),
        axis.title.x=element_blank(),
        axis.text.y.right=element_text(size=23,face="plain",color="black",hjust=0.2),
        axis.text.x=element_text(size=20,face="plain",color="black",angle=60,vjust=1.05,hjust=1.1))

ggsave("Bar-Omics.png",p,width=3.8,height=6,dpi=600)

#GO-Figure 5E
GO<-matrix[matrix$change=="Class I",]

gene<-bitr(GO$SYMBOL,fromType="SYMBOL",toType="ENTREZID",OrgDb=org.Mm.eg.db)

BP.params<-enrichGO(gene=gene$ENTREZID,OrgDb=org.Mm.eg.db,ont="BP",pAdjustMethod="BH",pvalueCutoff=0.05,qvalueCutoff=0.05,readable=T)
BP.list<-setReadable(BP.params,org.Mm.eg.db,keyType="ENTREZID")
go<-as.data.frame(BP.params)

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

go<-go[1:20,]

p <- ggplot(data=go, aes(x=GeneRatio, y=Description)) +
  geom_point(aes(color=p.adjust,size=Count)) +
  theme_test() +
  scale_color_gradientn(colours=c("#E02020","#FDE59A","#07388F"))+
  scale_y_discrete(labels = labels)+
  scale_x_continuous(limits=c(0.04,0.18))+
  coord_fixed(ratio=0.022)+
  ylab("GO terms of biological process") +
  xlab("Gene ratio") +
  ggtitle("GO terms entichment of class I genes")+
  theme(panel.background=element_rect(fill="#FFFFFF",colour="black",size=0.25),
        panel.border=element_rect(fill=NA,color="black",size=2,linetype="solid"),
        plot.title=element_text(colour="black",size=20,hjust=0.9),
        axis.ticks=element_line(colour="black",size=0.8),
        axis.ticks.length=unit(0.4,'line'),
        axis.line=element_line(colour="black",size=0),
        axis.title.x=element_text(size=18,face="plain",color="black",vjust=0),
        axis.title.y=element_text(size=18,face="plain",color="black",vjust=2.0),
        axis.text.x=element_text(size=14,face="plain",color="black",vjust=0),
        axis.text.y=element_text(size=14,face="plain",color="black",lineheight=0.6),
        legend.title=element_text(size=14,face="plain",colour="black"),
        legend.text=element_text(size=14,face="plain",colour="black"),
        legend.background=element_rect(fill="#FFFFFF"))

ggsave("go_point.png",p,width=10.5,height=6.4,dpi=600)

