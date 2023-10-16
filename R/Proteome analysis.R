#Proteome analysis
windowsFonts(ARL = windowsFont("Arial"))

#set work path
setwd("D:/Research/R/Data")

#Venn-DIA vs DDA
library(ggplot2)
library(VennDiagram)
library(corrplot)
library(PerformanceAnalytics)

#DIA-Figure 2A
Data<-read.csv("DIA_ProQ.csv",sep=",",header=T)

Data[Data=='Filtered']<-NA
Data[Data=='NaN']<-NA

I<-data.frame(PG=Data[,1],I=Data[,2])
II<-data.frame(PG=Data[,1],II=Data[,3])
III<-data.frame(PG=Data[,1],III=Data[,4])
IV<-data.frame(PG=Data[,1],IV=Data[,5])
V<-data.frame(PG=Data[,1],V=Data[,6])

I.na=na.omit(I)
II.na=na.omit(II)
III.na=na.omit(III)
IV.na=na.omit(IV)
V.na=na.omit(V)

venn.diagram(list(I=I.na$PG,II=II.na$PG,III=III.na$PG,IV=IV.na$PG,V=V.na$PG),
             filename="VennDiagram-DIA.tif",height=3000,width=3000,resolution=500,
             imagetype="tiff", fontfamily="ARL", cat.fontfamily="ARL",
             col="white",lwd=2,
             fill=c("dodgerblue","goldenrod1","darkorange1","seagreen3","orchid3"),
             alpha=0.5,cex=1.5,cat.cex=2.5,cat.fontface="plain",margin=0.05)

#DDA-Figure S4
Data<-read.csv("DDA_ProQ.csv",sep=",",header=T)

Data[Data=='Filtered']<-NA
Data[Data=='NaN']<-NA

I<-data.frame(PG=Data[,1],I=Data[,2])
II<-data.frame(PG=Data[,1],II=Data[,3])
III<-data.frame(PG=Data[,1],III=Data[,4])
IV<-data.frame(PG=Data[,1],IV=Data[,5])
V<-data.frame(PG=Data[,1],V=Data[,6])

I.na=na.omit(I)
II.na=na.omit(II)
III.na=na.omit(III)
IV.na=na.omit(IV)
V.na=na.omit(V)

venn.diagram(list(I=I.na$PG,II=II.na$PG,III=III.na$PG,IV=IV.na$PG,V=V.na$PG),
             filename="VennDiagram-DDA.tif",height=3000,width=3000,resolution=500,
             imagetype="tiff", fontfamily="ARL", cat.fontfamily="ARL",
             col="white",lwd=2,
             fill=c("dodgerblue","goldenrod1","darkorange1","seagreen3","orchid3"),
             alpha=0.5,cex=1.5,cat.cex=2.5,cat.fontface="plain",margin=0.05)

#Corr-DIA vs DDA
library(ggplot2)
library(GGally)
library(corrplot)
library(Hmisc)
library(PerformanceAnalytics)

#DIA-Figure 2B
Data<-read.csv("DIA_ProQ.csv",sep=",",header=T)

Data[Data=='Filtered']<-NA
Data[Data=='NaN']<-NA

Data.na=na.omit(Data)

DIA<-data.frame(PRO1=as.numeric(Data.na$DIA_1),
                PRO2=as.numeric(Data.na$DIA_2),
                PRO3=as.numeric(Data.na$DIA_3),
                PRO4=as.numeric(Data.na$DIA_4),
                PRO5=as.numeric(Data.na$DIA_5))


cor_func <- function(data, mapping, method, symbol, ...){
  
  x <- eval_data_col(data, mapping$x)
  y <- eval_data_col(data, mapping$y)
  corr <- cor(x, y, method=method,use='complete.obs')
  
  ggally_text(mapping=aes(),label = paste(symbol, as.character(sprintf("%0.3f",corr))),size=10,xP=0.45,yP=0.5,color="black")+
    theme(panel.grid.major=element_line(colour=NA),panel.background = element_rect(fill = "white"))
}

p<-ggpairs(log2(DIA+1), showStrips=T,
           lower=list(continuous=wrap("points", pch=21,size=3,stroke=0.15,fill="#A00000",color="white",alpha=0.2)),
           diag=list(continuous=wrap("densityDiag",color="#D00000")),
           upper=list(continuous=wrap(cor_func,method="pearson",symbol=expression(' ***\n'))))

p[1,1]<-p[1,1]+scale_y_continuous(limits=c(-0.05,0.25))
p[2,2]<-p[2,2]+scale_y_continuous(limits=c(-0.05,0.25))
p[3,3]<-p[3,3]+scale_y_continuous(limits=c(-0.05,0.25))
p[4,4]<-p[4,4]+scale_y_continuous(limits=c(-0.05,0.25))
p[5,5]<-p[5,5]+scale_y_continuous(limits=c(-0.05,0.25))

p<-p+
  theme(panel.background=element_rect(fill="white",colour="black",size=0.3),
        panel.border=element_rect(fill=NA,color="black",size=2,linetype="solid"),
        axis.text=element_text(colour="black",size=18),
        strip.background=element_rect(fill="white"),
        strip.text=element_text(colour="black",size=20,face="plain"))

ggsave("Figure 2B-DIA.png",p,width=7,height=7,dpi=600)

#DDA-Figure S4
Data<-read.csv("DDA_ProQ.csv",sep=",",header=T)

Data[Data=='Filtered']<-NA
Data[Data=='NaN']<-NA

Data.na=na.omit(Data)

DDA<-data.frame(I=as.numeric(Data.na$DDA_1),
                II=as.numeric(Data.na$DDA_2),
                III=as.numeric(Data.na$DDA_3),
                IV=as.numeric(Data.na$DDA_4),
                V=as.numeric(Data.na$DDA_5))


cor_func <- function(data, mapping, method, symbol, ...){
  
  x <- eval_data_col(data, mapping$x)
  y <- eval_data_col(data, mapping$y)
  corr <- cor(x, y, method=method,use='complete.obs')
  
  ggally_text(mapping=aes(),label = paste(symbol, as.character(sprintf("%0.3f",corr))),size=10,xP=0.45,yP=0.5,color="black")+
    theme(panel.grid.major=element_line(colour=NA),panel.background = element_rect(fill = "white"))
}

p<-ggpairs(log2(DDA+1), showStrips=T,
           lower=list(continuous=wrap("points", pch=21,size=3,stroke=0.15,fill="#A00000",color="white",alpha=0.2)),
           diag=list(continuous=wrap("densityDiag",color="#D00000")),
           upper=list(continuous=wrap(cor_func,method="pearson",symbol=expression(' ***\n'))))

p[1,1]<-p[1,1]+scale_x_continuous(limits=c(6,20))+scale_y_continuous(limits=c(-0.05,0.25))
p[2,2]<-p[2,2]+scale_x_continuous(limits=c(6,20))+scale_y_continuous(limits=c(-0.05,0.25))
p[3,3]<-p[3,3]+scale_x_continuous(limits=c(6,20))+scale_y_continuous(limits=c(-0.05,0.25))
p[4,4]<-p[4,4]+scale_x_continuous(limits=c(6,20))+scale_y_continuous(limits=c(-0.05,0.25))
p[5,5]<-p[5,5]+scale_x_continuous(limits=c(6,20))+scale_y_continuous(limits=c(-0.05,0.25))

p[2,1]<-p[2,1]+scale_x_continuous(limits=c(6,20))+scale_y_continuous(limits=c(6,20))
p[3,1]<-p[3,1]+scale_x_continuous(limits=c(6,20))+scale_y_continuous(limits=c(6,20))
p[3,2]<-p[3,2]+scale_x_continuous(limits=c(6,20))+scale_y_continuous(limits=c(6,20))
p[4,1]<-p[4,1]+scale_x_continuous(limits=c(6,20))+scale_y_continuous(limits=c(6,20))
p[4,2]<-p[4,2]+scale_x_continuous(limits=c(6,20))+scale_y_continuous(limits=c(6,20))
p[4,3]<-p[4,3]+scale_x_continuous(limits=c(6,20))+scale_y_continuous(limits=c(6,20))
p[5,1]<-p[5,1]+scale_x_continuous(limits=c(6,20))+scale_y_continuous(limits=c(6,20))
p[5,2]<-p[5,2]+scale_x_continuous(limits=c(6,20))+scale_y_continuous(limits=c(6,20))
p[5,3]<-p[5,3]+scale_x_continuous(limits=c(6,20))+scale_y_continuous(limits=c(6,20))
p[5,4]<-p[5,4]+scale_x_continuous(limits=c(6,20))+scale_y_continuous(limits=c(6,20))

p<-p+
  theme(panel.background=element_rect(fill="white",colour="black",size=0.3),
        panel.border=element_rect(fill=NA,color="black",size=2,linetype="solid"),
        axis.text=element_text(colour="black",size=18),
        strip.background=element_rect(fill="white"),
        strip.text=element_text(colour="black",size=20,face="plain"))

ggsave("Figure S4-DDA.png",p,width=7,height=7,dpi=600)

#Violin-DIA vs DDA
library(ggplot2)

#Figure S5
DIA<-read.csv("DIA_CV.csv",sep=",",header=T)
DDA<-read.csv("DDA_CV.csv",sep=",",header=T)

DIA<-data.frame(type="DIA",CV=DIA$CV)
DDA<-data.frame(type="DDA",CV=DDA$CV)

P<-rbind(DDA,DIA)

p<-ggplot(P, aes(x = type, y = CV, fill = type)) +
  geom_violin(alpha = 0.4) +
  geom_boxplot(width = 0.2, alpha = 0.6, show.legend = FALSE) +
  scale_x_discrete(name =element_blank())+
  scale_y_continuous(name = "Coefficent of variation (CV)",limits=c(0,1.7)) +
  scale_fill_manual(values=c("#C00000","#205080"),labels=c("DDA","DIA")) +
  coord_fixed(ratio=1.18)+
  theme(panel.background=element_rect(fill="white",colour="black",size=0),
        panel.border=element_rect(fill=NA,color="black",size=2,linetype="solid"),
        axis.line=element_line(colour="black",size=0),
        axis.ticks=element_line(colour="black",size=0.85),
        axis.ticks.length=unit(0.5,'line'),
        axis.title.x=element_text(size=24,face="bold",color="black"),
        axis.title.y=element_text(size=30,face="plain",color="black",vjust=1.5),
        axis.text.x=element_text(size=28,face="plain",color="black",vjust=-0.7),
        axis.text.y=element_text(size=26,face="plain",color="black"),
        legend.text=element_text(size=24,face="plain",colour="black"),
        legend.position=c(0.73,0.92),
        legend.title=element_blank(),
        legend.margin=margin(0.1,0.1,0.1,0.1,'cm'),
        legend.key.size=unit(0.8,'cm'),
        legend.spacing.x=unit(0.4,'cm'))+
  guides(fill=guide_legend(nrow=1,byrow=T))

ggsave("Violin-CV.png",p,width=7,height=7,dpi=1200)

#GV vs MII
#Bar-Figure 2C
library(ggplot2)

P<-read.csv("GVvsMII.csv",sep=",",header=T)

p<-ggplot()+
  geom_col(data=P,aes(x=X,y=Number,fill=Group), width=0.8)+
  scale_y_continuous(limits=c(0,4000))+
  scale_x_discrete(labels=P$Group)+
  theme_classic()+
  labs(x='Single oocytes',y='Protein groups')+
  scale_fill_manual(values = c("#CE5F56","#00789E","#FFFFFF"),labels= c("GV","MII","H"))+
  theme(panel.background=element_rect(fill="white",colour="black",size=0.01),
        panel.border=element_rect(fill=NA,color="black",size=1.6,linetype="solid"),
        axis.line=element_line(colour="black",size=0),
        axis.ticks=element_line(colour="black",size=0.7),
        axis.ticks.length=unit(0.5,'line'),
        axis.title.x=element_text(size=24,face="plain",color="black",vjust=0),
        axis.title.y=element_text(size=24,face="plain",color="black",vjust=2.0),
        axis.text=element_text(size=20,face="plain",color="black"),
        legend.text=element_text(size=18,face="plain",colour="black"),
        legend.position=c(0.92,0.88),
        legend.title=element_blank(),
        legend.margin=margin(0.01,0.01,0.01,0.01,'cm'),
        legend.key.size=unit(0.7,'cm'))+
  guides(fill=guide_legend(nrow=2,byrow=T))

ggsave("Bar-GVvsMII.png",p,width=10,height=5,dpi=1200)

#Venn-Figure 2D
library(VennDiagram)

GV<-read.csv("GV_18.csv",sep=",",header=T)
MII<-read.csv("MII_18.csv",sep=",",header=T)

venn.diagram(list(MII=MII$PG.ProteinGroups,GV=GV$PG.ProteinGroups),
             filename="VennDiagram-GVvsMII.tif",height=3000,width=3000,resolution=500,
             imagetype="tiff", fontfamily="ARL", cat.fontfamily="ARL",
             col="white",lwd=2,
             fill=c("#013EDF","#FF4F4F"),
             alpha=0.3,cex=3,cat.cex=4,cat.fontface="plain",margin=0.05)

#tSNE-Figure 2E
library(dplyr)
library(Seurat)
library(SeuratData)
library(SeuratDisk)
library(cowplot)
library("Matrix")

Data<-read.csv("GVvsMII_M.csv",sep=",",header=T,row.names=1)
Data[Data=='Filtered']<-NA
Data[Data=='NaN']<-NA
Data<-na.omit(Data)
Data_2<-as(as.matrix(Data),"dgCMatrix")

mmoc<-CreateSeuratObject(counts=Data_2,assay="PRO")
mmoc<-NormalizeData(mmoc, normalization.method="LogNormalize")
mmoc<-FindVariableFeatures(mmoc,selection.method="vst",nfeatures = 2000)
all.genes<-rownames(mmoc)
mmoc<-ScaleData(mmoc,features=all.genes)
mmoc<-RunPCA(mmoc, npcs=17, verbose=F)
mmoc<-JackStraw(mmoc, num.replicate = 100)
mmoc<-ScoreJackStraw(mmoc, dims = 1:17)
plot1<-JackStrawPlot(mmoc, dims = 1:10)
plot2<-ElbowPlot(mmoc)
plot1+plot2
mmoc<- FindNeighbors(mmoc, dims=1:10)
mmoc<- FindClusters(mmoc, resolution = 0.8)
mmoc<- RunTSNE(mmoc, dims = 1:10, perplexity=10)

Idents(mmoc)
new.cluster.ids=c("GV","MII")
names(new.cluster.ids)<-levels(mmoc)
mmoc<-RenameIdents(mmoc,new.cluster.ids)

p<-ggplot(data.frame(mmoc@meta.data,mmoc@reductions$tsne@cell.embeddings))+
  geom_point(aes(tSNE_1,tSNE_2,fill=mmoc@active.ident,color=mmoc@active.ident),size=5,shape=21,stroke=0.5)+
  scale_y_continuous(limits=c(-50,50))+
  scale_x_continuous(limits=c(-50,50))+
  scale_fill_manual(values = c("#EE776E","#00A9C9"),labels= c("GV","MII"))+
  scale_color_manual(values = c("#5D110B","#004F64"),labels= c("GV","MII"))+
  coord_fixed(ratio=1)+
  theme(panel.background=element_rect(fill="white",colour="black",size=0.25),
        panel.border=element_rect(fill=NA,color="black",size=2,linetype="solid"),
        axis.line=element_line(colour="black",size=0),
        axis.ticks=element_line(colour="black",size=0.85),
        axis.ticks.length=unit(0.5,'line'),
        axis.title.x=element_text(size=28,face="plain",color="black",vjust=-0.7),
        axis.title.y=element_text(size=28,face="plain",color="black",vjust=1.8),
        axis.text.x=element_text(size=25,face="plain",color="black",vjust=-0.2),
        axis.text.y=element_text(size=25,face="plain",color="black"),
        legend.text=element_text(size=25,face="plain",colour="black"),
        legend.position=c(0.88,0.88),
        legend.title=element_blank(),
        legend.margin=margin(0.1,0.1,0.1,0.1,'cm'),
        legend.key.height=unit(1.1,'cm'),
        legend.key=element_rect(fill="white"))+
  guides(fill=guide_legend(nrow=2,byrow=T))

ggsave("tSNE-GVvsMII.png",p,width=7,height=7,dpi=1200)

#Box-Figure 2F
library(ggplot2)
library(ggsignif)
library(rstatix)
library(ggpubr)
library(ggtext)

P<-read.csv("box.csv",sep=",",header=T)
P$Type<-factor(P$Type,levels=c("ACTB","DPPA5a","RPS9","DPPA3","BCL2l10","BUB1b","EZHIP","PLAT"))

df_p_val1 <- P%>% 
  group_by(Type) %>% 
  t_test(formula = Intensity ~ Group) %>% 
  add_significance(p.col = 'p',cutpoints = c(0,0.001,0.01,0.05,1),symbols = c('***','**','*','n.s.')) %>% 
  add_xy_position(x='Type')

df_p_val1

p<-ggplot()+
  geom_boxplot(P,mapping=aes(x=Type,y=Intensity,fill=Group),width=0.8)+
  stat_pvalue_manual(df_p_val1,label="{p.signif}",tip.length=0.01,size=6)+
  xlab('Biomarkers')+
  ylab(expression(Log[2](Intensity)))+
  scale_y_continuous(limits=c(1,16),breaks=c(3,6,9,12,15))+
  scale_fill_manual(values = c("#EE776E","#00A9C9"),labels= c("GV","MII"))+
  coord_fixed(ratio=0.255)+
  theme(panel.background=element_rect(fill="white",colour="black",size=0.25),
        panel.border=element_rect(fill=NA,color="black",size=1.6,linetype="solid"),
        axis.line=element_line(colour="black",size=0),
        axis.ticks=element_line(colour="black",size=0.6),
        axis.ticks.length=unit(0.5,'line'),
        axis.title.x=element_text(size=22,face="plain",color="black",vjust=-0.7),
        axis.title.y=element_text(size=22,face="plain",color="black",vjust=2),
        axis.text.x=element_text(size=18,face="plain",color="black",vjust=-0.2),
        axis.text.y=element_text(size=18,face="plain",color="black"),
        legend.text=element_text(size=18,face="plain",colour="black"),
        legend.title=element_blank(),
        legend.margin=margin(0.1,0.1,0.1,0.1,'cm'),
        legend.key.size=unit(0.7,'cm'),
        legend.key=element_rect(fill="white"),
        legend.position=c(0.92,0.87))+
  guides(fill=guide_legend(nrow=2,byrow=T))

ggsave("box-GVvsMII.png",p,width=9,height=5,dpi=1200)

#Benchmark
library(ggplot2)

#Violin-Figure S6
P<-read.csv("benchmark-oocyte.csv",sep=",",header=T)

i=1
data<-as.data.frame(matrix(nrow=0,ncol=4))
for(i in 1:dim(P)[1]){
  data[i,1]<-P[i,1]
  data[i,2]<-mean(c(P[i,3],P[i,4],P[i,5],P[i,6]))
  data[i,3]<-mean(c(P[i,7],P[i,8],P[i,9],P[i,10]))
  data[i,4]<-mean(c(P[i,11],P[i,12],P[i,13],P[i,14]))
}

Ratio1<-data.frame(Type="1.5 ng/0.5 ng",ratio=data$V4/data$V2)
Ratio2<-data.frame(Type="1.5 ng/1.0 ng",ratio=data$V4/data$V3)
Ratio3<-data.frame(Type="1.0 ng/0.5 ng",ratio=data$V3/data$V2)

median(Ratio1$ratio)
#3.02
#t-ratio = 3

median(Ratio2$ratio)
#1.43
#t-ratio = 1.5

median(Ratio3$ratio)
#2.12
#t-ratio = 2


P<-rbind(Ratio1,Ratio2,Ratio3)
P$Type<-factor(P$Type,levels=c("1.5 ng/0.5 ng","1.0 ng/0.5 ng","1.5 ng/1.0 ng"))

p<-ggplot(P, aes(x = Type, y = ratio, fill = Type)) +
  geom_violin(alpha = 0.4) +
  annotate("text",x=1.32,y=3.02,label=" Median:\n3.02",size=6.5)+
  annotate("text",x=2.32,y=2.12,label=" Median:\n2.12",size=6.5)+
  annotate("text",x=3.32,y=1.43,label=" Median:\n1.43",size=6.5)+
  geom_boxplot(width = 0.2, alpha = 0.6, show.legend = FALSE) +
  scale_fill_manual(values = c("#82B0D2","#FA7F6F","#FFBE7A"))+
  scale_y_log10(name = "Protein abundance ratio",limits=c(0.2,10)) +
  coord_fixed(ratio=1.8)+
  ggtitle("Mouse oocyte")+
  theme(panel.background=element_rect(fill="white",colour="black",size=0.25),
        panel.border=element_rect(fill=NA,color="black",size=2,linetype="solid"),
        plot.title=element_text(hjust=0.5,vjust=1.5,size=35,color="black"),
        axis.line=element_line(colour="black",size=0),
        axis.ticks=element_line(colour="black",size=0.85),
        axis.ticks.length=unit(0.5,'line'),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=30,face="plain",color="black",vjust=2),
        axis.text.x=element_text(size=24,face="plain",color="black"),
        axis.text.y=element_text(size=22,face="plain",color="black"),
        legend.text=element_text(size=24,face="plain",colour="black"),
        legend.position=c(0.80,0.90),
        legend.title=element_blank(),
        legend.margin=margin(0.1,0.1,0.1,0.1,'cm'),
        legend.key.size=unit(0.3,'cm'))+
  guides(fill=guide_legend(nrow=3,byrow=T))

ggsave("Violin-benchmark-oocyte.png",p,width=9,height=9,dpi=1200)


P<-read.csv("benchmark-yeast.csv",sep=",",header=T)

i=1
data<-as.data.frame(matrix(nrow=0,ncol=4))
for(i in 1:dim(P)[1]){
  data[i,1]<-P[i,1]
  data[i,2]<-mean(c(P[i,3],P[i,4],P[i,5],P[i,6]))
  data[i,3]<-mean(c(P[i,7],P[i,8],P[i,9],P[i,10]))
  data[i,4]<-mean(c(P[i,11],P[i,12],P[i,13],P[i,14]))
}

Ratio1<-data.frame(Type="1.5 ng/0.5 ng",ratio=data$V2/data$V4)
Ratio2<-data.frame(Type="1.5 ng/1.0 ng",ratio=data$V2/data$V3)
Ratio3<-data.frame(Type="1.0 ng/0.5 ng",ratio=data$V3/data$V4)

median(Ratio1$ratio)
#2.70
#t-ratio = 3

median(Ratio2$ratio)
#1.44
#t-ratio = 1.5

median(Ratio3$ratio)
#1.90
#t-ratio = 2


P<-rbind(Ratio1,Ratio2,Ratio3)
P$Type<-factor(P$Type,levels=c("1.5 ng/0.5 ng","1.0 ng/0.5 ng","1.5 ng/1.0 ng"))


p<-ggplot(P, aes(x = Type, y = ratio, fill = Type)) +
  geom_violin(alpha = 0.4) +
  annotate("text",x=1.32,y=2.70,label=" Median:\n2.70",size=6.5)+
  annotate("text",x=2.32,y=1.90,label=" Median:\n1.90",size=6.5)+
  annotate("text",x=3.32,y=1.44,label=" Median:\n1.44",size=6.5)+
  geom_boxplot(width = 0.2, alpha = 0.6, show.legend = FALSE) +
  scale_fill_manual(values = c("#82B0D2","#FA7F6F","#FFBE7A"))+
  scale_y_log10(name = "Protein abundance ratio",limits=c(0.2,10)) +
  coord_fixed(ratio=1.8)+
  ggtitle("Yeast")+
  theme(panel.background=element_rect(fill="white",colour="black",size=0.25),
        panel.border=element_rect(fill=NA,color="black",size=2,linetype="solid"),
        plot.title=element_text(hjust=0.5,vjust=1.5,size=35,color="black"),
        axis.line=element_line(colour="black",size=0),
        axis.ticks=element_line(colour="black",size=0.85),
        axis.ticks.length=unit(0.5,'line'),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=30,face="plain",color="black",vjust=2),
        axis.text.x=element_text(size=24,face="plain",color="black"),
        axis.text.y=element_text(size=22,face="plain",color="black"),
        legend.text=element_text(size=24,face="plain",colour="black"),
        legend.position=c(0.80,0.90),
        legend.title=element_blank(),
        legend.margin=margin(0.1,0.1,0.1,0.1,'cm'),
        legend.key.size=unit(0.3,'cm'))+
  guides(fill=guide_legend(nrow=3,byrow=T))

ggsave("Violin-benchmark-Yeast.png",p,width=9,height=9,dpi=1200)



#Compared with Published data
library(AnnotationDbi)
library(clusterProfiler)
library(org.Mm.eg.db)
library(ggplot2)
library(dplyr)
library(psych)
library(corrplot)
library(PerformanceAnalytics)
library(VennDiagram)
library(stringr)
library(tidyverse)
library(GGally)

#Venn-Figure S7
#Ref.S1-GV
Data1<-read.csv("GV_symbol.csv",sep=",",header=T)
Data2<-read.csv("GV_18_F.csv",sep=",",header=T)

gene2<-bitr(Data2$PG.ProteinGroups,fromType="UNIPROT",toType="SYMBOL",OrgDb=org.Mm.eg.db)

venn.diagram(list(scSTAP=gene2$SYMBOL,Ref_S1=Data1$GV_symbol),
             filename="VennDiagram-GV-Ref1.tif",height=3000,width=3000,resolution=500,
             imagetype="tiff",
             main="GV", main.pos=c(0.5,0.9), main.cex=4, main.fontfamily="ARL",
             fontfamily="ARL", cat.fontfamily="ARL",cat.pos = c(-40,40),
             col="white",lwd=2,
             fill=c("#013EDF","#FF4F4F"),
             alpha=0.3,cex=3,cat.cex=3.2,cat.fontface="plain",margin=0.2)

#Ref.S1-MII
Data1<-read.csv("MII_symbol.csv",sep=",",header=T)
Data2<-read.csv("MII_18_F.csv",sep=",",header=T)

gene2<-bitr(Data2$PG.ProteinGroups,fromType="UNIPROT",toType="SYMBOL",OrgDb=org.Mm.eg.db)

venn.diagram(list(scSTAP=gene2$SYMBOL,Ref_S1=Data1$MII_symbol),
             filename="VennDiagram-MII-Ref1.tif",height=3000,width=3000,resolution=500,
             imagetype="tiff",
             main="MII", main.pos=c(0.5,0.9), main.cex=4, main.fontfamily="ARL",
             fontfamily="ARL", cat.fontfamily="ARL",cat.pos = c(-40,40),
             col="white",lwd=2,
             fill=c("#013EDF","#FF4F4F"),
             alpha=0.3,cex=3,cat.cex=3.2,cat.fontface="plain",margin=0.2)

#Ref.S2-GV
Data1<-read.csv("TMT_group1_GV.csv",sep=",",header=T)
Data2<-read.csv("GV_DIA_18.csv",sep=",",header=T)

M<-strsplit(Data1$Accession,"\\|")

x=1
X=matrix(data=NA,nrow=dim(Data1)[1],ncol=1)
for(i in 1:dim(Data1)[1]){
  X[i,1]<-M[[i]][1]
}

Data1<-data.frame(PG.ProteinAccessions=as.character(X[,1]),GV=Data1$GV,
                  Rank=Data1$Rank_GV)

venn.diagram(list(scSTAP=Data2$PG.ProteinAccessions,Ref_S2=Data1$PG.ProteinAccessions),
             filename="VennDiagram-GV-Ref2.tif",height=3000,width=3000,resolution=500,
             imagetype="tiff",
             main="GV", main.pos=c(0.5,0.9), main.cex=4, main.fontfamily="ARL",
             fontfamily="ARL", cat.fontfamily="ARL",
             col="white",lwd=2,
             fill=c("#013EDF","#FF4F4F"),
             alpha=0.3,cex=3,cat.cex=3.2,cat.fontface="plain",margin=0.2)

#Ref.S2-MII
Data1<-read.csv("TMT_group1_MII.csv",sep=",",header=T)
Data2<-read.csv("MII_DIA_18.csv",sep=",",header=T)

M<-strsplit(Data1$Accession,"\\|")

x=1
X=matrix(data=NA,nrow=dim(Data1)[1],ncol=1)
for(i in 1:dim(Data1)[1]){
  X[i,1]<-M[[i]][1]
}

Data1<-data.frame(PG.ProteinAccessions=as.character(X[,1]),MII=Data1$MII,
                  Rank=Data1$Rank_MII)

venn.diagram(list(scSTAP=Data2$PG.ProteinAccessions,Ref_S2=Data1$PG.ProteinAccessions),
             filename="VennDiagram-MII-Ref2.tif",height=3000,width=3000,resolution=500,
             imagetype="tiff",
             main="MII", main.pos=c(0.5,0.9), main.cex=4, main.fontfamily="ARL",
             fontfamily="ARL", cat.fontfamily="ARL",
             col="white",lwd=2,
             fill=c("#013EDF","#FF4F4F"),
             alpha=0.3,cex=3,cat.cex=3.2,cat.fontface="plain",margin=0.2)

#Rank-Figure S7
#Ref.S2-GV
Data1<-read.csv("TMT_group1_GV.csv",sep=",",header=T)
Data2<-read.csv("GV_DIA_18.csv",sep=",",header=T)

M<-strsplit(Data1$Accession,"\\|")

x=1
X=matrix(data=NA,nrow=dim(Data1)[1],ncol=1)
for(i in 1:dim(Data1)[1]){
  X[i,1]<-M[[i]][1]
}

Data1<-data.frame(PG.ProteinAccessions=as.character(X[,1]),GV=Data1$GV,
                  Rank=Data1$Rank_GV)

Data3<-merge(Data2,Data1,by.x="PG.ProteinAccessions",by.y="PG.ProteinAccessions")
Data11<-anti_join(Data1,Data3,by="PG.ProteinAccessions")
Data21<-anti_join(Data2,Data3,by="PG.ProteinAccessions")

Data_Our_N<-data.frame(Accession=Data3$PG.ProteinAccessions,
                       Area=Data3$GV.x,Rank=Data3$Rank.x,
                       Group="Our",Sample="Identified in both two datasets")
Data_Paper_N<-data.frame(Accession=Data3$PG.ProteinAccessions,
                         Area=Data3$GV.y,Rank=Data3$Rank.y,
                         Group="Paper",Sample="Identified in both two datasets")


Data_Our_L<-data.frame(Accession=Data21$PG.ProteinAccessions,
                       Area=Data21$GV,Rank=Data21$Rank,
                       Group="Our",Sample="Identified in either one dataset")
Data_Paper_L<-data.frame(Accession=Data11$PG.ProteinAccessions,
                         Area=Data11$GV,Rank=Data11$Rank,
                         Group="Paper",Sample="Identified in either one dataset")

Data_Total<-rbind(Data_Our_L,Data_Our_N,Data_Paper_L,Data_Paper_N)
Data_Total<-arrange(Data_Total,Accession)

p<-ggplot(Data_Total)+
  geom_point(aes(x=Rank,y=log2(Area),fill=Sample),alpha=0.5,size=3,shape=21,stroke=0.2,color="white")+
  ggtitle("GV")+
  scale_fill_manual(values = c("#00A9C9","#EE776E"),labels= c("Identified in both\ntwo datasets","Identified in either\none dataset"))+
  ylab(expression(Log[2](area)))+
  theme_bw()+
  theme(panel.border=element_rect(fill=NA,color="black",size=1.6,linetype="solid"),
        plot.title=element_text(hjust=0.5,vjust=2,size=25,color="black"),
        axis.line=element_line(colour="black",size=0),
        axis.ticks=element_line(colour="black",size=0.7),
        axis.ticks.length=unit(0.3,'line'),
        axis.title.x=element_text(size=17,face="plain",color="black",vjust=0),
        axis.title.y=element_text(size=17,face="plain",color="black",vjust=2),
        axis.text.x=element_text(size=15,face="plain",color="black"),
        axis.text.y=element_text(size=15,face="plain",color="black"),
        legend.text=element_text(size=13,face="plain",colour="black",lineheight=1),
        legend.title=element_text(size=14,face="plain",colour="black"),
        legend.key.height=unit(1.5,"cm"))+
  labs(fill="Group")

ggsave("plot-Rank-GV.png",p,width=6,height=3,dpi=1200)

#Ref.S2-MII
Data1<-read.csv("TMT_group1_MII.csv",sep=",",header=T)
Data2<-read.csv("MII_DIA_18.csv",sep=",",header=T)

M<-strsplit(Data1$Accession,"\\|")

x=1
X=matrix(data=NA,nrow=dim(Data1)[1],ncol=1)
for(i in 1:dim(Data1)[1]){
  X[i,1]<-M[[i]][1]
}

Data1<-data.frame(PG.ProteinAccessions=as.character(X[,1]),MII=Data1$MII,
                  Rank=Data1$Rank_MII)

Data3<-merge(Data2,Data1,by.x="PG.ProteinAccessions",by.y="PG.ProteinAccessions")
Data11<-anti_join(Data1,Data3,by="PG.ProteinAccessions")
Data21<-anti_join(Data2,Data3,by="PG.ProteinAccessions")


Data_Our_N<-data.frame(Accession=Data3$PG.ProteinAccessions,
                       Area=Data3$MII.x,Rank=Data3$Rank.x,
                       Group="Our",Sample="Identified in both two datasets")
Data_Paper_N<-data.frame(Accession=Data3$PG.ProteinAccessions,
                         Area=Data3$MII.y,Rank=Data3$Rank.y,
                         Group="Paper",Sample="Identified in both two datasets")


Data_Our_L<-data.frame(Accession=Data21$PG.ProteinAccessions,
                       Area=Data21$MII,Rank=Data21$Rank,
                       Group="Our",Sample="Identified in either one dataset")
Data_Paper_L<-data.frame(Accession=Data11$PG.ProteinAccessions,
                         Area=Data11$MII,Rank=Data11$Rank,
                         Group="Paper",Sample="Identified in either one dataset")

Data_Total<-rbind(Data_Our_L,Data_Our_N,Data_Paper_L,Data_Paper_N)
Data_Total<-arrange(Data_Total,Accession)

p<-ggplot(Data_Total)+
  geom_point(aes(x=Rank,y=log2(Area),fill=Sample),alpha=0.5,size=3,shape=21,stroke=0.2,color="white")+
  ggtitle("MII")+
  scale_fill_manual(values = c("#00A9C9","#EE776E"),labels= c("Identified in both\ntwo datasets","Identified in either\none dataset"))+
  ylab(expression(Log[2](area)))+
  theme_bw()+
  theme(panel.border=element_rect(fill=NA,color="black",size=1.6,linetype="solid"),
        plot.title=element_text(hjust=0.5,vjust=2,size=25,color="black"),
        axis.line=element_line(colour="black",size=0),
        axis.ticks=element_line(colour="black",size=0.7),
        axis.ticks.length=unit(0.3,'line'),
        axis.title.x=element_text(size=17,face="plain",color="black",vjust=0),
        axis.title.y=element_text(size=17,face="plain",color="black",vjust=2),
        axis.text.x=element_text(size=15,face="plain",color="black"),
        axis.text.y=element_text(size=15,face="plain",color="black"),
        legend.text=element_text(size=13,face="plain",colour="black",lineheight=1),
        legend.title=element_text(size=14,face="plain",colour="black"),
        legend.key.height=unit(1.5,"cm"))+
  labs(fill="Group")

ggsave("plot-Rank-MII.png",p,width=6,height=3,dpi=1200)

#Corr-Figure S7
#Ref.S2-GV
Data1<-read.csv("TMT_group1_GV.csv",sep=",",header=T)
Data2<-read.csv("GV_DIA_18.csv",sep=",",header=T)

M<-strsplit(Data1$Accession,"\\|")

x=1
X=matrix(data=NA,nrow=dim(Data1)[1],ncol=1)
for(i in 1:dim(Data1)[1]){
  X[i,1]<-M[[i]][1]
}

Data1<-data.frame(PG.ProteinAccessions=as.character(X[,1]),GV=Data1$GV,
                  Rank=Data1$Rank_GV)

Data3<-merge(Data2,Data1,by.x="PG.ProteinAccessions",by.y="PG.ProteinAccessions")
Data3<-Data3[-which(Data3$Rank.x>2500),]
Data3<-Data3[-which(Data3$Rank.y>2500),]

Corr<-data.frame(Ref.S2=Data3$GV.y,scSTAP=Data3$GV.x)
cor(log2(Corr))

cor_func <- function(data, mapping, method, symbol, ...){
  
  x <- eval_data_col(data, mapping$x)
  y <- eval_data_col(data, mapping$y)
  corr <- cor(x, y, method=method,use='complete.obs')
  
  ggally_text(mapping=aes(),label = paste(symbol, as.character(sprintf("%0.3f",corr))),size=10,xP=0.45,yP=0.5,color="black")+
    theme(panel.grid.major=element_line(colour=NA),panel.background = element_rect(fill = "white"))
}

p<-ggpairs(log2(Corr), showStrips=T,
           lower=list(continuous=wrap("points", pch=21,size=3,stroke=0.15,fill="#A00000",color="white",alpha=0.2)),
           diag=list(continuous=wrap("densityDiag",color="#D00000")),
           upper=list(continuous=wrap(cor_func,method="pearson",symbol=expression(' ***\n'))))

p[1,1]<-p[1,1]+scale_x_continuous(limits=c(12.5,26))
p[2,1]<-p[2,1]+scale_x_continuous(limits=c(12.5,26))+scale_y_continuous(limits=c(4,16))
p[2,2]<-p[2,2]+scale_x_continuous(limits=c(4,16))

p<-p+
  theme(panel.background=element_rect(fill="white",colour="black",size=0.3),
        panel.border=element_rect(fill=NA,color="black",size=2,linetype="solid"),
        axis.text=element_text(colour="black",size=18),
        strip.background=element_rect(fill="white"),
        strip.text=element_text(colour="black",size=20,face="plain"))

ggsave("Figure S7-Corr-GV.png",p,width=5,height=5,dpi=600)


#Ref.S2-MII
Data1<-read.csv("TMT_group1_MII.csv",sep=",",header=T)
Data2<-read.csv("MII_DIA_18.csv",sep=",",header=T)

M<-strsplit(Data1$Accession,"\\|")

x=1
X=matrix(data=NA,nrow=dim(Data1)[1],ncol=1)
for(i in 1:dim(Data1)[1]){
  X[i,1]<-M[[i]][1]
}

Data1<-data.frame(PG.ProteinAccessions=as.character(X[,1]),MII=Data1$MII,
                  Rank=Data1$Rank_MII)

Data3<-merge(Data2,Data1,by.x="PG.ProteinAccessions",by.y="PG.ProteinAccessions")
Data3<-Data3[-which(Data3$Rank.x>2500),]
Data3<-Data3[-which(Data3$Rank.y>2500),]

Corr<-data.frame(Ref.S2=Data3$MII.y,scSTAP=Data3$MII.x)
cor(log2(Corr))

cor_func <- function(data, mapping, method, symbol, ...){
  
  x <- eval_data_col(data, mapping$x)
  y <- eval_data_col(data, mapping$y)
  corr <- cor(x, y, method=method,use='complete.obs')
  
  ggally_text(mapping=aes(),label = paste(symbol, as.character(sprintf("%0.3f",corr))),size=10,xP=0.45,yP=0.5,color="black")+
    theme(panel.grid.major=element_line(colour=NA),panel.background = element_rect(fill = "white"))
}

p<-ggpairs(log2(Corr), showStrips=T,
           lower=list(continuous=wrap("points", pch=21,size=3,stroke=0.15,fill="#A00000",color="white",alpha=0.2)),
           diag=list(continuous=wrap("densityDiag",color="#D00000")),
           upper=list(continuous=wrap(cor_func,method="pearson",symbol=expression(' ***\n'))))

p[1,1]<-p[1,1]+scale_x_continuous(limits=c(12.5,26))
p[2,1]<-p[2,1]+scale_x_continuous(limits=c(12.5,26))+scale_y_continuous(limits=c(4,16))
p[2,2]<-p[2,2]+scale_x_continuous(limits=c(4,16))

p<-p+
  theme(panel.background=element_rect(fill="white",colour="black",size=0.3),
        panel.border=element_rect(fill=NA,color="black",size=2,linetype="solid"),
        axis.text=element_text(colour="black",size=18),
        strip.background=element_rect(fill="white"),
        strip.text=element_text(colour="black",size=20,face="plain"))

ggsave("Figure S7-Corr-MII.png",p,width=5,height=5,dpi=600)


#Blank analysis
library(VennDiagram)
library(AnnotationDbi)
library(clusterProfiler)
library(org.Mm.eg.db)
library(dplyr)
library(grid)
library(stringr)
library(plyr)
library(Hmisc)
library(ggplot2)

#Venn-Figure S8
Blank1<-read.csv("Blank_1.csv",sep=",",header=T)
Blank2<-read.csv("Blank_2.csv",sep=",",header=T)
Blank3<-read.csv("Blank_3.csv",sep=",",header=T)
Blank4<-read.csv("Blank_4.csv",sep=",",header=T)

venn.diagram(list(Blank1=Blank1$PG.ProteinGroups,Blank2=Blank2$PG.ProteinGroups,Blank3=Blank3$PG.ProteinGroups,Blank4=Blank4$PG.ProteinGroups),
             filename="VennDiagram-Blank.png",height=3000,width=3900,resolution=500,
             imagetype="png",
             fontfamily="ARL", cat.fontfamily="ARL",
             col="white",lwd=2,
             fill=c("#8ECFC9","#FFBE7A","#FA7F6F","#82B0D2"),
             alpha=0.3,cex=2,cat.cex=2,cat.fontface="plain",margin=0.01)

#GO plot-Figure S8
Blank<-merge(Blank1,Blank2,by.x="PG.ProteinGroups",by.y="PG.ProteinGroups")
Blank<-merge(Blank,Blank3,by.x="PG.ProteinGroups",by.y="PG.ProteinGroups")
Blank<-merge(Blank,Blank4,by.x="PG.ProteinGroups",by.y="PG.ProteinGroups")

gene<-bitr(Blank$PG.ProteinGroups,fromType="UNIPROT",toType="ENTREZID",OrgDb=org.Mm.eg.db)

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

p <- ggplot(data=go, aes(x=Description, y=GeneRatio)) +
  geom_point(aes(color=p.adjust,size=Count)) + coord_flip() +
  theme_test() +
  scale_color_gradientn(colours=c("#E02020","#FDE59A","#07388F"))+
  scale_x_discrete(labels = labels)+
  scale_y_continuous(limits=c(0.07,0.4))+
  xlab("GO terms") +
  ylab("Gene ratio") +
  ggtitle("Biological process")+
  theme(panel.background=element_rect(fill="#FFFFFF",colour="black",size=0.25),
        panel.border=element_rect(fill=NA,color="black",size=2,linetype="solid"),
        plot.title=element_text(colour="black",size=20,hjust=0.4),
        axis.ticks=element_line(colour="black",size=0.8),
        axis.ticks.length=unit(0.4,'line'),
        axis.line=element_line(colour="black",size=0),
        axis.title.x=element_text(size=20,face="plain",color="black",vjust=0),
        axis.title.y=element_text(size=20,face="plain",color="black",vjust=2.0),
        axis.text.x=element_text(size=16,face="plain",color="black",vjust=0),
        axis.text.y=element_text(size=16,face="plain",color="black",lineheight=0.6),
        legend.title=element_text(size=14,face="plain",colour="black"),
        legend.text=element_text(size=14,face="plain",colour="black"),
        legend.background=element_rect(fill="#FFFFFF"))+
  guides(color=guide_colorbar(order = 1),size=guide_legend(order = 2))

ggsave("go_point-blank.png",p,width=8,height=6,dpi=600)
