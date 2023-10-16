#Consistency analysis of splitting
windowsFonts(ARL = windowsFont("Arial"))

#set work path
setwd("D:/Research/R/Data")

#Figure 1B-Violin
library(ggplot2)
library(pacman)
library(introdataviz)
library(tidyverse)

#Transcriptome
Data<-read.csv("Cell 1.csv",sep=",",header=T)

A1<-Data$GV1
A2<-Data$GV2

t.test(A1,A2,paired=T)

#data:  A1 and A2
#t = -0.19286, df = 14973, p-value = 0.8471
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -1.103644  0.905915
#sample estimates:
#  mean of the differences 
#-0.0988643 

Data<-read.csv("Cell 2.csv",sep=",",header=T)

A3<-Data$GV1
A4<-Data$GV2

t.test(A3,A4,paired=T)

#data:  A3 and A4
#t = -0.27018, df = 15368, p-value = 0.787
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -1.0078421  0.7636619
#sample estimates:
#  mean of the differences 
#-0.1220901

Data<-read.csv("Cell 3.csv",sep=",",header=T)

A5<-Data$GV2_1
A6<-Data$GV2_2

t.test(A5,A6,paired=T)

#data:  A5 and A6
#t = -0.029797, df = 15258, p-value = 0.9762
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -2.083669  2.021268
#sample estimates:
#  mean of the differences 
#-0.03120039

A1D<-data.frame(Intensity=A1,type="Part A",group="Cell1")
A2D<-data.frame(Intensity=A2,type="Part B",group="Cell1")

A3D<-data.frame(Intensity=A3,type="Part A",group="Cell2")
A4D<-data.frame(Intensity=A4,type="Part B",group="Cell2")

A5D<-data.frame(Intensity=A5,type="Part A",group="Cell3")
A6D<-data.frame(Intensity=A6,type="Part B",group="Cell3")

P<-rbind(A1D,A2D,A3D,A4D,A5D,A6D)

p<-ggplot(P, aes(x = group, y = Intensity, fill = type)) +
  introdataviz::geom_split_violin(alpha = 0.4, show.legend = FALSE) +
  geom_boxplot(width = 0.2, alpha = 0.6, show.legend = FALSE) +
  geom_text(x=1,y=-2.5,aes(label="italic(p)==0.8471"),parse=T,size=7)+
  geom_text(x=2,y=-2.5,aes(label=sprintf("italic(p)=='%0.4f'",0.7870)),parse=T,size=7)+
  geom_text(x=3,y=-2.5,aes(label="italic(p)==0.9762"),parse=T,size=7)+
  scale_x_discrete(name =element_blank(), labels = c("Cell 1#", "Cell 2#","Cell 3#")) +
  scale_y_log10(name = "Gene counts",limits=c(0.003,300000)) +
  scale_fill_manual(values=c("#C00000","#205080")) +
  coord_fixed(ratio=0.363)+
  theme(panel.background=element_rect(fill="white",colour="black",size=0.25),
        panel.border=element_rect(fill=NA,color="black",size=2,linetype="solid"),
        axis.line=element_line(colour="black",size=0),
        axis.ticks=element_line(colour="black",size=0.85),
        axis.ticks.length=unit(0.5,'line'),
        axis.title.x=element_text(size=24,face="bold",color="black"),
        axis.title.y=element_text(size=30,face="plain",color="black",vjust=1.8),
        axis.text.x=element_text(size=28,face="plain",color="black",vjust=-0.7),
        axis.text.y=element_text(size=22,face="plain",color="black"),
        legend.text=element_text(size=20,face="plain",colour="black"),
        legend.position=c(0.7,0.92),
        legend.title=element_blank(),
        legend.margin=margin(0.1,0.1,0.1,0.1,'cm'),
        legend.key.size=unit(0.8,'cm'),
        legend.spacing.x=unit(0.4,'cm'))+
  guides(fill=guide_legend(nrow=1,byrow=T))

ggsave("Figure 1B-T.png",p,width=7,height=7,dpi=1200)

#Proteome
Data<-read.csv("Cell 4.csv",sep=",",header=T)

A1<-Data$X44_1
A2<-Data$X44_2

t.test(A1,A2,paired=T)

#data:  A1 and A2
#t = 0.47529, df = 2431, p-value = 0.6346
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -30.55210  50.10074
#sample estimates:
#  mean of the differences 
#9.774321 

Data<-read.csv("Cell 5.csv",sep=",",header=T)

A3<-Data$X45_1
A4<-Data$X45_2

t.test(A3,A4,paired=T)

#data:  A3 and A4
#t = -0.97392, df = 2639, p-value = 0.3302
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -44.03259  14.80788
#sample estimates:
#  mean of the differences 
#-14.61236 

Data<-read.csv("Cell 6.csv",sep=",",header=T)

A5<-Data$X48_1
A6<-Data$X48_2

t.test(A5,A6,paired=T)

#data:  A5 and A6
#t = -1.2765, df = 2735, p-value = 0.2019
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -24.467890   5.172491
#sample estimates:
#  mean of the differences 
#-9.6477 

A1D<-data.frame(Intensity=A1,type="Part A",group="Cell4")
A2D<-data.frame(Intensity=A2,type="Part B",group="Cell4")

A3D<-data.frame(Intensity=A3,type="Part A",group="Cell5")
A4D<-data.frame(Intensity=A4,type="Part B",group="Cell5")

A5D<-data.frame(Intensity=A5,type="Part A",group="Cell6")
A6D<-data.frame(Intensity=A6,type="Part B",group="Cell6")

P<-rbind(A1D,A2D,A3D,A4D,A5D,A6D)

p<-ggplot(P, aes(x = group, y = Intensity, fill = type)) +
  introdataviz::geom_split_violin(alpha = 0.4, show.legend = FALSE) +
  geom_boxplot(width = 0.2, alpha = 0.6, show.legend = FALSE) +
  geom_text(x=1,y=-0.5,aes(label="italic(p)==0.6346"),parse=T,size=7)+
  geom_text(x=2,y=-0.5,aes(label="italic(p)==0.3302"),parse=T,size=7)+
  geom_text(x=3,y=-0.5,aes(label="italic(p)==0.2019"),parse=T,size=7)+
  scale_x_discrete(name =element_blank(), labels = c("Cell 4#", "Cell 5#","Cell 6#")) +
  scale_y_log10(name = "Protein intensity",limits=c(0.3,300000)) +
  scale_fill_manual(values=c("#C00000","#205080")) +
  coord_fixed(ratio=0.485)+
  theme(panel.background=element_rect(fill="white",colour="black",size=0.25),
        panel.border=element_rect(fill=NA,color="black",size=2,linetype="solid"),
        axis.line=element_line(colour="black",size=0),
        axis.ticks=element_line(colour="black",size=0.85),
        axis.ticks.length=unit(0.5,'line'),
        axis.title.x=element_text(size=24,face="bold",color="black"),
        axis.title.y=element_text(size=30,face="plain",color="black",vjust=1.8),
        axis.text.x=element_text(size=28,face="plain",color="black",vjust=-0.7),
        axis.text.y=element_text(size=22,face="plain",color="black"),
        legend.text=element_text(size=20,face="plain",colour="black"),
        legend.position=c(0.7,0.92),
        legend.title=element_blank(),
        legend.margin=margin(0.1,0.1,0.1,0.1,'cm'),
        legend.key.size=unit(0.8,'cm'),
        legend.spacing.x=unit(0.4,'cm'))+
  guides(fill=guide_legend(nrow=1,byrow=T))

ggsave("Figure 1B-P.png",p,width=7,height=7,dpi=1200)


#Figure 1C-Venn
library(ggplot2)
library(VennDiagram)
library(corrplot)
library(PerformanceAnalytics)

#Transcriptome
Intact<-read.csv("Intact_gene.csv",sep=",",header=T)
Half<-read.csv("Half_gene.csv",sep=",",header=T)

venn.diagram(list(Intact=Intact$gene_id,Half=Half$gene_id),
             filename="VennDiagram-T.tif",height=3000,width=3000,resolution=500,
             imagetype="tiff", fontfamily="ARL", cat.fontfamily="ARL",
             col="white",lwd=2,
             fill=c("#013EDF","#FF4F4F"),
             alpha=0.3,cex=3,cat.cex=4,cat.fontface="plain",margin=0.05)

#Proteome
Intact<-read.csv("Intact_Pro.csv",sep=",",header=T)
Half<-read.csv("Half_Pro.csv",sep=",",header=T)

venn.diagram(list(Intact=Intact$PG.ProteinGroups,Half=Half$PG.ProteinGroups),
             filename="VennDiagram-P.tif",height=3000,width=3000,resolution=500,
             imagetype="tiff", fontfamily="ARL", cat.fontfamily="ARL",
             col="white",lwd=2,
             fill=c("#013EDF","#FF4F4F"),
             alpha=0.3,cex=3,cat.cex=4,cat.fontface="plain",margin=0.05)

#Figure S2-Corr
library(ggplot2)

#Transcriptome
#Cell1
Data<-read.csv("Cell 1.csv",sep=",",header=T)

A1<-Data$GV1
A2<-Data$GV2

P<-data.frame(PA=A1,PB=A2)
cor(log2(P))

#PA        PB
#PA 1.0000000 0.9152723
#PB 0.9152723 1.0000000

b <- ggplot(P, aes(x = PA, y = PB))+
  geom_point(color="#A00000")+
  geom_text(x=2.1,y=-1.7,aes(label="italic(R)==0.915"),parse=T,size=16)+
  geom_text(x=-0.3,y=3.95,aes(label="Cell 1#"),parse=F,size=16)+
  scale_x_log10(limits=c(0.01,20000))+scale_y_log10(limits=c(0.01,20000))+
  labs(x='Part A',y='Part B')+
  theme_classic(base_family="ARL")+
  coord_fixed(ratio=1)+
  theme(panel.background=element_rect(fill="white",colour="black",size=0.25),
        panel.border=element_rect(fill=NA,color="black",size=3,linetype="solid"),
        axis.line=element_line(colour="black",size=0),
        axis.ticks=element_line(colour="black",size=1.3),
        axis.ticks.length=unit(0.5,'line'),
        axis.title.x=element_text(size=32,face="plain",color="black",vjust=-0.5),
        axis.title.y=element_text(size=32,face="plain",color="black",vjust=1.5),
        axis.text.x=element_text(size=25,face="plain",color="black"),
        axis.text.y=element_text(size=25,face="plain",color="black"),
        legend.text=element_text(size=10,face="plain",colour="black"),
        legend.position=c(0.85,0.92),
        legend.title=element_blank(),
        legend.margin=margin(0.05,0.05,0.05,0.05,'cm'),
        legend.key.size=unit(0.5,'cm'))

ggsave("Figure S2-T1.png",b,width=6,height=6,dpi=1200)

#Cell2
Data<-read.csv("Cell 2.csv",sep=",",header=T)

A3<-Data$GV1
A4<-Data$GV2

P<-data.frame(PA=A3,PB=A4)
cor(log2(P))

#PA        PB
#PA 1.0000000 0.9173341
#PB 0.9173341 1.0000000

b <- ggplot(P, aes(x = PA, y = PB))+
  geom_point(color="#A00000")+
  geom_text(x=2.1,y=-1.7,aes(label="italic(R)==0.917"),parse=T,size=16)+
  geom_text(x=-0.3,y=3.95,aes(label="Cell 2#"),parse=F,size=16)+
  scale_x_log10(limits=c(0.01,20000))+scale_y_log10(limits=c(0.01,20000))+
  labs(x='Part A',y='Part B')+
  theme_classic(base_family="ARL")+
  coord_fixed(ratio=1)+
  theme(panel.background=element_rect(fill="white",colour="black",size=0.25),
        panel.border=element_rect(fill=NA,color="black",size=3,linetype="solid"),
        axis.line=element_line(colour="black",size=0),
        axis.ticks=element_line(colour="black",size=1.3),
        axis.ticks.length=unit(0.5,'line'),
        axis.title.x=element_text(size=32,face="plain",color="black",vjust=-0.5),
        axis.title.y=element_text(size=32,face="plain",color="black",vjust=1.5),
        axis.text.x=element_text(size=25,face="plain",color="black"),
        axis.text.y=element_text(size=25,face="plain",color="black"),
        legend.text=element_text(size=10,face="plain",colour="black"),
        legend.position=c(0.85,0.92),
        legend.title=element_blank(),
        legend.margin=margin(0.05,0.05,0.05,0.05,'cm'),
        legend.key.size=unit(0.5,'cm'))

ggsave("Figure S2-T2.png",b,width=6,height=6,dpi=1200)

#Cell3
Data<-read.csv("Cell 3.csv",sep=",",header=T)

A5<-Data$GV2_1
A6<-Data$GV2_2

P<-data.frame(PA=A5,PB=A6)
cor(log2(P))

#PA        PB
#PA 1.0000000 0.9190987
#PB 0.9190987 1.0000000

b <- ggplot(P, aes(x = PA, y = PB))+
  geom_point(color="#A00000")+
  geom_text(x=2.1,y=-1.7,aes(label="italic(R)==0.919"),parse=T,size=16)+
  geom_text(x=-0.3,y=3.95,aes(label="Cell 3#"),parse=F,size=16)+
  scale_x_log10(limits=c(0.01,20000))+scale_y_log10(limits=c(0.01,20000))+
  labs(x='Part A',y='Part B')+
  theme_classic(base_family="ARL")+
  coord_fixed(ratio=1)+
  theme(panel.background=element_rect(fill="white",colour="black",size=0.25),
        panel.border=element_rect(fill=NA,color="black",size=3,linetype="solid"),
        axis.line=element_line(colour="black",size=0),
        axis.ticks=element_line(colour="black",size=1.3),
        axis.ticks.length=unit(0.5,'line'),
        axis.title.x=element_text(size=32,face="plain",color="black",vjust=-0.5),
        axis.title.y=element_text(size=32,face="plain",color="black",vjust=1.5),
        axis.text.x=element_text(size=25,face="plain",color="black"),
        axis.text.y=element_text(size=25,face="plain",color="black"),
        legend.text=element_text(size=10,face="plain",colour="black"),
        legend.position=c(0.85,0.92),
        legend.title=element_blank(),
        legend.margin=margin(0.05,0.05,0.05,0.05,'cm'),
        legend.key.size=unit(0.5,'cm'))

ggsave("Figure S2-T3.png",b,width=6,height=6,dpi=1200)

#Proteome
#Cell4
Data<-read.csv("Cell 4.csv",sep=",",header=T)

A1<-Data$X44_1
A2<-Data$X44_2

P<-data.frame(PA=A1,PB=A2)
cor(log2(P))

#PA        PB
#PA 1.0000000 0.9106018
#PB 0.9106018 1.0000000

b <- ggplot(P, aes(x = PA, y = PB))+
  geom_point(color="#A00000")+
  geom_text(x=3.5,y=0.3,aes(label="italic(R)==0.911"),parse=T,size=16)+
  geom_text(x=1.4,y=4.95,aes(label="Cell 4#"),parse=F,size=16)+
  scale_x_log10(limits=c(1,200000))+scale_y_log10(limits=c(1,200000))+
  labs(x='Part A',y='Part B')+
  theme_classic(base_family="ARL")+
  coord_fixed(ratio=1)+
  theme(panel.background=element_rect(fill="white",colour="black",size=0.25),
        panel.border=element_rect(fill=NA,color="black",size=3,linetype="solid"),
        axis.line=element_line(colour="black",size=0),
        axis.ticks=element_line(colour="black",size=1.3),
        axis.ticks.length=unit(0.5,'line'),
        axis.title.x=element_text(size=32,face="plain",color="black",vjust=-0.5),
        axis.title.y=element_text(size=32,face="plain",color="black",vjust=1.5),
        axis.text.x=element_text(size=25,face="plain",color="black"),
        axis.text.y=element_text(size=25,face="plain",color="black"),
        legend.text=element_text(size=10,face="plain",colour="black"),
        legend.position=c(0.85,0.92),
        legend.title=element_blank(),
        legend.margin=margin(0.05,0.05,0.05,0.05,'cm'),
        legend.key.size=unit(0.5,'cm'))

ggsave("Figure S2-P1.png",b,width=6,height=6,dpi=1200)

#Cell5
Data<-read.csv("Cell 5.csv",sep=",",header=T)

A3<-Data$X45_1
A4<-Data$X45_2

P<-data.frame(PA=A3,PB=A4)
cor(log2(P))

#PA        PB
#PA 1.0000000 0.9184416
#PB 0.9184416 1.0000000

b <- ggplot(P, aes(x = PA, y = PB))+
  geom_point(color="#A00000")+
  geom_text(x=3.5,y=0.3,aes(label="italic(R)==0.918"),parse=T,size=16)+
  geom_text(x=1.4,y=4.95,aes(label="Cell 5#"),parse=F,size=16)+
  scale_x_log10(limits=c(1,200000))+scale_y_log10(limits=c(1,200000))+
  labs(x='Part A',y='Part B')+
  theme_classic(base_family="ARL")+
  coord_fixed(ratio=1)+
  theme(panel.background=element_rect(fill="white",colour="black",size=0.25),
        panel.border=element_rect(fill=NA,color="black",size=3,linetype="solid"),
        axis.line=element_line(colour="black",size=0),
        axis.ticks=element_line(colour="black",size=1.3),
        axis.ticks.length=unit(0.5,'line'),
        axis.title.x=element_text(size=32,face="plain",color="black",vjust=-0.5),
        axis.title.y=element_text(size=32,face="plain",color="black",vjust=1.5),
        axis.text.x=element_text(size=25,face="plain",color="black"),
        axis.text.y=element_text(size=25,face="plain",color="black"),
        legend.text=element_text(size=10,face="plain",colour="black"),
        legend.position=c(0.85,0.92),
        legend.title=element_blank(),
        legend.margin=margin(0.05,0.05,0.05,0.05,'cm'),
        legend.key.size=unit(0.5,'cm'))

ggsave("Figure S2-P2.png",b,width=6,height=6,dpi=1200)

#Cell6
Data<-read.csv("Cell 6.csv",sep=",",header=T)

A5<-Data$X48_1
A6<-Data$X48_2

P<-data.frame(PA=A5,PB=A6)
cor(log2(P))

#PA        PB
#PA 1.0000000 0.9400282
#PB 0.9400282 1.0000000

b <- ggplot(P, aes(x = PA, y = PB))+
  geom_point(color="#A00000")+
  geom_text(x=3.5,y=0.3,aes(label="italic(R)==0.939"),parse=T,size=16)+
  geom_text(x=1.4,y=4.95,aes(label="Cell 6#"),parse=F,size=16)+
  scale_x_log10(limits=c(1,200000))+scale_y_log10(limits=c(1,200000))+
  labs(x='Part A',y='Part B')+
  theme_classic(base_family="ARL")+
  coord_fixed(ratio=1)+
  theme(panel.background=element_rect(fill="white",colour="black",size=0.25),
        panel.border=element_rect(fill=NA,color="black",size=3,linetype="solid"),
        axis.line=element_line(colour="black",size=0),
        axis.ticks=element_line(colour="black",size=1.3),
        axis.ticks.length=unit(0.5,'line'),
        axis.title.x=element_text(size=32,face="plain",color="black",vjust=-0.5),
        axis.title.y=element_text(size=32,face="plain",color="black",vjust=1.5),
        axis.text.x=element_text(size=25,face="plain",color="black"),
        axis.text.y=element_text(size=25,face="plain",color="black"),
        legend.text=element_text(size=10,face="plain",colour="black"),
        legend.position=c(0.85,0.92),
        legend.title=element_blank(),
        legend.margin=margin(0.05,0.05,0.05,0.05,'cm'),
        legend.key.size=unit(0.5,'cm'))

ggsave("Figure S2-P3.png",b,width=6,height=6,dpi=1200)
