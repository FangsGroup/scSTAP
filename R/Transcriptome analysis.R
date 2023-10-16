#Transcriptome analysis
windowsFonts(ARL = windowsFont("Arial"))

#set work path
setwd("D:/Research/R/Data")

#Matq-seq vs Smart-seq
library(ggplot2)
library(ggsignif)
library(rstatix)
library(ggpubr)
library(ggtext)

#Bar-Figure 3A
P<-read.csv("bar-MvsS.csv",sep=",",header=T)

p<-ggplot()+
  geom_col(data=P,aes(x=X,y=Number,fill=Group), width=0.8)+
  scale_y_continuous(limits=c(0,30000))+
  scale_x_discrete(labels=P$Group)+
  coord_fixed(ratio=0.000425)+
  labs(x='Single oocytes',y='Genes')+
  scale_fill_manual(values = c("#CE5B52","#18609C","#FFFFFF"),labels= c("MATQ-seq","Smart-seq2","H"))+
  theme(panel.background=element_rect(fill="white",colour="black",size=0.01),
        panel.border=element_rect(fill=NA,color="black",size=1.6,linetype="solid"),
        axis.line=element_line(colour="black",size=0),
        axis.ticks=element_line(colour="black",size=0.7),
        axis.ticks.length=unit(0.5,'line'),
        axis.title.x=element_text(size=22,face="plain",color="black",vjust=-0.2),
        axis.title.y=element_text(size=22,face="plain",color="black",vjust=2.0),
        axis.text=element_text(size=18,face="plain",color="black"),
        legend.text=element_text(size=18,face="plain",colour="black"),
        legend.position=c(0.75,0.88),
        legend.title=element_blank(),
        legend.margin=margin(0.01,0.01,0.01,0.01,'cm'),
        legend.key.size=unit(0.7,'cm'))+
  guides(fill=guide_legend(nrow=2,byrow=T))

ggsave("Bar-MvsS.png",p,width=5,height=5,dpi=1200)

#Box-Figure 3B
P<-read.csv("box-MvsS.csv",sep=",",header=T)

df_p_val1 <- P%>% 
  group_by(Group) %>% 
  wilcox_test(formula = Number ~ Type) %>% 
  add_significance(p.col = 'p',cutpoints = c(0,0.001,0.01,0.05,1),symbols = c('***','**','*','n.s.')) %>% 
  add_xy_position(x='Group')

df_p_val1

p<-ggplot()+
  geom_boxplot(P,mapping=aes(x=Group,y=Number, fill=Type),width=0.8)+
  scale_fill_manual(values=c("#CE5B52","#18609C"))+
  stat_pvalue_manual(df_p_val1,label="{p.signif}",tip.length=0.02,size=6)+
  labs(x=element_blank(),y='Genes')+
  scale_y_continuous(limits=c(0,18000))+
  coord_fixed(ratio=0.0001615)+
  theme(panel.background=element_rect(fill="white",colour="black",size=0.25),
        panel.border=element_rect(fill=NA,color="black",size=1.6,linetype="solid"),
        axis.line=element_line(colour="black",size=0),
        axis.ticks=element_line(colour="black",size=0.6),
        axis.ticks.length=unit(0.5,'line'),
        axis.title.x=element_text(size=22,face="plain",color="black",vjust=-0.7),
        axis.title.y=element_text(size=22,face="plain",color="black",vjust=2),
        axis.text.x=element_text(size=18,face="plain",color="black",vjust=0.9,angle=15,hjust=0.85),
        axis.text.y=element_text(size=18,face="plain",color="black"),
        legend.text=element_text(size=18,face="plain",colour="black"),
        legend.title=element_blank(),
        legend.margin=margin(0.1,0.1,0.1,0.1,'cm'),
        legend.key.size=unit(0.7,'cm'),
        legend.key=element_rect(fill="white"),
        legend.position=c(0.72,0.87))+
  guides(fill=guide_legend(nrow=2,byrow=T))

ggsave("box-MvsS.png",p,width=5,height=5,dpi=1200)

#Matq-seq
library(ggplot2)
library(VennDiagram)
library(corrplot)
library(PerformanceAnalytics)
library(GGally)
library(ggplot2)

#Line-Figure 3C
P<-read.csv("coverage.csv",sep=",",header=T)

p<-ggplot()+
  geom_line(data=P,aes(x=Percentile,y=Average),size=0.7,color="black")+
  geom_linerange(data=P,aes(x=Percentile,ymin=Average-SD,ymax=Average+SD),size=3, colour="#C03030",alpha=0.3)+
  scale_y_continuous(limits=c(0,1.5))+
  scale_x_continuous(limits=c(0,105))+
  coord_fixed(ratio=70)+
  labs(x="Gene body percentile (5' to 3')",y="Mapped reads %")+
  theme(panel.background=element_rect(fill="white",colour="black",size=0.01),
        panel.border=element_rect(fill=NA,color="black",size=1.7,linetype="solid"),
        axis.line=element_line(colour="black",size=0),
        axis.ticks=element_line(colour="black",size=1),
        axis.ticks.length=unit(0.5,'line'),
        axis.title.x=element_text(size=22,face="plain",color="black",vjust=-0.5,hjust=0.9),
        axis.title.y=element_text(size=22,face="plain",color="black",vjust=2.0),
        axis.text=element_text(size=20,face="plain",color="black"),
        legend.text=element_text(size=16,face="plain",colour="black"),
        legend.position=c(0.10,0.85),
        legend.title=element_blank(),
        legend.margin=margin(0.01,0.01,0.01,0.01,'cm'),
        legend.key.size=unit(0.7,'cm'))+
  guides(fill=guide_legend(nrow=2,byrow=T))

ggsave("Line-matq.png",p,width=5,height=5,dpi=1200)

#Corr-Figure 3D
RNA1<-read.csv("E1_exon.csv",sep=",",header=T)
RNA2<-read.csv("E2_exon.csv",sep=",",header=T)
RNA3<-read.csv("E3_exon.csv",sep=",",header=T)
RNA4<-read.csv("Y6_exon.csv",sep=",",header=T)
RNA5<-read.csv("Y9_exon.csv",sep=",",header=T)
RNA6<-read.csv("Y10_exon.csv",sep=",",header=T)

RNA1<-data.frame(gene_id=RNA1$gene_id,RNA1=RNA1$E1_CPM)
RNA2<-data.frame(gene_id=RNA2$gene_id,RNA2=RNA2$E2_CPM)
RNA3<-data.frame(gene_id=RNA3$gene_id,RNA3=RNA3$E3_CPM)
RNA4<-data.frame(gene_id=RNA4$gene_id,RNA4=RNA4$Y6_CPM)
RNA5<-data.frame(gene_id=RNA5$gene_id,RNA5=RNA5$Y9_CPM)
RNA6<-data.frame(gene_id=RNA6$gene_id,RNA6=RNA6$Y10_CPM)

RNA<-merge(RNA1,RNA2,by.x='gene_id',by.y='gene_id')
RNA<-merge(RNA,RNA3,by.x='gene_id',by.y='gene_id')
RNA<-merge(RNA,RNA4,by.x='gene_id',by.y='gene_id')
RNA<-merge(RNA,RNA5,by.x='gene_id',by.y='gene_id')
RNA<-merge(RNA,RNA6,by.x='gene_id',by.y='gene_id')

RNA<-data.frame(RNA1=RNA$RNA1,
                RNA2=RNA$RNA2,
                RNA3=RNA$RNA3,
                RNA4=RNA$RNA4,
                RNA5=RNA$RNA5,
                RNA6=RNA$RNA6)

cor_func <- function(data, mapping, method, symbol, ...){
  
  x <- eval_data_col(data, mapping$x)
  y <- eval_data_col(data, mapping$y)
  corr <- cor(x, y, method=method,use='complete.obs')
  
  ggally_text(mapping=aes(),label = paste(symbol, as.character(sprintf("%0.3f",corr))),size=8,xP=0.45,yP=0.5,color="black")+
    theme(panel.grid.major=element_line(colour=NA),panel.background = element_rect(fill = "white"))
}

p<-ggpairs(log2(RNA+1), showStrips=T,
           lower=list(continuous=wrap("points", pch=21,size=3,stroke=0.15,fill="#A00000",color="white",alpha=0.2)),
           diag=list(continuous=wrap("densityDiag",color="#D00000")),
           upper=list(continuous=wrap(cor_func,method="pearson",symbol=expression(' ***\n'))))

p[6,1]<-p[6,1]+scale_y_continuous(limits=c(0,13))
p[6,2]<-p[6,2]+scale_y_continuous(limits=c(0,13))
p[6,3]<-p[6,3]+scale_y_continuous(limits=c(0,13))
p[6,4]<-p[6,4]+scale_y_continuous(limits=c(0,13))
p[6,5]<-p[6,5]+scale_y_continuous(limits=c(0,13))
p[6,6]<-p[6,6]+scale_x_continuous(limits=c(0,13))

p<-p+
  theme(panel.background=element_rect(fill="white",colour="black",size=0.5),
        panel.border=element_rect(fill=NA,color="black",size=2,linetype="solid"),
        axis.text=element_text(colour="black",size=18),
        strip.background=element_rect(fill="white"),
        strip.text=element_text(colour="black",size=20,face="plain"))

ggsave("corr-matq.png",p,width=7,height=7,dpi=600)
