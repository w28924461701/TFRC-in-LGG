#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("reshape2")
#install.packages("ggpubr")


#引用包
library(limma)
library(reshape2)
library(ggpubr)

#读取表达数据文件
exp<-read.csv('tcgaclinic.csv')
exp$TFRC_group<-ifelse(exp$TFRC<median(exp$TFRC),'low','high')

#读取肿瘤微环境打分文件，并对数据进行整理
score=read.table(scoreFile, header=T, sep="\t", check.names=F, row.names=1)
score=score[,1:3]

#样品取交集
sameSample=intersect(row.names(exp), row.names(score))
exp=exp[sameSample,"TFRC_group",drop=F]
score=score[sameSample,,drop=F]
rt=cbind(score, exp)
rt$TFRC_group=factor(rt$TFRC_group, levels=c("low", "high"))

#将合并后的数据转换为ggplot2的输入文件
data=melt(rt, id.vars=c("TFRC_group"))
colnames(data)=c("Type", "scoreType", "Score")


#输出图

p1<-ggplot(data,aes(x=scoreType,y=Score,fill=Type))+geom_boxplot()+
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  theme_classic()+ theme(legend.position="top")+
  geom_jitter(shape=10, position=position_jitter(0.15))+
  stat_compare_means(aes(group=Type), method="wilcox.test",symnum.args=
   list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
       symbols = c("***", "wilcox.test,p<0.01", "*", " ")),label = "p.signif")+
  theme(
    axis.text =element_text(size =10,color = 'black'), 
    axis.line.x=element_line(size = 0.8),
    axis.line.y=element_line(size = 0.8))
       
pdf(file="ggplot.pdf", width=6, height=5)
print(p1)
dev.off()
  
