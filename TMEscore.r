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
  

##########绘制相关性散点图##########
outTab=data.frame()
for(i in colnames(rt)[1:(ncol(rt)-1)]){
  x=as.numeric(rt[,'TFRC'])
  y=as.numeric(rt[,i])
  if(sd(y)==0){y[1]=0.00001}
  cor=cor.test(x, y, method="spearman")
  outVector=cbind(Cell=i, cor=cor$estimate, pvalue=cor$p.value)
  outTab=rbind(outTab,outVector)
  if(cor$p.value<0.05){
    outFile=paste0("cor.", i, ".pdf")
    df1=as.data.frame(cbind(x,y))
    p1=ggplot(df1, aes(x, y)) + 
      xlab(paste0('TFRC', " expression")) + ylab(i)+
      geom_point() + geom_smooth(method="lm",formula = y ~ x) + theme_bw()+
      stat_cor(method = 'spearman', aes(x =x, y =y))
    p2=ggMarginal(p1, type="density", xparams=list(fill = "orange"), yparams=list(fill = "blue"))
    #相关性图形
    pdf(file=outFile, width=5.2, height=5)
    print(p2)
    dev.off()
  }
}
