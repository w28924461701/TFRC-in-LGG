
#引用包
library(limma)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(corrplot)

pFilter=0.001             #相关性检验pvalue的过滤条件
geneName="TFRC"           #目标基因名字
expFile="symbol.LGG.txt"      #表达数据文件
geneFile="gene.txt"       #免疫检查点的基因列表文件
setwd("D:/广东医大创课项目/何昊洋/immunity")     #设置工作目录

#读取基因表达文件,并对数据进行处理
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)

#读取基因文件,获取免疫检查点相关基因的表达量
gene=read.table(geneFile, header=F, sep="\t", check.names=F)
sameGene=intersect(row.names(data), as.vector(gene[,1]))
data=t(data[c(geneName, sameGene),])
data=log2(data+1)

#删除正常样品
group=sapply(strsplit(row.names(data),"\\-"),"[",4)
group=sapply(strsplit(group,""),"[",1)
group=gsub("2","1",group)
data=data[group==0,]
row.names(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", row.names(data))
data=t(avereps(data))

#对免疫检查点基因进行循环，找出与目标基因具有相关的免疫检查点
x=as.numeric(data[geneName,])
outTab=data.frame()
for(i in sameGene){
  if(i==geneName){next}
  y=as.numeric(data[i,])
  corT=cor.test(x, y, method = 'pearson')
  cor=corT$estimate
  pvalue=corT$p.value
  if(pvalue<pFilter){
    outTab=rbind(outTab, cbind(Query=geneName, Gene=i, cor, pvalue))
  }
}
#输出相关性结果文件
write.table(file="corResult.txt", outTab, sep="\t", quote=F, row.names=F)

#相关性矩阵
data=t(data[c(geneName, as.vector(outTab[,2])),])
M=cor(data)



#绘制相关性图形
pdf(file="corpotCHECKPOINT.pdf",width=7,height=7)
corrplot(M,
         method = "circle",
         order = "original",
         type = "upper",
         col=colorRampPalette(c("green", "white", "red"))(50)
)
dev.off()

pdf(file="corpotCHECKPOINT2.pdf",width=8,height=8)
corrplot(M,
         order="original",
         method = "color",
         number.cex = 0.7,
         addCoef.col = "black",
         diag = TRUE,
         tl.col="black",
         col=colorRampPalette(c("blue", "white", "red"))(50))
dev.off()

data<-as.data.frame(data)

data$Group_TFRC<-ifelse(data$TFRC<median(data$TFRC),'Low','High')

data$Group_TFRC<-as.factor(data$Group_TFRC)
data$TFRC<-NULL
data1=melt(data, id.vars=c("Group_TFRC"))
data<-data1
p1<-ggplot(data,aes(x=variable,y=value,fill=Group_TFRC))+geom_boxplot()+
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  theme_classic()+ theme(legend.position="top")+
  geom_jitter(shape=10, position=position_jitter(0.15))+
  stat_compare_means(aes(group=Group_TFRC), method="wilcox.test",symnum.args=
                       list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                            symbols = c("***", "**", "*", " ")),label = "p.signif")+
  theme(
    axis.text =element_text(size =10,color = 'black'), 
    axis.line.x=element_line(size = 0.8),
    axis.line.y=element_line(size = 0.8))

pdf(file="ggplotcheckpoint.pdf", width=6, height=5)
print(p1)
dev.off()
