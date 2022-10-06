
library("limma")  
setwd('D:\\广东医大创课项目\\何昊洋\\immunity')
expFile="symbol.LGG.txt"  
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]

#数据转换
v=voom(data, plot=F, save.plot=F)
out=v$E
out=rbind(ID=colnames(out), out)
write.table(out,file="uniq.symbol.txt",sep="\t",quote=F,col.names=F)     

#运行CIBERSORT，得到免疫细胞浸润的结果
source("Gene25.CIBERSORT.R")
results=CIBERSORT("ref.txt", "uniq.symbol.txt", perm=1000, QN=TRUE)


#引用包
library(limma)
library(reshape2)
library(ggpubr)
library(vioplot)
library(ggExtra)
gene='TFRC'
pFilter=0.05            #免疫细胞浸润结果的过滤条件


data=read.csv('tcgaclinic.csv')
data$X<-NULL
rownames(data)<-data$ID
data$ID<-NULL
#根据目标基因表达量对样品进行分组

data$TFRC_group=ifelse(data[,gene]>median(data[,gene]), "High", "Low")

#读取免疫细胞结果文件，并对数据进行整理

immune=results[results[,"P-value"]<pFilter,]
immune=as.matrix(immune[,1:(ncol(immune)-3)])

#删除正常样品
group=sapply(strsplit(row.names(immune),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
immune=immune[group==0,]
row.names(immune)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", row.names(immune))
immune=avereps(immune)

#数据合并
sameSample=intersect(row.names(immune), row.names(data))
rt=cbind(immune[sameSample,,drop=F], data[sameSample,,drop=F])


##################绘制箱线图##################
#把数据转换成ggplot2输入文件
data=rt[,-(ncol(rt)-1)]
data<-data[,-c(23:29)]
data=melt(data,id.vars=c("TFRC_group"))
colnames(data)=c("gene", "Immune", "Expression")
#绘制箱线图
group=levels(factor(data$gene))
data$gene=factor(data$gene, levels=c("Low","High"))
bioCol=c("#0066FF","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length(group)]
boxplot=ggboxplot(data, x="Immune", y="Expression", fill="gene",
                  xlab="",
                  ylab="Fraction",
                  legend.title=gene,
                  width=0.8,
                  palette=bioCol)+
  rotate_x_text(50)+
  stat_compare_means(aes(group=gene),symnum.args=list(cutpoints=c(0, 0.001, 0.01, 0.05, 1), symbols=c("***", "**", "*", "")), label="p.signif")
#输出图片
pdf(file="immune.diff.pdf", width=7, height=6)
print(boxplot)
dev.off()


##########绘制相关性散点图##########
outTab=data.frame()
for(i in colnames(rt)[1:(ncol(rt)-2)]){
  x=as.numeric(rt[,gene])
  y=as.numeric(rt[,i])
  if(sd(y)==0){y[1]=0.00001}
  cor=cor.test(x, y, method="spearman")
  outVector=cbind(Cell=i, cor=cor$estimate, pvalue=cor$p.value)
  outTab=rbind(outTab,outVector)
  if(cor$p.value<0.05){
    outFile=paste0("cor.", i, ".pdf")
    df1=as.data.frame(cbind(x,y))
    p1=ggplot(df1, aes(x, y)) + 
      xlab(paste0(gene, " expression")) + ylab(i)+
      geom_point() + geom_smooth(method="lm",formula = y ~ x) + theme_bw()+
      stat_cor(method = 'spearman', aes(x =x, y =y))
    p2=ggMarginal(p1, type="density", xparams=list(fill = "orange"), yparams=list(fill = "blue"))
    #相关性图形
    pdf(file=outFile, width=5.2, height=5)
    print(p2)
    dev.off()
  }
}
#输出相关性的结果文件
write.table(outTab,file="cor.result.txt",sep="\t",row.names=F,quote=F)
