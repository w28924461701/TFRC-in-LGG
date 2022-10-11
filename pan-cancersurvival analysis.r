

#install.packages("survival")
#install.packages("survminer")

library(survival)
library(survminer)
library(dplyr)
setwd("D:\\panCancer\\08.survival")                     #???ù???Ŀ¼
rt=read.table("expTime.txt",header=T,sep="\t",check.names=F,row.names=1)       #??ȡ?????ļ?
rt$futime=rt$futime/365
gene=colnames(rt)[3]
pFilter=0.05            #km????pvalue????????

#?????????ͽ???ѭ??
for(i in levels(as.factor(rt$CancerType))) {
  rt1=rt[(rt[,"CancerType"]==i),]
 # group=ifelse(rt1[,gene]>median(rt1[,gene]),"high","low")
  
  # rt1$group=case_when(rt1[,gene]<=unname(quantile(rt1[,gene],0.25))~ 'low',  
  #               rt1[,gene]>unname(quantile(rt1[,gene],0.25))&rt1[,gene]<unname(quantile(rt1[,gene],0.75))~'median',
              #  rt1[,gene]>=unname(quantile(rt1[,gene],0.75))~'high')
    #    rt1<-filter(rt1,!(group=='median'))
  
#  rt1<-filter(rt1,futime<10)
  
 res.cut <- surv_cutpoint(rt1, time = "futime", event = "fustat",  variables = c("TFRC"))  #???????о??Ļ???

  group=ifelse(rt1[,gene]> res.cut[["cutpoint"]][["cutpoint"]],"high","low")

  
  diff=survdiff(Surv(futime, fustat) ~group,data = rt1)
  pValue=1-pchisq(diff$chisq,df=1)
  if(pValue<pFilter){
    if(pValue<0.001){
      pValue="p<0.001"
    }else{
      pValue=paste0("p=",sprintf("%.03f",pValue))
    }
    fit <- survfit(Surv(futime, fustat) ~ group, data = rt1)
    #????????????
    surPlot=ggsurvplot(fit, 
                       data=rt1,
                       title=paste0("Cancer: ",i),
                       pval=pValue,
                       pval.size=6,
                       legend.labs=c("high","low"),
                       legend.title=paste0(gene," levels"),
                       font.legend=12,
                       xlab="Time(years)",
                       ylab="Overall survival",
                       break.time.by = 1,
                       palette=c("red","blue"),
                       conf.int=F,
                       fontsize=4,
                       risk.table=TRUE,
                       risk.table.title="",
                       risk.table.height=.25)
    pdf(file=paste0("survivalq.",i,".pdf"),onefile = FALSE,
        width = 6,             #ͼƬ?Ŀ???
        height =5)             #ͼƬ?ĸ߶?
    print(surPlot)
    dev.off()
  }
}


