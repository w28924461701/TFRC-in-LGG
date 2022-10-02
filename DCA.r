library(survival)
rt<-na.omit(rt)
Srv = Surv(rt$futime, rt$fustat)
coxmod1 = coxph(Srv ~ TFRC_group+Diagnosis_Age+Histologic_Grade+Subtype, data=rt)
rt$model1 = c(1 - (summary(survfit(coxmod1,newdata=rt), times=5)$surv))

coxmod2 = coxph(Srv ~ TFRC_group, data=rt)
rt$model2 = c(1 - (summary(survfit(coxmod2,newdata=rt), times=5)$surv))

coxmod3 = coxph(Srv ~TFRC_group+Histologic_Grade+Subtype, data=rt)
rt$model3 = c(1 - (summary(survfit(coxmod3,newdata=rt), times=5)$surv))

coxmod4 = coxph(Srv ~TFRC_group+Subtype, data=rt)
rt$model4 = c(1 - (summary(survfit(coxmod4,newdata=rt), times=5)$surv))

head(rt)
source("stdca.R") #stdca.R文件位于当前文件夹
mod1<-stdca(data=rt, outcome="fustat", ttoutcome="futime", timepoint=5, 
            predictors="model1", cmprsk=TRUE, smooth=TRUE, xstop=0.5,intervention="FALSE")

pdf("net_benefit.pdf",width = 6,height = 6)
stdca(data=rt, outcome="fustat", ttoutcome="futime", timepoint=5,
      predictors=c("model1","model2","model3","model4"), 
      cmprsk=TRUE, smooth=TRUE, 
      xstop=0.5,intervention="FALSE")
dev.off()

mod1<-stdca(data=rt,outcome="fustat", ttoutcome="futime", timepoint=5,
            predictors="model1", cmprsk=TRUE, smooth=TRUE, xstop=0.5,intervention="TRUE")
pdf("net_reduction1.pdf",width = 6,height = 6)
stdca(data=rt,outcome="fustat", ttoutcome="futime", timepoint=5,
      predictors=c("model1","model2","model3","model4"), 
      cmprsk=TRUE, smooth=TRUE, 
      xstop=0.5,intervention="TRUE")
dev.off()
library(ggplotify)
library(magick)
library(cowplot)
fnames<-Sys.glob("net_*.pdf")

p<-lapply(fnames,function(i){
  pn<-as.ggplot(image_read_pdf(i))
})
pdf("ne.pdf",width = 6,height = 8)
plot_grid(plotlist = p, ncol=2,
          #labels = c("(A)","(B)"),
          label_size = 25, #A和B字体大小
          label_y = 0.75 #A和B的位置，默认值为1，太高
)
dev.off()
