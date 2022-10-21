library(ggplot2)
library(ggpubr)
rt<-read.csv('clipboard',sep = '\t',check.names = F)
p <- ggplot(rt1, aes(x=Group, y=`TFRC-IOD/Sum`,fill=Group)) + 
  geom_boxplot()
p +scale_fill_manual(values=c("#E69F00", "#56B4E9"))+ theme_classic()+
  geom_dotplot(binaxis='y', stackdir='center',
               position=position_dodge(0.2))+ theme(legend.position="top")+
  theme(axis.text=element_text(size=14)
        ,axis.title.x=element_text(size=16),axis.title.y=element_text(size=14))+
  stat_compare_means()

