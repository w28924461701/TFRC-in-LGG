cli<-read.csv('clipboard',sep = '\t')
library(ggplot2)
library(tidyverse)
library(gapminder)
library(reshape2)
short2long = melt(cli)
short2long<-na.omit(short2long)
bp <- ggplot(short2long, aes(x=variable, y=value, fill=variable)) + 
  geom_boxplot()+labs(x="TFRC methylation",
       y =rt1<-rt1[,c(1:4,6,9,13,16)] "Methylation level(beta value)")+
  geom_point(position=position_jitterdodge(),alpha=0.13)+theme_bw(base_size = 13)+
  theme(axis.text.x = element_blank())


r<-read.csv('tcgaclinic.csv')
rt<-merge(r,cli)
library(GGally)
library(ggplot2)
library(wooldridge)
ggcorr(rt1,palette = "RdBu")

ggpairs(rt1, title="correlogram with ggpairs()") 
