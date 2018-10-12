require(ggplot2)
require(reshape2)
setwd('/Users/erfan/Main/Library/ASTRAL-multiind/data/distique-results')
d<-read.csv("results.D1.txt",sep=" ",header=F)
ggplot(data=d,aes(x=as.factor(V3),y=V6))+geom_boxplot()+theme_bw()

d2<-read.csv("comparison.D2.txt",sep=" ",header=F)
d2$V1<-factor(d2$V1,levels=c("model.200-5.500000.0.000001","model.200-5.1000000.0.000001","model.200-5.2000000.0.000001"))
ggplot(data=d2,aes(x=as.factor(V4),y=V7,fill=as.factor(V3)))+geom_boxplot()+facet_wrap(~V1)+theme_bw()+
  scale_fill_brewer(name="",palette = "Dark2")+theme(legend.position = "bottom")

