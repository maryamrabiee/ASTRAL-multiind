require(ggplot2)
require(reshape2)
setwd('/Users/erfan/Main/Library/ASTRAL-multiind')


d<-read.csv('data/newdataset/semi.internal.multi.stat',sep=" ",header=T)
ds<-melt(data=d,id.vars = c("Replicate"))
ds$version<-"astral-multi"
ds[ds$variable %in% c("internal","semiterminal"),]$version <- "astral-II"

p<-read.csv('data/parameter.log.info1',sep=" ",header=T)
m<-merge(x=ds,y=p,by.x="Replicate",by.y="Replicate")
a=10
m$generation<-cut(ecdf(m$Generations)(m$Generations),breaks = 0:a*1/a)
ggplot(data=m,
       aes(x=generation,y=value,color=variable,group=variable))+stat_summary(fun.y="mean",geom="line")+
  stat_summary(fun.y="mean",geom="point")+
  stat_summary(fun.ymax = function(x){mean(x)+sd(x)/sqrt(length(x))},
               fun.ymin=function(x){mean(x)-sd(x)/sqrt(length(x))},geom="errorbar",width=0.3)+theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        legend.position = c(0.8,0.8))+scale_color_brewer(name="",palette = "Set1",labels=c("ASTRAL-III (internal)","ASTRAL-III (semiinternal)","ASTRAL-multiind"))

