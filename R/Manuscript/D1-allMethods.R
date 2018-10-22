setwd("/Users/maryamrabiee/Documents/Research/ASTRAL-multiind/results/ASTRAL-multiind/ASTRAL-multi paper")
d<-read.csv("all.stat",sep=" ",header=F)
head(d)
require(ggplot2)
d$V3 <- as.factor(d$V3)
d$V2 <- as.factor(d$V2)
levels(d$V3)
levels(d$V2)
levels(d$V4)
#levels(d$V3) <- list("50 genes" = "50","200 genes" = "200","500 genes" = "500", "1000 genes" = "1000")
ggplot(data=d,aes(x=as.factor(V3), y= V8,color=factor(interaction(V5,V4),c("ASTRAL.original","ASTRAL.contracted","ASTRAL.half","Njst.original","Njst.half"))))+xlab("Number of genes")+ylab("FN rate")+
  stat_summary(aes(group=as.factor(interaction(V4,V5))),geom="line",linetype=1)+
  stat_summary(aes(group=as.factor(interaction(V4,V5))),linetype=1,geom = "errorbar",width=0.2)+
  scale_color_manual(values = c("#045a8d", "#43a2ca", "#7bccc4","#ffeda0","#feb24c"),labels=c("ASTRAL-multi","ASTRAL-multi-5%","ASTRAL-half","NJst","NJst-half"))+
  theme_bw()+facet_wrap(~V2,scales = "free_y")+theme(legend.title=element_blank())

ggsave("figures/Allmethods_D1_FN.pdf",width = 12.2 ,height = 4.8)
#  theme(legend.position=c(0.9,0.5))


#---------------------------------------Figure 5-----------------------------

d$method<- paste(d$V4,d$V5,sep="-")
head(d)
dcast(V2~V3,data=d[d$method=="original-ASTRAL",],fun.aggregate = function(x) sum(x!=0))
fd=d
fd$V2 <-as.factor(fd$V2)
levels(fd$V2)
levels(fd$V2) <-list("0.5M"="500000","1M"="1000000","2M"="2000000")
ggplot(data=fd[fd$method%in% c("original-Njst","original-ASTRAL","contracted-ASTRAL"),],
       aes(x=as.factor(V3), y= V8,color=method))+
  xlab("Number of genes")+ylab("Species tree error (RF distance)")+
  stat_summary(aes(group=as.factor(interaction(V4,V5))),geom="line",linetype=1)+
  stat_summary(aes(group=as.factor(interaction(V4,V5))),linetype=1,geom = "errorbar",width=0.2)+
  scale_color_manual(values = c("#045a8d", "#43a2ca","#feb24c"),labels=c("ASTRAL-multi-5%","ASTRAL-multi","NJst original"))+
  theme_bw()+facet_wrap(~V2,scales = "free_y")+theme(legend.title=element_blank(),legend.margin = margin(0,0,0,0))

ggsave("figures/Allmethods_D1_main_FN.pdf",width = 8 ,height = 2.8)
# stat_summary(geom="errorbar",fun.ymin=function(x) {mean(x)-sd(x)/sqrt(length(x))},
#               fun.ymax = function(x) {mean(x)+sd(x)/sqrt(length(x))},width=0.2,linetype=2)+
  
#-----------------------------------------------------------------------------------------------------
  
  ggplot(data=d,aes(x=as.factor(V3), y= V8,fill=factor(interaction(V5,V4),c("ASTRAL.original","ASTRAL.contracted","ASTRAL.half","Njst.original","Njst.half"))))+xlab("Number of genes")+ylab("FN rate")+
    geom_boxplot(outlier.alpha = 0.3)+theme_bw()+facet_wrap(~V2,scales = "free_y")+
    scale_fill_manual(values = c("#045a8d", "#43a2ca", "#7bccc4","#fff7bc","#fec44f"),labels=c("ASTRAL-multi","ASTRAL-multi-5%","ASTRAL-half","NJst","NJst-half"))+
    theme_bw()+facet_wrap(~V2,scales = "free_y")+theme(legend.title=element_blank())
  
  ggsave("figures/Allmethods_D1_bar_FN.pdf",width = 12.2 ,height = 4.8)  
  
  head(d)
  nrow(d)
  nrow(d[d$V2 == 500000,])
  summary(aov(V8~V5*(V2+V3),data=d[d$V2 == 500000  & (d$V4 == "contracted" | (d$V5 == "Njst" & d$V4 == "original")),]))
  summary(aov(V8~V5*(V2+V3),data=d[d$V2 != 500000  & (d$V4 == "contracted" | (d$V5 == "Njst" & d$V4 == "original")),]))
 
  summary(aov(V8~V5*V3,data=d[d$V2 == 500000 & d$V3 != 1000 & (d$V4 == "contracted" | (d$V5 == "Njst" & d$V4 == "original")),]))   
  
  summary(aov(V8~V5,data=d[d$V2 == 500000 & d$V3 != 1000 & (d$V4 == "contracted" | (d$V5 == "Njst" & d$V4 == "original")),]))   
  
  summary(aov(V8~V5*V3,data=d[d$V2 != 500000 & d$V3 != 1000 & (d$V4 == "contracted" | (d$V5 == "Njst" & d$V4 == "original")),]))
  
  summary(aov(V8~V5*V3,data=d[d$V2 == 2000000 & d$V3 != 1000 & (d$V4 == "contracted" | (d$V5 == "Njst" & d$V4 == "original")),]))
  t.test(d[d$V2 == 500000 & d$V3 != 1000 & (d$V4 == "contracted"),"V8"]-d[d$V2 == 500000 & d$V3 != 1000 & ((d$V5 == "Njst" & d$V4 == "original")),"V8"])
  t.test(d[d$V2 != 500000 & d$V3 != 1000 & (d$V4 == "contracted"),"V8"]-d[d$V2 != 500000 & d$V3 != 1000 & ((d$V5 == "Njst" & d$V4 == "original")),"V8"])
  