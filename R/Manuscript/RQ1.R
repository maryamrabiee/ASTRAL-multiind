require(ggplot2)
require(reshape2)
setwd('/Users/maryamrabiee/Documents/Research/ASTRAL-multiind/results/ASTRAL-multiind/ASTRAL-multi paper')

#------------------------------------------------Table 2------------------------------------------------

# d<-read.csv('semi.internal.multi.stat',sep=" ",header=T)
df<-read.csv('semi.internal.multi.stat.final',sep=" ",header=T)
# ds<-melt(data=d,id.vars = c("Replicate"))
dfm<-melt(data=df,id.vars = c("Replicate"))
dfm$version<-"astral-multi"
dfm[dfm$variable %in% c("internal","semiterminal"),]$version <- "astral-II"
# ds$version<-"astral-multi"
# ds[ds$variable %in% c("internal","semiterminal"),]$version <- "astral-II"
head(d)
median(df$astral.multi)
median(df$semiterminal)
median(df$internal)
mean(df$astral.multi)
mean(df$semiterminal)
mean(df$internal)
sd(df$astral.multi)
sd(df$semiterminal)
sd(df$internal)
t.test(df$internal,df$astral.multi,paired = T)

#------------------------------------------------------------------------------------------------------------
p<-read.csv('../data/parameter.log.info1',sep=" ",header=T)
m<-merge(x=dfm,y=p,by.x="Replicate",by.y="Replicate")
a=10
m$generation<-cut(ecdf(m$Generations)(m$Generations),
                  breaks = 0:a*1/a,labels=1:a)
m$ILS<-cut(ecdf(m$TrueQscore)(m$TrueQscore),
                  breaks = 0:a*1/a,labels=1:a)

#----------------------------------------------Figure 4-A-----------------------------------------
nrow(m)
head(m)
ggplot(data=m,
       aes(x=ILS,y=value,color=variable,group=variable))+stat_summary(fun.y=mean,geom="line")+
  stat_summary(fun.y="mean",geom="point",size=1)+
  stat_summary(fun.ymax = function(x){mean(x)+sd(x)/sqrt(length(x))},
               fun.ymin=function(x){mean(x)-sd(x)/sqrt(length(x))},
               geom="errorbar",width=0.1,size=.3)+theme_bw()+
  theme(axis.text.x = element_text(angle=45, hjust = 1),
        legend.position = c(0.7,0.75),legend.title = element_blank(), legend.margin = margin(0,0,0,0))+
  scale_x_discrete(labels=function(x) {paste("<",round(100*quantile(m$TrueQscore,0:a*1/a)[as.numeric(x)+1])/100)})+
      #  scale_color_brewer(name="",palette = "Set2",labels=c("ASTRAL-III (internal)","ASTRAL-III (semiinternal)","ASTRAL-multiind"))+
        scale_color_manual(values = c("#74a9cf", "#2b8cbe","#feb24c"),labels=c("ASTRAL-ind (internal)","ASTRAL-ind (semi-terminal)","ASTRAL-multi"))+
        xlab("Quartet score")+ylab(" Species tree error (FN %)")+scale_y_continuous(labels = percent)+
  ggsave("RQ1-ILS.pdf",height = 3.2,width = 4)

#-------------------------------------------------------------------------------------------------
head(d)
ggplot(data=m,aes(x=generation,y=value,color=variable,group=variable))+stat_summary(fun.y=mean,geom="line")+
  stat_summary(fun.y="mean",geom="point")+
  stat_summary(fun.ymax = function(x){mean(x)+sd(x)/sqrt(length(x))},
               fun.ymin=function(x){mean(x)-sd(x)/sqrt(length(x))},
               geom="errorbar",width=0.1,size=.3)+theme_bw()+
  theme(axis.text.x = element_text(angle = 20, hjust = 1),legend.position = c(0.8,0.8),legend.title = element_blank())+
  scale_x_discrete(labels=function(x) {paste("<",round(quantile(m$Generations,0:a*1/a)[as.numeric(x)]))})+
  
  #  scale_color_brewer(name="",palette = "Set2",labels=c("ASTRAL-III (internal)","ASTRAL-III (semiinternal)","ASTRAL-multiind"))+
  scale_color_manual(values = c("#74a9cf", "#2b8cbe","#feb24c"),labels=c("ASTRAL-ind (internal)","ASTRAL- (semi-terminal)","ASTRAL-multi"))+
  xlab("Number of Generations")+ylab("Incorrect branches of species tree (%)")+

ggsave("RQ1-gen.pdf",height = 4.6,width = 7.4)


t.test(md$internal,md$astral.multi,paired = T)



#----------------------------------------------Figure 4-B-----------------------------------------
setwd("/Users/maryamrabiee/Documents/Research/ASTRAL-multiind/results/ASTRAL-multiind/ASTRAL-multi paper")
require(reshape2)
library("scales")
d=read.csv("monophyly.all.stat",sep=" ",head=F)
d2=read.csv("monophyly.all.stat.final",sep=" ",head=F)
unique(d$V4)
# d0=read.csv("0.contraction.monophyly.stat",sep=" ",head=F)
# nrow(d)
# dd=rbind(d0,d)
head(d)

levels(as.factor(d$V2))
xscl = function(x){-sqrt(1-x)}
require(scales)
ggplot(aes(x=as.factor(V2),y=V4/V3,group=1),data = d2)+ 
  stat_summary(geom="line")+ xlab("Support contraction level")+
  stat_summary(geom="point")+ 
  ylab("Incompatible with monophyly (%)")+ theme_bw()+
  stat_summary(fun.ymax = function(x){mean(x)+sd(x)/sqrt(length(x))},
               fun.ymin=function(x){mean(x)-sd(x)/sqrt(length(x))},
               geom="errorbar",width=0.25,size=.3)+
  scale_y_continuous(labels=percent)+
  #scale_x_continuous(breaks=sapply(c(0.25,0.5,0.75,0.9,0.95,0.99),xscl),labels = c(0.25,0.5,0.75,0.9,0.95,0.99))+
  theme(axis.title.y = element_text(size = rel(0.9), angle = 90), plot.margin = margin(5, 5, 19 , 5))+
  ggsave("monophyly.pdf",height = 3.2,width = 4)

# unique(d$V2)
# ggplot(aes(x=-log10(1-V2,y=V4/V3),data = d)+ stat_summary(geom="line")+ xlab("Contraction level (log scale)")+
#   ylab("Species incompatible with monophyly (%)")+ theme_bw()
# 


reverselog_trans <- function(base = exp(1)) {
  trans <- function(x) -log(x, base)
  inv <- function(x) base^(-x)
  trans_new(paste0("reverselog-", format(base)), trans, inv, 
            log_breaks(base = base), 
            domain = c(1e-100, Inf))
}

# ggplot(aes(x=(1-V2),y=V4/V3),data = d)+ stat_summary(geom="line")+ xlab("Contraction level (log scale)")+
#   ylab("Portion of species incompatible with monophyly")+ theme_bw()+
#   stat_summary(fun.ymax = function(x){mean(x)+sd(x)/sqrt(length(x))},
#                fun.ymin=function(x){mean(x)-sd(x)/sqrt(length(x))},
#                geom="errorbar",width=0.025,size=.3)+
#                 scale_x_log10(breaks=c(0.25,0.5,0.75,0.9,0.95,0.99),trans=reverselog_trans(10))
#   scale_x_continuous(breaks=c(0.25,0.5,0.75,0.9,0.95,0.99),labels = c(0.25,0.5,0.75,0.9,0.95,0.99))+
#   ggsave("monophyley.pdf",height = 4.6,width = 7.4)
  
 
