require(ggplot2)
require(reshape2)
require(cowplot)
d<-read.csv('/Users/erfan/Main/Library/ASTRAL-multiind/data/newdataset/branchlengths_semiterminal.txt',header=F,sep=" ")
p<-read.csv('/Users/erfan/Main/Library/ASTRAL-multiind/data/parameter.log.info1',sep=" ",header=T)
k<-merge(y=p,x=d,by.x="V1",by.y="Replicate")
k$br<-k$V4
k[k$V2=="trueBr",]$br<-k[k$V2=="trueBr",]$V4/k[k$V2=="trueBr",]$Haploid_efective_population_size
kt<-k[,c(1,2,3,6,8:20)]

kt1<-kt[kt$V2=="estimatedBr",]
kt2<-kt[kt$V2=="trueBr",]
nk<-names(kt1)
f=merge(x=kt1,y=kt2,by.x=nk[c(1,3,4,5:16)],by.y=nk[c(1,3,4,5:16)])
f$V7<-"terminal branches"
f$V3<-factor(f$V3,levels=c("truegt","estgt"))
levels(f$V3)<-list("true gene trees"="truegt","estimated gene trees"="estgt")
a=10

l$V9<-cut(f$br.y ,quantile(ecdf(round(f$br.y *10^8)/10^8),c(0:a*1/a)))

plot.semi<-ggplot(data=f,aes(x=(br.y),y=(br.x+10^(-6))))+geom_point(alpha=0.05,color="#084594")+
  geom_abline()+theme_bw()+facet_grid(V7~V3)+
  scale_x_log10()+scale_y_log10(breaks=c(10^-6,10^-4,10^-1,10^2),labels=function(x) {c("0",10^-4,10^-1,10^2)})+ 
  coord_cartesian(xlim = c(10^-6,10^1.7),ylim=c(10^-6,10^1.7))+geom_hline(yintercept = 10^-6,linetype=2)+
  xlab("true branch lengths (cu)")+ylab("estimated branch lengths (cu)")
#ggsave('/Users/erfan/Main/Library/ASTRAL-multiind/figures/branchlengths-semiinternal.pdf',width=7.08 ,height=3.71)

plot.line<-ggplot(data=d[d$V2=="estimatedBr" ,],aes(x=(V5),group=V3,color=V3))+stat_ecdf(size=1.5)+theme_bw()+
  scale_color_manual(name="",labels=c("estimated gt","true gt"),values=c("#1f78b4","#33a02c"))+
  xlab('LocalPP values for terminal branches')+ylab('ecdf')+
  theme(legend.position = c(0.2,0.7))
  
#ggsave('/Users/erfan/Main/Library/ASTRAL-multiind/figures/localpp-semiinternal-ecdf.pdf',width=7.08,height=7.48)

d2<-read.csv('/Users/erfan/Main/Library/ASTRAL-multiind/data/newdataset/branchlengths_iterminal.txt',header=F,sep=" ")
dp<-merge(x=d2,y=p,by.x="V1",by.y="Replicate")
dp$V3<-dp$V3/dp$Haploid_efective_population_size
dp$V2<-factor(dp$V2,levels=c("truegt","estgt"))
levels(dp$V2)<-list("true gene trees"="truegt","estimated gene trees"="estgt")
dp$V6<-"Internal branches"
plot.internal<-ggplot(data=dp,aes(x=(V3),y=(V4)+10^-8))+geom_point(alpha=0.05,color="#084594")+
  geom_abline()+theme_bw()+facet_grid(V6~V2)+
  scale_x_log10()+scale_y_log10(breaks=c(10^-8,10^-6,10^-4,10^-1,10^2),labels=function(x) {c("0",10^-6,10^-4,10^-1,10^2)})+ 
  coord_cartesian(xlim = c(10^-8,10^1.7),ylim=c(10^-8,10^1.7))+geom_hline(yintercept = 10^-8,linetype=2)+
  xlab("true branch lengths (cu)")+ylab("estimated branch lengths (cu)")
#ggsave('/Users/erfan/Main/Library/ASTRAL-multiind/figures/branchlengths-internal.pdf',width=7.08 ,height=3.71)
plot_grid(plot_grid(plot.semi,plot.internal,nrow=2,labels=c("a)","b)")),plot.line,labels=c("","c)"))
ggsave('/Users/erfan/Main/Library/ASTRAL-multiind/figures/branchlengths-total.pdf',width=9.19 , height=5.27)



l<-read.csv('/Users/erfan/Main/Library/ASTRAL-multiind/data/newdataset/branchlengths-D2.txt.1',sep=" ",header=F)
l$V6<-as.numeric(as.character(l$V6))
l$V7<-as.numeric(as.character(l$V7))
l$V8<-as.numeric(as.character(l$V8))
l$V6<-l$V6/200000
a=100

l$V9<-cut(l$V6,quantile(ecdf(round(l$V6*10^8)/10^8),c(0:a*1/a)))

l$V5<-factor(l$V5,levels = c("truegt","estgt"))
l$V1<-factor(l$V1,levels = c("model.200-5.500000.0.000001","model.200-5.1000000.0.000001","model.200-5.2000000.0.000001"))

levels(l$V1)<-list("0.5M"="model.200-5.500000.0.000001","1M"="model.200-5.1000000.0.000001","2M"="model.200-5.2000000.0.000001")
levels(l$V5)<-list("true gene trees"="truegt","estimated gene trees"="estgt")
l$V11<-as.factor(paste(l$V3,l$V4,sep="-"))
levels(l$V11)<-list("1ind-1000gt"="1ind-1000","5ind-1000gt"="multiind-1000","5ind-200gt"="multiind-200")
l$V11<-factor(l$V11,levels=c("5ind-1000gt","1ind-1000gt","5ind-200gt"))

l$V10<-l$V7/l$V6

l$V10<-l$V10+(10^(-6))
ggplot(data=l,aes(x=V9,y=V10,group=V11,color=V11))+
  facet_grid(V1~V5,scales = "free_y")+stat_summary(fun.y="mean",geom="line",size=0.6)+
  stat_summary(fun.ymax=function(x){mean(x)+sd(x)/sqrt(length(x))},
               fun.ymin=function(x){mean(x)-sd(x)/sqrt(length(x))},geom="errorbar",
               size=0.6,width=0.1)+scale_y_log10()+
  theme_bw()+scale_color_brewer(name="",palette = "Set1")+
  geom_hline(yintercept = 1,color="red",linetype=2)+
  theme(legend.position = "bottom",axis.text.x = element_text(size=8,hjust = 1,angle = 30))+
  xlab("true internal branch lenghts (cu.)")+ylab("(estimated bl (cu.))/(true bl (cu.)) in log scale ")
ggsave('/Users/erfan/Main/Library/ASTRAL-multiind/figures/lineplot-internal-branchLengths-D2.pdf',width=9.03, height=7.12)


ggplot(data=l,aes(x=as.factor(V1),y=V10,color=V11))+
  facet_grid(~V5)+geom_violin()+
  theme_bw()+scale_fill_brewer(name="",palette = "Paired")+
  theme(legend.position = "bottom",axis.text.x = element_text(size=12,hjust = 1,angle = 90))+
  xlab("model condition")+
  ylab("estimated br/(true br)(cu.)")+scale_y_log10()
ggsave('/Users/erfan/Main/Library/ASTRAL-multiind/figures/boxplot-ratio-internal-branchLengths-D2.pdf',width=9.03,height=9)


ggplot(data=l,aes(x=V11,y=V10,group=V11,color=V11))+
  facet_grid(V1~V5)+geom_boxplot()+
  theme_bw()+scale_fill_brewer(name="",palette = "Paired")+
  theme(legend.position = "bottom",axis.text.x = element_text(size=12,hjust = 1,angle = 90))+
  xlab("model condition")+coord_cartesian(ylim=c(10^-3,10^3))+
  ylab("estimated br/(true br)(cu.)")+scale_y_log10()

