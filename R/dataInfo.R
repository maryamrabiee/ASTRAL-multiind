require(ggplot2)
require(reshape2)
setwd('/Users/erfan/Main/Library/ASTRAL-multiind')
rd1 = "#225ea8"
rd2 = '#41b6c4'
rd3 = '#a1dab4'
bl1 = '#ffffcc'
d<-read.csv('data/gtError-simphy-multiind.csv',sep=" ",header=F)
d<-d[,c(1:7)]
h<-d[!(d$V3 %in% c("-")),]

h$V2<-as.numeric(as.character(factor(h$V2)))
h$V3<-as.numeric(as.character(factor(h$V3)))
h$V4<-as.numeric(as.character(factor(h$V4)))
h$V5<-as.numeric(as.character(factor(h$V5)))
h$V6<-as.numeric(as.character(factor(h$V6)))
h$V7<-as.numeric(as.character(factor(h$V7)))

h$rf<-(h$V3+h$V6)/(h$V2+h$V5)
names(h)<-c("V2","V3","V4","V5","V6","V7","V8","rf")
h$V1<-"S330"
h<-h[,c(9,1:8)]

d2<-read.csv('data/newdataset/gterror.csv',sep=" ",header=F)
d2$V3<-as.numeric(as.character(factor(d2$V3)))
d2$V4<-as.numeric(as.character(factor(d2$V4)))
d2$V5<-as.numeric(as.character(factor(d2$V5)))
d2$V6<-as.numeric(as.character(factor(d2$V6)))
d2$V7<-as.numeric(as.character(factor(d2$V7)))
d2$rf<-(d2$V4+d2$V7)/(d2$V3+d2$V6)

htotal<-rbind(h,d2)
htotal$V1<-as.factor(htotal$V1)
levels(htotal$V1)<-list("D2-1M"="model.200-5.1000000.0.000001","D2-2M"="model.200-5.2000000.0.000001","D2-0.5M"="model.200-5.500000.0.000001","D1"="S330")
htotal$V1<-factor(htotal$V1,levels=c("D2-0.5M","D2-1M","D2-2M","D1"))
ggplot(data=htotal,aes(x=rf,group=V1,fill=V1))+geom_density(alpha=0.8,adjust=1.5)+
  theme_bw()+xlab('Gt Error (NRF ration)')+ylab('Density')+
  theme(legend.position = c(0.8,0.8),axis.text=element_text(size=12,color="black"),
        axis.title=element_text(size=12,color="black"))+
  scale_fill_manual(name="",values=c(rd1,rd2,rd3,bl1))
ggsave('figures/gtError.pdf')


d<-read.csv('data/newdataset/score-speciesTrees.score',sep=" ",header=F)
d<-d[,c(1,2,3)]

d2<-read.csv('data/quartetScores.csv',sep=" ",header=T)
d2<-d2[1:nrow(d2),c(1,2)]
d2$V1<-"S330"
d2<-d2[,c(3,1,2)]
names(d)<-c("V1","Rep","TrueQscore")
finalsc<-rbind(d,d2)
levels(finalsc$V1)<-list("D2-1M"="model.200-5.1000000.0.000001","D2-2M"="model.200-5.2000000.0.000001","D2-0.5M"="model.200-5.500000.0.000001","D1"="S330")
finalsc$V1<-factor(finalsc$V1,levels=c("D2-0.5M","D2-1M","D2-2M","D1"))
ggplot(data=finalsc,aes(x=TrueQscore,group=V1,fill=V1))+geom_density(alpha=0.8,adjust=1.5)+
  theme_bw()+xlab('Normalized ASTRAL quartet score')+ylab('Density')+
  theme(legend.position = c(0.8,0.8),axis.text=element_text(size=12,color="black"),
        axis.title=element_text(size=12,color="black"))+
  scale_fill_manual(name="",values=c(rd1,rd2,rd3,bl1))+scale_x_continuous(breaks=c(0.33,0.4, 0.6,0.8,1))
ggsave('figures/quartetScore.pdf')

#d<-read.csv('data/result2.csv',sep=" ",header=F)
#qplot(data=d,x=reorder(as.factor(V1),V5),y=V5,color=as.factor(V2),geom='jitter')+theme_bw()+xlab('replicate')+ylab('FN error')+theme(axis.ticks = element_blank(), axis.text.x = element_text(angle = 90, hjust = 0, colour = "grey50"))+scale_color_brewer(name="Number of rounds",palette="Dark2")
#ggsave('EstimatedGeneTrees.pdf')
p<-read.csv('data/parameter.log.info1',sep="\t",header=T)

h<-read.csv('data/newdataset/score-speciesTrees.score',sep=" ",header=F)

dcast(data=h,V1~.,fun.aggregate = mean)

y<-read.csv('data/newdataset/seqLength.txt',header=F,sep=" ")
y$V3<-as.numeric(as.character(y$V3))
y$V2<-as.factor(y$V2)
g<-dcast(data=y,V1+V2~.,fun.aggregate=median,value.var="V3")


q<-read.csv('data/newdataset/gterror.csv',sep=" ",header=F)
q$V3<-as.numeric(as.character(q$V3))
q$V4<-as.numeric(as.character(q$V4))
q$V6<-as.numeric(as.character(q$V6))
q$V7<-as.numeric(as.character(q$V7))
q$rf<-(q$V4+q$V7)/(q$V3+q$V6)
dcast(data=q,V1~.,fun.aggregate = mean,value.var = "rf")
