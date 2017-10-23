require(ggplot2)
require(reshape2)

d=read.csv("data/accuracy.stat.truegenes",sep=" ",head=F)

h=read.csv('data/gtError-simphy-multiind.csv',sep=" ",header=F)
h<-h[,c(1,2,3,4,5,6,7)]
h<-h[!(h$V3 %in% c("-")),]

h$V2<-as.numeric(as.character(factor(h$V2)))
h$V3<-as.numeric(as.character(factor(h$V3)))
h$V4<-as.numeric(as.character(factor(h$V4)))
h$V5<-as.numeric(as.character(factor(h$V5)))
h$V6<-as.numeric(as.character(factor(h$V6)))
h$V7<-as.numeric(as.character(factor(h$V7)))

h$rf<-(h$V3+h$V6)/(h$V2+h$V5)
agt=dcast(V1~., data=h[,c(1,9)],fun.aggregate = function(x) mean (x,na.rm = T))
names(agt)[2] = "RF"

a = merge(merge(d,agt,by.x="V1",by.y="V1"),read.csv('data/parameter.log.info1.final',header = T,sep= ' '),by.x="V1", by.y="Replicate")
names(a)[5] = "stError"
names(a)[1] = "Replicate"

ggplot(aes(x=Number_of_Genes,y=stError, color=log(3/2*(TrueQscore-1/3))), data=a)+
  geom_point()+theme_bw()+geom_smooth()+theme(legend.position=c(.6,.86),legend.direction = 1)+
  scale_color_continuous(name="Quartet score",breaks=c(-1,-2,-3),labels=format(2/3*exp(-1:-3)+1/3,digits = 2))+
  xlab("# of genes")+ylab("Species tree error (RF distance)")

ggsave("st-truegt.pdf",width = 4.5,height = 3.5)

