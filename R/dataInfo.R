require(ggplot2)
require(reshape2)

d<-read.csv('../data/gtError-simphy-multiind.csv',sep=" ",header=F)

h<-d[!(d$V3 %in% c("-")),]

h$V2<-as.numeric(as.character(factor(h$V2)))
h$V3<-as.numeric(as.character(factor(h$V3)))
h$V4<-as.numeric(as.character(factor(h$V4)))
h$V5<-as.numeric(as.character(factor(h$V5)))
h$V6<-as.numeric(as.character(factor(h$V6)))
h$V7<-as.numeric(as.character(factor(h$V7)))

h$rf<-(h$V3+h$V6)/(h$V2+h$V5)
ggplot(data=h,aes(x=rf))+geom_density(fill='red',alpha=0.5,adjust=1.5)+theme_bw()+xlab('RF Distance (True vs Estimated)')+ylab('Density')

d<-read.csv('../data/result2.csv',sep=" ",header=F)
qplot(data=d,x=reorder(as.factor(V1),V5),y=V5,color=as.factor(V2),geom='jitter')+theme_bw()+xlab('replicate')+ylab('FN error')+theme(axis.ticks = element_blank(), axis.text.x = element_text(angle = 90, hjust = 0, colour = "grey50"))+scale_color_brewer(name="Number of rounds",palette="Dark2")
ggsave('EstimatedGeneTrees.pdf')
