setwd("/Users/maryamrabiee/Documents/Research/ASTRAL-multiind/results/ASTRAL-multiind/ASTRAL-multi paper")
require(ggplot2)
a=read.csv('clusters.stat',sep=' ',header=T)
nrow(a)
b=read.csv('../data/allReplicates.parameter.info',sep=' ',header=T)
b[b$Replicate==323,]
ss=merge(a,b)
cols =subset(ss, , -c(Replicate, clusters, qScore,Number_of_leaves,Number_of_Genes))
write.csv(ss, file ="cluster.stat.csv")
ggplot(aes(x=Number_of_leaves,y=clusters,fill="#9102bd"),data=cols)+geom_point()+theme_bw()+theme(legend.position=c(.4,.8))+xlab("# of Generations")+ylab("Species tree error (RF distance)")
ggplot(aes(x=Number_of_leaves,y=clusters,color=log(3/2*(qScore-1/3))),data=cols)+geom_point()+theme_bw()+xlab("# of Leaves")+ylab("Number of Clusters")+geom_smooth()
ggplot(aes(x=Number_of_Genes,y=clusters,color=log(3/2*(qScore-1/3))),data=cols)+geom_point()+theme_bw()+xlab("# of Genes")+ylab("Number of Clusters")+geom_smooth()
ggplot(aes(x=qScore,y=clusters,color="#2b8cbe"),data=cols)+geom_point(color="#2b8cbe")+theme_bw()+xlab("Quartet Score")+ylab("Number of Clusters")+geom_smooth()
qplot(log(Number_of_leaves),log(clusters),data=d)+stat_smooth(method="lm")
qplot(log(Number_of_leaves*Number_of_Genes),log(clusters),data=d)+stat_smooth(method="lm")+geom_abline(intercept=0,color="red",linetype=2)+theme_classic()

head(ss)
s <- ss[log(ss$clusters)> 8.5,]
nrow(s)
qplot(log(Number_of_leaves*Number_of_Genes),log(clusters),color=log(3/2*(EstQscore-1/3)),data=ss)+
  stat_smooth(method="lm",se=F,color="black",size=1.3,data = s)+
  geom_abline(intercept=1,color="red",linetype=2,size=1.5)+theme_classic()+theme(legend.position=c(.15,.8))+scale_color_continuous(name="Quartet score")+
  xlab("Number of leaves x Number of Genes")+ylab("Search space size")
ggsave("D1-clusters.pdf",height = 4 , width = 5.5)

