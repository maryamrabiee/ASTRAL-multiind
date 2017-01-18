require(ggplot2)
require(reshape2)
require(reshape)

d<-read.csv('../data/gradient_comp.csv',header=F,sep=" ")
names(d)<-c("method_astralII","method_greedy","replicate","numBranch","missingBranch","missRate")
d$replicate <- as.factor(d$replicate)
d$method<-as.factor(paste(d$method_astralII,d$method_greedy,sep="-"))
ggplot(data=d, aes(method, missRate,color="red")) + 
  stat_summary(fun.y ="mean",geom = "point") + theme_bw() + 
  xlab('Method')+
  theme(text = element_text(size=12),
        axis.text.x = element_text(size=12,angle = 0),
        axis.text.y = element_text(size=12,angle = 0),
        legend.position="none")+
  ylab('RF distance')
ggsave('../figures/Accuracy-trueSpecies.pdf')
