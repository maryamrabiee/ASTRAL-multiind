setwd("/Users/maryamrabiee/Documents/Research/ASTRAL-multiind/results/ASTRAL-multiind/ASTRAL-multi paper/data")
require(ggplot2)

d=read.csv("RT.all",sep=" ",head=T)


head(d)
a = merge(d,read.csv('../../data/parameter.log.info1',header = T,sep= ' '))
head(a)
a$Input_Size <- a$Number_of_leaves*a$Number_of_Genes
max(a$Input_Size)
a$InputSQ=cut(ecdf(a$Input_Size)(a$Input_Size), breaks=c(0,0.25,0.5,0.75,1),labels=c("Low\n[1886,22032)","Medium\n[22032,45864)", "High\n[45864,76209)","Very High\n[76209,162837]"))

#a =subset(a, , c(Replicate, AST, AST.CONTR , NJST, TrueQscore,Number_of_leaves,Number_of_Genes,InputSQ))

m = melt(a,measure.vars = c("AST","AST.CONTR","NJST"))
head(m)

quantile(a$Input_Size,c(0,0.25,0.5,0.75,1))
ggplot(aes(x=as.factor(InputSQ) ,y=value, fill=as.factor(variable)), data=m)+ geom_boxplot(data=m)+theme_bw()+theme(legend.position=c(.2,.86),legend.direction = 1) + 
         xlab("Input size (# genes * # species)")+ylab("Running Time")+
  #  scale_fill_brewer(palette="Paired",labels=c( "ASTRAL original","ASTRAL-5% contraction","Njst"))+
  theme(legend.title=element_blank())+ scale_fill_manual(values = c("#2b8cbe","#2171b5","#bae4bc"),labels=c(  "ASTRAL-multi","ASTRAL-multi-5%","NJst"))
ggsave('RT-Inputsize.pdf',width =6,height =4.35)                                                                   

m$ILS=cut(ecdf(m$TrueQscore)(m$TrueQscore), breaks=c(0,0.25,0.5,0.75,1),labels=rev(c("Low\n[0.55,0.98]","Medium\n[0.46,0.55)", "High\n[0.411,0.46)", "Very High\n[0.355,0.41]")))
quantile(m$TrueQscore,c(0,0.25,0.5,0.75,1))
ggplot(aes(x=as.factor(ILS) ,y=value, fill=as.factor(variable)), data=m)+ geom_boxplot(data=m,outlier.alpha = 0.8)+theme_bw()+theme(legend.position=c(.89,.87),legend.direction = 1) + 
  xlab("ILS")+ylab("Running Time")+
 # scale_fill_brewer(palette="Paired",labels=c( "ASTRAL original","ASTRAL-5% contraction","Njst"))+
  theme(legend.title=element_blank())+ scale_fill_manual(values = c("#2b8cbe","#2171b5","#bae4bc"),labels=c( "ASTRAL-multi","ASTRAL-multi-5%","NJst"))+
  theme(legend.title=element_blank())

ggsave('RT-ILS.pdf',width =6,height =4.35)

