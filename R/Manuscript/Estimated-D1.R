setwd("/Users/maryamrabiee/Documents/Research/ASTRAL-multiind/results/ASTRAL-multiind/ASTRAL-multi paper/")
require(ggplot2)
require(reshape2)
require(segmented)
am<-data.frame()
d=read.csv("true.est.stat",sep=" ",head=T)
#d=read.csv("truegt.stat",sep=" ",head=F)
head(d)
nrow(d)
names(d)[1] = "Replicate"
#dd=read.csv("estgt.stat",sep=" ",head=F)

h=read.csv('../data/gtError-simphy-multiind.csv',sep=" ",header=F)
h<-h[h$V3!="-",]
h$V2<-as.numeric(as.character(factor(h$V2)))
h$V3<-as.numeric(as.character(factor(h$V3)))
h$V4<-as.numeric(as.character(factor(h$V4)))
h$V5<-as.numeric(as.character(factor(h$V5)))
h$V6<-as.numeric(as.character(factor(h$V6)))
h$V7<-as.numeric(as.character(factor(h$V7)))

h$agt<-(h$V3+h$V6)/(h$V2+h$V5)
agt=dcast(V1~., data=h[,c(1,9)],fun.aggregate = function(x) mean (x,na.rm = T))
names(agt)[2] = "agt"
names(agt)[1] = "Replicate"

a = merge(d,read.csv('../data/parameter.log.info1',header = T,sep= ' '))
head(agt)
a = merge(a,agt)
# names(b)[7] = "True Gene Trees"
# names(b)[4] = "Estimated Gene Trees"
# names(b)[1] = "Replicate"
#m = melt(b,measure.vars = c("True Gene Trees","Estimated Gene Trees"))

a$AGTE=cut(ecdf(a$agt)(a$agt), breaks=c(0,0.25,0.5,0.75,1),labels=c("[0, 0.25)","[0.25, 0.5)", "[0.5, 0.75)","[0.75, 1]"))
head(a)
m = melt(a,measure.vars = c("true","est"))
head(m)
a$ILS=cut(ecdf(a$EstQscore)(a$EstQscore), breaks=c(0,0.33,0.66,1),labels=c("High","Medium", "Low"))

#ggplot(aes(x=ILS,y=value, fill=variable), data=m)+ geom_boxplot()+theme_bw()+theme(legend.position=c(.6,.86),legend.direction = 1)+
#  scale_fill_brewer(palette = "Blues")+
#  xlab("ILS")+ylab("Species tree error (RF distance)")

#ggsave("st-truegt.pdf",width = 4.5,height = 3.5)


head(a)
names(a)[2]="True Gene trees"
names(a)[3]="Estimated Gene trees"
#am$Genes <- as.factor(am$Genes)
am$AGTE <- as.factor(am$AGTE)
am[am$value > 0.3,]
am = melt(a, measure.vars = c("True Gene trees","Estimated Gene trees"))
am$Genes=cut(ecdf(am$Number_of_Genes)(am$Number_of_Genes), breaks=c(0,0.25,0.5,0.75,1),labels=c("[0, 0.25)","[0.25, 0.5)", "[0.5, 0.75)","[0.75, 1]"))
head(am)

lmm = lm(value~agt,data=am[am$variable=="Estimated Gene trees",]); 
segmented(lmm,seg.Z = ~agt)

am$AGTE[am$variable == "True Gene trees"] = "[0, 0.25)"
head(am)
g11 = ggplot(aes(x=agt,y=value, color=log(3/2*(TrueQscore-1/3))), data=am[am$variable=="Estimated Gene trees",])+ geom_point(aes(shape=Genes))+ theme_bw()+
  scale_shape_manual(values=c(16, 17, 15,18))+geom_smooth(aes(group=agt<0.2954),method="glm",color="red",se =F)+
  scale_color_continuous(name="Quartet score",breaks=c(-0.6,-1.8,-3.1),labels=function(x) format(2/3*exp(x)+1/3,digits = 2))+ 
  xlab("Average gene tree error (AGTE)")+ylab("Species tree error (RF distance)")+
  theme(legend.margin = margin(0,0,0,0))
  #theme(axis.text=element_text(size=12),axis.title=element_text(size=14),strip.text = element_text(size=15))   
ggsave("estimated/RF-AGT-D3.pdf",width = 4.95,height = 3.5)



g12 = ggplot(aes(x=Number_of_Genes,y=value, color=log(3/2*(TrueQscore-1/3))), data=am[am$variable=="Estimated Gene trees",])+ 
  geom_point(aes(shape=AGTE))+theme_bw()+
  scale_shape_manual(values=c(16, 17, 15,18))+
  geom_smooth(aes(group=Number_of_Genes<403),method="glm",color="red",se =F)+
  scale_color_continuous(name="Quartet score",breaks=c(-0.6,-1.8,-3.1),labels=function(x) format(2/3*exp(x)+1/3,digits = 2))+
  xlab("Number of genes")+ylab("Species tree error (RF distance)")+
  theme(legend.margin = margin(0,0,0,0))
 # theme(axis.text=element_text(size=12),axis.title=element_text(size=14),strip.text = element_text(size=15)) +  
ggsave("estimated/RF-Genes-D3-est.pdf",width = 5,height = 3.5)


g21=ggplot(aes(x=Number_of_leaves,y=value, color=log(3/2*(TrueQscore-1/3))), data=am[am$variable=="Estimated Gene trees",])+ 
  geom_point()+theme_bw()+geom_smooth(method="glm",color="red",se =F)+
  scale_color_continuous(name="Quartet score",breaks=c(-0.6,-1.8,-3.1),labels=function(x) format(2/3*exp(x)+1/3,digits = 2))+
  xlab("Number of species")+ylab("Species tree error (RF distance)")+
  theme(legend.margin = margin(0,0,0,0))
  #theme(axis.text=element_text(size=12),axis.title=element_text(size=14),strip.text = element_text(size=15))   
ggsave("estimated/RF-species-D3.pdf",width = 5,height = 3.5)

lmm = lm(value~genlog,data=am[am$variable=="Estimated Gene trees",]); 
segmented(lmm,seg.Z = ~genlog)

g31 = ggplot(aes(x=log10(Generations),y=value, color=log(3/2*(TrueQscore-1/3))), data=am[am$variable=="Estimated Gene trees",])+ geom_point()+theme_bw()+
  geom_smooth(method="glm",color="red",se =F)+
  scale_color_continuous(name="Quartet score",breaks=c(-0.6,-1.8,-3.1),labels=function(x) format(2/3*exp(x)+1/3,digits = 2))+
  xlab("Log of number of generations")+ylab("Species tree error (RF distance)")+
  theme(legend.margin = margin(0,0,0,0))
  #theme(axis.text=element_text(size=12),axis.title=element_text(size=14),strip.text = element_text(size=15))   
ggsave("estimated/RF-generations-D3.pdf",width = 5,height = 3.5)


lmm = lm(value~Haploid_efective_population_size,data=am[am$variable=="Estimated Gene trees",]); 
segmented(lmm,seg.Z = ~Haploid_efective_population_size)

g32= ggplot(aes(x=log10(Haploid_efective_population_size),y=value, color=log(3/2*(TrueQscore-1/3))), data=am[am$variable=="Estimated Gene trees",])+
  geom_point(aes(shape=AGTE))+theme_bw()+
  scale_color_continuous(name="Quartet score",breaks=c(-0.6,-1.8,-3.1),labels=function(x) format(2/3*exp(x)+1/3,digits = 2))+
  geom_smooth(method="glm",color="red",se =F)+
  scale_shape_manual(values=c(16, 17, 15,18))+
  xlab("Haploid efective population size")+ylab("Species tree error (RF distance)")+
  theme(legend.margin = margin(0,0,0,0))
  #theme(axis.text=element_text(size=12),axis.title=element_text(size=14),strip.text = element_text(size=15))   
ggsave("estimated/RF-population-D3.pdf",width = 5,height = 3.5)

am$qt = log(3/2*(am$TrueQscore-1/3))
            
lmm = lm(value~qt,data=am[am$variable=="Estimated Gene trees",]); 
segmented(lmm,seg.Z = ~qt)

qtf =  function(x) format(2/3*exp(x)+1/3)

head(am)
am$Replicate[is.na(as.numeric(as.character(am$value)))]
g22=ggplot( aes(x=log(3/2*(TrueQscore-1/3)),y=value ),data=am[am$variable=="Estimated Gene trees",])+ 
  geom_point(aes(color=Number_of_Genes,shape=AGTE))+
  geom_smooth(color="red",se=F,method="lm")+ theme_bw()+
  scale_color_continuous(name="Genes")+
  scale_shape_manual(values=c(16, 17, 15,18))+
  xlab("True quartet score (log)")+ylab("Species tree error (RF distance)")+
  theme(legend.margin = margin(0,10,0,0))+
 # theme(axis.text=element_text(size=12), axis.title=element_text(size=14),strip.text = element_text(size=15))+
  scale_x_continuous(labels = ,breaks=c(-.6,-1.8,-3.1))
ggsave("estimated/RF-trueQscore-D3.pdf",width = 6,height = 4.3)


am[am$agt<0.34 & am$Number_of_Genes < 130 & am$TrueQscore < 0.4,]
lmm = lm(value~Number_of_Genes,data=am[am$variable=="Estimated Gene trees",]); 
segmented(lmm,seg.Z = ~Number_of_Genes)

lmm = lm(value~TrueQscore,data=am[am$variable=="Estimated Gene trees",]); 
segmented(lmm,seg.Z = ~TrueQscore, psi=c(0.42))

first <- am[am$Number_of_Genes < 403.8,]
second <- am[am$Number_of_Genes > 403.8,]
cor(first$Number_of_Genes,first$value)
cor(second$Number_of_Genes,second$value)



summary(aov(formula=value~Number_of_Genes+qt+agt, data=am[am$variable=="Estimated Gene trees",]))
summary(aov(formula=value~ngt+qt+agt, data=am[am$variable=="Estimated Gene trees",]))
lmm = lm(value~Haploid_efective_population_size,data=am[am$variable=="Estimated Gene trees",]);
segmented(lmm,seg.Z = ~Haploid_efective_population_size)

am$height=am$Generations/am$Haploid_efective_population_size
am$inputsize=am$Number_of_leaves*am$Number_of_Genes
am$ngt = log(am$Number_of_Genes)
am$birthdeath = log(am$Speciation_rate)

b<-summary(aov(formula=value~qt+ngt+agt, data=am[am$variable=="Estimated Gene trees",]))
sum((b[[1]]$'Sum Sq')[1:3])/sum((b[[1]]$'Sum Sq'))

library(cowplot)

allp<-plot_grid(plot_grid(g11,g12,labels=c("A)", "B)")),plot_grid(g21,g22,labels=c("C)", "D)")),plot_grid(g31,g32,labels=c("E)", "F)")),nrow = 3)
save_plot("all.pdf",base_width = 9,base_height = 9.2, allp)
#                      labels=c("A", "B", "C", "D","E","F"),nrow = 3)

