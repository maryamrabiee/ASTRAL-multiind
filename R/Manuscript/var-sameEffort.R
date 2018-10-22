setwd("/Users/maryamrabiee/Documents/Research/ASTRAL-multiind/results/ASTRAL-multiind/ASTRAL-multi paper/var-fixed-effort")

est=read.csv('est.all.stat.RQ1',sep= ' ',h=F)
true=read.csv('true.all.stat.RQ1',sep= ' ',h=F)
head(est)
head(true)

est2 <- est
true2 <- true
tail(est2)
est2$V4 <-as.numeric(est2$V4)
est2$V2 <- as.factor(est2$V2)
levels(est2$V2) <-list("1"="1ind","2"="2ind","5"="5ind")
est2$V2 <-as.numeric(as.character(est2$V2))
true2$V4 <-as.numeric(true2$V4)
true2$V2 <- as.factor(true2$V2)
levels(true2$V2) <-list("1"="1ind","2"="2ind","5"="5ind")
est2$V2 <-as.numeric(as.character(est2$V2))

summary(aov(V7~V2*V4,data=est2[est2$V3=="variableEffort",]))
summary(aov(V7~V2*V4,data=est2[true2$V3=="variableEffort",]))

means <- dcast(V4~V2,value.var = "V7",data=est[est$V3=="variableEffort",],fun.aggregate = mean)
means[,4]

x = (mean(a[a$V3=="1IND" & a$V4=="0.5M",c("V2")])+ mean(a[a$V3=="1IND" & a$V4=="1M",c("V2")])+mean(a[a$V3=="1IND" & a$V4=="2M",c("V2")]))/3
x
y = (mean(a[a$V3=="5IND" & a$V4=="0.5M",c("V2")])+ mean(a[a$V3=="5IND" & a$V4=="1M",c("V2")])+mean(a[a$V3=="5IND" & a$V4=="2M",c("V2")]))/3
y
(mean(a[a$V3=="1IND" & a$V4=="0.5M",c("V2")])-mean(a[a$V3=="5IND" & a$V4=="0.5M",c("V2")])) / mean(a[a$V3=="1IND" & a$V4=="0.5M",c("V2")])

summary(aov(V2~V4*V3,data=a))
is.numeric(is.character(est$V7))
levels(est$V4) <-list("2 M"="2000000","1 M"="1000000","0.5 M"="500000")

#--------------------------------EST GENE TREES------------------------ 
est$V4 <-as.factor(est$V4)
levels(est$V4)
levels(est$V4) <-list("0.5M"="500000","1M"="1000000","2M"="2000000")
ggplot(aes(x=as.factor(V4),y=V7,fill=as.factor(V2)),data=est[est$V3=="variableEffort",])+geom_boxplot(outlier.alpha=0.8,width=.6)+
  xlab("Model condition (tree height)")+ylab("Species tree error (RF distance)")+
  scale_fill_brewer(palette=1,labels=c("1 individual x 500 genes","2 individuals x 500 genes","5 individuals x 500 genes"),name="")+
  scale_x_discrete()+theme_classic()+coord_cartesian(ylim=c(0,0.2)) +
#  theme(legend.position=c(.85,.71),legend.margin=margin(0, 0, 0, 0),legend.title = element_blank(), plot.margin=margin(5, 25, 5, 5))+
  guides(fill=FALSE)+
  
ggsave("variableEffort.pdf",width =4.4,height =3)


ggplot(aes(x=as.factor(V4),y=V7,fill=as.factor(V2)),data=est[est$V3=="fixedEffort",])+geom_boxplot(outlier.alpha=0.8,width=.6)+
  xlab("Model condition (tree height)")+ylab("Species tree error (RF distance)")+
  scale_fill_brewer(palette=1,labels=c("1 individual x 1000 genes","2 individuals x 500 genes","5 individuals x 200 genes"),name="")+
  scale_x_discrete()+theme_classic()+coord_cartesian(ylim=c(0,0.2))+
  #theme(legend.position=c(.85,.70),legend.margin=margin(0, 0, 0, 0),legend.title = element_blank(),plot.margin=margin(5, 25, 5, 5))+
  guides(fill=FALSE)+
ggsave("fixedEffort.pdf",width =4.4,height =3)

b=read.csv('125ind.D1.stat',sep= ' ',h=F)
summary(aov(V2~V4*V3,data=b))
head(b)
#b$V4 <- as.factor(b$V4)
#levels(b$V4) <-list("2 M"="2M","1 M"="1M","0.5 M"="0.5M")
# ggplot(aes(x=as.factor(V4),y=V2,fill=as.factor(V3)),data=b)+geom_boxplot(outlier.alpha=0.8,width=.6)+
#   xlab("Model condition (tree height)")+ylab("Species tree error (RF distance)")+
#   scale_fill_brewer(palette=1,labels=c("1 ind x 1000 genes","2 ind x 500 genes","5 ind x 200 genes"),name="")+
#   scale_x_discrete()+theme_classic()+theme(legend.position=c(.89,.69))+coord_cartesian(ylim=c(0,0.2))
# ggsave("fixedEffort.pdf",width =5.47,height =3.99)

head(est)
dcast(V2~V4,data=est[true$V3=="fixedEffort",],fun.aggregate = function(x) sum(x!=0))

#d$V4 <- as.factor(V4)
#levels(d$V4) <-list("2 M"="2000000","1 M"="1000000","0.5 M"="500000")
#levels(d$V3)

#--------------------------------TRUE GENE TREES------------------------
d = true
head(d)
unique(d$V4)
d$V4 <-as.factor(d$V4)
levels(d$V4)
levels(d$V4) <-list("0.5M"="500000","1M"="1000000","2M"="2000000")
ggplot(aes(x=as.factor(V4),y=V7,fill=as.factor(V2)),data = d[ d$V3 == "fixedEffort", ])+geom_boxplot(outlier.alpha=0.8,width=.6)+
  xlab("Model condition (tree height)")+ylab("Species tree error (RF distance)")+
  scale_fill_brewer(palette=1,labels=c("1 individual x 1000 genes","2 individuals x 500 genes","5 individuals x 200 genes"),name="")+
  #scale_x_discrete()+theme_classic()+theme(legend.position=c(.85,.65),legend.title = element_blank(),plot.margin=margin(5, 25, 5, 5))+coord_cartesian(ylim=c(0,0.2))
  scale_x_discrete()+theme_classic()+theme(legend.position=c(.8,.65),legend.title = element_blank())+coord_cartesian(ylim=c(0,0.2))

ggsave("true-fixedEffort.pdf",width =4.4,height =3)

ggplot(aes(x=as.factor(V4),y=V7,fill=as.factor(V2)),data = d[ d$V3 == "variableEffort", ])+geom_boxplot(outlier.alpha=0.8,width=.6)+
  xlab("Model condition (tree height)")+ylab("Species tree error (RF distance)")+
  scale_fill_brewer(palette=1,labels=c("1 individual x 500 genes","2 individuals x 500 genes","5 individuals x 500 genes"),name="")+
 # scale_x_discrete()+theme_classic()+ theme(legend.position=c(.85,.65),legend.title = element_blank(),plot.margin=margin(5, 25, 5, 5))+
  scale_x_discrete()+theme_classic()+ theme(legend.position=c(.8,.65),legend.title = element_blank())+
  coord_cartesian(ylim=c(0,0.2))
ggsave("true-varEffort.pdf",width =4.4,height =3)
