require(reshape2)
require(ggplot2)
require(Hmisc)

d <- read.csv('../data/seqLength.csv',header=F,sep=" ")

d$V2 <- as.numeric(as.character(d$V2))

h <- dcast(d, V1~.,value.var="V2",fun.aggregate=median)

names(h) <- c("V1","medianSeqLen")
q<-read.csv('../data/parameter.log.info',header=T,sep=" ")
print(h[h$medianSeqLen<300,1])
print(h[h$medianSeqLen>2000,1])
stat<-merge(x=h,y=q,by.x="V1",by.y="Replicate")


f<-read.csv('../data/estimatedBrLengths.csv',sep=" ",header=F)
br<-merge(f,q,by.x="V1",by.y="Replicate",all=TRUE)
br$true<-br$V2/(br$Haploid_efective_population_size)
br$logl <- abs(log10(br[,3])-log10(br[,14]))
br$mrs <- (br[,14]-br[,3])^2
ggplot(data=br[br$V3>0,],aes(x=log10(true),y=log10(V3),color=Number_of_Genes))+
  geom_point(alpha=0.5,size=1.5)+theme_bw()+
  theme(legend.position="bottom")+
  xlab("true branch length (log scale)")+
  geom_abline(linetype=2)+ ylab("Estimated branch length (log scale)")+
  coord_cartesian(ylim=c(-5.3,2.3),xlim=c(-5.3,2.3))+geom_smooth()+scale_color_continuous( low='red',high="white")

ggsave('../figures/BL-alpha0.5.pdf',width=8,height=6.2)

br$labelG <- "H"
br[br$Number_of_Genes>770,]$labelG <- "High NumGenes"
#br[br$Number_of_Genes<=750 & br$Number_of_Genes>500,]$labelG <-"Med NumGenes"
br[br$Number_of_Genes<=770,]$labelG <- "Low NumGenes"
br$labelG <- factor(br$labelG, levels = c("Low NumGenes", "Med NumGenes","High NumGenes" ))

ggplot(data=br[br$V3>0,],aes(x=log10(true),y=log10(V3)))+facet_wrap(~labelG)+
  geom_point(alpha=0.1,size=1.5)+theme_bw()+
  theme(legend.position="bottom")+
  xlab("true branch length (log scale)")+
  geom_abline(linetype=2)+ ylab("Estimated branch length (log scale)")+
  coord_cartesian(ylim=c(-5.3,2.3),xlim=c(-5.3,2.3))+geom_smooth()
ggsave('../figures/BL-diffLevelsNumGenes.pdf',height=9.5, width=5)


bla=recast(melt(data=br[br$V3!=0,c(1,15,16,17)],
                id.vars=c("V1","labelG")),.~variable,
           measure.var = "value",
           fun.aggregate=function(x)mean(x[!is.infinite(x)]))
bla[,2:3] = sqrt(bla[,2:3])
l=latex(format(bla,digits=3),file="../figures/bl.tex",rowname = NULL)

