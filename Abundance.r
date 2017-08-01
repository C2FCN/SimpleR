library(ggplot2)
library(reshape2)
library(Rmisc) 
df<-read.table("Order_Genus_Unclassified_new.xls",row.names = 1,header = T,check.names = F,sep = "\t",stringsAsFactors = F)
all<-read.table("all.xls",row.names = 1,header = T,check.names = F,sep = "\t",stringsAsFactors = F)
map<-read.table("map.txt",header = T,sep = "\t",stringsAsFactors = F)
#df<-rbind(df[1:4,],others)
#rownames(df)[5]<-"Others"

group<-unique(map[,2])

t<-c()
for (i in 1:length(group)){
  tmp<-df[,as.character(map[map[,2]==group[i],1])]
  tt<-rowSums(tmp)
  t<-rbind(t,tt)
}
groupdf<-t(t)
colnames(groupdf)<-group

all<-groupdf[drop=F,nrow(groupdf),]
groupdf<-groupdf[-nrow(groupdf),]
groupdf<-prop.table(as.matrix(groupdf),2) #
groupdf[,1]<-groupdf[,1]*all[,1]
groupdf[,2]<-groupdf[,2]*all[,2]
groupdf<-melt(groupdf) #
SE<-summarySE(groupdf,measurevar = "value",groupvars = "Var2")

mytheme<-theme(legend.position="right",
                 panel.grid.major.x=element_blank(),
                 panel.grid.minor.x=element_blank(),
                 panel.grid.major.y=element_blank(),
                 panel.grid.minor.y=element_blank(),
                 panel.background=element_rect(fill="white"),
                 axis.text.x = element_text(angle=0, hjust = 0.5,vjust=0.5))

mycol<-c("#FF2600","#0522FF","#E98005","#008F00","#919191","#E2EDF8","#0093C6","#FFFB00",
         "#5AB9D1","#672D00","#0D0D00","#4180FF","#C063BC","#EAC282")


SE<-SE[,-c(2,4,6)]
SE<-as.data.frame(SE)
groupdf[groupdf[,2]=="hindgut",4]<-SE[2,3]
groupdf[groupdf[,2]=="foregut",4]<-SE[1,3]
groupdf[groupdf[,2]=="hindgut",5]<-SE[2,2]
groupdf[groupdf[,2]=="foregut",5]<-SE[1,2]
colnames(groupdf)[4:5]<-c("SEse","SEvalue")
groupdf[,c(4,5)]<-groupdf[,c(4,5)]*10
pdf("Order_Genus_Unclassified_new.pdf",width = 5)
p<-ggplot(groupdf)+geom_bar(aes(x=Var2,y=value,fill=Var1),stat = "identity")+
  xlab(NULL)+ylab("Relative abundance(%)")+mytheme+
  scale_fill_manual(values = mycol)+
  guides(fill=guide_legend(title = NULL,ncol=1,size=5))
p+geom_errorbar(aes(x=Var2,y=value,ymin=SEvalue-SEse,ymax=SEvalue+SEse),width=.2,position=position_dodge(.9))
dev.off()
