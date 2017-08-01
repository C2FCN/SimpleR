library(ggplot2)
library(reshape2)
mytheme<-theme(legend.position="right",
               panel.grid.major.x=element_blank(),
               panel.grid.minor.x=element_blank(),
               panel.grid.major.y=element_blank(),
               panel.grid.minor.y=element_blank(),
               panel.background=element_rect(fill="white"),
               axis.text.x = element_text(angle=0, hjust = 0.5,vjust=0.5))

mycol<-c("#FFFB00","#5AB9D1","#ff7c00","#FF2600","#0522FF","#E98005","#008F00","#919191","#E2EDF8","#0093C6",
         "#672D00","#0D0D00","#4180FF","#C063BC","#EAC282")
df<-read.table("Phylum_ex_1_percent.xls",row.names = 1,header = T,check.names = F,sep = "\t",stringsAsFactors = F)
map<-read.table("map.txt",header = T,sep = "\t",stringsAsFactors = F)
group<-unique(map[,2])
t<-c()
for (i in 1:length(group)){
  tmp<-df[,as.character(map[map[,2]==group[i],1])]
  tt<-rowSums(tmp)
  t<-rbind(t,tt)
}
groupdf<-t(t)
colnames(groupdf)<-group
groupdf<-prop.table(as.matrix(groupdf),2)
zero<-rep(0,4)

groupdf<-cbind(groupdf,zero)
groupdf<-melt(groupdf)
groupdf[,2]<-factor(groupdf[,2],levels=c("zero","hindgut","foregut"))

pdf("Phylum_ex_1percent_Pie.pdf")
ggplot(groupdf)+geom_bar(aes(x=Var2,y=value,fill=Var1),stat = "identity",width =1,color=c("black"),size=1.3)+
  xlab(NULL)+ylab("Relative abundance(%)")+
  scale_fill_manual(values = mycol)+mytheme+
  guides(fill=guide_legend(title = NULL,ncol=1,size=5))+coord_polar(theta = "y")
dev.off()
