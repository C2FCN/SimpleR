library(ggplot2)
library(grid)
otu<-read.table("table.xls",check.names = F,row.names = 1,header = T,sep = "\t",stringsAsFactors = F)
map<-read.table("group.txt",check.names = F,header = T,sep = "\t",stringsAsFactors = F)
otu<-otu[,as.character(map[,1])]
otu<-otu[apply(otu,1,sum)!=0,]
pca<-prcomp(t(otu),scal=TRUE)
pc123 <-pca$x[,1:3]
pc123<-as.data.frame(pc123)
pc123[,4]<-map[,2]
colnames(pc123)[4]<-"Group"
pc <-summary(pca)$importance[2,]*100
p1<-ggplot(pc123,aes(x=PC1,y=PC2,colour=Group))+geom_point()+stat_ellipse(lwd=1)+
    theme_grey()+xlab(paste("\nPC1-Percent variant explained ",round(pc[1],2),"%",sep="") )+
    ylab(paste("\nPC2-Percent variant explained ",round(pc[2],2),"%\n",sep=""))+
    ggtitle("PCA-PC1 vs PC2\n")
p2<-ggplot(pc123[,c(1,4)],aes(x=Group,y=PC1,fill=Group))+geom_boxplot()+theme_grey()+theme(legend.position="none")
p3<-ggplot(pc123[,c(2,4)],aes(x=Group,y=PC2,fill=Group))+geom_boxplot()+theme_grey()+coord_flip()+theme(legend.position="none")
write.table(pc123[,1:3],file ="pca_data.xls",quote = F,sep = "\t",row.names = T,col.names = T )
pdf("PCA.pdf")
grid.newpage()  
pushViewport(viewport(layout = grid.layout(4,4))) 
vplayout <- function(x,y){
  viewport(layout.pos.row = x, layout.pos.col = y)
}
print(p1, vp = vplayout(1:3,1:3))  
print(p2, vp = vplayout(1:3,4))   
print(p3, vp = vplayout(4,1:1:4)) 
dev.off()
