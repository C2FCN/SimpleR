library(ggplot2)
library(gtable)
library(grid)
#box data
df<-read.table("JTYP2.subreads.len",header = F)
df<-as.data.frame(as.numeric(sort(df[,1])))
colnames(df)<-"sublen"
#line
totalbase<-sum(df[,1])
msublen<-c()
for (i in 0:ceiling(max(df[,1]/100))){
  tmp<-c(totalbase-sum(df[df[,1]<=i*100,1]),i*100)
  msublen<-rbind(msublen,tmp)
}
msublen[,1]<-msublen[,1]/1000000
colnames(msublen)<-c("value","binnum")
msublen<-as.data.frame(msublen)

# two plots
p1 <- ggplot()+
  geom_histogram(data=df,aes(x=sublen,y = ..count..),colour = "#9BD280", fill = "#51B121", binwidth = 150) + 
  theme_classic()+xlab('Subread Length')+ylab('Subreads')+xlim(0,25000)
p2 <-ggplot()+geom_line(data=msublen,aes(x=binnum,y=value),stat = 'identity', colour="black") + xlim(0,25000)+
  theme_classic() %+replace% theme(panel.background = element_rect(fill = NA) )+
  ylab('Mb > Subread Length')

# extract gtable
g1 <- ggplot_gtable(ggplot_build(p1))
g2 <- ggplot_gtable(ggplot_build(p2))

# overlap the panel of 2nd plot on that of 1st plot
pp <- c(subset(g1$layout, name == "panel", se = t:r))
g <- gtable_add_grob(g1, g2$grobs[[which(g2$layout$name == "panel")]], pp$t, 
                     pp$l, pp$b, pp$l)

# axis tweaks
ia <- which(g2$layout$name == "axis-l")
ga <- g2$grobs[[ia]]
ax <- ga$children[[2]]
ax$widths <- rev(ax$widths)
ax$grobs <- rev(ax$grobs)
ax$grobs[[1]]$x <- ax$grobs[[1]]$x - unit(1, "npc") + unit(0.15, "cm")
g <- gtable_add_cols(g, g2$widths[g2$layout[ia, ]$l], length(g$widths) - 1)
g <- gtable_add_grob(g, ax, pp$t, length(g$widths) - 1, pp$b)


# draw it
#pdf("123.pdf")
grid.draw(g)
#dev.off()
