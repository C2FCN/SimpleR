#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
my %opts;
my $VERSION = "XXXX";
GetOptions( \%opts,"i=s","o=s","w=f","h=f");
my $usage = <<"USAGE";
       Program : $0  
       Discription:Ternary Plot
       Version : $VERSION
       Contact : XXXXX
       Usage :perl $0 [options]
                        -i      * Input otu table file (e.g.: relative abundance otu table) .
                        -o      * Output file name (e.g.:Ternary.pdf)
			-w      Default:12 .
                        -h      Default:6 .
Example:perl plot-ternary.pl -i otu_table.xls -o Ternary.pdf
USAGE

die $usage if (!(defined $opts{i} && $opts{o}));
$opts{w}=defined $opts{w}?$opts{w}:10;
$opts{h}=defined $opts{h}?$opts{h}:6;
open CMD,">cmd.r";
print CMD "
library(ggtern)
df<-read.table(\"$opts{i}\",header = T,sep = \"\\t\",row.names = 1,stringsAsFactors = F,check.names = F)
df[,4]<-apply(df,1,mean)
sample<-colnames(df)
colnames(df)[4]<-\"Abundance\"
ggtern(data=df,aes(x=sample[1],y=sample[2], z=sample[3]))+
  geom_point(aes(fill = rownames(df), size = Abundance),colour=\"transparent\", shape = 21)+
  labs(fill = \"Type\")+scale_size_continuous(guide=F)+
  theme(text=element_text(size=16),
        tern.axis.line.T=element_line(size=0.5,colour=\"black\"),
        tern.axis.line.L=element_line(size=0.5,colour=\"black\"),
        tern.axis.line.R=element_line(size=0.5,colour=\"black\"),
        axis.title=element_text(size=10,face=\"bold\"))
ggsave(\"$opts{o}\",width = \"$opts{w}\",height = \"$opts{h}\")
";
`R --restore --no-save < cmd.r`;
#`rm cmd.r`;
