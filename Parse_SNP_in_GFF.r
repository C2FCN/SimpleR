library(rtracklayer)
df<-import("name.gff")
input<-read.table("name.snps.xls",header = T,sep = "\t",stringsAsFactors = F,check.names = F)
df<-df[-1]
type<-unique(df@elementMetadata["type"])[-3,]
input[,c(8,9)]<-NA

pastenew<-function(df){
  num<-length(df)
  tmp<-c()
  for (i in 1:num){
  tmp<-paste(df[i],tmp,sep = "|")
  }
  return(tmp)
}

for (i in 1:nrow(input)){
  tmp<-df[(input[i,3]>=df@ranges@start) & (input[i,3]<=(df@ranges@start+df@ranges@width-1))]
  tmptype<-tmp@elementMetadata["type"]
  if(length(tmp)==0){
      input[i,8]<-"Intergenic"
      input[i,9]<-"-"
  }else if(length(unique(tmp@ranges))==1){  
          if(sum(tmptype %in% type) >= 2 ){
              input[i,9]<-tmp[2]@elementMetadata["product"][[1]]
              input[i,8]<-as.vector(tmptype[2,])
          }else{
              input[i,9]<-tmp@elementMetadata["product"][[1]]
              input[i,8]<-as.vector(tmptype[[1]])
          }
  }else if(length(unique(tmp@ranges))==2){
    input[i,8]<-"CDS"
    input[i,9]<-pastenew(unique(tmp[(tmp@elementMetadata["type"] %in% "CDS")[[1]]]@elementMetadata["product"])[[1]])
  }
}
colnames(input)[8:9]<-c("Func","Description")
write.table(input,"name.snps.anno.xls",quote = F,sep = "\t",row.names = F,col.names = T)

