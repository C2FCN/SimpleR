library(seqinr)
#Input data
df<-read.table("CDS_clusters_sta.xls",header = F,sep="\t",row.names = 1,stringsAsFactors = F,fill = T)
fadf<-read.fasta(file = "all.cds.proteins.sorted.fasta",seqtype = "AA",as.string = T)
nafadf<-read.fasta(file = "all_nucleic_acid_sequences.fasta",seqtype = "DNA",as.string = T)

###rename geneid function
renamefun<-function(gta){
  #gta => gene tree fasta;agn=> all gene names from gta;dpgn=>duplicated gene names
  #!!!!!!!You should give correct Regular Expression to extract species name,here we use '_WP_[0-9.]*$|_orf[0-9_]*$'
  #!!!!!!!You should give correct Regular Expression to extract species name,here we use '_WP_[0-9.]*$|_orf[0-9_]*$'
  #!!!!!!!You should give correct Regular Expression to extract species name,here we use '_WP_[0-9.]*$|_orf[0-9_]*$'
  #!!!!!!!You should give correct Regular Expression to extract speciescd name,here we use '_WP_[0-9.]*$|_orf[0-9_]*$'
  #!!!!!!!You should give correct Regular Expression to extract species name,here we use '_WP_[0-9.]*$|_orf[0-9_]*$'
  dpindex<-which(duplicated(gsub("_WP_[0-9.]*$|_orf[0-9_]*$","",names(gta))))
  agn<-gsub("_WP_[0-9.]*$|_orf[0-9_]*$","",names(gta))
  if(length(dpindex)==0){
    names(gta)<-agn
    return(gta) 
  }else{
    dpgn<-agn[duplicated(agn)]
    for (n in 1:length(dpgn)){
      m<-length(agn[agn==dpgn[n]])
      agn[agn==dpgn[n]]<-paste(dpgn[n],seq(1,m),sep = "_")
    }
    names(gta)<-agn
    return(gta)
  }
}
#renamefun(gsub("_WP_[0-9.]*$|_orf[0-9_]*$","",names(genetreefa)))

#Initializtion 
count<-0
snum<-26
j=1
spe<-c()
single_CCC<-c()
genetreeid<-c()
core_gene_names<-c()
specific_record<-c()
specific_gene_df<-c()
#we have 26 species
spedf<-as.data.frame(matrix(0,ncol=2,nrow=snum+1))
clustermatrix<-as.data.frame(matrix(0,ncol=nrow(df),nrow=snum))
dir.create("Single_copy_core_genes")
dir.create("More_copy_core_genes")
dir.create("Morethan2lessthan900genes")
dir.create("More_copy_genes")
dir.create("Nucleic_acid_core_genes")
dir.create("Morethan900genes_AA")
dir.create("Morethan900genes_NA")
###define row and col names of clustermatrix
###check total number of species,use df[,1] or df[,2] or check yourself
tmpname<-unique(gsub("_WP_[0-9.]*$|_orf[0-9_]*$","",c(df[,1],unlist(strsplit(df[,1]," ")))))

if (length(tmpname)< snum){
  tmpname<-unique(gsub("_WP_[0-9.]*$|_orf[0-9_]*$","",c(df[,2],unlist(strsplit(df[,2]," ")))))
  try(if(length(tmpname) < snum) stop("please check the number of species!"))
}else if (length(tmpname) == snum){
  rownames(clustermatrix)<-tmpname
}
colnames(clustermatrix)<-df[,1]

#Parse core genes,specific genes and single copy core genes
for (i in 1:nrow(df)){
  allgenes<-gsub("_WP_[0-9.]*$|_orf[0-9_]*$","",c(df[i,1],unlist(strsplit(df[i,2]," "))))
  deduplication<-unique(gsub("_WP_[0-9.]*$|_orf[0-9_]*$","",c(df[i,1],unlist(strsplit(df[i,2]," ")))))
  
  if(length(deduplication)==snum){
    #record all core gene names for output of core gene fasta file(include singlecopy and multicopy)
    core_gene_names<-c(core_gene_names,df[i,1])
    #1:share,0:none
    clustermatrix[,i]<-1
    #Core genes
    #Single copy core genes can not be affected by deduplication.
    if(length(allgenes)==length(deduplication)){
      single_CCC<-c(single_CCC,df[i,1])
      ###generate gene tree(one cluster containg 26 genes) fasta file 
      ###Nucleic_acid_core_genes
      nagenetreefa<-nafadf[c(df[i,1],unlist(strsplit(df[i,2]," ")))]
      write.fasta(nagenetreefa,names = names(nagenetreefa),as.string = T,file.out =paste("./Nucleic_acid_core_genes/",paste(df[i,1],".fasta",sep=""),sep = "") ,open="w") 
      
      genetreefa<-fadf[c(df[i,1],unlist(strsplit(df[i,2]," ")))]
      genetreefa<-renamefun(genetreefa)    #rename fasta gene name
      genetreefa<-genetreefa[sort(names(genetreefa))] #sort name
      write.fasta(genetreefa,names = names(genetreefa),as.string = T,file.out =paste("./Single_copy_core_genes/",paste(df[i,1],".fasta",sep=""),sep = "") ,open="w") 
      genetreefa<-renamefun(genetreefa)
      write.fasta(genetreefa,names = names(genetreefa),as.string = T,file.out =paste("./Morethan2lessthan900genes/",paste(df[i,1],".fasta",sep=""),sep = "") ,open="w") 
    }else if (length(allgenes)<=900){
      ###generate gene tree(one cluster containg more than 26 genes) fasta file 
      ###Nucleic_acid_core_genes
      nagenetreefa<-nafadf[c(df[i,1],unlist(strsplit(df[i,2]," ")))]
      write.fasta(nagenetreefa,names = names(nagenetreefa),as.string = T,file.out =paste("./Nucleic_acid_core_genes/",paste(df[i,1],".fasta",sep=""),sep = "") ,open="w") 
      
      genetreefa<-fadf[c(df[i,1],unlist(strsplit(df[i,2]," ")))]
      write.fasta(genetreefa,names = names(genetreefa),as.string = T,file.out =paste("./More_copy_core_genes/",paste(df[i,1],".fasta",sep=""),sep = "") ,open="w") 
      #genetreefa<-renamefun(genetreefa)
      write.fasta(genetreefa,names = names(genetreefa),as.string = T,file.out =paste("./Morethan2lessthan900genes/",paste(df[i,1],".fasta",sep=""),sep = "") ,open="w") 
    }else if (length(allgenes)>900){
      ###Nucleic_acid_core_genes
      nagenetreefa<-nafadf[c(df[i,1],unlist(strsplit(df[i,2]," ")))]
      write.fasta(nagenetreefa,names = names(nagenetreefa),as.string = T,file.out =paste("./Morethan900genes_NA/",paste(df[i,1],".fasta",sep=""),sep = "") ,open="w") 
      
      genetreefa<-fadf[c(df[i,1],unlist(strsplit(df[i,2]," ")))]
      write.fasta(genetreefa,names = names(genetreefa),as.string = T,file.out =paste("./Morethan900genes_AA/",paste(df[i,1],".fasta",sep=""),sep = "") ,open="w") 
      }
    count<-count+1
  }else if(length(deduplication)==1){
    if (deduplication %in% specific_record){
      specific_gene_df[specific_gene_df[,1]==deduplication,2]<-paste(specific_gene_df[specific_gene_df[,1]==deduplication,2],df[i,1])
    }else{
      specific_record<-c(specific_record,deduplication)
      tmp_specific_gene_df<-c(deduplication,df[i,1])
      specific_gene_df<-rbind(specific_gene_df,tmp_specific_gene_df)
    }
    #2 or more genes from the same genome(one species) won't be condidered
    clustermatrix[which(tmpname==deduplication),i]<-1
    #Specific genes   
    if (deduplication %in% spe){
      spedf[which(spedf[,1]==deduplication),2]<-as.numeric(spedf[which(spedf[,1]==deduplication),2])+1
    }else{
      j=j+1
      spedf[j,]<-c(deduplication,1)
      spe<-c(spe,deduplication)
    }
  }else{
    clustermatrix[(tmpname %in% deduplication),i]<-1
    if ((length(allgenes)>=3) & (length(allgenes)<=900)){
      ###generate gene tree(one cluster containg 3 or more genes) fasta file 
      genetreefa<-fadf[c(df[i,1],unlist(strsplit(df[i,2]," ")))]
      write.fasta(genetreefa,names = names(genetreefa),as.string = T,file.out =paste("./More_copy_genes/",paste(df[i,1],".fasta",sep=""),sep = "") ,open="w")  
      #genetreefa<-renamefun(genetreefa)
      write.fasta(genetreefa,names = names(genetreefa),as.string = T,file.out =paste("./Morethan2lessthan900genes/",paste(df[i,1],".fasta",sep=""),sep = "") ,open="w") 
    }
    
  }
}

#colnames(genetreeid)<-c("Cluster_name","Num")
dim(single_CCC)<-c(length(single_CCC),1)
colnames(single_CCC)<-c("Single_Copy_Genes")
spedf[1,]<-c("Core_Genes",count)
colnames(spedf)<-c("Type","count")
tclustermatrix<-t(clustermatrix)
tclustermatrix<-rbind(rownames(clustermatrix),tclustermatrix)
rownames(tclustermatrix)[1]<-"Cluster"
clustermatrix<-rbind(colnames(clustermatrix),clustermatrix)
rownames(clustermatrix)[1]<-"Species"
single_CCC_fasta<-fadf[single_CCC]
core_gene_fasta<-fadf[core_gene_names]
colnames(specific_gene_df)<-c("species","specific_genes")
rownames(specific_gene_df)<-specific_gene_df[,1]


#Output
for (k in 1:nrow(specific_gene_df)){
  specific_gene_fastadf<-fadf[as.matrix(as.data.frame(strsplit(specific_gene_df[k,2]," ")))]
  write.fasta(specific_gene_fastadf,names = names(specific_gene_fastadf),as.string = T,file.out =paste0(specific_gene_df[k,1],"_specific.fasta") ,open="w") 
}
write.fasta(core_gene_fasta,names = names(core_gene_fasta),as.string = T,file.out ="Core_genes.fasta" ,open="w") 
write.fasta(single_CCC_fasta,names = names(single_CCC_fasta),as.string = T,file.out ="Single_Copy_Genes.fasta" ,open="w") 
write.table(single_CCC,"Single_Copy_Genes.xls",quote = F,sep = "\t",row.names = F,col.names = T)
write.table(spedf,"Core_and_specific_genes.xls",quote = F,sep = "\t",row.names = F,col.names = T)
write.table(clustermatrix,"cluster_gain_loss.xls",quote = F,sep = "\t",row.names = T,col.names = F)
write.table(tclustermatrix,"tcluster_gain_loss.xls",quote = F,sep = "\t",row.names = T,col.names = F)
