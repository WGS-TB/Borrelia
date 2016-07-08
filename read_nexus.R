read.nexus.data("DatasetW1.nex.txt")
for (i in seq(1,32,by=1) ){
  s<-x[i]
  name<-paste("ref",i)
 write.dna(s,name,format ="fasta",colsep ="",append=FALSE )
  
}

