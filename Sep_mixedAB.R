#5:9
#i:10-i
library("pracma")
for(i in 5:9){
fileName0<-paste("AB.",i,".",10-i,"_R",1,".fq",sep="")
fileName<-paste("/home/ali/Desktop/new.simulated data/",fileName0,sep="")
conn <- file(fileName,open="r")
linn <-readLines(conn)
close(conn)

j<-1
index=list()
index[j]<-1
s1<-substring(linn[1],1,2)
for (k in seq(1,length(linn),4)) {
  s2<-substring(linn[k],1,2)
  if (!strcmp(s1,s2))
    {
    j<-j+1
    index[j]<-k
    s1<-s2
  }
  print(substring(s2,1,2))
  
}
Outname0<-paste("/home/ali/Desktop/new.simulated data/seperated/",fileName0,sep="")
for(t in 1:(length(index)-1)){

  Outname<-paste(Outname0,t)
  fileConn<-file(Outname)
  writeLines( linn[ index[[t]]:index[[t+1]]-1 ],fileConn )
  close(fileConn)
  
}
t<-t+1
Outname<-paste(Outname0,t)
fileConn<-file(Outname)
writeLines( linn[ index[[t]]:length(linn) ],fileConn )
close(fileConn)
}
