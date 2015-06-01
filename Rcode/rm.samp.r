rm.samp<-function(comi,rm.samplist.file)
{
  rm.samplist1=read.table(file=paste("input/",rm.samplist.file,sep=""),header=T,sep=",",row.names=1,stringsAsFactors=FALSE)
  jd=match(rownames(comi),rm.samplist1[,1])
  comm=comi[is.na(jd),]
  comm
}