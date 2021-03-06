minread<-function(comi,samplist=NA,prefix,read.limit,write.otu=c(TRUE,FALSE),code.wd)
{
# from the raw OTU table to new OTU table with correct names and also return min reads for resampling

comi=data.frame(comi)
if(is.null(nrow(samplist)))
{
comc=comi
}else{
source(file=paste(code.wd,"/rename.samp.r",sep=""))
comc=rename.samp(comi,samplist)
}
reads=rowSums(comc)
if(is.na(read.limit))
{
  com.ok=comc
  minread=min(reads)
  samp.not=NA
}else{
  id.reads.ok=which(reads>=read.limit)
  id.reads.not=which(reads<read.limit)
  com.ok=comc[id.reads.ok,]
  minread=min(reads[id.reads.ok])
  samp.not=reads[id.reads.not]
}
if(write.otu)
{
  write.table(cbind(OTU=colnames(com.ok),t(com.ok)),file=paste("output/",prefix,".00.rawOTU.txt",sep=""),sep="\t",row.names=FALSE,quote=FALSE)
  write.csv(data.frame(samp.not),file=paste("output/",prefix,".00.SampNotInclude.csv",sep=""))
  write(paste("minimum reads number is ",minread,sep=""),file=paste("output/",prefix,".00.MinReads.txt",sep=""))
}
output=list(sample.not=samp.not,min.read=minread,read.samp=reads,com.all=comc,com.okay=com.ok)
output
}