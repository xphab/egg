ieg.upload<-function(comi,treat,prefix,category=NA)
{
  comi=comi[,colSums(comi)>0]
  sp.name=colnames(comi)
  samp.name=rownames(comi)
  treat=treat[match(samp.name,rownames(treat)),]
  if(is.na(category))
  {
    category=sp.name
  }
  comi[comi==0]=""
  comm=data.frame(cate=category,t(comi))
  write.table(data.frame(GeneID=rownames(comm),comm),file=paste("output/",prefix,".ieg.comm.txt",sep=""),quote=FALSE,sep="\t",row.names=FALSE)
  tr.lev=levels(as.factor(treat[,1]))
  treatment=data.frame(treat=matrix(,nrow=length(tr.lev),ncol=1))
  for(i in 1:length(tr.lev))
  {
    samp=paste(rownames(treat)[treat[,1]==tr.lev[i]],collapse=",")
    treatment[i,1]=paste("T",tr.lev[i],":",samp,sep="")
  }
  write.table(treatment,file=paste("output/",prefix,".ieg.treat.txt",sep=""),quote=FALSE,row.names=FALSE,col.names=FALSE)
  output=list(com=comm,treatment=treatment)
  output
}