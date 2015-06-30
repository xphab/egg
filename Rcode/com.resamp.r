com.resamp<-function(comi,reads.limit=10000,nworker=1)
{
  # resampling
  comt=t(comi)
  old.reads=colSums(comt)
  id.yes=(old.reads>reads.limit)
  if(sum(id.yes)==0)
  {
    message("The reads in all samples are less than reads.limit. Resampling is not applicable.")
    comi
  }else{
    com.res=data.frame(comt)
    comc=comt[,id.yes,drop=FALSE]
    
    resamp<-function(i,comc,reads.limit)
    {
      comv=comc[,i]
      id=(comv>0)
      comv.cum=cumsum(comv[id])
      seq.num=sum(comv)
      resamp=sample(1:seq.num,reads.limit)
      h=hist(resamp,breaks=c(0,comv.cum),plot = FALSE)
      comv[id]=h$counts
      comv
    }
    
    if(ncol(comc)<nworker){nworker=ncol(comc)}
    
    if(nworker==1)
    {
      res=sapply(1:ncol(comc),resamp,comc=comc,reads.limit=reads.limit)
      com.res[,id.yes]=data.frame(res)
    }else{
      library(parallel)
      c1<-makeCluster(nworker,type="PSOCK")
      res=parSapply(c1,1:ncol(comc),resamp,comc=comc,reads.limit=reads.limit)
      stopCluster(c1)
      gc()
      com.res[,id.yes]=data.frame(res)
    }
    t(com.res)
  }
}