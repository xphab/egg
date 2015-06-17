bmntd<-function(comm, pd, abundance.weighted = TRUE, exclude.conspecifics = FALSE,time.output=FALSE)
{
  if(sum(colnames(comm)!=colnames(pd))>0)
  {
    sp.name=colnames(pd)
    comm=comm[,match(sp.name,colnames(comm))]
  }
  comt=comm
  comt[comt>0]=1
  if(!abundance.weighted){com.10=comt}
  N=nrow(comm)
  time1=Sys.time()
  min.d=comm[1,]
  if(exclude.conspecifics)
  {
    pd[pd==0]=NA
    for(i in 1:N)
    {
      id=(comm[i,]==0)
      min.d=apply(pd[!id,,drop=FALSE],2,min,na.rm=TRUE)
      comt[i,]=min.d
    }
  }else{
    for(i in 1:N)
    {
      id=(comm[i,]==0)
      min.d[!id]=0
      min.d[id]=apply(pd[!id,id,drop=FALSE],2,min)
      comt[i,]=min.d
    }
  }
  time2=Sys.time()
  if(abundance.weighted)
  {
    comm.p=comm/rowSums(comm)
    time3=Sys.time()
    res=comt %*% (t(comm.p))
    time4=Sys.time()
    res=(res+t(res))/2
  }else{
    res=comt %*% (t(com.10))
    time3=Sys.time()
    samp.n=rowSums(com.10)
    com.n=matrix(samp.n,nrow = N,ncol = N)
    com.n=com.n+t(com.n)
    time4=Sys.time()
    res=(res+t(res))/com.n
  }
  res=as.dist(res)
  time5=Sys.time()
  if(time.output)
  {
    time=c(time5,time4,time3,time2)-c(time4,time3,time2,time1)
    output=list(result=res,time=time)
  }else{
    output=res
  }
  output
}