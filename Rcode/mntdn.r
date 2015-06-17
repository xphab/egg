mntdn<-function(comm, pd, abundance.weighted = TRUE)
{
  if(sum(colnames(comm)!=colnames(pd))>0)
  {
    sp.name=colnames(pd)
    comm=comm[,match(sp.name,colnames(comm))]
  }
  N=nrow(comm)
  id=(comm>0)
  diag(pd)=NA
  if(abundance.weighted)
  {
    min.d=matrix(0,nrow = N,ncol = ncol(comm))
    for(i in 1:N)
    {
      pdx=pd[id[i,],id[i,],drop=FALSE]
      min.d[i,id[i,]]=apply(pdx,2,min,na.rm=TRUE)
    }
    comm.p=comm/rowSums(comm)
    res=min.d * comm.p
    res=rowSums(res)
    res
  }else{
    res=comm[,1]
    for(i in 1:N)
    {
      pdx=pd[id[i,],id[i,],drop=FALSE]
      res[i]=mean(apply(pdx,2,min,na.rm=TRUE))
    }
    res
  }
}