bmpd<-function(comm, pd, abundance.weighted = TRUE, time.output=FALSE)
{
  if(sum(colnames(comm)!=rownames(pd))>0)
  {
    sp.name=intersect(colnames(comm),rownames(pd))
    comm=comm[,match(sp.name,colnames(comm))]
    pd=pd[match(sp.name,rownames(pd)),match(sp.name,rownames(pd))]
  }
  comt=comm
  if(!abundance.weighted){comt[comt>0]=1}
  N=nrow(comm)
  time1=Sys.time()
  comt=comt/rowSums(comt)
  comt=as.matrix(comt)
  pd=as.matrix(pd)
  time2=Sys.time()
  comd=comt %*% pd
  time3=Sys.time()
  res=comd %*% t(comt)
  time4=Sys.time()
  res=(res+t(res))/2
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