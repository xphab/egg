com.simp<-function(comi,new.limit=1000,resamp.meth=c("sequence","taxa"),check.cor=TRUE,times=10,out=c("max.r","mean"),nworker=1,code.wd)
{
  # to generate simplified community matrix by randomly resample sequence or taxa
  library(vegan)
  if(resamp.meth[1]=="sequence")
  {
    source(file = paste(code.wd,"/com.resamp.r",sep = ""))
    res.array=replicate(times,com.resamp(comi = comi,reads.limit = new.limit,nworker = nworker))
  }else{
    sp.name=colnames(comi)
    res.array=replicate(times,comi[,match(sample(sp.name,new.limit),sp.name)])
  }
  
  if(out[1]=="max.r")
  {
    old.dis=vegdist(comi)
    res.list=lapply(1:dim(res.array)[3],function(i,arr){arr[,,i]},arr=res.array)
    new.dis.list=lapply(res.list,vegdist)
    cor.res=sapply(new.dis.list,function(xx,yy){mantel(xx,yy)$statistic},yy=old.dis)
    max.cor=max(cor.res)
    message("-------------max Mantel r=",max.cor,".-----------------",date())
    res=res.array[,,match(max.cor,cor.res)]
  }else{
    res=apply(res.array,c(1,2),mean)
    if(check.cor)
    {
      old.dis=vegdist(comi)
      new.dis=vegdist(res)
      cor.res=mantel(old.dis,new.dis)$statistic
      message("Output is mean, the Mantel r=",cor.res,". ",date())
    }
  }
  res
}