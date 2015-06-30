beta.dis<-function(comi,method=c("euclidean","manhattan","jaccard","bray","binomial"),binary=c(FALSE,TRUE,"Both"),out.format=c("3col","dist"))
{
  library(vegan)
  beta.in<-function(meth,comi,bin,out.f)
  {
    dis=vegdist(comi,method=meth,binary = bin)
    if(out.f=="3col")
    {
      res=as.vector(dis)
    }else{
      res=dis
    }
    res
  }
  samp.num=nrow(comi)
  rown=matrix(1:samp.num,nrow=samp.num,ncol = samp.num)
  rown.v=as.vector(as.dist(rown))
  coln.v=as.vector(as.dist(t(rown)))
  name1=rownames(comi)[rown.v]
  name2=rownames(comi)[coln.v]
  
  if(binary[1]=="Both")
  {
    beta1=lapply(method, beta.in, comi=comi,bin=FALSE,out.f=out.format[1])
    beta2=lapply(method, beta.in, comi=comi,bin=TRUE,out.f=out.format[1])
    if(out.format[1]=="3col")
    {
      res=data.frame(matrix(c(unlist(beta1),unlist(beta2)),nrow=length(beta1[[1]])))
      colnames(res)=c(paste(method,"abun.weight",sep = "."),paste(method,"not.weight",sep = "."))
      output=data.frame(sample1=name1,sample2=name2,res)
    }else{
      output=c(beta1,beta2)
      names(output)=c(paste(method,"abun.weight",sep = "."),paste(method,"not.weight",sep = "."))
    }
  }else{
    beta1=lapply(method, beta.in, comi=comi,bin=binary[1],out.f=out.format[1])
    if(binary[1]){weight="not.weight"}else{weight="abun.weight"}
    if(out.format[1]=="3col")
    {
      res=data.frame(matrix(unlist(beta1),nrow=length(beta1[[1]])))
      colnames(res)=paste(method,weight,sep = ".")
      output=data.frame(sample1=name1,sample2=name2,res)
    }else{
      output=beta1
      names(output)=paste(method,weight,sep = ".")
    }
  }
  output
}