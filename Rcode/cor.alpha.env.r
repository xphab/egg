cor.alpha.env<-function(comi,env,cor.meth=c("pearson","kendall","spearman"),alpha.div=NA,index=NA,code.wd=NA)
{
  library(vegan)
  library(data.table)
  if(is.null(nrow(comi))&is.null(nrow(alpha.div)))
  {
    message("Neither com or alpha is available. return NA")
    output=NA
  }else{
    if(is.null(nrow(alpha.div))|sum(rownames(alpha.div)!=rownames(env))>0)
    {
      source(file=paste(code.wd,"/alpha.r",sep=""))
      if(sum(rownames(comi)!=rownames(env))>0)
      {
        samp.name=intersect(rownames(comi),rownames(env))
        comi=comi[match(samp.name,rownames(comi)),]
        env=env[match(samp.name,rownames(env)),]
        message("Some samples are not match in com and env matrix, so removed.")
      }
      alpha.div=alpha(comi)
    }
  }
  if(is.na(index)[1]){index=c(1:5,7)}
  
  alpha.list=lapply(index, function(i,xxx){xxx[,i]},xxx=alpha.div)
  names(alpha.list)=colnames(alpha.div)[index]
  env.list=lapply(1:ncol(env), function(i,env){env[,i]}, env=env)
  names(env.list)=colnames(env)
  
  cor.a.env<-function(i,env.list,alpha.list,cor.meth)
  {
    cor.a.a<-function(alphai,envi,cor.meth)
    {
      cor.a.m<-function(meth,xx,yy)
      {
        cor11<-cor.test(xx,yy,method=meth)
        out=c(cor11$estimate[[1]],cor11$p.value)
        out
      }
      coraa=unlist(lapply(cor.meth,cor.a.m,xx=alphai,yy=envi))
      coraa
    }
    envi=env.list[[i]]
    envi.name=names(env.list)[i]
    corae<-lapply(alpha.list,cor.a.a,envi=envi,cor.meth=cor.meth)
    res=data.frame(t(matrix(unlist(corae),nrow=length(corae[[1]]))))
    out=data.frame(Env=rep(envi.name,nrow(res)),Alpha=names(alpha.list),res)
    out
  }
  
  res<-lapply(1:ncol(env),cor.a.env,env.list=env.list,alpha.list=alpha.list,cor.meth=cor.meth)
  res=data.frame(rbindlist(res))
  colnames(res)=c("Factor","Alpha.Index",as.vector(sapply(cor.meth,paste,c("r","p"),sep=".")))
  res
}