cor.beta.env<-function(comi,env,cor.meth=c("pearson","kendall","spearman"),code.wd)
{
  library(vegan)
  library(data.table)
  if(sum(rownames(comi)!=rownames(env))>0)
  {
    samp.name=intersect(rownames(comi),rownames(env))
    comi=comi[match(samp.name,rownames(comi)),]
    env=env[match(samp.name,rownames(env)),]
    message("Some samples are not match in com and env matrix, so removed.")
  }
  
  source(file = paste(code.wd,"/beta.dis.r",sep = ""))
  beta.div=beta.dis(comi,binary = "Both",out.format = "dist")
  
  cor.b.e<-function(i,env.list,beta.list,cor.meth)
  {
    cor.b.b<-function(betai,envi,cor.meth)
    {
      cor.b.m<-function(meth,betai,envi)
      {
        envi.dis=dist(envi)
        corbm=mantel(betai,envi.dis,method = meth,na.rm = TRUE)
        res=c(corbm$statistic,corbm$signif)
      }
      corbb=unlist(lapply(cor.meth, cor.b.m,betai=betai,envi=envi))
      corbb
    }
    envi=env.list[[i]]
    envi.name=names(env.list)[i]
    message("Now mantel test between beta distances and ", envi.name,". env i=",i," in ", length(env.list),". ",date())
    corbe=lapply(beta.list, cor.b.b,envi=envi,cor.meth=cor.meth)
    res=data.frame(t(matrix(unlist(corbe),nrow=length(corbe[[1]]))))
    out=data.frame(Env=rep(envi.name,nrow(res)),Beta=names(beta.list),res)
    out
  }
  
  env.list=lapply(1:ncol(env), function(i,env){env[,i]}, env=env)
  names(env.list)=colnames(env)
  res<-lapply(1:ncol(env),cor.b.e,env.list=env.list,beta.list=beta.div,cor.meth=cor.meth)
  res=data.frame(rbindlist(res))
  res.name=c("Factor","beta.index",as.vector(sapply(paste("Mantel",cor.meth,sep = "."), paste, c("r","p"),sep=".")))
  setnames(res,res.name)
  res
}