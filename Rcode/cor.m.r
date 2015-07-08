cor.m<-function(dis,env,nworker=4,method=c("pearson","kendall","spearman"))
{
  # correlation between dis and each env or all by mantel test
  dis=as.matrix(dis)
  if(sum(rownames(dis)!=rownames(env))>0)
  {
    samp.name=intersect(rownames(dis),rownames(env))
    dis=dis[match(samp.name,rownames(dis)),match(samp.name,rownames(dis))]
    env=env[match(samp.name,rownames(env)),]
  }
  library(vegan)
  env.num=ncol(env)
  envdis<-function(i,env){dist(env[,i])}
  env.d=lapply(1:env.num, envdis, env=env)
  env.d[[length(env.d)+1]]=dist(env)
  mtest<-function(i,X.l,Y,meth,nworker)
  {
    X=X.l[[i]]
    message("now mantel test, i=",i," in ",length(X.l),". ",date())
    test=mantel(X,Y,method = meth,na.rm = TRUE,parallel = nworker)
    (res=c(test$statistic,test$signif))
  }
  
  res=data.frame(matrix(nrow = env.num+1,ncol = 2*length(method)))
  for(i in 1:length(method))
  {
    test=lapply(1:length(env.d), mtest,X.l=env.d,Y=dis,meth=method[i],nworker=nworker)
    res[,(2*i-1):(2*i)]=t(matrix(unlist(test),nrow = 2))
  }
  rownames(res)=c(colnames(env),"All")
  colnames(res)[(2*(1:length(method)))]=paste("p.mantel.",method,sep = "")
  colnames(res)[(2*(1:length(method))-1)]=paste("r.mantel.",method,sep = "")
  res
}