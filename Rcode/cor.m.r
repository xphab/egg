cor.m<-function(dis,env,method=c("pearson","kendall","spearman"))
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
  mtest<-function(X,Y,meth)
  {
    test=mantel(X,Y,method = meth)
    (res=c(test$statistic,test$signif))
  }
  
  res=data.frame(matrix(nrow = env.num+1,ncol = 2*length(method)))
  for(i in 1:length(method))
  {
    test=lapply(env.d, mtest,Y=dis,meth=method[i])
    res[,(2*i-1):(2*i)]=t(matrix(unlist(test),nrow = 2))
  }
  rownames(res)=c(colnames(env),"All")
  colnames(res)[(2*(1:length(method)))]=paste("p.mantel.",method,sep = "")
  colnames(res)[(2*(1:length(method))-1)]=paste("r.mantel.",method,sep = "")
  res
}