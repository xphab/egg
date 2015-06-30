cor.v<-function(vect,env,method=c("pearson","kendall","spearman"))
{
  # correlation between vector and each env or all by cor test
  library(vegan)
  vect=as.vector(as.matrix(vect))
  env.num=ncol(env)
  env.l=lapply(1:env.num, function(i,env){env[,i]},env=env)
  vtest<-function(X,Y,meth)
  {
    test=cor.test(X,Y,method = meth)
    (res=c(test$estimate[[1]],test$p.value))
  }
  res=data.frame(matrix(nrow = env.num,ncol = 2*length(method)))
  for(i in 1:length(method))
  {
    test=lapply(env.l, vtest,Y=vect,meth=method[i])
    res[,(2*i-1):(2*i)]=t(matrix(unlist(test),nrow = 2))
  }
  rownames(res)=c(colnames(env))
  colnames(res)[(2*(1:length(method)))]=paste("p.",method,sep = "")
  colnames(res)[(2*(1:length(method))-1)]=paste("r.",method,sep = "")
  res
}