cca.test<-function(com.test,env,cases=NA,code.wd)
{
  if(sum(rownames(com.test)!=rownames(env))>0)
  {
    samp.name=intersect(rownames(com.test),rownames(env))
    com.test=com.test[match(samp.name,rownames(com.test)),]
    env=env[match(samp.name,rownames(env)),]
    message("some sample names are not the same in com.test and env, thus discarded")
  }
  library(vegan)
  env.num=ncol(env)
  env.name=colnames(env)
  if(is.na(cases)[1])
  {
    case.num=2^env.num-1
    source(file = paste(code.wd,"/intTobit.r",sep = ""))
    cases=lapply(1:case.num, intTobit,len=env.num)
  }else{
    if(ncol(cases)!=env.num)
    {
      message("the element number in cases is not the same as env number.")
      break
    }else{
      case.num=nrow(cases)
      cases=lapply(1:case.num,function(i,x){x[i,]},x=as.matrix(cases))
    }
  }

  envs<-function(c,env){env[c==1]}
  env.case=lapply(cases,envs,env=env)

  cca.t<-function(X,Y)
  {
    library(vegan)
    A<-cca(Y~.,data=X)
    prop=sum(A$CCA$eig)/A$tot.chi
    p=anova(A)$'Pr(>F)'[1]
    mv=max(vif.cca(A))
    res=c(mv,p,prop)
    res
  }
  
  cca.index=lapply(env.case,cca.t,Y=com.test)
  cca.index=t(matrix(unlist(cca.index),nrow=3))
  case.l=t(matrix(unlist(cases),nrow=env.num))
  res=data.frame(cbind(case.l,cca.index))
  colnames(res)=c(env.name,"max_Vif","p_value","proportion")
  res
}