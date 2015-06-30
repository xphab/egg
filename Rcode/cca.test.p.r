cca.test.p<-function(com.test,env,cases=NA,summ=FALSE,nworker=4,memory.G=30,list.limit=8000,code.wd)
{
  if(sum(rownames(com.test)!=rownames(env))>0)
  {
    samp.name=intersect(rownames(com.test),rownames(env))
    com.test=com.test[match(samp.name,rownames(com.test)),]
    env=env[match(samp.name,rownames(env)),]
    message("some sample names are not the same in com.test and env, thus discarded")
  }
  memory.limit(size = memory.G*1024)
  library(vegan)
  library(parallel)
  env.num=ncol(env)
  env.name=colnames(env)
  samp.num=nrow(com.test)
  if(is.na(cases)[1])
  {
    f.num=min(env.num,samp.num-1)
    case.num=2^env.num-1
    source(file = paste(code.wd,"/intTobit.r",sep = ""))
    cases=lapply(1:case.num, intTobit,len=env.num)
    casum=unlist(lapply(cases,sum))
    cases=cases[casum<=f.num]
  }else{
    if(sum(!(colnames(cases) %in% env.name))>0)
    {
      message("the factor names in cases do not match env name.")
      break
    }else{
      if(sum(colnames(cases)!=env.name)>0)
      {
        cases=cases[,match(env.name,colnames(cases))]
        cases[is.na(cases)]=0
      }
      case.num=nrow(cases)
      cases=lapply(1:case.num,function(i,x){x[i,]},x=as.matrix(cases))
    }
  }
  message("Total case numer is ", length(cases), ". ",date())
  cca.t<-function(i,com,env,cases,summ)
  {
    library(vegan)
    X=env[cases[[i]]==1]
    A<-cca(com~.,data=X)
    prop=sum(A$CCA$eig)/A$tot.chi
    p=anova(A)$'Pr(>F)'[1]
    mv=max(vif.cca(A))
    res=c(cases[[i]],mv,p,prop)
    if(summ)
    {
      AS=summary(A)
      res=list(res=res,sp=AS$species,st=AS$sites,en=AS$biplot,pr=t(AS$cont[[1]]))
    }
    res
  }
  
  # time estimate
  message("test time cost of one run. begin at ", date())
  t1=Sys.time()
  cca.index=cca.t(i=1,com=com.test,env,cases,summ)
  t2=Sys.time()
  message("test time cost of one run. end at ", date())
  message("time difference of ", difftime(t2,t1,units = "secs"), " secs.")
  gc()
  
  par.time=ceiling(case.num/list.limit)
  case.num.i=matrix(nrow = par.time,ncol = 2)
  case.num.i[,1]=((0:(par.time-1))*list.limit+1)
  case.num.i[,2]=(case.num.i[,1]+list.limit-1)
  case.num.i[par.time,2]=case.num
  res=list()
  
  for(i in 1:par.time)
  {
    message("Now test CCA of cases ",case.num.i[i,1]," to ",case.num.i[i,2],". ",date())
    nw=min(nworker,(case.num.i[i,2]-case.num.i[i,1]+1))
    c1<-makeCluster(nw,type="PSOCK")
    res.i<-parLapply(c1,case.num.i[i,1]:case.num.i[i,2],cca.t,com=com.test,env=env,cases=cases,summ=summ)
    stopCluster(c1)
    gc()
    res[case.num.i[i,1]:case.num.i[i,2]]=res.i
  }
  
  if(summ)
  {
    res1=t(sapply(1:case.num,function(i,L){L[[i]][[1]]},L=res))
    res1=data.frame(res1)
    colnames(res1)=c(env.name,"max_Vif","p_value","proportion")
    output=list(index=res1,detail=res)
  }else{
    res=t(matrix(unlist(res),nrow=(env.num+3)))
    output=data.frame(res)
    colnames(output)=c(env.name,"max_Vif","p_value","proportion")
  }
  output
}