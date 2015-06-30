cor.cate.env<-function(comi,category,env,nworker=4)
{
  library(vegan)
  library(data.table)
  # match names
  if(sum(colnames(comi)!=rownames(category))>0)
  {
    sp.name=intersect(colnames(comi),rownames(category))
    comi=comi[,match(sp.name,colnames(comi))]
    category=category[match(sp.name,rownames(category)),]
  }
  if(sum(rownames(comi)!=rownames(env))>0)
  {
    samp.name=intersect(rownames(comi),rownames(env))
    comi=comi[match(samp.name,rownames(comi)),]
    env=env[match(samp.name,rownames(env)),]
  }
  
  # function
  cor.cate<-function(i,env.list,comm,catelist,nworker)
  {
    cor.cateone<-function(j,catesamp.list,comm,env1)
    {
      cate1.samp=catesamp.list[[j]]
      message("j=",j)
      library(vegan)
      com.cate=comm[,match(cate1.samp,colnames(comm)),drop=FALSE]
      com.sum=rowSums(com.cate)
      cor1=cor.test(com.sum,env1,method = "pearson")
      cor2=cor.test(com.sum,env1,method = "kendall")
      cor3=cor.test(com.sum,env1,method = "spearman")
      
      com.dis.eu=dist(com.cate)
      com.dis.bray=vegdist(com.cate)
      env1.dis=dist(env1)
      cor4=mantel(com.dis.eu,env1.dis,method = "pearson")
      cor5=mantel(com.dis.eu,env1.dis,method = "kendall")
      cor6=mantel(com.dis.eu,env1.dis,method = "spearman")
      cor7=mantel(com.dis.bray,env1.dis,method = "pearson",na.rm = TRUE)
      cor8=mantel(com.dis.bray,env1.dis,method = "kendall",na.rm = TRUE)
      cor9=mantel(com.dis.bray,env1.dis,method = "spearman",na.rm = TRUE)
      
      if(sum(com.cate)==0)
      {
        cca.prop=NA
        cca.p=NA
      }else{
        id.c=(rowSums(com.cate)>0)
        com.cate.c=com.cate[id.c,]
        env1.c=env1[id.c]
        cca.v=cca(com.cate.c~env1.c)
        cca.prop=sum(cca.v$CCA$eig)/cca.v$tot.chi
        cca.p=anova(cca.v)$'Pr(>F)'[1]
      }
      
      
      res=c(cor1$estimate[[1]],cor1$p.value, cor2$estimate[[1]],cor2$p.value, cor3$estimate[[1]],cor3$p.value,
            cor4$statistic,cor4$signif, cor5$statistic,cor5$signif, cor6$statistic,cor6$signif, 
            cor7$statistic,cor7$signif, cor8$statistic,cor8$signif, cor9$statistic,cor9$signif,
            cca.prop,cca.p)
      res
    }
    
    env1=env.list[[i]]
    env1.name=names(env.list)[i]
    message("now test correlation between each category and ", env1.name,". i=", i," in ",length(env.list),". ",date())
    library(parallel)
    gc()
    c1<-makeCluster(nworker,type="PSOCK")
    cor.res<-parLapply(c1,1:length(catelist),cor.cateone,catesamp.list=catelist,comm=comm,env1=env1)
    stopCluster(c1)
    gc()
    
    cor.res<-data.frame(t(matrix(unlist(cor.res),nrow = length(cor.res[[1]]))))
    cor.out=data.frame(Env=rep(env1.name,nrow(cor.res)),Category=names(catelist),cor.res)
    cor.out
  }
  
  # calculation
  env.list=lapply(1:ncol(env), function(i,env){env[,i]}, env=env)
  names(env.list)=colnames(env)
  cate.list<-function(cate1,category){rownames(category)[category[,1]==cate1]}
  cate.lev=levels(as.factor(as.matrix(category[,1])))
  catelist<-lapply(cate.lev,cate.list,category)
  names(catelist)=cate.lev
  if(nworker>length(cate.lev)){nworker=length(cate.lev)}
  
  result=lapply(1:ncol(env), cor.cate, env.list=env.list,comm=comi,catelist=catelist,nworker=nworker)
  result=rbindlist(result)
  colnames(result)=c("Factor","Taxa.category","sum.Pearson.r","sum.Pearson.p","sum.Kendall.r","sum.Kendall.p","sum.Spearman.r","sum.Spearman.p",
                     "Mantel.Eucl.Pearson.r","Mantel.Eucl.Pearson.p","Mantel.Eucl.Kendall.r","Mantel.Eucl.Kendall.p","Mantel.Eucl.Spearman.r","Mantel.Eucl.Spearman.p",
                     "Mantel.Bray.Pearson.r","Mantel.Bray.Pearson.p","Mantel.Bray.Kendall.r","Mantel.Bray.Kendall.p","Mantel.Bray.Spearman.r","Mantel.Bray.Spearman.p",
                     "CCA.proportion","CCA.p")
  result
}