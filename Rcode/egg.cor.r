egg.cor<-function(comi,category,env,alpha.env.yn=1,beta.env.yn=1,cate.env.yn=1,nworker=4,write.output=TRUE,prefix,code.wd)
{
  library(vegan)
  library(data.table)
  # match names
  comi=comi[rowSums(comi)>0,colSums(comi)>0]
  if(!is.null(nrow(category)))
  {
    if(sum(colnames(comi)!=rownames(category))>0)
    {
      sp.name=intersect(colnames(comi),rownames(category))
      comi=comi[,match(sp.name,colnames(comi))]
      category=category[match(sp.name,rownames(category)),]
    }
  }
  if(sum(rownames(comi)!=rownames(env))>0)
  {
    samp.name=intersect(rownames(comi),rownames(env))
    comi=comi[match(samp.name,rownames(comi)),]
    env=env[match(samp.name,rownames(env)),]
  }
  
  # Pearson, Kendall, Spearman, Mantel, CCA between whole community and each env factor
  if(alpha.env.yn==0)
  {
    cor.alpha=NA
  }else{
    message("Now correlation test between alpha indexes and each env factor. ",date())
    source(file = paste(code.wd,"/cor.alpha.env.r",sep=""))
    cor.alpha=cor.alpha.env(comi = comi,env = env,code.wd = code.wd)
    if(write.output){write.csv(cor.alpha,paste("output/",prefix,".cor.01.env_VS_alpha.csv",sep=""))}
  }
  
  if(beta.env.yn==0)
  {
    cor.beta=NA
  }else{
    message("Now correlation test between beta indexes and each env factor. ",date())
    source(file = paste(code.wd,"/cor.beta.env.r",sep=""))
    cor.beta=cor.beta.env(comi=comi,env = env,code.wd = code.wd)
    if(write.output){write.csv(cor.beta,paste("output/",prefix,".cor.02.env_VS_beta.csv",sep=""))}
  }
  
  # Pearson, Kendall, Spearman, Mantel, CCA between each member of each category and each env factor
  if(cate.env.yn==0|is.null(nrow(category)))
  {
    cor.cate=NA
  }else{
    source(file = paste(code.wd,"/cor.cate.env.r",sep = ""))
    cor.cate=list()
    cate.name=colnames(category)
    for(i in 1:ncol(category))
    {
      message("Now correlation test between each env and each catogory of---",cate.name[i],"----. cate i=",i," in ", ncol(category),". ", date())
      cor.cate[[i]]<-cor.cate.env(comi=comi,category=category[,i,drop=FALSE],env=env,nworker=nworker)
      if(write.output){write.csv(cor.cate[[i]],paste("output/",prefix,".cor.03.",i,".env_VS_",cate.name[i],".csv",sep=""))}
    }
    names(cor.cate)=cate.name
  }
  
  # output
  output=list(cor.alpha=cor.alpha,cor.beta=cor.beta,cor.cate=cor.cate)
  output
}