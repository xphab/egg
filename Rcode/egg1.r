egg1<-function(comi,treat,com.raw,classif,level=5,env,prefix,write.output=TRUE,code.wd,statement.yn=1,alpha.yn=1,beta.yn=1,DCA.yn=1,Dissim.yn=1,taxa.yn=1,cateDCA.g=NA,cateSum.g=NA,cateSum.DCA=0,cateSum.cateDCA=0)
{
# Package "egg" (environment genomic gadgets) initiated by Daliang Ning (ningdaliang@gmail.com)
# includes all commonly used methods for microbial community data
# comi: a row is a sample, a column is a species. so rownames are sample names, and colnames are species names.

library(vegan)

# remove zero-read species
if(is.null(nrow(comi)))
{
  warning("The final community table has not been input. Analysis is stopped.")
}else{
  comm=comi[,colSums(comi)>0]
  comm.raw=com.raw[,colSums(com.raw)>0]
  # taxa number
  message("now calculating total number of taxa, sequences and samples. ",date())
  sp.name=colnames(comm)
  sp.num=ncol(comm)
  samp.name=rownames(comm)
  reads.t=sum(comm)
  samp.num=nrow(comm)
  resamp.read=rowSums(comm)
  if(statement.yn!=0)
  {
    sp.raw.num=ncol(comm.raw)
    reads.raw.t=sum(comm.raw)
    samp.raw.num=nrow(comm.raw)
    if(min(resamp.read)!=max(resamp.read)){message("error! reads after resampling are not equal");break}
    state.1=paste("A total of ",sp.num," taxa were detected from ",reads.t," sequences in ",samp.num," samples after resampled to ",resamp.read[1]," sequences per sample. ",
                  "A total of ", sp.raw.num, " taxa were detected from ", reads.raw.t, " sequences in ", samp.raw.num, " samples before resampled.", sep="")
    if(write.output){write(state.1,paste("output/",prefix,".01.statement.txt",sep=""))}
  }
  
  # match
  if(!is.null(nrow(treat)))
  {
    treat.name=colnames(treat)
    treat=data.frame(treat[match(samp.name,rownames(treat)),])
    colnames(treat)=treat.name
  }
  if(!is.null(nrow(classif)))
  {
    classif=classif[match(sp.name,rownames(classif)),]
  }
  if(!is.null(nrow(env))){env=env[match(samp.name,rownames(env)),]}
  if(!is.null(nrow(com.raw))){com.raw=com.raw[match(samp.name,rownames(com.raw)),]}
  if(!is.null(nrow(cateDCA.g)))
  {
    cateDCA.g.name=colnames(cateDCA.g);rname=rownames(cateDCA.g)
    cateDCA.g=data.frame(cateDCA.g[match(sp.name,rownames(cateDCA.g)),])
    colnames(cateDCA.g)=cateDCA.g.name;rownames(cateDCA.g)=rname
    
  }
  
  if(!is.null(nrow(cateSum.g)))
  {
    cateSum.g.name=colnames(cateSum.g);rname=rownames(cateSum.g)
    cateSum.g=data.frame(cateSum.g[match(sp.name,rownames(cateSum.g)),])
    colnames(cateSum.g)=cateSum.g.name;rownames(cateSum.g)=rname
  }
  
  
  # alpha diversity
  if(alpha.yn!=0)
  {
    message("now calculating alpha diversity indexes. ",date())
    source(file=paste(code.wd,"/alpha.r",sep=""))
    alpha.div=alpha(comm,com.raw)
    if(write.output){write.csv(alpha.div,paste("output/",prefix,".02.alpha.csv",sep=""))}
  }else{
    alpha.div=NA
  }

  # beta diversity
  ## distance
  if(beta.yn!=0)
  {
    message("now calculating beta diversity indexes. ",date())
    source(file = paste(code.wd,"/beta.dis.r",sep=""))
    beta.div=beta.dis(comm,binary = "Both")
    if(write.output){write.csv(beta.div,paste("output/",prefix,".05.1.beta.csv",sep=""))}
  }else{
    beta.div=NA
  }
  
  ## DCA
  if(DCA.yn!=0)
  {
    message("now calculating DCA. ",date())
    dca.res=decorana(comm)
    dca.sum=summary(dca.res)
    if(write.output){write.csv(dca.sum$site.scores,paste("output/",prefix,".03.DCAsite.csv",sep=""))}
    if(write.output){write.csv(dca.sum$spec.scores,paste("output/",prefix,".03.DCAspecies.csv",sep=""))}
  }else{
    dca.res=NA
  }
  #### DCA in each category
  if(is.null(nrow(cateDCA.g)))
  {
    cateDCA.res=NA
  }else{
    message("now calculating DCA for each cate. ",date())
    source(file=paste(code.wd,"/cateDCA.r",sep=""))
    cateDCA.res=cateDCA(comm,cateDCA.g,prefix.cDCA=paste("output/",prefix,".06.cateDCA.",sep = ""),write.output=write.output)
  }

    
  ## Dissimilarity test between teatments
  if(Dissim.yn!=0&!is.null(nrow(treat)))
  {
    message("now doing dissimilarity test between treatments. ",date())
    source(file=paste(code.wd,"/dissim.r",sep=""))
    dis.test=list()
    for(i in 1:ncol(treat))
    {
      dis.test[[i]]=dissim(comm,treat[,i],dist.method=c("euclidean","jaccard","bray"))
      write.csv(dis.test[[i]],paste("output/",prefix,".05.2.DissimiTest.",colnames(treat)[i],".csv",sep=""))
    }
    names(dis.test)=colnames(treat)
  }else{
    dis.test=NA
  }
  
  # correlation
  ## env vs diverity linear model
  
  
  # taxonomic composition
  ## percentage
  if(taxa.yn!=0&!is.null(nrow(classif)))
  {
    message("now calculating percentage of each taxon. ",date())
    source(file=paste(code.wd,"/samp.taxa.r",sep=""))
    taxa.samp=samp.taxa(comm,classif,level=level,code.wd)
    if(write.output)
    {
      write.csv(taxa.samp$abun,paste("output/",prefix,".04.TaxaComp.Abun.csv",sep=""))
      write.csv(taxa.samp$percent,paste("output/",prefix,".04.TaxaComp.percent.csv",sep=""))
    }
  }else{
    taxa.samp=NA
  }
  
  # category sum
  ## community table of category
  if(!is.null(nrow(cateSum.g)))
  {
    message("now sum each category to get new community table. ", date())
    com.cateS=list()
    pbar <- txtProgressBar(min = 0, max = 20, style = 3)
    for(i in 1:ncol(cateSum.g))
    {
      cate.lev=levels(as.factor(as.vector(cateSum.g[,i])))
      com.cateS[[i]]=data.frame(matrix(NA,nrow = nrow(comm),ncol = length(cate.lev)))
      for(j in 1:length(cate.lev))
      {
        com.cateS[[i]][,j]=rowSums(data.frame(comm[,cateSum.g[,i]==cate.lev[j]]))
        setTxtProgressBar(pbar, round((20*(((i-1)/ncol(cateSum.g))+(j/length(cate.lev)/ncol(cateSum.g)))),1))
      }
      colnames(com.cateS[[i]])=cate.lev;rownames(com.cateS[[i]])=rownames(comm)
      if(write.output)
      {
        write.csv(t(com.cateS[[i]]),file = paste("output/",prefix,".07.cateSum.",i,".",colnames(cateSum.g)[i],".csv",sep = ""))
      }
    }
    close(pbar)
    names(com.cateS)=colnames(cateSum.g)
  }else{
    com.cateS=NA
  }
  
  ### do DCA for new community table
  if(cateSum.DCA!=0)
  {
    dca.cS.res=list()
    for(m in 1:length(com.cateS))
    {
      message("now calculating DCA for composition of ",names(com.cateS)[m],". ",date())
      dca.cS.res[[m]]=decorana(comm)
      dca.cS.sum=summary(dca.cS.res[[m]])
      if(write.output){write.csv(dca.cS.sum$site.scores,paste("output/",prefix,".08.DCAcate.",m,".",names(com.cateS)[m],".All.Site.csv",sep=""))}
      if(write.output){write.csv(dca.cS.sum$spec.scores,paste("output/",prefix,".08.DCAcate.",m,".",names(com.cateS)[m],".All.species.csv",sep=""))}
    }
    names(dca.cS.res)=names(com.cateS)
  }else{
    dca.cS.res=NA
  }
  
  #### do DCA for new community table, category by category
  if(cateSum.cateDCA!=0)
  {
    source(file=paste(code.wd,"/cateMatch.r",sep=""))
    cate.list=cateMatch(from.list=cateDCA.g,to.list=cateSum.g)
    source(file=paste(code.wd,"/cateDCA.r",sep=""))
    cs.cDCA.res=list()
    for(i in 1:length(cate.list))
    {
      cs.cDCA.res[[i]]=cateDCA(com.cateS[[i]],cate.list[[i]],prefix.cDCA=paste("output/",prefix,".09.DCAcate.",i,".",names(com.cateS)[i],".",sep = ""),write.output=write.output)
    }
    names(cs.cDCA.res)=names(cate.list)
  }else{
    cs.cDCA.res=NA
  }
 
  # output
  output=list(sample.num=samp.num,species.num=sp.num,sequences=reads.t,resample.reads=resamp.read,
              alpha=alpha.div,beta=beta.div,dca=dca.res,taxa=taxa.samp, dissimilary=dis.test,categoryDCA=cateDCA.res,
              categorySum=com.cateS,cateSum.DCA=dca.cS.res,cateSum.catDCA=cs.cDCA.res)
  output
}
}