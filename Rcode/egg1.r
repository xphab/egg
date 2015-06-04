egg1<-function(comi,treat,com.raw,classif,level=5,env,prefix,write.output=TRUE,code.wd,statement.yn=1,alpha.yn=1,DCA.yn=1,Dissim.yn=1,taxa.yn=1)
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
  
  ## Dissimilarity test between teatments
  if(Dissim.yn!=0&!is.null(nrow(treat)))
  {
    message("now doing dissimilarity test between treatments. ",date())
    source(file=paste(code.wd,"/dissim.r",sep=""))
    dis.test=list()
    for(i in 1:ncol(treat))
    {
      dis.test[[i]]=dissim(comm,treat[,i],dist.method=c("euclidean","jaccard","bray"))
      write.csv(dis.test[[i]],paste("output/",prefix,".05.DissimiTest.",colnames(treat)[i],".csv",sep=""))
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
  # output
  output=list(sample.num=samp.num,species.num=sp.num,sequences=reads.t,resample.reads=resamp.read,
              alpha=alpha.div,dca=dca.res,taxa=taxa.samp, dissimilary=dis.test)
  output
}
}