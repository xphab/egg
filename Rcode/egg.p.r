egg.p<-function(comi,tree=tree,pd=pd,nworker=8,memory.G=50,prefix,PD.yn=0,MPD.yn=0,NRI.yn=0,MNTD.yn=0,NTI.yn=0,bMPD.yn=0,bNRI.yn=0,bMNTD.yn=0,bNTI.yn=0,ab.weight=TRUE,exclude.consp=FALSE,rand.times=1000,code.wd)
{
  comm=comi[,colSums(comi)>0]
  library(picante)
  # calculate phylogenetic distance matrix if not available yet
  if(is.null(nrow(pd)))
  {
    if(is.na(tree)[[1]])
    {
      message("neither pd nor tree is available. stop.")
      break
    }else{
      message("Now calculating phylogenetic distances. ",date())
      source(file = paste(code.wd,"/pdist.p.r",sep = ""))
      pd=pdist.p(tree,nworker=nworker,memory.G=memory.G)
      gc()
      write.csv(pd,file=paste("output/",prefix,".p01.phylodist.csv",sep=""))
    }
  }
  
  # match species names
  sp.name=intersect(colnames(comm),rownames(pd))
  comm=comm[,match(sp.name,colnames(comm))]
  pd=pd[match(sp.name,rownames(pd)),match(sp.name,rownames(pd))]
  gc()
  lost=data.frame(Discarded=c("com",setdiff(colnames(comm),rownames(pd)),"","pd",setdiff(rownames(pd),colnames(comm))))
  write.csv(lost,file=paste("output/",prefix,".p02.lostSpecies.csv",sep = ""))
  samp.name=rownames(comm)
  samp.num=nrow(comm)
  
  # phylogenetic diversity
  
  gc()
  # MPD and NRI
  if(MPD.yn!=0)
  {
    message("Now calculating MPD. ",date())
    source(file = paste(code.wd,"/mpdn.r",sep = ""))
    mpdist=mpdn(comm = comm,pd = pd,abundance.weighted = ab.weight)
    mpdist=data.frame(mpdist)
    write.csv(mpdist,file = paste("output/",prefix,".p03.MPD.csv",sep = ""))
    gc()
  }else{
    mpdist=NA
  }
  if(NRI.yn!=0)
  {
    message("Now calculating NRI. ",date())
    source(file = paste(code.wd,"/NRI.p.r",sep = ""))
    NRI.res=NRI.p(comm,pd,nworker = nworker,memo.size.GB = memory.G,samp.group = NA,weighted = ab.weight,grouping = FALSE,rand = rand.times,output.MPD = FALSE,code.wd=code.wd)
    write.csv(NRI.res,file=paste("output/",prefix,".p04.NRI.csv",sep = ""))
    gc()
  }else{
    NRI.res=NA
  }
  
  # MNTD and NTI
  if(MNTD.yn!=0)
  {
    message("Now calculating MNTD. ",date())
    source(file = paste(code.wd,"/mntdn.r",sep = ""))
    mntdist=mntdn(comm,pd,abundance.weighted = ab.weight)
    mntdist=data.frame(mntdist)
    write.csv(mntdist,file = paste("output/",prefix,".p05.MNTD.csv",sep = ""))
    gc()
  }else{
    mntdist=NA
  }
  if(NTI.yn!=0)
  {
    message("Now calculating NTI. ",date())
    source(file = paste(code.wd,"/NTI.p.r",sep = ""))
    NTI.res=NTI.p(comm=comm, dis=pd, nworker=nworker, memo.size.GB=memory.G, samp.group=NA, weighted=ab.weight,grouping=FALSE,rand=rand.times,output.MNTD=FALSE,code.wd=code.wd)
    #sesmntd=ses.mntd(samp = comm,dis = pd,null.model = "taxa.labels",abundance.weighted = ab.weight,runs = rand.times)
    #sesmntd[,9]=-sesmntd$mntd.obs.z
    #colnames(sesmntd)[9]="NTI"
    write.csv(NTI.res,file=paste("output/",prefix,".p06.NTI.csv",sep = ""))
    gc()
  }else{
    sesmntd=NA
  }
  
  # betaMPD and betaNRI
  if(bMPD.yn!=0)
  {
    message("Now calculating betaMPD. ",date())
    source(file = paste(code.wd,"/bmpd.r",sep = ""))
    bMPD=as.matrix(bmpd(comm,pd,abundance.weighted = ab.weight))
    write.csv(bMPD,file=paste("output/",prefix,".p07.betaMPD.csv",sep = ""))
    gc()
  }else{
    bMPD=NA
  }
  if(bNRI.yn!=0)
  {
    message("Now calculating betaNRI. ",date())
    source(file = paste(code.wd,"/bNRI.p.r",sep = ""))
    bNRI=bNRI.p(comm,pd,nworker = nworker,memo.size.GB = memory.G,samp.group = NA,weighted = ab.weight,grouping = FALSE,rand = rand.times,output.bMPD = FALSE,code.wd=code.wd)
    write.csv(bNRI,file=paste("output/",prefix,".p08.betaNRI.csv",sep = ""))
    gc()
  }else{
    bNRI=NA
  }
  
  # betaMNTD and betaNTI
  if(bMNTD.yn!=0)
  {
    message("Now calculating betaMNTD. ",date())
    source(file = paste(code.wd,"/bmntd.r",sep = ""))
    bMNTD=as.matrix(bmntd(comm, pd, abundance.weighted = ab.weight, exclude.conspecifics = exclude.consp))
    write.csv(bMNTD,file=paste("output/",prefix,".p09.betaMNTD.csv",sep = ""))
    gc()
  }else{
    bMNTD=NA
  }
  if(bNTI.yn!=0)
  {
    message("Now calculating betaNTI. ",date())
    source(file = paste(code.wd,"/bNTI3.p.r",sep = ""))
    bNTI=bNTI3.p(comm, pd, nworker=nworker, memo.size.GB=memory.G, samp.group=NA, weighted=ab.weight,grouping=FALSE,exclude.consp=exclude.consp,rand=rand.times,output.bMNTD=FALSE,code.wd=code.wd)
    write.csv(bNTI,file=paste("output/",prefix,".p10.betaNTI.csv",sep = ""))
    gc()
  }else{
    bNTI=NA
  }
  
  # altogether
  alt.alpha=data.frame(MPD=mpdist,NRI=NRI.res,MNTD=mntdist,NTI=NTI.res)
  m=matrix(nrow=samp.num,ncol=samp.num)
  alt.beta=data.frame(sample1=samp.name[as.vector(as.dist(row(m)))],sample2=samp.name[as.vector(as.dist(col(m)))],
                      betaMPD=as.vector(as.dist(bMPD)),betaNRI=as.vector(as.dist(bNRI)),
                      betaMNTD=as.vector(as.dist(bMNTD)),betaNTI=as.vector(as.dist(bNTI)))
  write.csv(alt.alpha,file = paste("output/",prefix,".p11.phylo_alpha_all.csv",sep = ""))
  write.csv(alt.beta,file = paste("output/",prefix,".p12.phylo_beta_all.csv",sep = ""))
  # output
  output=list(dist=pd,MPD=mpdist,NRI=NRI.res,MNTD=mntdist,NTI=NTI.res,
              betaMPD=bMPD,betaNRI=bNRI,betaMNTD=bMNTD,betaNTI=bNTI)
}