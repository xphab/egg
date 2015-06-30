egg.p<-function(comi,tree=tree,pd=pd,env=env,nworker=8,memory.G=50,prefix,PD.yn=0,MPD.yn=0,NRI.yn=0,MNTD.yn=0,NTI.yn=0,bMPD.yn=0,bNRI.yn=0,bMNTD.yn=0,bNTI.yn=0,cor.MPD.env.yn=0, cor.NRI.env.yn=0, cor.MNTD.env.yn=0, cor.NTI.env.yn=0, cor.bMPD.env.yn=0, cor.bNRI.env.yn=0, cor.bMNTD.env.yn=0, cor.bNTI.env.yn=0, ab.weight=TRUE,exclude.consp=FALSE,rand.times=1000,code.wd)
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
  
  samp.name=intersect(rownames(comm),rownames(env))
  comm=comm[match(samp.name,rownames(comm)),]
  env=env[match(samp.name,rownames(env)),]
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
  
  ## MPD and NRI # Correlation with env factors
  if(cor.MPD.env.yn!=0)
  {
    if((sum(!is.na(mpdist))==0)&(!file.exists(paste("output/",prefix,".p03.MPD.csv",sep = ""))))
    {
      message("MPD is not available for correlation test.")
      mpd.cor=NA
    }else{
      source(file = paste(code.wd,"/cor.v.r",sep = ""))
      message("now correlation test between env and MPD. ",date())
      if(sum(!is.na(mpdist))==0)
      {
        mpdist=read.table(file = paste("output/",prefix,".p03.MPD.csv",sep = ""),header = T,sep = ",",row.names = 1)
      }
      mpd.cor=cor.v(vect = mpdist,env = env)
      write.csv(mpd.cor,file=paste("output/",prefix,".p03.2.MPDvsEnv.csv",sep = ""))
    }
  }else{
    mpd.cor=NA
  }
  
  if(cor.NRI.env.yn!=0)
  {
    if((sum(!is.na(NRI.res))==0)&(!file.exists(paste("output/",prefix,".p04.NRI.csv",sep = ""))))
    {
      message("NRI is not available for correlation test.")
      nri.cor=NA
    }else{
      source(file = paste(code.wd,"/cor.v.r",sep = ""))
      message("now correlation test between env and NRI. ",date())
      if(sum(!is.na(NRI.res))==0)
      {
        NRI.res=read.table(file = paste("output/",prefix,".p04.NRI.csv",sep = ""),header = T,sep = ",",row.names = 1)
      }
      nri.cor=cor.v(vect = NRI.res,env = env)
      write.csv(nri.cor,file=paste("output/",prefix,".p04.2.NRIvsEnv.csv",sep = ""))
    }
  }else{
    nri.cor=NA
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
    NTI.res=NA
  }
  
  ## MNTD and NTI # correlation with env factor
  if(cor.MNTD.env.yn!=0)
  {
    if((sum(!is.na(mntdist))==0)&(!file.exists(paste("output/",prefix,".p05.MNTD.csv",sep = ""))))
    {
      message("MNTD is not available for correlation test.")
      mntd.cor=NA
    }else{
      source(file = paste(code.wd,"/cor.v.r",sep = ""))
      message("now correlation test between env and MNTD. ",date())
      if(sum(!is.na(mntdist))==0)
      {
        mntdist=read.table(file = paste("output/",prefix,".p05.MNTD.csv",sep = ""),header = T,sep = ",",row.names = 1)
      }
      mntd.cor=cor.v(vect = mntdist,env = env)
      write.csv(mntd.cor,file=paste("output/",prefix,".p05.2.MNTDvsEnv.csv",sep = ""))
    }
  }else{
    mntd.cor=NA
  }
  
  if(cor.NTI.env.yn!=0)
  {
    if((sum(!is.na(NTI.res))==0)&(!file.exists(paste("output/",prefix,".p06.NTI.csv",sep = ""))))
    {
      message("NTI is not available for correlation test.")
      nti.cor=NA
    }else{
      source(file = paste(code.wd,"/cor.v.r",sep = ""))
      message("now correlation test between env and NTI. ",date())
      if(sum(!is.na(NTI.res))==0)
      {
        NTI.res=read.table(file = paste("output/",prefix,".p06.NTI.csv",sep = ""),header = T,sep = ",",row.names = 1)
      }
      nti.cor=cor.v(vect = NTI.res,env = env)
      write.csv(nti.cor,file=paste("output/",prefix,".p06.2.NTIvsEnv.csv",sep = ""))
    }
  }else{
    nti.cor=NA
  }
  gc()
  
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
  
  ## betaMPD and betaNRI # correlation with env factor
  #### bMPD vs ENV
  if(cor.bMPD.env.yn!=0)
  {
    if((sum(!is.na(bMPD))==0)&(!file.exists(paste("output/",prefix,".p07.betaMPD.csv",sep = ""))))
    {
      message("bMPD is not available for correlation test.")
      bmpd.cor=NA
    }else{
      source(file = paste(code.wd,"/cor.m.r",sep = ""))
      message("now correlation test between env and betaMPD. ",date())
      if(sum(!is.na(bMPD))==0)
      {
        bMPD=read.table(file = paste("output/",prefix,".p07.betaMPD.csv",sep = ""),header = T,sep = ",",row.names = 1)
      }
      bmpd.cor=cor.m(dis = bMPD,env = env)
      write.csv(bmpd.cor,file=paste("output/",prefix,".p07.2.bMPDvsEnv.csv",sep = ""))
    }
  }else{
    bmpd.cor=NA
  }
  
  ### bNRI vs ENV
  if(cor.bNRI.env.yn!=0)
  {
    if((sum(!is.na(bNRI))==0)&(!file.exists(paste("output/",prefix,".p08.betaNRI.csv",sep = ""))))
    {
      message("bNRI is not available for correlation test.")
      bNRI.cor=NA
    }else{
      source(file = paste(code.wd,"/cor.m.r",sep = ""))
      message("now correlation test between env and betaNRI. ",date())
      if(sum(!is.na(bNRI))==0)
      {
        bNRI=read.table(file = paste("output/",prefix,".p08.betaNRI.csv",sep = ""),header = T,sep = ",",row.names = 1)
      }
      bNRI.cor=cor.m(dis = bNRI,env = env)
      write.csv(bNRI.cor,file=paste("output/",prefix,".p08.2.bNRIvsEnv.csv",sep = ""))
    }
  }else{
    bNRI.cor=NA
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
  
  ## betaMNTD and betaNTI # correlation with env factor
  #### bMNTD vs ENV 
  if(cor.bMNTD.env.yn!=0)
  {
    if((sum(!is.na(bMNTD))==0)&(!file.exists(paste("output/",prefix,".p09.betaMNTD.csv",sep = ""))))
    {
      message("bMNTD is not available for correlation test.")
      bMNTD.cor=NA
    }else{
      source(file = paste(code.wd,"/cor.m.r",sep = ""))
      message("now correlation test between env and betaMNTD. ",date())
      if(sum(!is.na(bMNTD))==0)
      {
        bMNTD=read.table(file = paste("output/",prefix,".p09.betaMNTD.csv",sep = ""),header = T,sep = ",",row.names = 1)
      }
      bMNTD.cor=cor.m(dis = bMNTD,env = env)
      write.csv(bMNTD.cor,file=paste("output/",prefix,".p09.2.bMNTDvsEnv.csv",sep = ""))
    }
  }else{
    bMNTD.cor=NA
  }
  
  #### bNTI vs ENV
  if(cor.bNTI.env.yn!=0)
  {
    if((sum(!is.na(bNTI))==0)&(!file.exists(paste("output/",prefix,".p10.betaNTI.csv",sep = ""))))
    {
      message("bNTI is not available for correlation test.")
      bNTI.cor=NA
    }else{
      source(file = paste(code.wd,"/cor.m.r",sep = ""))
      message("now correlation test between env and betaNTI. ",date())
      if(sum(!is.na(bNTI))==0)
      {
        bNTI=read.table(file = paste("output/",prefix,".p10.betaNTI.csv",sep = ""),header = T,sep = ",",row.names = 1)
      }
      bNTI.cor=cor.m(dis = bNTI,env = env)
      write.csv(bNTI.cor,file=paste("output/",prefix,".p10.2.bNTIvsEnv.csv",sep = ""))
    }
  }else{
    bNTI.cor=NA
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