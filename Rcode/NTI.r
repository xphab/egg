NTI<-function(comm, dis, memo.size.GB=50, samp.group=NA, weighted=c(TRUE,FALSE),grouping=c(FALSE,TRUE),rand=1000,output.MNTD=c(FALSE,TRUE),code.wd)
{
  # calculate NTI based on MNTD, need package "picante" by parallel compute#
  ## written by Daliang Ning (ningdaliang@gmail.com) ##
  # version beta p1.0: 2015.6.13
  ## cite this version of R script as personal communication from Daliang Ning in University of Oklahoma ##
  ## grouping: If grouping=TRUE, randomization will perform within group. If group="N", randomization will be across all samples. Default is FALSE.##
  
  #load package
  library(picante)
  source(file = paste(code.wd,"/mntdn.r",sep = ""))
  memory.limit(size=memo.size.GB*1024)
  weighted=weighted[1]
  grouping=grouping[1]
  output.MNTD=output.MNTD[1]
  
  # match
  if(sum(colnames(comm)!=rownames(dis))>0)
  {
    sp.name=intersect(colnames(comm),rownames(dis))
    comm=comm[,match(sp.name,colnames(comm))]
    dis=dis[match(sp.name,rownames(dis)),match(sp.name,rownames(dis))]
  }
  
  if(grouping)
  {
    samp.name=rownames(comm)
    result=data.frame(matrix(NA,nrow=length(samp.name),ncol = 1))
    colnames(result)="NTI"
    rownames(result)=samp.name
    result.mntd=result
    
    # calculate MNTD and NTI within group #
    group=levels(as.factor(samp.group[,2]))
    group.n=length(group)
    MNTD=list()
    NTI=list()
    for(m in 1:group.n)
    {
      samp.namex=samp.group[samp.group[,2]==group[m],1]#choose group
      comx=comm[match(samp.namex,rownames(comm)),]#remove others
      comx=comx[,colSums(comx)!=0]#remove undetected species in this group
      spnamex=colnames(comx)
      disx=dis[match(spnamex,rownames(dis)),match(spnamex,colnames(dis))]#remove undetected species
      message("Now calculating observed MNTD. Begin at ", date(),". Please wait...")
      MNTD.obs<-mntdn(comx, disx, abundance.weighted = weighted)# calculate observed MNTD.
      MNTD[[m]]=as.matrix(MNTD.obs)
      message("Now randomizing. m=",m," in ",group.n,". begin at ", date(),". Please wait...")
      MNTD.rand<-t(replicate(rand,mntdn(comx,taxaShuffle(disx),abundance.weighted = weighted)))
      MNTD.rand.mean<-apply(MNTD.rand,2,mean,na.rm=TRUE)
      MNTD.rand.sd<-apply(MNTD.rand,2,sd,na.rm=TRUE)
      NTI[[m]]=t((MNTD.rand.mean-t(MNTD[[m]]))/MNTD.rand.sd)
      result[match(rownames(NTI[[m]]),samp.name),]=NTI[[m]]
      result.mntd[match(samp.namex,samp.name),]=MNTD[[m]]
    }
  }else{
    # calculate across all samples #
    message("Now calculating observed MNTD. Begin at ", date(),". Please wait...")
    MNTD.obs<-as.matrix(mntdn(comm, dis, abundance.weighted = weighted)) # calculate observed MNTD.
    spname=colnames(comm)
    message("Now randomizing. Begin at ", date(),". Please wait...")
    MNTD.rand<-t(replicate(rand,mntdn(comm,taxaShuffle(dis),abundance.weighted = weighted)))
    MNTD.rand.mean<-apply(MNTD.rand,2,mean,na.rm=TRUE)
    MNTD.rand.sd<-apply(MNTD.rand,2,sd,na.rm=TRUE)
    result=t((MNTD.rand.mean-t(MNTD.obs))/MNTD.rand.sd)
    result.mntd=MNTD.obs
  }
  if(output.MNTD)
  {
    output=list(NTI=result,MNTD=result.mntd)
    output
  }else{
    result
  }
}