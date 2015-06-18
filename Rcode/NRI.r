NRI<-function(comm, dis, memo.size.GB=50, samp.group=NA, weighted=c(TRUE,FALSE),grouping=c(FALSE,TRUE),rand=1000,output.MPD=c(FALSE,TRUE),code.wd)
{
  # calculate NRI based on MPD#
  ## written by Daliang Ning (ningdaliang@gmail.com) ##
  # version beta p1.0: 2015.6.13
  ## cite this version of R script as personal communication from Daliang Ning in University of Oklahoma ##
  ## grouping: If grouping=TRUE, randomization will perform within group. If group="N", randomization will be across all samples. Default is FALSE.##
  
  #load package
  library(picante)
  source(file = paste(code.wd,"/mpdn.r",sep = ""))
  memory.limit(size=memo.size.GB*1024)
  weighted=weighted[1]
  grouping=grouping[1]
  output.MPD=output.MPD[1]
  
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
    colnames(result)="NRI"
    rownames(result)=samp.name
    result.mntd=result
    
    # calculate MPD and NRI within group #
    group=levels(as.factor(samp.group[,2]))
    group.n=length(group)
    MPD=list()
    NRI=list()
    for(m in 1:group.n)
    {
      samp.namex=samp.group[samp.group[,2]==group[m],1]#choose group
      comx=comm[match(samp.namex,rownames(comm)),]#remove others
      comx=comx[,colSums(comx)!=0]#remove undetected species in this group
      spnamex=colnames(comx)
      disx=dis[match(spnamex,rownames(dis)),match(spnamex,colnames(dis))]#remove undetected species
      message("Now calculating observed MPD. Begin at ", date(),". Please wait...")
      MPD.obs<-mpdn(comx, disx, abundance.weighted = weighted)# calculate observed MPD.
      MPD[[m]]=as.matrix(MPD.obs)
      message("Now randomizing. m=",m," in ",group.n,". begin at ", date(),". Please wait...")
      MPD.rand<-t(replicate(rand,mpdn(comx,taxaShuffle(disx),abundance.weighted = weighted)))
      MPD.rand.mean<-apply(MPD.rand,2,mean,na.rm=TRUE)
      MPD.rand.sd<-apply(MPD.rand,2,sd,na.rm=TRUE)
      NRI[[m]]=t((MPD.rand.mean-t(MPD[[m]]))/MPD.rand.sd)
      result[match(rownames(NRI[[m]]),samp.name),]=NRI[[m]]
      result.mntd[match(samp.namex,samp.name),]=MPD[[m]]
    }
  }else{
    # calculate across all samples #
    message("Now calculating observed MPD. Begin at ", date(),". Please wait...")
    MPD.obs<-as.matrix(mpdn(comm, dis, abundance.weighted = weighted)) # calculate observed MPD.
    spname=colnames(comm)
    message("Now randomizing. Begin at ", date(),". Please wait...")
    MPD.rand<-t(replicate(rand,mpdn(comm,taxaShuffle(dis),abundance.weighted = weighted)))
    MPD.rand.mean<-apply(MPD.rand,2,mean,na.rm=TRUE)
    MPD.rand.sd<-apply(MPD.rand,2,sd,na.rm=TRUE)
    result=t((MPD.rand.mean-t(MPD.obs))/MPD.rand.sd)
    result.mntd=MPD.obs
  }
  if(output.MPD)
  {
    output=list(NRI=result,MPD=result.mntd)
    output
  }else{
    result
  }
}