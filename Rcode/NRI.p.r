NRI.p<-function(comm, dis, nworker=4, memo.size.GB=50, samp.group=NA, weighted=c(TRUE,FALSE),grouping=c(FALSE,TRUE),rand=1000,output.MPD=c(FALSE,TRUE),code.wd)
{
  # calculate NRI based on MPD#
  ## written by Daliang Ning (ningdaliang@gmail.com) ##
  # version beta p1.0: 2015.6.13
  ## cite this version of R script as personal communication from Daliang Ning in University of Oklahoma ##
  ## grouping: If grouping=TRUE, randomization will perform within group. If group="N", randomization will be across all samples. Default is FALSE.##
  
  #load package
  #library(picante)
  library(parallel)
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
  samp.name=rownames(comm)
  result=data.frame(matrix(NA,nrow=length(samp.name),ncol = 1))
  colnames(result)="NRI"
  rownames(result)=samp.name
  result.mpd=result
  
  ## randomization function ##
  MPD.random<-function(i,diss,com,weighted,code.wd)
  {
    source(file=paste(code.wd,"/mpdn.r",sep=""))
    diss.rand=diss
    rand.name=sample(colnames(diss))
    colnames(diss.rand)=rand.name
    rownames(diss.rand)=rand.name
    gc()
    MPD.rand<-as.matrix(mpdn(com, diss.rand, abundance.weighted = weighted))
    MPD.rand
  }
  
  if(grouping)
  {
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
      gc()
      message("Now calculating observed MPD. Begin at ", date(),". Please wait...")
      MPD.obs<-mpdn(comx, disx, abundance.weighted = weighted)# calculate observed MPD.
      MPD[[m]]=as.matrix(MPD.obs)
      c1<-makeCluster(nworker,type="PSOCK")
      message("Now randomizing by parallel computing. m=",m," in ",group.n,". begin at ", date(),". Please wait...")
      MPD.rand<-parLapply(c1,1:rand,MPD.random,diss=disx,com=comx,weighted=weighted,code.wd=code.wd)
      stopCluster(c1)
      gc()
      MPD.rand<-array(unlist(MPD.rand),dim=c(nrow(MPD.rand[[1]]),ncol(MPD.rand[[1]]),length(MPD.rand)))
      NRI[[m]]=(apply(MPD.rand,c(1,2),mean)-MPD[[m]])/(apply(MPD.rand,c(1,2),sd))
      
      result[match(rownames(NRI[[m]]),samp.name),]=NRI[[m]]
      result.mpd[match(samp.namex,samp.name),]=MPD[[m]]
      gc()
    }
  }else{
    # calculate across all samples #
    message("Now calculating observed MPD. Begin at ", date(),". Please wait...")
    gc()
    MPD.obs<-as.matrix(mpdn(comm, dis, abundance.weighted = weighted)) # calculate observed MPD.
    spname=colnames(comm)
    gc()
    c1<-makeCluster(nworker,type="PSOCK")
    message("Now randomizing by parallel computing. Begin at ", date(),". Please wait...")
    MPD.rand<-parLapply(c1,1:rand,MPD.random,diss=dis,com=comm,weighted=weighted,code.wd=code.wd)
    stopCluster(c1)
    gc()
    MPD.rand<-array(unlist(MPD.rand),dim=c(nrow(MPD.rand[[1]]),ncol(MPD.rand[[1]]),length(MPD.rand)))
    NRI=(apply(MPD.rand,c(1,2),mean)-MPD.obs)/(apply(MPD.rand,c(1,2),sd))
    result=NRI
    result.mpd=MPD.obs
    gc()
  }
  if(output.MPD)
  {
    output=list(NRI=result,MPD=result.mntd)
  }else{
    output=data.frame(NRI=result)}
  output
}