NTI.p<-function(comm, dis, nworker=4, memo.size.GB=50, samp.group=NA, weighted=c(TRUE,FALSE),grouping=c(FALSE,TRUE),rand=1000,output.MNTD=c(FALSE,TRUE),code.wd)
{
  # calculate NTI based on MNTD#
  ## written by Daliang Ning (ningdaliang@gmail.com) ##
  # version beta p1.0: 2015.6.13
  ## cite this version of R script as personal communication from Daliang Ning in University of Oklahoma ##
  ## grouping: If grouping=TRUE, randomization will perform within group. If group="N", randomization will be across all samples. Default is FALSE.##
  
  #load package
  #library(picante)
  library(parallel)
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
  samp.name=rownames(comm)
  result=data.frame(matrix(NA,nrow=length(samp.name),ncol = 1))
  colnames(result)="NTI"
  rownames(result)=samp.name
  result.mntd=result
  
  ## randomization function ##
  MNTD.random<-function(i,diss,com,weighted,code.wd)
  {
    source(file=paste(code.wd,"/mntdn.r",sep=""))
    diss.rand=diss
    rand.name=sample(colnames(diss))
    colnames(diss.rand)=rand.name
    rownames(diss.rand)=rand.name
    gc()
    MNTD.rand<-as.matrix(mntdn(com, diss.rand, abundance.weighted = weighted))
    MNTD.rand
  }
  
  if(grouping)
  {
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
      gc()
      message("Now calculating observed MNTD. Begin at ", date(),". Please wait...")
      MNTD.obs<-mntdn(comx, disx, abundance.weighted = weighted)# calculate observed MNTD.
      MNTD[[m]]=as.matrix(MNTD.obs)
      c1<-makeCluster(nworker,type="PSOCK")
      message("Now randomizing by parallel computing. m=",m," in ",group.n,". begin at ", date(),". Please wait...")
      MNTD.rand<-parLapply(c1,1:rand,MNTD.random,diss=disx,com=comx,weighted=weighted,code.wd=code.wd)
      stopCluster(c1)
      gc()
      MNTD.rand<-array(unlist(MNTD.rand),dim=c(nrow(MNTD.rand[[1]]),ncol(MNTD.rand[[1]]),length(MNTD.rand)))
      NTI[[m]]=(apply(MNTD.rand,c(1,2),mean)-MNTD[[m]])/(apply(MNTD.rand,c(1,2),sd))
      
      result[match(rownames(NTI[[m]]),samp.name),]=NTI[[m]]
      result.mntd[match(samp.namex,samp.name),]=MNTD[[m]]
      gc()
    }
  }else{
    # calculate across all samples #
    message("Now calculating observed MNTD. Begin at ", date(),". Please wait...")
    gc()
    MNTD.obs<-as.matrix(mntdn(comm, dis, abundance.weighted = weighted)) # calculate observed MNTD.
    spname=colnames(comm)
    gc()
    c1<-makeCluster(nworker,type="PSOCK")
    message("Now randomizing by parallel computing. Begin at ", date(),". Please wait...")
    MNTD.rand<-parLapply(c1,1:rand,MNTD.random,diss=dis,com=comm,weighted=weighted,code.wd=code.wd)
    stopCluster(c1)
    gc()
    MNTD.rand<-array(unlist(MNTD.rand),dim=c(nrow(MNTD.rand[[1]]),ncol(MNTD.rand[[1]]),length(MNTD.rand)))
    NTI=(apply(MNTD.rand,c(1,2),mean)-MNTD.obs)/(apply(MNTD.rand,c(1,2),sd))
    result=NTI
    result.mntd=MNTD.obs
    gc()
  }
  if(output.MNTD)
  {
    output=list(NTI=result,MNTD=result.mntd)
  }else{
    output=data.frame(NTI=result)}
  output
}