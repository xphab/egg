bNRI.p<-function(comm, dis, nworker=4, memo.size.GB=50, samp.group=NA, weighted=c(TRUE,FALSE),grouping=c(FALSE,TRUE),rand=1000,output.bMPD=c(FALSE,TRUE),code.wd)
{
# calculate betaNRI based on betaMPD#
## written by Daliang Ning (ningdaliang@gmail.com) ##
# version beta p1.0: 2015.6.13
## cite this version of R script as personal communication from Daliang Ning in University of Oklahoma ##
## grouping: If grouping=TRUE, randomization will perform within group. If group="N", randomization will be across all samples. Default is FALSE.##

#load package
library(parallel)
source(file = paste(code.wd,"/bmpd.r",sep = ""))
memory.limit(size=memo.size.GB*1024)
weighted=weighted[1]
grouping=grouping[1]
output.bMPD=output.bMPD[1]

# match
if(sum(colnames(comm)!=rownames(dis))>0)
{
  sp.name=intersect(colnames(comm),rownames(dis))
  comm=comm[,match(sp.name,colnames(comm))]
  dis=dis[match(sp.name,rownames(dis)),match(sp.name,rownames(dis))]
}
samp.name=rownames(comm)
result=data.frame(matrix(NA,length(samp.name),length(samp.name)))
colnames(result)=samp.name
rownames(result)=samp.name
result.mpd=result

## randomization function ##
bMPD.random<-function(i,diss,com,weighted,code.wd)
{
  source(file=paste(code.wd,"/bmpd.r",sep=""))
  diss.rand=diss
  rand.name=sample(colnames(diss))
  colnames(diss.rand)=rand.name
  rownames(diss.rand)=rand.name
  gc()
  bMPD.rand<-as.matrix(bmpd(com, diss.rand, abundance.weighted = weighted))
  bMPD.rand
}

if(grouping)
{
# calculate betaMPD and betaNRI within group #
group=levels(as.factor(samp.group[,2]))
group.n=length(group)
bMPD=list()
bNRI=list()
	for(m in 1:group.n)
	{
	samp.namex=samp.group[samp.group[,2]==group[m],1]#choose group
	if(length(samp.namex)<2){bNRI[[m]]=NA;bMPD[[m]]=NA}else{
	comx=comm[match(samp.namex,rownames(comm)),]#remove others
	comx=comx[,colSums(comx)!=0]#remove undetected species in this group
	spnamex=colnames(comx)
	disx=dis[match(spnamex,rownames(dis)),match(spnamex,colnames(dis))]#remove undetected species
	gc()
	message("Now calculating observed betaMPD. Begin at ", date(),". Please wait...")
	bMPD.obs<-bmpd(comx, disx, abundance.weighted = weighted)# calculate observed betaMPD.
	bMPD[[m]]=as.matrix(bMPD.obs)
	c1<-makeCluster(nworker,type="PSOCK")
	message("Now randomizing by parallel computing. m=",m," in ",group.n,". begin at ", date(),". Please wait...")
	bMPD.rand<-parLapply(c1,1:rand,bMPD.random,diss=disx,com=comx,weighted=weighted,code.wd=code.wd)
	stopCluster(c1)
	gc()
	bMPD.rand<-array(unlist(bMPD.rand),dim=c(nrow(bMPD.rand[[1]]),ncol(bMPD.rand[[1]]),length(bMPD.rand)))
	bNRI[[m]]=(bMPD[[m]]-apply(bMPD.rand,c(1,2),mean))/(apply(bMPD.rand,c(1,2),sd))
									}
	result[match(rownames(bNRI[[m]]),samp.name),match(colnames(bNRI[[m]]),samp.name)]=bNRI[[m]]
	result.mpd[match(samp.namex,samp.name),match(samp.namex,samp.name)]=bMPD[[m]]
	gc()
	}
}else{
# calculate across all samples #
message("Now calculating observed betaMPD. Begin at ", date(),". Please wait...")
gc()
bMPD.obs<-as.matrix(bmpd(comm, dis, abundance.weighted = weighted)) # calculate observed betaMPD.
spname=colnames(comm)
gc()
c1<-makeCluster(nworker,type="PSOCK")
message("Now randomizing by parallel computing. Begin at ", date(),". Please wait...")
bMPD.rand<-parLapply(c1,1:rand,bMPD.random,diss=dis,com=comm,weighted=weighted,code.wd=code.wd)
stopCluster(c1)
gc()
bMPD.rand<-array(unlist(bMPD.rand),dim=c(nrow(bMPD.rand[[1]]),ncol(bMPD.rand[[1]]),length(bMPD.rand)))
bNRI=(bMPD.obs-apply(bMPD.rand,c(1,2),mean))/(apply(bMPD.rand,c(1,2),sd))
result=bNRI
result.mpd=bMPD.obs
gc()
}
if(output.bMPD)
{
output=list(betaNRI=result,betaMPD=result.mpd)
}else{
output=result}
output
}
