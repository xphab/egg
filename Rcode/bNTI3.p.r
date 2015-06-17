bNTI3.p<-function(comm, dis, nworker=4, memo.size.GB=50, samp.group=NA, weighted=c(TRUE,FALSE),grouping=c(FALSE,TRUE),exclude.consp=FALSE,rand=1000,output.bMNTD=c(FALSE,TRUE),code.wd)
{
# calculate betaNTI based on betaMNTD, need package "picante" by parallel compute#
## written by Daliang Ning (ningdaliang@gmail.com) ##
# version beta p3.0: 2015.4.1# add setting for exclude.consp.
## Improvement from Stegen's paper: you can calculate within group.
## cite the original reference: Stegen JC, Lin X, Fredrickson JK, Chen X, Kennedy DW, Murray CJ et al. (2013). Quantifying community assembly processes and identifying features that impose them. Isme Journal 7: 2069-2079.##
## cite this version of R script as personal communication from Daliang Ning in University of Oklahoma ##
## grouping: If grouping=TRUE, randomization will perform within group. If group="N", randomization will be across all samples. Default is FALSE.##

#load package
library(picante)
library(parallel)
source(file =paste(code.wd,"/bmntd.r",sep = ""))
memory.limit(size=memo.size.GB*1024)

samp.name=rownames(comm)
result=data.frame(matrix(NA,length(samp.name),length(samp.name)))
colnames(result)=samp.name
rownames(result)=samp.name
result.mntd=result

## randomization function ##
bMNTD.random<-function(i,diss,com,weighted,exclude.consp,code.wd)
{
  source(file=paste(code.wd,"/bmntd.r",sep=""))
  diss.rand=diss
  rand.name=sample(colnames(diss))
  colnames(diss.rand)=rand.name
  rownames(diss.rand)=rand.name
  gc()
  bMNTD.rand<-as.matrix(bmntd(com, diss.rand, abundance.weighted = weighted, exclude.conspecifics = exclude.consp))
  bMNTD.rand
}

if(grouping)
{
# calculate betaMNTD and betaNTI within group #
group=levels(as.factor(samp.group[,2]))
group.n=length(group)
bMNTD=list()
bNTI=list()
	for(m in 1:group.n)
	{
	samp.namex=samp.group[samp.group[,2]==group[m],1]#choose group
	if(length(samp.namex)<2){bNTI[[m]]=NA;bMNTD[[m]]=NA}else{
	comx=comm[match(samp.namex,rownames(comm)),]#remove others
	comx=comx[,colSums(comx)!=0]#remove undetected species in this group
	spnamex=colnames(comx)
	disx=dis[match(spnamex,rownames(dis)),match(spnamex,colnames(dis))]#remove undetected species
	gc()
	message("Now calculating observed betaMNTD. Begin at ", date(),". Please wait...")
	bMNTD.obs<-bmntd(comx, disx, abundance.weighted = weighted, exclude.conspecifics = exclude.consp)# calculate observed betaMNTD.
	bMNTD[[m]]=as.matrix(bMNTD.obs)
	c1<-makeCluster(nworker,type="PSOCK")
	message("Now randomizing by parallel computing. m=",m," in ",group.n,". begin at ", date(),". Please wait...")
	bMNTD.rand<-parLapply(c1,1:rand,bMNTD.random,diss=disx,com=comx,weighted=weighted,exclude.consp=exclude.consp,code.wd=code.wd)
	stopCluster(c1)
	gc()
	bMNTD.rand<-array(unlist(bMNTD.rand),dim=c(nrow(bMNTD.rand[[1]]),ncol(bMNTD.rand[[1]]),length(bMNTD.rand)))
	bNTI[[m]]=(bMNTD[[m]]-apply(bMNTD.rand,c(1,2),mean))/(apply(bMNTD.rand,c(1,2),sd))
									}
	result[match(rownames(bNTI[[m]]),samp.name),match(colnames(bNTI[[m]]),samp.name)]=bNTI[[m]]
	result.mntd[match(samp.namex,samp.name),match(samp.namex,samp.name)]=bMNTD[[m]]
	gc()
	}
}else{
# calculate across all samples #
message("Now calculating observed betaMNTD. Begin at ", date(),". Please wait...")
gc()
bMNTD.obs<-as.matrix(bmntd(comm, dis, abundance.weighted = weighted, exclude.conspecifics = exclude.consp)) # calculate observed betaMNTD.
spname=colnames(comm)
gc()
c1<-makeCluster(nworker,type="PSOCK")
message("Now randomizing by parallel computing. Begin at ", date(),". Please wait...")
bMNTD.rand<-parLapply(c1,1:rand,bMNTD.random,diss=dis,com=comm,weighted=weighted,exclude.consp=exclude.consp,code.wd=code.wd)
stopCluster(c1)
gc()
bMNTD.rand<-array(unlist(bMNTD.rand),dim=c(nrow(bMNTD.rand[[1]]),ncol(bMNTD.rand[[1]]),length(bMNTD.rand)))
bNTI=(bMNTD.obs-apply(bMNTD.rand,c(1,2),mean))/(apply(bMNTD.rand,c(1,2),sd))
result=bNTI
result.mntd=bMNTD.obs
gc()
}
if(output.bMNTD)
{
output=list(betaNTI=result,betaMNTD=result.mntd)
}else{
output=result}
output
}
