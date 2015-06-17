RC.p<-function(comm,method="bray",rand=1000,portion=c(TRUE,FALSE),nworker=4,memory.G=50)
{
# by Daliang Ning (ningdaliang@gmail.com) on 2015.2.12 #
library(vegan)
library(parallel)
memory.limit(size=memory.G*1024)
com<-comm[,colSums(comm)>0]
BC.obs<-as.matrix(vegdist(com,method=method))
com.rd0=com
com.rd0[]=0
id<-(1:ncol(com))
prob.sp<-colSums(com>0)
prob.ab<-colSums(com)
Si<-rowSums(com>0)
Ni<-rowSums(com)
samp.num=nrow(com)

BC.rand<-function(j,com.rd0,samp.num,id,prob.sp,prob.ab,Si,Ni,method)
{
library(vegan)
com.rd=com.rd0
for(i in 1:samp.num)
{
id.sp<-sample(id,Si[i],replace=FALSE,prob=prob.sp)
if(length(id.sp)==1){count=rep(id.sp,Ni[i])}else{
count<-sample(id.sp,Ni[i],replace=TRUE,prob=prob.ab[id.sp])
								}
table<-table(count)
com.rd[i,as.numeric(names(table))]=as.vector(table)
}
BCrand=as.matrix(vegdist(com.rd,method=method))
BCrand
}

c1<-makeCluster(nworker,type="PSOCK")
message("Now parallel computing. begin at ", date(),". Please wait...")
BC.rd<-parLapply(c1,1:rand,BC.rand,com.rd0=com.rd0,samp.num=samp.num,id=id,prob.sp=prob.sp,prob.ab=prob.ab,Si=Si,Ni=Ni,method=method)
stopCluster(c1)
BC.rd=array(unlist(BC.rd),dim=c(nrow(BC.rd[[1]]),ncol(BC.rd[[1]]),length(BC.rd)))
gc()
comp<-function(x,c){(x<c)+0.5*(x==c)}
message("----now calculating rc at",date(),"----")
alpha=matrix(rowSums(apply(BC.rd,3,comp,c=BC.obs)),nrow=nrow(BC.obs))/rand
rc=(alpha-0.5)*2
rownames(rc)=rownames(BC.obs)
colnames(rc)=colnames(BC.obs)
if(portion)
{
message("----now calculating portion and distribution at",date(),"----")
rc.v=as.vector(as.dist(rc))
rc.higher=sum(rc.v>0.95)/length(rc.v)
rc.lower=sum(rc.v<(-0.95))/length(rc.v)
rc.drift=1-rc.higher-rc.lower
portion.rc=data.frame(rc.lower,rc.drift,rc.higher,num=length(rc.v))
Kernal=density(rc.v,from=-1,to=1)
distrib=data.frame(x.axis=Kernal$x,y.axis=Kernal$y)
res=list(rc=rc,proportion=portion.rc,distribution=distrib)
res
}else{rc}
}
