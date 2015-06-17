pdist.p<-function(tree,nworker=4,memory.G=50)
{
# by Daliang Ning (ningdaliang@gmail.com) on 2015.2.23
# because cophenetic function in package "ape" often returns error when dealing with large tree. so i wrote this function.
memory.limit(size=memory.G*1024)
library(picante)
library(parallel)
edge<-tree$edge
edge.len<-tree$edge.length
tip.num=length(tree$tip)

pathf<-function(i,edge,edge.len)
{
pathf<-list()
j=1
pathf[[1]]=edge[edge[,2]==i,1]
pathf[[2]]=edge.len[edge[,2]==i]
while(length(edge[edge[,2]==pathf[[1]][j],])!=0){
pathf[[1]][j+1]=edge[edge[,2]==pathf[[1]][j],1]
pathf[[2]][j+1]=edge.len[edge[,2]==pathf[[1]][j]]
j=j+1
}
pathf
}

c1<-makeCluster(nworker,type="PSOCK")
message("Now computing path. begin at ", date(),". Please wait...")
path<-parLapply(c1,1:tip.num,pathf,edge=edge,edge.len=edge.len)
stopCluster(c1)
gc()

pdf<-function(m, tip.num, path)
{
pdf=rep(0,tip.num)
for(n in (m+1):tip.num)
{
a=path[[m]][[1]]
b=path[[n]][[1]]
al=path[[m]][[2]]
bl=path[[n]][[2]]
pdf[n]=sum(al[c(TRUE,is.na(match(a,b)[1:(length(a)-1)]))])+sum(bl[c(TRUE,is.na(match(b,a)[1:(length(b)-1)]))])
}
pdf
}

c2<-makeCluster(nworker,type="PSOCK")
message("Now computing pdist. begin at ", date(),". Please wait...")
pdist.l<-parLapply(c2,1:(tip.num-1),pdf,tip.num=tip.num,path=path)
stopCluster(c2)
gc()
pdist=matrix(unlist(pdist.l),nrow=tip.num)
pdist=data.frame(pdist,rep(0,tip.num))
pdist=pdist+t(pdist)
rownames(pdist)=tree$tip.label
colnames(pdist)=tree$tip.label
pdist
}