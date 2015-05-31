alpha<-function(comm,com.raw=NA)
{
# by Daliang Ning 2015.5.23
library(vegan)
sp.pool=t(estimateR(comm))
if(!is.null(nrow(com.raw))){sp.pool.b=t(estimateR(com.raw));colnames(sp.pool.b)=paste(colnames(sp.pool.b),".before.resamp",sep="");sp.pool=cbind(sp.pool,sp.pool.b)}
alpha.method=c("shannon","simpson","invsimpson")
samp.num=nrow(comm)
samp.name=rownames(comm)
alpha=data.frame(matrix(,nrow=samp.num,ncol=length(alpha.method)+ncol(sp.pool)))
colnames(alpha)=c(alpha.method,colnames(sp.pool))
rownames(alpha)=samp.name
for(i in 1:length(alpha.method))
{
alpha[i]=diversity(comm,index=alpha.method[i])
}
alpha[,(length(alpha.method)+1):ncol(alpha)]=sp.pool
alpha
}