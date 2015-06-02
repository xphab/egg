dissim<-function(comm,treat,dist.method=NA)
{
# by Daliang Ning(ningdaliang@gmail.com) on 2014.4.8
library(vegan)
if(is.na(dist.method)[1])
{
dist.method=c("euclidean","manhattan","jaccard","bray","binomial")
}
treat.level=levels(as.factor(as.vector(treat)))
treat.num=length(treat.level)
dis.test=data.frame(matrix(,nrow=(treat.num*(treat.num-1)/2),ncol=(2+(length(dist.method)*3*2))))
colnames(dis.test)[1:2]=c("group1","group2")
colnames(dis.test)[3:(length(dist.method)*2+2)]=as.vector(rbind(paste("mrpp.delta",dist.method,sep="."),paste("mrpp.p",dist.method,sep=".")))
colnames(dis.test)[(length(dist.method)*2+3):(length(dist.method)*4+2)]=as.vector(rbind(paste("anosim.R",dist.method,sep="."),paste("anosim.p",dist.method,sep=".")))
colnames(dis.test)[(length(dist.method)*4+3):(length(dist.method)*6+2)]=as.vector(rbind(paste("adonis.R2",dist.method,sep="."),paste("adonis.p",dist.method,sep=".")))
k=1
for(i in 1:(treat.num-1))
{
for(j in (i+1):treat.num)
{
dis.test[k,1:2]=c(treat.level[i],treat.level[j])
id=(treat==treat.level[i]|treat==treat.level[j])
grouping=as.vector(treat[id])
data=comm[id,]
for(m in 1:length(dist.method))
	{
	testm=mrpp(data,grouping,distance=dist.method[m])
	dis.test[k,((m-1)*2+3)]=testm$delta
	dis.test[k,((m-1)*2+4)]=testm$Pvalue
	dist.res=vegdist(data,method=dist.method[m])
	testas=anosim(data,grouping,distance=dist.method[m])
	dis.test[k,((m-1)*2+(length(dist.method)*2)+3)]=testas$statistic
	dis.test[k,((m-1)*2+(length(dist.method)*2)+4)]=testas$signif
	testad=adonis(data~grouping,method=dist.method[m])
	dis.test[k,((m-1)*2+(length(dist.method)*4)+3)]=testad$aov.tab$R2[1]
	dis.test[k,((m-1)*2+(length(dist.method)*4)+4)]=testad$aov.tab$Pr[1]
	message("now i=",i,"/",(treat.num-1)," j=",j,"/",treat.num," k=",k," m=",m,"/",length(dist.method))
	}
k=k+1
}
}
dis.test
}