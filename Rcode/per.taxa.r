per.taxa<-function(classif,abun=NA,level=4)
{
# percentage of each taxon, by incidence or abundance

# match sample names
if(is.na(sum(abun)))
{
res=data.frame(Freq=-1,percent=-1)
i=1
rownames(res)=colnames(classif)[i]
tab=as.matrix(table(classif[,i]))
per=data.frame(Freq=tab,percent=tab[,1]/sum(tab[,1]))
res=rbind(res,per)
for(i in 2:level)
{
classif[,i]=paste(classif[,i-1],classif[,i],sep=".")
res[nrow(res)+2,]=c(-1,-1);rownames(res)[nrow(res)]=colnames(classif)[i]
tab=as.matrix(table(classif[,i]))
per=data.frame(Freq=tab,percent=tab[,1]/sum(tab[,1]))
res=rbind(res,per)
}
}else{
res=data.frame(Abun=-1,Percent=-1)
abun=data.frame(abun)
samp.name=rownames(abun)
classif=classif[match(samp.name,rownames(classif)),]
i=1
rownames(res)=colnames(classif)[i]
tax.lev=levels(as.factor(classif[,i]))
per=data.frame(matrix(,nrow=length(tax.lev),ncol=2))
rownames(per)=tax.lev;colnames(per)=c("Abun","Percent")
for(j in 1:length(tax.lev))
	{
	per[j,1]=sum(abun[which(classif[,i]==tax.lev[j]),1])
	}
per[,2]=per[,1]/sum(per[,1])
res=rbind(res,per)
for(i in 2:level)
	{
	classif[,i]=paste(classif[,i-1],classif[,i],sep=".")
	res[nrow(res)+2,]=c(-1,-1);rownames(res)[nrow(res)]=colnames(classif)[i]
	tax.lev=levels(as.factor(classif[,i]))
	per=data.frame(matrix(,nrow=length(tax.lev),ncol=2))
	rownames(per)=tax.lev;colnames(per)=c("Abun","Percent")
	for(j in 1:length(tax.lev))
		{
		per[j,1]=sum(abun[which(classif[,i]==tax.lev[j]),1])
		}
	per[,2]=per[,1]/sum(per[,1])
	res=rbind(res,per)
	}
}
res[is.na(res)]=""
res
}