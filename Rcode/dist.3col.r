dist.3col<-function(dist)
{
# by Daliang Ning (ningdaliang@gmail.com) on 2015.5.17
# transfrom distance matrix to a 3-column table
# dist muast be a matrix
dist=as.matrix(dist)
rowname=rownames(dist)
colname=colnames(dist)
rown=row(dist)
coln=col(dist)
dist.v=as.vector(as.dist(dist))
rown.v=as.vector(as.dist(rown))
coln.v=as.vector(as.dist(coln))
res=data.frame(name1=rowname[rown.v],name2=colname[coln.v],dis=dist.v)
res
}