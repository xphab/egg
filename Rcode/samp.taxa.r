samp.taxa<-function(comi,classif,level=4,code.wd)
{
  source(file=paste(code.wd,"/per.taxa.r",sep=""))
  sp.name=colnames(comi)
  classif=classif[match(sp.name,rownames(classif)),]
  abun.all=colSums(comi)
  taxa.p=per.taxa(classif,abun=abun.all,level=level)
  st.ab=data.frame(matrix(,nrow=nrow(taxa.p),ncol=(nrow(comi)+1)))
  rownames(st.ab)=rownames(taxa.p)
  colnames(st.ab)=c("All",rownames(comi))
  st.p=st.ab
  st.ab[,1]=taxa.p[,1]
  st.p[,1]=taxa.p[,2]
  for(i in 1:nrow(comi))
  {
    taxa.p=per.taxa(classif,abun=comi[i,],level=level)
    st.ab[,(i+1)]=taxa.p[,1]
    st.p[,(i+1)]=taxa.p[,2]
  }
  res=list(abun=st.ab,percent=st.p)
  res
}