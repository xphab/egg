cateDCA<-function(comi,cate.list,prefix.cDCA,write.output=TRUE)
{
  cateDCA=list()
  # match
  sp.name=colnames(comi)
  cate.list.name=colnames(cate.list)
  cate.list=data.frame(cate.list[match(sp.name,rownames(cate.list)),])
  colnames(cate.list)=cate.list.name
  cate.list[is.na(cate.list)]="MissingCategory"  
  for (i in 1:ncol(cate.list))
  {
    cateDCA[[i]]=list()
    cate.lev=levels(as.factor(as.vector(cate.list[,i])))
    for (j in 1:length(cate.lev))
    {
      message("category DCA: i=",i,"/",ncol(cate.list)," j=",j,"/",length(cate.lev),". ",date())
      com.cate=comi[,cate.list[,i]==cate.lev[j]]
      cateDCA[[i]][[j]]=decorana(com.cate)
      dca.sum.cate=summary(cateDCA[[i]][[j]])
      if(write.output){write.csv(dca.sum.cate$site.scores,paste(prefix.cDCA,i,".",j,".",colnames(cateDCA.g)[i],".",cate.lev[j],".site.csv",sep=""))}
      if(write.output){write.csv(dca.sum.cate$spec.scores,paste(prefix.cDCA,i,".",j,".",colnames(cateDCA.g)[i],".",cate.lev[j],".species.csv",sep=""))}
    }
    names(cateDCA[[i]])=cate.lev
  }
  names(cateDCA)=colnames(cate.list)
  cateDCA
}