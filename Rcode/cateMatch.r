cateMatch<-function(from.list,to.list)
{
  # to match two category lists, generate new category list using to.list's levels as rownames.
  sp.name=levels(as.factor(c(as.vector(rownames(from.list)),as.vector(rownames(to.list)))))
  c1name=colnames(from.list);c2name=colnames(to.list);r1name=rownames(from.list);r2name=rownames(to.list)
  list1=data.frame(from.list[match(sp.name,rownames(from.list)),])
  list1[is.na(list1)]="MissingCategory"
  list2=data.frame(to.list[match(sp.name,rownames(to.list)),])
  list2[is.na(list2)]="MissingCategory"
  colnames(list1)=c1name;colnames(list2)=c2name;rownames(list1)=r1name;rownames(list2)=r2name
  res=list()
  for(i in 1:ncol(list2))
  {
    lev=levels(as.factor(as.vector(list2[,i])))
    res[[i]]=data.frame(matrix(NA,nrow=length(lev),ncol = ncol(list1)))
    for(j in 1:length(lev))
    {
      id=(list2[,i]==lev[j])
      for(m in 1:ncol(list1))
      {
        tb=table(as.vector(list1[id,m]))
        if(length(tb)>1){warning(colnames(list2)[i]," ",lev[j]," VS ",colnames(list1)[m]," has duplicates.")}
        res[[i]][j,m]=names(tb)[tb==max(tb)][1]
      }
    }
    rownames(res[[i]])=lev
    colnames(res[[i]])=colnames(list1)
  }
  names(res)=colnames(list2)
  res
}