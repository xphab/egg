wd="C:/Users/Daliang/Dropbox/ToolDevelop/github/cut/input"
files=c(
  "test1.csv",
  "test2.csv",
  "test3.csv"
)
output.file="test.combine.csv"
########################
setwd(wd)
comb=list()
comb[[1]]=read.table(file = files[1],header = TRUE,sep = ",",row.names = 1)
col.name=colnames(comb[[1]])

for(i in 2:length(files))
{
  filei=read.table(file = files[i],header = TRUE,sep = ",",row.names = 1)
  if(sum(colnames(filei)!=col.name)>0)
  {
    message(files[i]," has different header from the first file. Some columns may be removed.")
    filei=filei[,match(col.name,colnames(filei))]
  }
  comb[[i]]=filei
}
library(data.table)
comb=data.frame(rbindlist(comb))
rownames(comb)=1:nrow(comb)
write.csv(comb,file = output.file)
# V1.0 Daliang Ning (ningdaliang@gmail.com) on 2015.7.8