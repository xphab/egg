
## input ##
# install.packages("vegan")
wd="C:/Users/Daliang/Dropbox/ToolDevelop/github/egg" # work directory name
prefix="test" # project name
read.limit=100 # min reads per sample
gene="ITS" #ITS or 16S
its.conf=0.8 # ITS classification confidence threshold
memory.G=50 # memory limitation

# check following file names if you used other names.
com.file="otu_table.txt" # file name of otu table before resample
treat.file="treat.csv" # file name of treatment information
samplist.file="samplist.csv" # file for sample name correction
comr.file="otu_resampled.txt" # file name of otu table after resample
classif.file="classifier" # classification file

## loading files ##
library(vegan)
memory.limit(size=memory.G*1024)

setwd(wd)
com.raw=read.table(file=paste("input/",com.file,sep=""),header=T,sep="\t",row.names=1);com.raw=t(com.raw)
if(grepl("tabular",com.file))
{
rownames(com.raw)=sub("abel.","",rownames(com.raw))
rownames(com.raw)=sub("_",".",rownames(com.raw))
}
dim(com.raw)
write.csv(data.frame(sample=rownames(com.raw)),file=paste("output/",prefix,".00.oldSampName.csv",sep=""))
samplist=read.table(file=paste("input/",samplist.file,sep=""),header=T,sep=",",row.names=1)
## 0.1 ## prepare otu table for resample
source("Rcode/minread.r");source("Rcode/rename.samp.r")
com.read=minread(com.raw,samplist,prefix,read.limit=read.limit,write.otu=TRUE)
com.b=com.read$com.okay
com.read$min.read
dim(com.b)

## loading file ##
com.a=read.table(file=paste("input/",comr.file,sep=""),header=T,sep="\t",row.names=1);com.a=t(com.a)
dim(com.a)
treat=read.table(file=paste("input/",treat.file,sep=""),header=T,sep=",",row.names=1)
if(gene=="ITS")
{
source("Rcode/tran.ITSrdp.r")
classif.all=tran.ITSrdp(file=paste("input/",classif.file,".csv",sep=""),conf=its.conf)
classif=classif.all$taxa
}else{
classif=read.table(file=paste("input/",classif.file,".txt",sep=""),header=T,sep="\t",row.names=1)
}

## 1.1 ## basic diversity analysis, alpha, DCA, taxa overall composition
source("Rcode/egg1.r")
com.egg=egg1(comi=com.a,treat=treat,com.raw=com.b,classif=classif,level=5,env=NA,prefix)

 