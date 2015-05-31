
## input ##
wd="C:/Users/Daliang/Documents/R/GWMC/16S-Merge-0527" # work directory name
prefix="GW150530M-16SUP" # project name
read.limit=27000 # min reads per sample
gene="16S" #ITS or 16S

# check following file names if you used other names.
com.file="otu_table_UP.tabular" # file name of otu table before resample
treat.file="treat.csv" # file name of treatment information
samplist.file="samplist.csv" # file for sample name correction
comr.file="otu_resampled.txt" # file name of otu table after resample
classif.file="classifier_16S_UP" # classification file


## loading files ##
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
classif.all=tran.ITSrdp(file=paste("input/",classif.file,".csv",sep=""))
classif=classif.all$taxa
}else{
classif=read.table(file=paste("input/",classif.file,".txt",sep=""),header=T,sep="\t",row.names=1)
}

## 1.1 ## basic diversity analysis
source("Rcode/egg1.r")
com.egg=egg1(comi=com.a,treat=treat,com.raw=com.b,classif=classif,env=NA,prefix)

## 