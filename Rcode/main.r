
setwd("C:/Users/Daliang/Dropbox/GWMC-analysis/150523-ITS-merge/UPARSE")
# install.packages("vegan")

#########################################
## input ##
inp=read.table(file="input/1.Input.csv",header=T,sep=",",row.names=1,as.is=TRUE)
wd=inp[1,1]# work directory name
code.wd=inp[2,1] # the directory you save the Rcode
prefix=inp[3,1] # project name
read.limit=as.numeric(inp[4,1]) # min reads per sample
gene=inp[5,1] #ITS or 16S
its.conf=as.numeric(inp[6,1]) # ITS classification confidence threshold
memory.G=as.numeric(inp[7,1])# memory limitation
com.file=inp[8,1] # file name of otu table before resample
treat.file=inp[9,1] # file name of treatment information
samplist.file=inp[10,1] # file for sample name correction
comr.file=inp[11,1] # file name of otu table after resample
classif.file=inp[12,1] # classification file
rm.samplist=inp[13,1] # samples that need be removed. name as "remove.samp.csv"

## loading files ##
library(vegan)
memory.limit(size=memory.G*1024)
com.raw=read.table(file=paste("input/",com.file,sep=""),header=T,sep="\t",row.names=1);com.raw=t(com.raw)
if(grepl("tabular",com.file))
{
rownames(com.raw)=sub("abel.","",rownames(com.raw))
rownames(com.raw)=sub("_",".",rownames(com.raw))
}
dim(com.raw)
write.csv(data.frame(sample=rownames(com.raw)),file=paste("output/",prefix,".00.oldSampName.csv",sep=""))
samplist=read.table(file=paste("input/",samplist.file,sep=""),header=T,sep=",",row.names=1)
sum(is.na(match(samplist[,1],rownames(com.raw))))
## 0.1 ## prepare otu table for resample
source(file=paste(code.wd,"/minread.r",sep=""))
source(file=paste(code.wd,"/rename.samp.r",sep=""))
com.read=minread(com.raw,samplist,prefix,read.limit=read.limit,write.otu=TRUE,code.wd=code.wd)
com.b=com.read$com.okay
com.read$min.read
dim(com.b)

## loading file ##
com.a=read.table(file=paste("input/",comr.file,sep=""),header=T,sep="\t",row.names=1);com.a=t(com.a)
dim(com.a)
if(!is.na(rm.samplist))
{
  source(file=paste(code.wd,"/rm.samp.r",sep=""))
  com.a=rm.samp(com.a,rm.samplist)
  com.b=rm.samp(com.b,rm.samplist)
}

treat=read.table(file=paste("input/",treat.file,sep=""),header=T,sep=",",row.names=1)
if(gene=="ITS")
{
source(file=paste(code.wd,"/tran.ITSrdp.r",sep=""))
classif.all=tran.ITSrdp(file=paste("input/",classif.file,sep=""),conf=its.conf)
classif=classif.all$taxa
}else{
classif=read.table(file=paste("input/",classif.file,sep=""),header=T,sep="\t",row.names=1)
}

## 1.1 ## basic diversity analysis, alpha, DCA, taxa overall composition
source(file=paste(code.wd,"/egg1.r",sep=""))
com.egg=egg1(comi=com.a,treat=treat,com.raw=com.b,classif=classif,level=5,env=NA,prefix,write.output=TRUE,code.wd=code.wd)

# 2 # generate community file and treatment file for ieg pipeline
source(file=paste(code.wd,"/ieg.upload.r",sep=""))
ieg.up=ieg.upload(com.a,treat,prefix,category=NA)

# please feel free to contact Daliang Ning (ningdaliang@gmail.com)
# If you use it, you may cite this version as
# Daliang Ning. 2015. Egg. Retrived Jun 1, 2015, from https://github.com/DaliangNing/egg
##### end ####