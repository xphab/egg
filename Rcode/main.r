
setwd("C:/Users/Daliang/Dropbox/NSF-Macrosystem/GiganteITS_analyze")

# setwd("C:/Users/Daliang/Dropbox/ToolDevelop/github/egg")
# install.packages("vegan")
# install.packages("picante")
# install.packages("data.table")

#########################################
## input ##
inp=read.table(file="input/1.Input.csv",header=T,sep=",",row.names=1,as.is=TRUE)
wd=inp["wd",1]
code.wd=inp["code.wd",1]
prefix=inp["prefix",1]
read.limit=as.numeric(inp["read.limit",1])
gene=inp["gene",1]
its.conf=inp["its.conf",1]
memory.G=as.numeric(inp["memory.G",1])
com.file=inp["com.file",1]
treat.file=inp["treat.file",1]
samplist.file=inp["samplist.file",1]
comr.file=inp["comr.file",1]
classif.file=inp["classif.file",1]
rm.samplist=inp["rm.samplist",1]
prep.resamp=inp["prep.resamp",1]
statement.yn=inp["statement.yn",1]
alpha.yn=inp["alpha.yn",1]
beta.yn=inp["beta.yn",1]
DCA.yn=inp["DCA.yn",1]
Dissim.yn=inp["Dissim.yn",1]
taxa.yn=inp["taxa.yn",1]
ieg.yn=inp["ieg.yn",1]
cateDCA.list=inp["cateDCA.list",1]
cateSum.list=inp["cateSum.list",1]
cateSum.DCA=inp["cateSum.DCA",1]
cateSum.cateDCA=inp["cateSum.cateDCA",1]
cateCor.file=inp["cateCor.file",1]
alpha.env.yn=inp["alpha.env.yn",1]
beta.env.yn=inp["beta.env.yn",1]
cate.env.yn=inp["cate.env.yn",1]
pd.file=inp["pd.file",1]
tree.file=inp["tree.file",1]
nworker=as.numeric(inp["nworker",1])
phylo.yn=inp["phylo.yn",1]
PD.yn=inp["PD.yn",1]
MPD.yn=inp["MPD.yn",1]
NRI.yn=inp["NRI.yn",1]
MNTD.yn=inp["MNTD.yn",1]
NTI.yn=inp["NTI.yn",1]
bMPD.yn=inp["bMPD.yn",1]
bNRI.yn=inp["bNRI.yn",1]
bMNTD.yn=inp["bMNTD.yn",1]
bNTI.yn=inp["bNTI.yn",1]
ab.weight=as.logical(inp["ab.weight",1])
exclude.consp=as.logical(inp["exclude.consp",1])
rand.times=as.numeric(inp["rand.times",1])
env.file=inp["env.file",1]
cor.MPD.env.yn=inp["cor.MPD.env.yn",1]
cor.NRI.env.yn=inp["cor.NRI.env.yn",1]
cor.MNTD.env.yn=inp["cor.MNTD.env.yn",1]
cor.NTI.env.yn=inp["cor.NTI.env.yn",1]
cor.bMPD.env.yn=inp["cor.bMPD.env.yn",1]
cor.bNRI.env.yn=inp["cor.bNRI.env.yn",1]
cor.bMNTD.env.yn=inp["cor.bMNTD.env.yn",1]
cor.bNTI.env.yn=inp["cor.bNTI.env.yn",1]

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
if(!is.na(samplist.file))
{
  write.csv(data.frame(sample=rownames(com.raw)),file=paste("output/",prefix,".00.oldSampName.csv",sep=""))
  samplist=read.table(file=paste("input/",samplist.file,sep=""),header=T,sep=",",row.names=1)
  sum(is.na(match(samplist[,1],rownames(com.raw))))
}
## 0.1 ## prepare otu table for resample
if(file.exists(paste("output/",prefix,".00.rawOTU.txt",sep="")))
{
  com.b=read.table(file=paste("output/",prefix,".00.rawOTU.txt",sep=""),header=T,sep="\t",row.names = 1)
  com.b=t(com.b)
}else{
  if(prep.resamp!=0)
  {
    source(file=paste(code.wd,"/minread.r",sep=""))
    com.read=minread(com.raw,samplist,prefix,read.limit=read.limit,write.otu=TRUE,code.wd=code.wd)
    com.b=com.read$com.okay
    com.read$min.read
  }else{
    if(!is.na(samplist.file))
    {
      source(file=paste(code.wd,"/rename.samp.r",sep=""))
      com.b=rename.samp(com.raw,samplist)
      write.table(cbind(OTU=colnames(com.b),t(com.b)),file=paste("output/",prefix,".00.rawOTU.txt",sep=""),sep="\t",row.names=FALSE,quote=FALSE)
    }else{
      com.b=com.raw
    }
  }
}
dim(com.b)

## loading file ##
### resampled otu table
if(file.exists(paste("input/",comr.file,sep="")))
{
  com.a=read.table(file=paste("input/",comr.file,sep=""),header=T,sep="\t",row.names=1);com.a=t(com.a)
}else{
  if(prep.resamp==0)
  {
    com.a=com.b
    warning("using the raw community table instead of resampled table.")
  }else{
    com.a=NA
    warning("you are preparing or need to prepare resampled table.")
  }
}
dim(com.a)

### remove sample list
if(!is.na(rm.samplist))
{
  source(file=paste(code.wd,"/rm.samp.r",sep=""))
  com.a=rm.samp(com.a,rm.samplist)
  com.b=rm.samp(com.b,rm.samplist)
}

### treatment
if(file.exists(paste("input/",treat.file,sep="")))
{
  treat=read.table(file=paste("input/",treat.file,sep=""),header=T,sep=",",row.names=1)
}else{
  treat=NA
}
dim(treat)

### env
if(file.exists(paste("input/",env.file,sep="")))
{
  env=read.table(file=paste("input/",env.file,sep=""),header=T,sep=",",row.names=1)
}else{
  env=NA
}

### classification information
if(file.exists(paste("input/",classif.file,sep="")))
{
  if(gene=="ITS")
  {
    source(file=paste(code.wd,"/tran.ITSrdp.r",sep=""))
    classif.all=tran.ITSrdp(file=paste("input/",classif.file,sep=""),conf=its.conf)
    classif=classif.all$taxa
  }else{
    classif=read.table(file=paste("input/",classif.file,sep=""),header=T,sep="\t",row.names=1)
  }
}else{
  classif=NA
}
dim(classif)

### category DCA list
if(file.exists(paste("input/",cateDCA.list,sep="")))
{
  cateDCA.g=read.table(file=paste("input/",cateDCA.list,sep=""),header=T,sep=",",row.names=1)
}else{
  cateDCA.g=NA
}

#### category sum list
if(file.exists(paste("input/",cateSum.list,sep="")))
{
  cateSum.g=read.table(file=paste("input/",cateSum.list,sep=""),header=T,sep=",",row.names=1)
}else{
  cateSum.g=NA
}

## 1 ## basic diversity analysis, alpha, DCA, taxa overall composition
source(file=paste(code.wd,"/egg1.r",sep=""))
com.egg=egg1(comi=com.a,treat=treat,com.raw=com.b,classif=classif,level=5,env=NA,prefix,
             write.output=TRUE,code.wd=code.wd,statement.yn=statement.yn,alpha.yn=alpha.yn,beta.yn=beta.yn,
             DCA.yn=DCA.yn,Dissim.yn=Dissim.yn,taxa.yn=taxa.yn,cateDCA.g=cateDCA.g,cateSum.g=cateSum.g,
             cateSum.DCA=cateSum.DCA,cateSum.cateDCA=cateSum.cateDCA)

## 2 ## correlation test
source(file = paste(code.wd,"/egg.cor.r",sep=""))
if(file.exists(paste("input/",cateCor.file,sep="")))
{
  cateCor.g=read.table(file=paste("input/",cateCor.file,sep=""),header=T,sep=",",row.names=1)
}else{
  cateCor.g=NA
}
cor.egg=egg.cor(comi = com.a,category = cateCor.g,env = env,
                alpha.env.yn = alpha.env.yn,beta.env.yn = beta.env.yn,cate.env.yn = cate.env.yn,
                nworker=nworker,write.output = TRUE,prefix = prefix,code.wd = code.wd)

## 3 ## generate community file and treatment file for ieg pipeline
if(ieg.yn!=0)
{
  source(file=paste(code.wd,"/ieg.upload.r",sep=""))
  ieg.up=ieg.upload(com.a,treat,prefix,category=NA)
}

# 3 # Phylogenetic analysis
if(phylo.yn!=0)
{
  library(picante)
  # loading files
  if(file.exists(paste("input/",pd.file,sep="")))
  {
    pd=read.table(file=paste("input/",pd.file,sep=""),header=T,sep=",",row.names=1)
    if(sum(colnames(pd)!=rownames(pd))>0)
    {
      colnames(pd)=rownames(pd)
      message("force the colnames of pd to be the same as rownames")
    }
  }else{
    pd=NA
  }
  
  if(file.exists(paste("input/",tree.file,sep="")))
  {
    tree=read.tree(file=paste("input/",tree.file,sep=""))
  }else{
    tree=NA
  }
  
  source(file=paste(code.wd,"/egg.p.r",sep=""))
  phylo.egg=egg.p(comi=com.a,tree=tree,pd=pd,env=env,nworker=nworker,memory.G=memory.G,prefix=prefix,
                  PD.yn=PD.yn,MPD.yn=MPD.yn,NRI.yn=NRI.yn,MNTD.yn=MNTD.yn,NTI.yn=NTI.yn,
                  bMPD.yn=bMPD.yn,bNRI.yn=bNRI.yn,bMNTD.yn=bMNTD.yn,bNTI.yn=bNTI.yn,
                  cor.MPD.env.yn=cor.MPD.env.yn, cor.NRI.env.yn=cor.NRI.env.yn,
                  cor.MNTD.env.yn=cor.MNTD.env.yn, cor.NTI.env.yn=cor.NTI.env.yn,
                  cor.bMPD.env.yn=cor.bMPD.env.yn, cor.bNRI.env.yn=cor.bNRI.env.yn,
                  cor.bMNTD.env.yn=cor.bMNTD.env.yn, cor.bNTI.env.yn=cor.bNTI.env.yn,
                  ab.weight=ab.weight,exclude.consp=exclude.consp,rand.times=rand.times,code.wd=code.wd)
  
}

# END # save work space
save.image(file=paste("output/",prefix,".",format(Sys.time(),"%Y%b%d"),".RData",sep = ""))

# please feel free to contact Daliang Ning (ningdaliang@gmail.com)
# If you use it, you may cite this version as
# Daliang Ning. 2015. Egg. Retrived Jun 26, 2015, from https://github.com/DaliangNing/egg
##### end ####