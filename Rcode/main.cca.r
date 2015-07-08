
# Step 1# input #
wd="C:/Users/Daliang/Dropbox/forCongWangCCA"
code.wd="C:/Users/Daliang/Dropbox/ToolDevelop/github/egg/Rcode"

com.file="OTU.txt"
env.file="env.csv"
top.case.file="case.top.csv"
prefix="cca.test"

com.simple.yn=FALSE # simplify community matrix?
cca.test.yn=TRUE # test CCA using simplified commmunity?
factor.num.up=9
factor.num.down=7
cca.topcases.yn=FALSE # you have got top cases and want to do CCA using full community matrix?
cca.top1.case=FALSE # input the case number(s) if you have decided which case you want to use and want to output data for plot.
cca.one.yn=FALSE # do you want to do CCA for each env factor with full community matrix?

new.limit=1000
resample.method="taxa" # "taxa" for geochip, "sequence" for sequencing results.
nworker=4
resample.times=20
memory.Gb=50
list.limit=5000

####### once you set all parameters above, run all scripts ##########################
#####################################################################################
####### you do not need to change anything below unless you mean it #################
# Step 2# match env and OTU sample names #
setwd(wd)
library(vegan)
com=read.table(file = paste("input/",com.file,sep = ""),header = T,sep = "\t",row.names = 1)
comm=t(com)
dim(comm)
env=read.table(file = paste("input/",env.file,sep = ""),header = T,sep = ",",row.names = 1)
dim(env)
if(sum(rownames(comm)!=rownames(env))>0)
{
  samp.name=intersect(rownames(comm),rownames(env))
  comm=coms[match(samp.name,rownames(comm)),]
  env=env[match(samp.name,rownames(env)),]
  message("some sample names are not the same between community matrix and env. removed")
}

# Step 3 # resample community matrix to get a simplified matrix # CCA test using simplified community matrix
coms.log=FALSE
if(com.simple.yn)
{
  
  source(file = paste(code.wd,"/com.simp.r",sep = ""))
  coms=com.simp(comi = comm,new.limit = new.limit,resamp.meth = resample.method, check.cor = TRUE,times = resample.times,out = "max.r",nworker = nworker,code.wd = code.wd)
  dim(coms)
  coms.log=TRUE
  write.csv(t(coms),file = paste("output/",prefix,".com.simple.csv",sep = ""))
}

if(cca.test.yn)
{
  if(!coms.log)
  {
    if(!file.exists(paste("output/",prefix,".com.simple.csv",sep = "")))
    {
      message("you have not simplified the community matrix")
      break
    }else{
      coms=t(read.table(file = paste("output/",prefix,".com.simple.csv",sep = ""),header = T,sep = ",",row.names = 1))
    }
  }
  env.num=ncol(env)
  env.name=colnames(env)
  samp.num=nrow(comm)
  f.num=min(env.num,samp.num-1)
  case.num=2^env.num-1
  source(file = paste(code.wd,"/intTobit.r",sep = ""))
  cases=lapply(1:case.num, intTobit,len=env.num)
  cases=data.frame(t(matrix(unlist(cases),nrow = length(cases[[1]]))))
  colnames(cases)=env.name
  casum=rowSums(cases)
  cases.select=cases[casum<=f.num & casum<=factor.num.up & casum>=factor.num.down,]
  
  source(file = paste(code.wd,"/cca.test.p.r",sep = ""))
  cca.test=cca.test.p(com.test = coms,env = env,cases = cases.select, summ = FALSE,nworker = nworker,memory.G = memory.Gb,list.limit = list.limit,code.wd = code.wd)
  write.csv(cca.test,file = paste("output/",prefix,".cca.simple.test.",factor.num.down,"to",factor.num.up,".csv",sep = ""))
}

# Step 4 # CCA top cases you picked up #
if(cca.topcases.yn)
{
  if(!file.exists(paste("input/",top.case.file,sep = "")))
  {
    message("you have not generate any file indicating the top cases.")
  }else{
    topcase=read.table(file = paste("input/",top.case.file,sep = ""),header = T,sep = ",",row.names = 1)
    source(file = paste(code.wd,"/cca.test.p.r",sep = ""))
    cca.all=cca.test.p(com.test = comm,env = env,cases = topcase,summ = TRUE,nworker = nworker,memory.G = memory.Gb,list.limit = list.limit,code.wd = code.wd)
    topcases=TRUE
    write.csv(cca.all$index,file = paste("output/",prefix,".cca.topcases.index.csv",sep = ""))
    save(cca.all,file = paste("output/",prefix,".cca.topcases.RData",sep = ""))
  }
}

if(is.numeric(cca.top1.case))
{
  if(!topcases)
  {
    if(!file.exists(paste("output/",prefix,".cca.topcases.RData",sep = "")))
    {
      message("you have not finish CCA for top cases using full community matrix.")
      break
    }else{
      load(file = paste("output/",prefix,".cca.topcases.RData",sep = ""))
    }
  }
  for(i in 1:length(cca.top1.case))
  {
    write.csv(cca.all$detail[[cca.top1.case[i]]]$sp,file = paste("output/",prefix,".cca.topcase",cca.top1.case[i],".sp.csv",sep = ""))
    write.csv(cca.all$detail[[cca.top1.case[i]]]$st,file = paste("output/",prefix,".cca.topcase",cca.top1.case[i],".site.csv",sep = ""))
    write.csv(cca.all$detail[[cca.top1.case[i]]]$en,file = paste("output/",prefix,".cca.topcase",cca.top1.case[i],".env.csv",sep = ""))
    write.csv(cca.all$detail[[cca.top1.case[i]]]$pr,file = paste("output/",prefix,".cca.topcase",cca.top1.case[i],".prop.csv",sep = ""))
  }
}


if(cca.one.yn)
{
  env.num=ncol(env)
  case1=matrix(0,nrow=env.num,ncol=env.num)
  diag(case1)=1
  case1=data.frame(case1)
  rownames(case1)=1:env.num
  colnames(case1)=colnames(env)
  source(file = paste(code.wd,"/cca.test.p.r",sep = ""))
  cca.one=cca.test.p(com.test = comm,env = env,cases = case1,summ = FALSE,nworker = nworker,memory.G = memory.Gb,list.limit = list.limit,code.wd = code.wd)
  write.csv(cca.one,file = paste("output/",prefix,"cca.eachEnv.csv",sep = ""))
}
