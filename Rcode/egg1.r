egg1<-function(comi,treat,com.raw,classif,level=5,env,prefix,write.output=TRUE,code.wd)
{
# Package "egg" (environment genomic gadgets) initiated by Daliang Ning (ningdaliang@gmail.com)
# includes all commonly used methods for microbial community data
# comi: a row is a sample, a column is a species. so rownames are sample names, and colnames are species names.

library(vegan)

# remove zero-read species
comm=comi[,colSums(comi)>0]
comm.raw=com.raw[,colSums(com.raw)>0]
# taxa number
message("now calculating total number of taxa, sequences and samples. ",date())
sp.name=colnames(comm)
sp.num=ncol(comm);sp.raw.num=ncol(comm.raw)
reads.t=sum(comm);reads.raw.t=sum(comm.raw)
samp.name=rownames(comm)
samp.num=nrow(comm);samp.raw.num=nrow(comm.raw)
resamp.read=rowSums(comm)
if(min(resamp.read)!=max(resamp.read)){message("error! reads after resampling are not equal");break}
state.1=paste("A total of ",sp.num," taxa were detected from ",reads.t," sequences in ",samp.num," samples after resampled to ",resamp.read[1]," sequences per sample. ",
"A total of ", sp.raw.num, " taxa were detected from ", reads.raw.t, " sequences in ", samp.raw.num, " samples before resampled.", sep="")
if(write.output){write(state.1,paste("output/",prefix,".01.statement.txt",sep=""))}
# match
treat=treat[match(samp.name,rownames(treat)),]
classif=classif[match(sp.name,rownames(classif)),]
if(!is.null(nrow(env))){env=env[match(samp.name,rownames(env)),]}
if(!is.null(nrow(com.raw))){com.raw=com.raw[match(samp.name,rownames(com.raw)),]}

# alpha diversity
message("now calculating alpha diversity indexes. ",date())
source(file=paste(code.wd,"/alpha.r",sep=""))
alpha.div=alpha(comm,com.raw)
if(write.output){write.csv(alpha.div,paste("output/",prefix,".02.alpha.csv",sep=""))}

# beta diversity
## DCA
message("now calculating DCA. ",date())
dca.res=decorana(comm)
dca.sum=summary(dca.res)
if(write.output){write.csv(dca.sum$site.scores,paste("output/",prefix,".03.DCAsite.csv",sep=""))}
if(write.output){write.csv(dca.sum$spec.scores,paste("output/",prefix,".03.DCAspecies.csv",sep=""))}

# correlation
## env vs diverity linear model


# taxonomic composition
## percentage
message("now calculating percentage of each taxon. ",date())
source(file=paste(code.wd,"/samp.taxa.r",sep=""))
taxa.samp=samp.taxa(comm,classif,level=level,code.wd)
if(write.output)
{
  write.csv(taxa.samp$abun,paste("output/",prefix,".04.TaxaComp.Abun.csv",sep=""))
  write.csv(taxa.samp$percent,paste("output/",prefix,".04.TaxaComp.percent.csv",sep=""))
}

# output
output=list(sample.num=samp.num,species.num=sp.num,sequences=reads.t,resample.reads=resamp.read,
alpha=alpha.div,dca=dca.res,taxa=taxa.samp)
output
}