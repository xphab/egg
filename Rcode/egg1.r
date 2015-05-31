egg1<-function(comi,treat,com.raw,classif,env,prefix,write.output=TRUE)
{
# Package "egg" (environment genomic gadgets) initiated by Daliang Ning (ningdaliang@gmail.com)
# includes all commonly used methods for microbial community data
# comi: a row is a sample, a column is a species. so rownames are sample names, and colnames are species names.

library(vegan)

# remove zero-read species
comm=comi[,colSums(comi)>0]

# taxa number
message("now calculating total number of taxa, sequences and samples. ",date())
sp.name=colnames(comm)
sp.num=ncol(comm)
reads.t=sum(comm)
samp.name=rownames(comm)
samp.num=nrow(comm)
resamp.read=rowSums(comm)
if(min(resamp.read)!=max(resamp.read)){message("error! reads after resampling are not equal");break}
state.1=paste("A total of ",sp.num," taxa were detected from ",reads.t," sequences in ",samp.num," samples after resampled to ",resamp.read[1]," sequences per sample.",sep="")
if(write.output){write(state.1,paste("output/",prefix,".01.statement.txt",sep=""))}
# match
treat=treat[match(samp.name,rownames(treat)),]
if(!is.null(nrow(env))){env=env[match(samp.name,rownames(env)),]}
if(!is.null(nrow(com.raw))){com.raw=com.raw[match(samp.name,rownames(com.raw)),]}

# alpha diversity
message("now calculating alpha diversity indexes. ",date())
source("Rcode/alpha.r")
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
source("Rcode/per.taxa.r")
abun=colSums(comm)
taxa.per=per.taxa(classif,abun=abun,level=2)
if(write.output){write.csv(taxa.per,paste("output/",prefix,".04.TaxaCom.csv",sep=""))}

# output
output=list(sample.num=samp.num,species.num=sp.num,sequences=reads.t,resample.reads=resamp.read,
alpha=alpha.div,dca=dca.res,taxa=taxa.per)
output
}