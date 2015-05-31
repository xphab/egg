rename.samp<-function(comi,samplist)
{
# to correct the sample names in community data and merge the data of the same sample.
# samplist is a two-column table, of which the first column is old names and the second is the correct names
# comi: a row is a sample.
comi=data.frame(comi)
old.name=rownames(comi)
new.name=as.vector(samplist[,2][match(old.name,samplist[,1])])
if(sum(is.na(new.name))>0)
{
warning("error! some old names are not in sample list, and were removed")
}
freq=data.frame(table(new.name))
new2.name=as.vector(freq[,1])
#new2.name=levels(as.factor(new.name))
comm=data.frame(matrix(,nrow=length(new2.name),ncol=ncol(comi)))
id.new.sig=(which(freq[,2]==1));id.new.nos=(which(freq[,2]!=1))
pbar <- txtProgressBar(min = 0, max = 20, style = 3)
comm[id.new.sig,]=comi[match(new2.name[id.new.sig],new.name),]
setTxtProgressBar(pbar, round((20*length(id.new.sig)/length(new2.name)),1))
for(i in 1:length(id.new.nos))
{
comm[id.new.nos[i],]=colSums(comi[which(new.name==new2.name[id.new.nos[i]]),])
setTxtProgressBar(pbar, round((20*(length(id.new.sig)+i)/length(new2.name)),1))
}
close(pbar)
rownames(comm)=as.vector(new2.name)
colnames(comm)=colnames(comi)
comm
}
 
