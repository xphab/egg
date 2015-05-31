tran.ITSrdp<-function(file,conf=0.8)
{
# to transform the output of rdp ITS classifer to csv file we need according to a certain confidence
oldfile=read.table(file,header=F,sep=",",row.names=NULL,colClasses="character",fill=TRUE)
dim(oldfile)
id=((1:7)*2-2+(which(oldfile[8,]=="Fungi")))
new.tax=oldfile[8:nrow(oldfile),id]
new.conf=oldfile[8:nrow(oldfile),id+1]
for(i in 1:ncol(new.conf))
{
new.conf[,i]=as.numeric(sub("%","",new.conf[,i]))/100
}
colnames(new.tax)=c("domain","phylum","class","order","family","genus","species")
colnames(new.conf)=paste(colnames(new.tax),"conf",sep=".")
rownames(new.tax)=oldfile[8:nrow(oldfile),1]
rownames(new.conf)=oldfile[8:nrow(oldfile),1]
tax.conf=cbind(new.tax,new.conf)
if(is.na(conf)|conf==0)
{
output=list(taxa=new.tax,confidence=tax.conf)
}else{
id.ok=(new.conf>=conf)
new.tax[!id.ok]<-"unclassified"
output=list(taxa=new.tax,confidence=tax.conf)
}
output
}
