intTobit<-function(x,len=NA,vect=TRUE,num=TRUE)
{
  # convert decimal to binary.
  x1=rev(as.numeric(intToBits(x)))
  if(is.na(len))
  {
    begin=match(1,x1,length(x1))
  }else{
    begin=length(x1)-len+1
  }
  res=x1[begin:length(x1)]
  if(!vect)
  {
    res=paste(res,collapse = "")
    if(num)
    {
      res=as.numeric(res)
    }
  }
  res
}