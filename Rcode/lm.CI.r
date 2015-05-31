n<-50
x<-sample(40:70,n,rep=T)
y<-.7*x+rnorm(n,sd=5)
plot(x,y,xlim=c(20,90),ylim=c(0,80))
mylm<-lm(y~x)
abline(mylm,col="red")
newx<-seq(20,90)
prd<-predict(mylm,newdata=data.frame(x=newx),interval = c("confidence"), 
level = 0.90,type="response")
lines(newx,prd[,2],col="red",lty=2)
lines(newx,prd[,3],col="red",lty=2)