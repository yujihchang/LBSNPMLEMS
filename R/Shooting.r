
rm(list=ls(all=TRUE))
setwd("C:/Users/h93h96/Dropbox")
source("JingQin.r")
setwd("C:/Users/h93h96/Dropbox/§E¤é¹ü thesis/NPMLE package")
source("NPMLE.r")
source("predict haz.r")

AAH=BBH=NULL
AAC=BBC=NULL
pall=seq(0,5,length.out=100)

i=1
data=read.table(paste("C:/Users/h93h96/Dropbox/sim n 400/30%/30%/r=0/data",i,".txt",sep=""),header=TRUE)
#data=data[1:5,]
Y=data[,1]
A=data[,2]
X=cbind(data[,c(3,4)],runif(400),rbinom(400,1,0.5),rnorm(400),rnorm(400),rnorm(400),
rnorm(400),rnorm(400),rbinom(400,1,0.5),rbinom(400,1,0.5),rbinom(400,1,0.5),rnorm(400),rnorm(400),rnorm(400),
rbinom(400,1,0.5),rbinom(400,1,0.5),rbinom(400,1,0.5),rbinom(400,1,0.5),rbinom(400,1,0.5),
rbinom(400,1,0.5),rbinom(400,1,0.5),rbinom(400,1,0.5),rbinom(400,1,0.5),rbinom(400,1,0.5),
rbinom(400,1,0.5),rbinom(400,1,0.5),rbinom(400,1,0.5),rbinom(400,1,0.5),rbinom(400,1,0.5),
rbinom(400,1,0.5),rbinom(400,1,0.5),rbinom(400,1,0.5),rbinom(400,1,0.5),rbinom(400,1,0.5),
rbinom(400,1,0.5),rbinom(400,1,0.5),rbinom(400,1,0.5),rbinom(400,1,0.5),rbinom(400,1,0.5),
rbinom(400,1,0.5),rbinom(400,1,0.5),rbinom(400,1,0.5),rbinom(400,1,0.5),rbinom(400,1,0.5),
rbinom(400,1,0.5),rbinom(400,1,0.5),rbinom(400,1,0.5),rbinom(400,1,0.5),rbinom(400,1,0.5))
D=data[,5]

BB=NPMLE(A,Y,D,as.matrix(X),0)

DL=-(BB$IF+t(BB$IF))/2
XL=chol(DL)
YL=XL%*%BB$coef


shoot=function(lambda){
BBnew=BB$coef
P=length(BBnew)
BBold=rep(Inf,P)
lambdaK=lambda/abs(BB$coef)
while(max(abs(BBnew-BBold))>10^-5){
BBold=BBnew
for(k in 1:P){
AL=sum((YL-XL[,-k]%*%BBnew[-k])*XL[,k])
Xk2=sum(XL[,k]^2)
if(AL<0&abs(AL)>lambdaK[k])BBnew[k]=(AL+lambdaK[k])/(Xk2)
if(AL>0&abs(AL)>lambdaK[k])BBnew[k]=(AL-lambdaK[k])/(Xk2)
if(abs(AL)<=lambdaK[k])BBnew[k]=0
}
}
BBnew
}

BCOLL=NULL
lan=1000
lambda=max(abs(t(YL)%*%XL)*abs(BB$coef))
parL=seq(0,lambda,length.out=lan)

for(i in 1:lan){
BCOLL=cbind(BCOLL,shoot(parL[i]))
}


plot(parL,BCOLL[1,],type="l",ylim=c(0,max(BCOLL)),col=1,xlim=c(0,lambda),
xlab=expression(lambda),ylab="Regression coefficient",
main="Coefficient path",lwd=2)
for(path in 2:NROW(BCOLL)){
points(parL,BCOLL[path ,],type="l",col=path,lwd=2 )
}
#r=1
nlambda=NCOL(BCOLL)
aic=bic=NULL
X=as.matrix(X)

R0=BB$bas.haz

setwd("C:\\Users\\h93h96\\Dropbox\\§E¤é¹ü thesis\\NPMLE package\\LASO")
source("AIC BIC.r")

for(ab in 1:nlambda){
ans=AIC.BIC(A,Y,D,X,1,BCOLL[,ab],R0)
aic=c(aic,ans$AIC)
bic=c(bic,ans$BIC)
R0=ans$Rout
print(c(parL[ab],ans$df,ans$like,ans$AIC,ans$BIC))
plot(parL[1:ab],aic,type="l",ylim=c(min(c(aic,bic)),max(c(aic,bic))))
points(parL[1:ab],bic,type="l",col=2)

}

choose.matrix=cbind(parL,aic,bic)
CAIC=which(aic==min(aic))
CBIC=which(bic==min(bic))



plot(aic,type="l")
points(bic,type="l",co=2)
BCOLL[,CAIC]
BCOLL[,CBIC]



