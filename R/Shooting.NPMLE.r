
Shooting.NPMLE=function(A,Y,D,X,r,nlan){
#r=0.4;nlan=100
#setwd("/home/Yujihzhang/NPMLE package/npmle and aLasso")
source("NPMLE.r")
source("predict haz.r")
source("AIC BIC.r")
source("estR.r")
source("tran.npmle.sd.r")

BB=NPMLE(A,Y,D,as.matrix(X),r)
#IF=(BB$IF2+t(BB$IF2))/2
IF=-(BB$IF+t(BB$IF))/2
#IF=-BB$IF

DL=IF
XL=chol(DL ,pivot = FALSE)
YL=XL%*%BB$coef

shoot=function(lambda){
rstep=0
BBnew=BB$coef
P=length(BBnew)
BBold=rep(Inf,P)
lambdaK=lambda/abs(BB$coef)
kk=0
while(max(abs(BBnew-BBold))>10^-5&rstep<500){
rstep=rstep+1
print(rstep)
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
lan=nlan-1
lambda=max(abs(t(YL)%*%XL)*abs(BB$coef))
parL=c(0,exp(seq(log(10^-8),log(lambda),length.out=lan)))

for(i in 1:(lan+1)){
BCOLL=cbind(BCOLL,shoot(parL[i]))
#cat("i=",i,"\n")
}

#plot(parL,BCOLL[1,],type="l",ylim=c(0,max(BCOLL)),col=1,xlim=c(0,lambda),
#xlab=expression(lambda),ylab="Regression coefficient",
#main="Coefficient path",lwd=2)
#for(path in 2:NROW(BCOLL)){
#points(parL,BCOLL[path ,],type="l",col=path,lwd=2 )
#}
#r=1
nlambda=NCOL(BCOLL)
aic=bic=NULL
X=as.matrix(X)

R0=BB$bas.haz


for(ab in 1:nlambda){
ans=AIC.BIC(A,Y,D,X,1,BCOLL[,ab],R0)
aic=c(aic,ans$AIC)
bic=c(bic,ans$BIC)
#print(c("ab=",ab))
R0=ans$Rout
#print(c(parL[ab],ans$df,ans$like,ans$AIC,ans$BIC))
#plot(parL[1:ab],aic,type="l",ylim=c(min(c(aic,bic)),max(c(aic,bic))))
#points(parL[1:ab],bic,type="l",col=2)

}

choose.matrix=cbind(parL,aic,bic)
CAIC=which(aic==min(aic))
CBIC=which(bic==min(bic))

#plot(aic,type="l")
#points(bic,type="l",co=2)

coef.bic=BCOLL[,CBIC]
Rbic=ESTR(coef.bic,BB$bas.haz,Y,D,as.matrix(X),r)
coef.aic=BCOLL[,CAIC]
Raic=ESTR(coef.aic,BB$bas.haz,Y,D,as.matrix(X),r)

Bic1=ifelse(coef.bic!=0,1,0)
LBic.beta1=which(Bic1==1) 
LBic.beta2=which(Bic1!=1)
LBic.beta=c(LBic.beta1,LBic.beta2)

Aic1=ifelse(coef.aic!=0,1,0)
LAic.beta1=which(Aic1==1)
LAic.beta2=which(Aic1!=1)
LAic.beta=c(LAic.beta1,LAic.beta2)


IFbic=asyvar(coef.bic[LBic.beta],Rbic,Y,D,as.matrix(X)[,LBic.beta],r)
IFbic=-(IFbic$IF+t(IFbic$IF))/2

IFaic=asyvar(coef.aic[LAic.beta],Raic,Y,D,as.matrix(X)[,LAic.beta],r)
IFaic=-(IFaic$IF+t(IFaic$IF))/2

P=length(BB$coef)

GBIC=IFbic
LB11=length(LBic.beta1)
GB11=GBIC[1:LB11,1:LB11]
GB21=GBIC[1:LB11,-(1:LB11)]
GB12=t(GB21)
GB22=GBIC[-(1:LB11),-(1:LB11)]
invGB11=solve(GB11)

PLbic=parL[CBIC]

if(LB11>0){
if(LB11>1){EB=GB22-GB12%*%invGB11%*%GB21}
if(LB11==1){GB12=matrix(GB12);GB21=matrix(GB21)
EB=GB22-GB12%*%invGB11%*%t(GB21)}


if(LB11>1) AB11=diag(coef.bic[LBic.beta1])^2
if(LB11==1) AB11=(coef.bic[LBic.beta1])^2


tilde.GB11=GB11-PLbic*AB11
inv.tilde.GB11=solve(tilde.GB11)
if(LB11>1&LB11<P){covBic=invGB11-(invGB11-inv.tilde.GB11)%*% GB21 %*%solve(EB) %*% GB12 %*% (invGB11-inv.tilde.GB11)}
if(LB11==P){covBic=invGB11}

if(LB11==1){covBic=invGB11-(invGB11-inv.tilde.GB11)%*% t(GB21) %*%solve(EB) %*% GB12 %*% (invGB11-inv.tilde.GB11)}

sd.coef.bic=sqrt(diag(covBic))
}
if(LB11==0)sd.coef.bic=NA

GAIC=IFaic
LA11=length(LAic.beta1)
GA11=GAIC[1:LA11,1:LA11]
GA21=GAIC[1:LA11,-(1:LA11)]
GA12=t(GA21)
GA22=GAIC[-(1:LA11),-(1:LA11)]
invGA11=solve(GA11)

PLaic=parL[CAIC]


#print(LA11)

if(LA11>0){ 
if(LA11>1){EA=GA22-GA12%*%invGA11%*%GA21}
if(LA11==1){GA12=matrix(GA21);GA12=matrix(GA21)
EA=GA22-GA12%*%invGA11%*%t(GA21)}

if(LA11>1) AA11=diag(coef.aic[LAic.beta1])^2
if(LA11==1) AA11=(coef.aic[LAic.beta1])^2

tilde.GA11=GA11-PLaic*AA11
inv.tilde.GA11=solve(tilde.GA11)
if(LA11>1&LA11<P){covAic=invGA11-(invGA11-inv.tilde.GA11)%*% GA21 %*%solve(EA) %*% GA12 %*% (invGA11-inv.tilde.GA11)}
if(LA11==P){covAic=invGA11}
if(LA11==1){covAic=invGA11-(invGA11-inv.tilde.GA11)%*% t(GA21) %*%solve(EA) %*% GA12 %*% (invGA11-inv.tilde.GA11)}

sd.coef.aic=sqrt(diag(covAic))
}
if(LA11==0)sd.coef.aic=NA



Ansi=list(
npmle=BB$coef,haz=BB$b,
npmle.sd=BB$sd,
XL=XL,YL=YL,allB=BCOLL,
X=X,
alllambda=parL,lambda.aic=PLaic,lambda.bic=PLbic,AIC=min(aic),BIC=min(aic),
coef.aic=coef.aic,
coef.aic.sd=sd.coef.aic,
coef.bic=coef.bic,
coef.bic.sd=sd.coef.bic,
choose.X.aic=LAic.beta1,
choose.X.bic=LBic.beta1
)

return(Ansi)

}






#ans=Shooting.NPMLE(A,Y,D,X,1,100)


