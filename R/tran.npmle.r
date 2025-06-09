library(compiler)
library(survival)
source("Ties-TTSMLike(beta-type).r")
source("predict haz.r")



tran.npmle=function(A,Y,D,X,r){
p=length(X)/length(Y)

lc=which(D==1)
Aa=c(A,(Y-A)[lc])
if(p==1)Xa=c(X,X[lc]) else Xa=rbind(X,X[lc,])
Ya=c(Y,Y[lc])
Da=c(D,D[lc])
DF=which(A==Y)
LDF=length(DF)
if(LDF>0){
if(p>1) {B=coxph(Surv(A[-DF],Y[-DF],D[-DF])~X[,-DF],method="breslow")$coef}  
if(p==1){B=coxph(Surv(A[-DF],Y[-DF],D[-DF])~X[-DF],method="breslow")$coef }
          }
if(LDF==0){ B=coxph(Surv(A,Y,D)~X,method="breslow")$coef}  


                                                 
int.ans=TTSMLike_beta(B,Aa,Ya,Da,Xa,r)          
#cat(r,"ok")
n=length(Y)                                      

R=predit.haz(sort(unique(Y)),int.ans$b,Y[D==1])+cumsum(runif(length(unique(Y)),0,10^-10))
#R= (1:length(unique(Y)))/length(unique(Y))
                                                   
B=c(int.ans$coef)                                   

                                  
Y=sort(Y,index.return=TRUE)
D=D[Y$ix]
                                                 
                                 #covariate matrix Xij_n*p ; n=sample size 
if(p==1){X=matrix(X[Y$ix],n,p)}                    
if(p>1){ X=matrix(X[Y$ix,],n,p)}
Y=Y$x

#--------------------------------------------------------------------------------------

Yp=unique(Y)
tiefun_1=function(cid){
sum(Y[D==1]==cid)
}
tiefun_2=function(cid){
sum(Y==cid)
}

nYjp=apply(matrix(Yp),1,tiefun_1)                 
Y.times=apply(matrix(Yp),1,tiefun_2)             
lc.jump=which(diff(c(0,Y))!=0)                   



np=length(Yp)                                    
#R=R+(1:np)/(np) #R=R[lc.jump]
#R=(1:np)/(np) #R=R[lc.jump] 

###########################  transformation function  ##################################
#--------------------------------------------------------------------------------------- 
# LA is cumulative based-line
# la is dLA
#---------------------------------------------------------------------------------------
G=function(t){
if(r==0){AA=t}
if(r>0){AA=log(1+r*t)/r }
AA}

g=function(t){
if(r==0){AA=1}
if(r>0){AA=1/(1+r*t) }
AA}

dg=function(t){
if(r==0){AA=0}
if(r>0){AA=-r/(1+r*t)^2 }
AA}

kpa=function(t){dg(t)/g(t)}
#############################################################################
# eat R
#############################################################################

estR=function(B,R){

dRk=diff(c(0,R))
if(p==1){BX= c(X)*c(B)}                             
if(p>1){BX=apply(matrix(t(t(X)*B),n,p),1,sum)}      

eBX=exp(BX)                                         
BXM=matrix(BX,n,np)
eBXM=exp(BXM)

RM=matrix(R,n,np,byrow=TRUE)
dtM=matrix(diff(c(0,Yp)) ,n,np,byrow=TRUE)
Sz=cbind(rep(1,n),exp(-G(RM*eBXM))[,-np])

muZM =matrix(apply(Sz*dtM,1,sum),n,np)
dtM.k=matrix(c(diff(Yp),0) ,n,np,byrow=TRUE)

Szgz=(exp(-G(RM*eBXM))*g(RM*eBXM)*eBXM*dtM.k)
muZk.1=t(apply(Szgz[,np:1],1,cumsum))[,np:1]
muZk=apply(muZk.1/muZM,2,sum)

up=(nYjp+dRk*muZk)
Rtise=rep(as.numeric(R),times=Y.times)
down=( cumsum((eBX*(g(Rtise*eBX)-kpa(Rtise*eBX)*D))[n:1])[n:1] )[lc.jump]
cumsum(up/down ) 


}


estB=function(B,R){

Rtise=rep(R,times=Y.times)

if(p==1){BX=X*c(B)}
if(p>1){BX=apply(matrix(t(t(X)*B),n,p),1,sum)}
eBX=exp(BX)
BXM=matrix(BX,n,np)
eBXM=exp(BXM)

  RM=matrix(R,n,np,byrow=TRUE)
  RM0=matrix(c(0,R[-np]),n,np,byrow=TRUE)

dtM=matrix(diff(c(0,Yp)) ,n,np,byrow=TRUE)

  Sz=exp(-G(RM0*eBXM))
  gz=g(RM0*eBXM)

muZ=apply(Sz*dtM,1,sum)

part3=apply(Sz*gz*eBXM*RM0*dtM,1,sum)/muZ
part2=-g(Rtise*eBX)*Rtise*eBX
part1=(kpa(Rtise*eBX)*Rtise*eBX+1)*D
partall=part1+part2+part3
equ=rep(NA,p)
for(h in 1:p){
Zi=X[,h]
equ[h]=sum(partall*Zi)
}
equ
}

#############################################################################
#J-matrix
##############################################################################

newton=function(B,R){
JM=matrix(NA,p,p)
eva=1/n
goal=estB(B,R)
for(j in 1:p){
ee=rep(0,p)
ee[j]<-eva

Be=B+ee
JM[j,]=(estB(Be,R)-goal)/eva
}
rse<-list(JM=JM,goal=goal)
return(rse)
}
###############################################################################
R0=0
#reptime=0
RRR=0
while(abs(sum(R0)-sum(R))>10^-3&RRR<50){
R0=R
R=estR(B,R0)
#reptime=reptime+1
#print(reptime)
RRR=RRR+1
}
goal=estB(B,R)

if(sum(goal^2)^0.5>10^(-3)){
BBnew=B
BBold=rep(Inf,p)
}
if(sum(goal^2)^0.5<10^(-3)){
BBold=BBnew=B
}

###############################################################################
kk=0
while((max(abs(BBold-BBnew))>10^-3|sum(goal^2)^0.5>10^(-6) )&kk<50){
BBold=BBnew
rrr=0
while(sum(goal^2)^0.5>10^(-3)&rrr<50){
resJM=newton(BBnew,R)
if(p==1){BBnew=BBnew-c(1/(resJM$JM)*resJM$goal)}
if(p>1){BBnew=BBnew-c(solve(resJM$JM)%*%resJM$goal)}
#BBnew=BBnew-c(solve(resJM$JM)%*%resJM$goal)
#cat(BBnew,resJM$goal,"\n")
goal=resJM$goal
rrr=rrr+1
}

R0=R
R=estR(BBnew,R0)

RRR=0
while(abs(sum(R0)-sum(R))>10^-3&RRR<10){
R0=R
R=estR(BBnew,R0)
RRR=RRR+1
#cat(RRR)
}





goal=estB(BBnew,R)
cat(BBnew," ",goal,kk,"\n")
kk=kk+1
}

#################################################################################
if(kk<600){
rse<-list(coef=BBnew,bas.haz=R,iter=kk)
return(rse)}

if(kk==600) cat("divergece","\n") 

}

###################################################################################

tran.npmle=cmpfun(tran.npmle)




