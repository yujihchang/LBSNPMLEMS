
TTSMLike_beta=function(B,A,Y,D,X,r){
n=length(Y)                  #sample size   
Y=sort(Y,index.return=TRUE)
D=D[Y$ix]
p=dim(X)[2]                   #number of covariate
if(sum(p)==0) p=1
X=matrix(X,n,p)              #covariate matrix Xij_n*p ; n=sample size 
X=X[Y$ix,]
if(sum(p)==0) p=1
X=matrix(X,n,p)          
A=A[Y$ix]

Y=Y$x

#---------------------------------------------------------------------------------------
Yp=unique(Y[D==1])
tiefun=function(cid){
sum(Y[D==1]==cid)
}
nYjp=apply(matrix(Yp),1,tiefun)
np=length(Yp)                    #number of failure
###########################  transformation function  ##################################
#--------------------------------------------------------------------------------------- 
# LA is cumulative based-line
# la is dLA
#---------------------------------------------------------------------------------------
LA=function(t){
if(r==0){AA=t}
if(r>0){AA=log(1+r*t)/r }
AA}

la=function(t){
if(r==0){AA=1}
if(r>0){AA=1/(1+r*t) }
AA}

dla=function(t){
if(r==0){AA=0}
if(r>0){AA=-r/(1+r*t)^2 }
AA}
kpa=function(t){dla(t)/la(t)}
###########################   estimating equation    ####################################
Chen2009=function(B,X,r){
BM=matrix(rep(B,each=n),n,p) # Bij=Bj for i=1~n  
SBX=apply(X*BM,1,sum)        # SBXi  i=i~n  is sum^p_j X_ij * B_j                         
R=rep(NA,np)
#---------------------------------------------------------------- calculate for based-line

R1=function(x){ sum(
ifelse(A<=Yp[1],1,0)*
ifelse(Yp[1]<=Y,1,0)*
LA(exp( SBX )*x) 
                    )-nYjp[1]}
R[1]=uniroot(R1,c(0,10))$root #-------------------------------------------------疊代的第一點

if(np>=2){
for(i in 2:np) {
R[i]=R[i-1]+nYjp[i]/sum(
                     ifelse(A<=Yp[i],1,0)*
                     ifelse(Yp[i]<=Y,1,0)*              
                     la(R[i-1]*exp( SBX ) )*
                               exp( SBX )  
                   )                         
                }
}
#################################### 第2部分 -- 建立jacobian matrix #########################

#根據 Chen 2009 的方式疊代出，dR(t,B)/dB  

RB=matrix(NA,p,np)
RB1=RB1D=ifelse( Yp[1]>=A & Y>=Yp[1] ,1,0)*la(exp(SBX)*R[1])*exp(SBX)

RB1U=matrix( RB1,p,n,byrow=TRUE)*t(X)*R[1]


RB1Dr=sum(RB1D)
RB[,1]=-apply(RB1U,1,sum)/rep(RB1Dr,p) #----------------------------------------疊代的第一點

if(np>=2){

for(i in 2:np) {
RBi=RBiD=ifelse( Yp[i]>=A & Y>=Yp[i] ,1,0)*la( exp(SBX)*R[i-1] )* exp(SBX)   

RBiU=matrix( RBi,p,n,byrow=TRUE)*
      (  t(X) +
                matrix( dla( exp(SBX)*R[i-1] )/ la( exp(SBX)*R[i-1] ), p , n , byrow=TRUE)*
              ( matrix( RB[,i-1],p,n) + matrix( R[i],p,n ) * t(X) )* 
                matrix( exp(SBX),p,n,byrow=TRUE) 
      )

RBiDr=sum(RBiD)

RB[,i]=as.vector(RB[,i-1]) - nYjp[i]*apply(RBiU,1,sum)/(rep(RBiDr,p)^2)

}
}
#-------------------------------------------------------------------------------------------------------

plug=function(t){

cid=sum(t>=Yp)
if(cid<1) ans=0
if(1<=cid & cid <=np) ans=R[cid]
if(np<cid) ans=R[np]
ans
}

Rin=apply(matrix(Y),1,plug)

plug.RB=function(t){

cid=sum(t>=Yp)
if(cid<1) ans=rep(0,p)
if(1<=cid & cid <=np) ans=RB[,cid]
if(np<cid) ans=RB[,np]
ans
}

RBin=matrix(NA,p,n)
for(i in 1:n){  RBin[,i]=plug.RB(Y[i])  }


#####################################  建立 S^(p)(t;B) ,p=0,1,2  ########################################
#建立                                                                                                   #
#    phi (Yp,BB)=  X  +      dg(Yp;X)/g(Yp;X)           *    (R(Yp)*X+dBR(Yp))    * e^(B*X)            #
#                -----   ----------------------------       ----------------      -------------        #   
#              crrm(X)+  dg( rrm(YpM)*exp(crrm(XBM)) )  *     crrm(RM)*crrm(X)   *exp(crrm(XBM))       #
#                        /g( rrm(YpM)*exp(crrm(XBM)) )       +crrm(dBRM)                               #
#                                                                                                      #
########################################################################################################
phi.t=function(i,j){
d=length(i)
SB=SBX[i]
if(d>1) esb=matrix(exp(SB),p,d,byrow=TRUE) else esb=exp(SB)
TX=t(X[i,])
as.vector(TX+dla( Rin[j]*esb )/la( Rin[j]*esb )*
(Rin[j]*TX+RBin[,j])*esb)
}


################################################################
S0j=function(j){
sum(ifelse(Y[j]>=A & Y>=Y[j],1,0)*la( Rin[j]*exp(SBX) )*exp(SBX))
                }

#----------------------------------------------------------------
#all S0(tj) at jump points
S0t=apply(matrix(which(D==1)),1,S0j)


#################################################################
S1j=function(j){
risk=which(Y[j]>=A & Y>=Y[j] )
d=length(risk)
SB=SBX[risk]
if(d>1) {esb=matrix(exp(SB),p,d,byrow=TRUE)
         Ans= apply(la(Rin[j]*esb)*esb*phi.t(risk,j),1,sum) } 
 
if(d==1){ esb=rep(exp(SB),each=p) 
          Ans=apply(matrix(la(Rin[j]*esb)*esb*phi.t(risk,j)),1,sum) 
         }
Ans
}
#----------------------------------------------------------------
#all S0(tj) at jump points
S1t=apply(matrix(which(D==1)),1,S1j)


#----------------------------------------------------------------
#all S1(tj)/S0(tj)  and {S1(tj)/S0(tj)}^2  at jump points

S10=S1t/matrix(S0t,p,sum(D),byrow=TRUE)
S10.2=S10%*%t(S10)

#################################################################


S2j=function(j){
risk=which(Y[j]>=A & Y>=Y[j] )
d=length(risk)
SB=SBX[risk]
esb=rep(exp(SB),each=p)
ge=(la(Rin[j]*esb)*esb)^(0.5)*phi.t(risk,j)
gm=matrix(ge,p,d)   
gm%*%t(gm)

}

#################################################################
#all S20=S2(tj)/S0(tj) at jump points
VS2=function(j){as.vector(S2j(j))}

S20=matrix(
apply(apply(matrix(which(D==1)),1,VS2)/matrix(S0t,p*p,sum(D),byrow=TRUE),1,sum)
,p,p)


#------------------------------------------------------------------------------
#建立 partial score

Ap=-S20+S10.2



###---------------------------------------------------------------------------------
#用計方程式決定收斂
jp=which(D==1)
Dphi.t=function(i){phi.t(jp[i],jp[i])}
if(p>1)eq=apply( mapply(Dphi.t,1:sum(D)),1,sum)-apply(S10,1,sum)
if(p==1)eq=sum(mapply(Dphi.t,1:sum(D)))-apply(S10,1,sum)

rse0<-list(equ=eq,J=Ap,R=R)
return(rse0)

###---------------------------------------------------------------------------------

                    }  #----end of Chen2009



k=0
path=NULL
BB=B
p=length(X)/n                #number of covariate
CH=Chen2009(BB,X,r)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
while(max(abs(CH$equ))>10^(-7)&k<500){

BB=BB-solve(CH$J)%*%CH$equ
CH=Chen2009(BB,X,r)
path=c(path,BB)
k=k+1
#cat(k,BB,"\t",max(abs(CH$equ)),"\n")

                                        }
##############################################################################
rse<-list(fal.time=Yp,tun.time=A,cen.index=D,coef=BB,base.haz=CH$R,k=k,equ=CH$equ)
return(rse)
}
##########################################################################################
##########################################################################################
##########################################################################################

#x is the timeing that you want to predit.
#R is the estimated jump function of hazard function (breslow-type). 
predit.haz=function(x,R,Y){
Ys=unique(sort(Y))
prefun=function(t){
cid=sum(t>=Ys)
if(cid<1) ans=0
if(1<=cid& cid<length(Ys)) ans=R[cid]
if(length(Ys)<=cid) ans=R[length(Ys)]
ans
}
apply(matrix(x), 1 ,prefun)
}
##########################################################################################
##########################################################################################
##########################################################################################

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
B=c(int.ans$coef)                                   
Y=sort(Y,index.return=TRUE)
D=D[Y$ix]
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
RRR=0
while(abs(sum(R0)-sum(R))>10^-3&RRR<50){
R0=R
R=estR(B,R0)
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
##########################################################################################
##########################################################################################
##########################################################################################

AIC.BIC=function(A,Y,D,X,r,B,R0){
p=NCOL(X)

lc=which(D==1)
Aa=c(A,(Y-A)[lc])
if(p==1)Xa=c(X,X[lc]) else Xa=rbind(X,X[lc,])
Ya=c(Y,Y[lc])
Da=c(D,D[lc])
DF=which(A==Y)
LDF=length(DF)
n=length(Y)                                       #sample size  
                                  
Y=sort(Y,index.return=TRUE)
D=D[Y$ix]
                                                 
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

nYjp=apply(matrix(Yp),1,tiefun_1)                  # 計算每個 failure 時間點 jump 的個數
Y.times=apply(matrix(Yp),1,tiefun_2)               # 計算每個時間點 jump 的個數  
lc.jump=which(diff(c(0,Y))!=0)                     # 排序後的 Y 在向量中哪個位置發勝 jump.

np=length(Yp)                                      #number of distinct obs's.

###########################  transformation function  ##################################
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
################################################

estR=function(B,R){

dRk=diff(c(0,R))
if(p==1){BX= c(X)*c(B)}                             # \beta^{'} X_i, i=1,2,3,...,n
if(p>1){BX=apply(matrix(t(t(X)*B),n,p),1,sum)}      # \beta^{'} X_i, i=1,2,3,...,n

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

Rold=rep(Inf,length(R0))
Rnew=R0
stop=0
while(max(abs(Rnew-Rold))>10^-3&stop<50){
Rold=Rnew
Rnew=estR(B,Rnew)
stop=stop+1
#print(c(max(abs(Rnew-Rold)),stop))
}

R=Rnew

###############################################
like=function(B,R){
if(p==1){BX=X*c(B)}
if(p>1){BX=apply(matrix(t(t(X)*B),n,p),1,sum)}
eBX=exp(BX)
BXM=matrix(BX,n,np)
eBXM=exp(BXM)

RM0=matrix(c(0,R[-np]),n,np,byrow=TRUE)
dtM=matrix(diff(c(0,Yp)) ,n,np,byrow=TRUE)
Sz=exp(-G(RM0*eBXM))
muZ=apply(Sz*dtM,1,sum)
part3=-log(muZ)

dR=diff(c(0,R))
part1=(BX+log(dR)+log(g(R*eBX) ))[D==1]
part2=-G(R*eBX) 

sum(part1)+sum(part2)+sum(part3)
}

k=sum(abs(B)>0)
likeBR=like(B,R)
AIC=2*k-2*likeBR
BIC=log(n)*k-2*likeBR

ans=list(AIC=AIC,BIC=BIC,Rout=R,like=likeBR,df=k)
return(ans)

}
##########################################################################################
##########################################################################################
###################################################################################

asyvar=function(B,R,Y,D,X,r){

B=c(B)

n=length(Y)                  #sample size   
Y=sort(Y,index.return=TRUE)
D=D[Y$ix]

p=length(B)                   #number of covariate

X=matrix(X,n,p)              #covariate matrix Xij_n*p ; n=sample size 

if(p>1){X=matrix(X[Y$ix,],n,p)}
if(p==1){X=matrix(X[Y$ix],n,p)}

Y=Y$x


#---------------------------------------------------------------------------------------
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

np=length(Yp)                    #number of distinct obs's.
#R=(1:np)/np #R=R[lc.jump]
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
# eat R                                                                     #
#############################################################################

estR=function(B,R){
dRk=diff(c(0,R))
if(p==1){BX=X*c(B)}
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
Rtise=rep(R,times=Y.times)
down=( cumsum((eBX*(g(Rtise*eBX)-kpa(Rtise*eBX)*D))[n:1])[n:1] )[lc.jump]
cumsum(up/down ) 
}
##########################################################################################
##########################################################################################
##########################################################################################

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
equi=NULL

for(h in 1:p){
Zi=X[,h]
equ[h]=sum(partall*Zi)
equi=rbind(equi,partall*Zi)
}
ansequ=list(equ=equ,equi=equi)
return(ansequ)
}

#############################################################################
#J-matrix
##############################################################################

JM=matrix(NA,p,p)
eva=1/(n)
goal0=estB(B,R)
goal=goal0$equ

for(j in 1:p){
ee=rep(0,p)
ee[j]<-eva

Be=B+ee

Re=estR(Be,R)
for(k in 1:5){Re=estR(Be,Re)}
JM[j,]=(estB(Be,Re)$equ-goal)/eva
}

JM2=0
if(p>1)goal2=goal0$equi
if(p==1)goal2=matrix(goal0$equi,p,n,byrow=TRUE)

if(p>1)for( ee in 1:n){JM2=JM2+goal2[,ee]%*%t(goal2[,ee]) }
if(p==1)for( ee in 1:n){JM2=JM2+goal2[,ee]*(goal2[,ee]) }

IF2=JM2/n
ans=list(sd=diag(solve(-JM))^0.5, IF=JM,sd2=diag(solve(JM2))^0.5,IF2=JM2)
return(ans)
}

###################################################################################


ESTR=function(B,R,Y,D,X,r){

B=c(B)

n=length(Y)                  #sample size   
Y=sort(Y,index.return=TRUE)
D=D[Y$ix]

p=length(B)                   #number of covariate

X=matrix(X,n,p)              #covariate matrix Xij_n*p ; n=sample size 

if(p>1){X=matrix(X[Y$ix,],n,p)}
if(p==1){X=matrix(X[Y$ix],n,p)}


Y=Y$x


#---------------------------------------------------------------------------------------
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



np=length(Yp)                    #number of distinct obs's.
#R=(1:np)/np #R=R[lc.jump]
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

if(p==1){BX=X*c(B)}
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

Rtise=rep(R,times=Y.times)

down=( cumsum((eBX*(g(Rtise*eBX)-kpa(Rtise*eBX)*D))[n:1])[n:1] )[lc.jump]

cumsum(up/down ) 

}

R=estR(B,R)
R01=0
R02=max(R)
stop=0
while(abs(R01-R02)>10^-6&stop<100){
stop=stop+1
R01=R02
R=estR(B,R)
R02=max(R)
#print(stop)
}

R
}

###################################################################################


NPMLE=function(A,Y,D,X,r){
ans=tran.npmle(A,Y,D,X,r)
B=ans$coef
Rb=ans$b
R=predit.haz(sort(Y),Rb,Y)
iter=ans$iter
SD=asyvar(B,Rb,Y,D,X,r) 

rse<-list(iter=iter,bas.haz=R,coef=B,sd=SD$sd,IF=SD$IF,sd2=SD$sd2,IF2=SD$IF2)
return(rse)
}



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


