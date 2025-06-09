
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