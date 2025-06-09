library(compiler)

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

###################################################################################
AIC.BIC=cmpfun(AIC.BIC)




