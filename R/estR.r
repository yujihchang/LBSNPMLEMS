
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




