
source("tran.npmle.r")
source("tran.npmle.sd.r")

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