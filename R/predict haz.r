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



