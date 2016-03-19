getwd()
setwd("/home/th4/gitrepo/pca")

usair.dat<-source("chap3usair.dat")$value

usair.dat


usair.dat=t(usair.dat)
usair.dat
dim(usair.dat) #p x n

#
#相関行列
usair.cor=cor(t(usair.dat))
usair.cor

X=usair.dat
X
dim(X)

Xvar=function(X){
    return((1/ncol(X))*rowSums(X))#各行の和
}

xvar=Xvar(X)
xvar#標本平均ベクトル

n=ncol(usair.dat)
n
p=nrow(usair.dat)
p

#標本分散共分散行列の初期化
S=matrix(0,nrow=p,ncol=p)
#S
#(X[,i]-xvar)%*%t(X[,i]-xvar)
for(i in 1:n){
    S=S+(X[,i]-xvar)%*%t(X[,i]-xvar)
}
S=(1/n)*S
S

#相関行列を求める
#標本分散共分散行列の初期化
R=matrix(0,nrow=p,ncol=p)

Z=matrix(0,nrow=p,ncol=n)
Z

#標本分散共分散行列の初期化
#基準化

for(i in 1:n){
    for(j in 1:p){
        Z[j,i]=(X[j,i]-xvar[j])/sqrt(S[j,j])
    }
}

Z


for(j in 1:p){
    for(k in 1:p){
        R[j,k]=S[j,k]/(sqrt(S[j,j])*sqrt(S[k,k]))
    }
}
R

Sys.sleep(15)


#Rの固有値
Revalues=eigen(R)$value
Revalues

#Rの固有ベクトル
Revectors=eigen(R)$vectors
Revectors[,1]=-Revectors[,1]
Revectors



#寄与率
#i番目の寄与率
for(i in 1:p){
    cat(Revalues[i]/sum(Revalues),",")
}


#累積寄与率(cumulative contribution ratio)
cum=0

for(j in 1:p){
    for(i in 1:j){
        cum=cum+Revalues[i]/sum(Revalues)
    }
    cat(cum,",")
    cum=0
}



X

y1=matrix(0,nrow=n,ncol=1)
for(i in 1:n){
    y1[i,]=t(Revectors[,1])%*%Z[,i]
    y2[i,]=t(Revectors[,2])%*%Z[,i]
    cat(y1[i,],y2[i,],"\n")
    
}
Z
y1

Revectors[,2]
y2=matrix(0,nrow=n,ncol=1)
for(i in 1:n){
    y2[i,]=t(Revectors[,2])%*%Z[,i]
}

y2

y3=matrix(0,nrow=n,ncol=1)
for(i in 1:n){
    y3[i,]=t(Revectors[,3])%*%Z[,i]
}

y3



#-------解釈をわかりやすくするために主成分にマイナスを掛けて変換-----------

#マイナスを掛ける
#第1主成分は、「住みやすさ」を表している。(製造業者数、SO2とPopが高く大気汚染に関係していそう)
#値が大きければ、大気汚染や雨の少ない住みやすい環境と解釈できる。

#マイナスを掛ける
#第2主成分は、「雨天」を表している。(PrecipとDaysの値が-0.62,-0.67に注目する)
#値が大きければ、雨の少ない環境と解釈できる。

#マイナスを掛ける
#第3主成分は、「気候のタイプ」を表している。(Neg.Temp:0.67713618,Precip:-0.49203630であるから、
#Neg.TempとPrecipが対比していることに注目する)
#値が大きくなるほど、気温が高い降水量の多い都市を表す。

y1=-y1;y2=-y2;y3=-y3

plot(y1,y2,main="USAir data PC1-PC2",
     xlim=c(min(y1),max(y1)),
     ylim=c(min(y2),max(y2)),
     xlab="第1主成分(住みやすさ)",
     ylab="第2主成分(雨天)",type="n",lwd=2)

text(y1,y2,
     labels=abbreviate(row.names(t(usair.dat))),cex=0.7,lwd=2)


plot(y1,y3,main="USAir data PC1-PC3",
     xlim=c(min(y1),max(y1)),
     ylim=c(min(y3),max(y3)),
     xlab="第1主成分(住みやすさ)",
     ylab="第3主成分(気候タイプ)",type="n",lwd=2)

text(y1,y3,
     labels=abbreviate(row.names(t(usair.dat))),cex=0.7,lwd=2)

plot(y2,y3,main="USAir data PC2-PC3",
     xlim=c(min(y2),max(y2)),
     ylim=c(min(y3),max(y3)),
     xlab="第2主成分(雨天)",
     ylab="第3主成分(気候タイプ)",type="n",lwd=2)

text(y2,y3,
     labels=abbreviate(row.names(t(usair.dat))),cex=0.7,lwd=2)

