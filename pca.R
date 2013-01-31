usair.dat<-source("chap3usair.dat")$value

usair.dat


usair.dat=t(usair.dat)
usair.dat

#
#���֍s��
usair.cor=cor(t(usair.dat))
usair.cor

X=usair.dat
X
Xvar=function(X){
	return((1/n)*rowSums(X))#�e�s�̘a
}

xvar=Xvar(X)
xvar#�W�{���σx�N�g��

n=ncol(usair.dat)
n
p=nrow(usair.dat)
p

#�W�{���U�����U�s��̏�����
S=matrix(0,nrow=p,ncol=p)
#S
#(X[,i]-xvar)%*%t(X[,i]-xvar)
for(i in 1:n){
	S=S+(X[,i]-xvar)%*%t(X[,i]-xvar)
}
S=(1/n)*S
S

#���֍s������߂�
#�W�{���U�����U�s��̏�����
R=matrix(0,nrow=p,ncol=p)

Z=matrix(0,nrow=p,ncol=n)
Z

#�W�{���U�����U�s��̏�����
#���

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



#R�̌ŗL�l
Revalues=eigen(R)$value
Revalues

#R�̌ŗL�x�N�g��
Revectors=eigen(R)$vectors
Revectors[,1]=-Revectors[,1]
Revectors



#��^��
#i�Ԗڂ̊�^��
for(i in 1:p){
	cat(Revalues[i]/sum(Revalues),",")
}


#�ݐϊ�^��(cumulative contribution ratio)
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



#-------���߂��킩��₷�����邽�߂Ɏ听���Ƀ}�C�i�X���|���ĕϊ�-----------

#�}�C�i�X���|����
#��1�听���́A�u�Z�݂₷���v��\���Ă���B(�����ƎҐ��ASO2��Pop��������C�����Ɋ֌W���Ă�����)
#�l���傫����΁A��C������J�̏��Ȃ��Z�݂₷�����Ɖ��߂ł���B

#�}�C�i�X���|����
#��2�听���́A�u�J�V�v��\���Ă���B(Precip��Days�̒l��-0.62,-0.67�ɒ��ڂ���)
#�l���傫����΁A�J�̏��Ȃ����Ɖ��߂ł���B

#�}�C�i�X���|����
#��3�听���́A�u�C��̃^�C�v�v��\���Ă���B(Neg.Temp:0.67713618,Precip:-0.49203630�ł��邩��A
#Neg.Temp��Precip���Δ䂵�Ă��邱�Ƃɒ��ڂ���)
#�l���傫���Ȃ�قǁA�C���������~���ʂ̑����s�s��\���B

y1=-y1;y2=-y2;y3=-y3

#------------------


plot(y1,y2,main="USAir data PC1-PC2",
	xlim=c(min(y1),max(y1)),
	ylim=c(min(y2),max(y2)),
	xlab="��1�听��(�Z�݂₷��)",
	ylab="��2�听��(�J�V)",type="n",lwd=2)

text(y1,y2,
labels=abbreviate(row.names(t(usair.dat))),cex=0.7,lwd=2)


plot(y1,y3,main="USAir data PC1-PC3",
	xlim=c(min(y1),max(y1)),
	ylim=c(min(y3),max(y3)),
	xlab="��1�听��(�Z�݂₷��)",
	ylab="��3�听��(�C��^�C�v)",type="n",lwd=2)

text(y1,y3,
labels=abbreviate(row.names(t(usair.dat))),cex=0.7,lwd=2)

plot(y2,y3,main="USAir data PC2-PC3",
	xlim=c(min(y2),max(y2)),
	ylim=c(min(y3),max(y3)),
	xlab="��2�听��(�J�V)",
	ylab="��3�听��(�C��^�C�v)",type="n",lwd=2)

text(y2,y3,
labels=abbreviate(row.names(t(usair.dat))),cex=0.7,lwd=2)