## R script, need evd package 

## 3 input arguments:
## 	1: symmetric matrix
## 	2: pvalue
## 	3: outputfile name

## assume for each matrix col an extreme value distribution (gev), fit gev 
## compute for each matrix row the elements which are above a significance level pval
## output: 0/1 matrix, 1 = p(score)>pval  

library(evd)

args <- commandArgs(trailingOnly=TRUE)

matr<-as.matrix(read.table(args[1]));
outMat=matr
#print(matr)

pval=as.numeric(args[2])
outfile=args[3];
#print(pval)

mdim=dim(matr);
#print(mdim)

for (i in c(1:mdim[2])){
	#print(i)
	x=matr[,i]
	x=x[x!=0]
	#print(x);
	#x=rgev(100, loc = 10, scale = 1.1, shape = 0.2)
	rowGEV=fgev(x,std.err = FALSE,method="BFGS",control=list(maxit=1000))
	#print(rowGEV)
	rowLoc=as.vector(rowGEV$estimate[1])
	rowScale=as.vector(rowGEV$estimate[2])
	rowShape=as.vector(rowGEV$estimate[3])
	if (pval<1 && pval>0){
		cutoff=qgev((1-pval),loc=rowLoc,scale=rowScale,shape=rowShape)
	} else if (pval == 1){
		cutoff = 0
	} else{
		cutoff = max(x);
	}
	print(c("col=",i,"cutoff=",cutoff),quote=FALSE)
	for (j in c(1:mdim[1])){
		#print(c("MD:",mdim[1]));
		#print(c("cutoff=",cutoff,outMat[j,i]))
		if (outMat[j,i]>=as.numeric(cutoff)){
		#	print(c(j,cutoff,as.numeric(outMat[i,j])),quote=FALSE)
			outMat[j,i]=1
		}else{
			outMat[j,i]=0
		}
		if ( i==j ) {
			outMat[j,i]=1
		}
	
	}
}
#print(outMat)
#print(outfile);
if (!is.na(outfile)){
	write.table(outMat,file=outfile,row.names=FALSE,col.names=FALSE)
}
#warnings()