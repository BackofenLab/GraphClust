#!/bin/bash

## script that computes the probability that a certain subset of ids
## occurs by chance, subset is drawn from complete set without replacement
## therefore the hypergeometric distribution applies for that model 

## 1) each id in each line in $svectorfile is mapped to its class given by classHash file,
##    saved in cluster.class.all  
## 2) first line of cluster.class.all is used to get the class sizes of the full set (svectorfile)
##    i.e. it is assumed that each line (in svectorfile) contain the same class sizes,
##    perl is used to create file "class.size", format: count class
## 3) for each line in svector files the classes of a cluster 1..S are determined, saved in cluster.class
## 4) for each line in cluster.class, the largest class is determined
##     ouput values: class_idx:class_idx O:maxclass_occurence C:size_maxclass S:size_cluster D:size_fullset
##     output file: cluster.counts
## 5) R is used for each line of cluster counts, 
##    determine probability that in a drawn set of S we see O instances of class_idx of size C out of D instances 
##      

# cat data.names | awk '{print $1,$6}' > class.hash
# cat data.names | awk '{print $6}' | sort -n | uniq -c > class.size

classHash=$1
svectorFile=$2
clustFile=$3
S=$4
outDir=$5
classSizeFile=$6
datamap=$7

mkdir -p $outDir

#svectorFile='/home/heyne/SCRATCH_heyne/DATA/Dros_clustering/root_Rfam_testset/SVECTOR/data.svector.fast_cluster'
#classHash="/home/heyne/SCRATCH_heyne/DATA/Dros_clustering/root_Rfam_testset/FASTA/class.hash"
#classSize="/home/heyne/SCRATCH_heyne/DATA/Dros_clustering/root_Rfam_testset/FASTA/class.size"


cat $svectorFile | awk -v HF=$classHash '
 BEGIN{ 
 	while (getline < HF ) {CH[$1]=$2} 
 	} 
 	{for(i=1;i<=NF;i++)
 		{printf("%s ",CH[$(i)])} 
 	printf("\n")  
 	}' > $outDir/cluster.class.all 

cat $clustFile | awk -v S=$S  -v HF=$datamap '
 BEGIN{ 
 	while (getline < HF ) {CH[$1]=$2} 
 	} 
 	{A=S; if (S==-1) A=NF;
 	for(i=1;i<=A;i++)
 		{printf("%s ",CH[$(i)])} 
 	printf("\n")  
 	}' > $outDir/cluster.class.frags 
 
head -n 1 $outDir/cluster.class.all | perl -ane 'my %sh=();map{$sh{$_}++} @F; map{print $sh{$_}." ".$_."\n";}sort keys %sh;' | sort -n > $outDir/class.size
#classSize="$outDir/class.size"
classSize="$classSizeFile"

 cat $clustFile | awk -v S=$S  -v HF=$classHash '
 BEGIN{ 
 	while (getline < HF ) {CH[$1]=$2} 
 	} 
 	{A=S; if (S==-1) A=NF;
 	for(i=1;i<=A;i++)
 		{printf("%s ",CH[$(i)])} 
 	printf("\n")  
 	}' > $outDir/cluster.class 
 	
cat $outDir/cluster.class | awk -v S=$S -v D=$(cat $classSize | awk '{s+=$1}END{print s}') -v HF2=$classSize '
 	BEGIN{ 
 		while(getline < HF2 ) {SH[$2]=$1}
 		}
 	{
 	delete H; maxv=0; maxi=0; 
 	for(i=1;i<=NF;i++){H[$(i)]++} 
 	for (v in H){
 		if(maxv<H[v]&&v!=0){maxv=H[v];maxi=v}
 	}
 	A=S;if (A==-1) A=NF;
 	printf("class_idx: %s O: %s C: %s S: %s D: %s\n",maxi,maxv,SH[maxi],A,D)
 	}' > $outDir/cluster.counts
 	
 cat $outDir/cluster.counts | awk '{id=$2;O=$4;C=$6;S=$8;D=$10; st="echo \"D=" D "; C=" C "; S=" S "; O=" O "; r=0; for (k in c(O:min(S,C))) r=r+( (choose(C,k)*choose(D-C,S-k))/(choose(D,S))); print(r);\" | R --vanilla --slave";
 						printf("center_id: %s %s ",NR,$0); system(st);}'
