#!/usr/bin/perl -w 

use strict;

## add class name to each row/col of given matrix, 
## matrix is ordered as given in names file
## class of each  id given with class hash file
## matrix output to sdtout, redirect to file as you want 

my $chf	= shift; ## class hash file
my $nf		= shift; ## names file
my $mf		= shift; ## matrix file;
my $prec   = shift; ## output number precicion

## read in class hash
my %ch;
my $eval_class = 0;
if ($chf ne "0"){
  open(HASH,$chf) ;
  while (my $line=<HASH>){
	 chomp($line);
	 my @ent=split(" ",$line);
	 $ch{$ent[0]}=$ent[1];
  } 
  close(HASH);
  $eval_class = 1;
}

## read in names, map idx, i.e. line-idx, to name
open(NAMES,$nf);
my @names;
my $idx=0;
while (my $line=<NAMES>){
	chomp($line);
	$names[$idx]=$line;
	$idx++;
} 
close(NAMES);

## create new matrix
my @newMat;

## col names
my $colNames1 = "# s ";
my $colNames2 = "s c ";
my $colNames3 = "ID ";
foreach my $name(@names){
	$colNames1 .= $name." " if ($eval_class);
	$colNames2 .= "c".$ch{$name}." " if ($eval_class);
  $colNames3 .= $name." ";
}

$newMat[0] = $colNames1."\n" if ($eval_class);
$newMat[1] = $colNames2."\n" if ($eval_class);
$newMat[0] = $colNames3."\n" if (!$eval_class);

## add class to each row
open(MAT,$mf);
$idx=0;
while (my $line=<MAT>){
	chomp $line;
	my @ent=split(" ",$line);
	foreach my $idxt (0..$#ent){
	  $ent[$idxt] = sprintf("%1.".$prec."f",$ent[$idxt]);
	}

	my $anot;
  $anot = "$names[$idx] c".$ch{$names[$idx]} if ($eval_class);
  $anot = "$names[$idx] " if (!$eval_class);
  
	push(@newMat,$anot." ".join(" ",@ent)."\n");
	$idx++;
}
close(MAT);

## output
print @newMat;
