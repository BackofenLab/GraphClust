## annotate newick tree with classes
## class from class hash file
## each id in tree is replaced by its class, new id: "class$classidx"

$chf	= shift; ## class hash file
$tf		= shift; ## tree file
#$of		= shift; ## out file

open(HASH,$chf);
my %ch;
while ($line=<HASH>){
	@ent=split(" ",$line);
	$ch{$ent[0]}=$ent[1];
} 
close(HASH);

open(TREE,$tf);
my $tree = <TREE>;
close(TREE);
chomp $tree;

foreach my $key (sort keys %ch){
	my $class = $ch{$key};
	if ($tree =~ /\($key\:/ ){
		$tree =~ s/\($key\:/\(class$class\:/;
	} 
	if ($tree =~ /\,$key\:/){
		$tree =~ s/\,$key\:/\,class$class\:/;
	} 
}

print $tree."\n";

# system("echo '$tree' > $of");

