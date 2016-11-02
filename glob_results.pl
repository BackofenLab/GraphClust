#!/usr/bin/perl
use warnings;
use strict;

use Getopt::Long;
use Cwd qw(abs_path getcwd);
use FindBin;
use lib "$FindBin::Bin";

use Array::Utils qw(:all);
use POSIX qw(ceil floor);
use List::Util qw/ min max /;

use GraphClust;
use vars qw(%CONFIG);
*CONFIG = \%GraphClust::CONFIG;

my $rootDir;
my $verbose;
my $merge_cluster_ol;
my $merge_overlap;
my $newInf = 0;

my $do_summary_only = 0;    ## do only f-measure summary based on last round
## write new soft partition (e.g. used for blacklist)
my $redo_all = 0;
my $round_last;

usage()
  unless GetOptions(
  "root=s"       => \$rootDir,
  "verbose"      => \$verbose,
  "cluster-ol=f" => \$merge_cluster_ol,
  "frag-ol=f"    => \$merge_overlap,
  "summary-only" => \$do_summary_only,
  "all"          => \$redo_all,
  "last-ci=i"    => \$round_last,
  "infernal_1_1" => \$newInf,
  );

my $bin_dir = abs_path($FindBin::Bin);

die "Please provide rootdir  with --root <DIR>! Exit...\n\n" if ( !$rootDir );

################################################################################
## Load options from root-dir/config file
%GraphClust::CONFIG = readConfigFile( "$rootDir/config", 1 );
SUBSECTION("options loaded for glob_results.pl");

#printConfig( \%CONFIG ) if ($verbose);
mkdir("$rootDir/RESULTS");

my $evaluate = $GraphClust::CONFIG{evaluate};
$merge_cluster_ol = $GraphClust::CONFIG{results_merge_cluster_ol} if ( !$merge_cluster_ol );
$merge_overlap = $GraphClust::CONFIG{results_merge_overlap} if ( !$merge_overlap );

my $min_cluster_size = $GraphClust::CONFIG{results_min_cluster_size};

my $cm_min_bitscore       = $CONFIG{cm_min_bitscore};
my $cm_bitscore_sig       = $CONFIG{cm_bitscore_sig};
my $cm_max_eval           = $CONFIG{cm_max_eval};
my $evaluate_min_overlap  = $CONFIG{evaluate_min_overlap};
my $evaluate_class0_as_fp = $CONFIG{evaluate_class0_as_fp};
my $do_randi              = 0;

print "evaluate        : $evaluate\n";
print "evaluate_min_ol : $evaluate_min_overlap\n" if ($evaluate);
print "evaluate class0 : $evaluate_class0_as_fp\n";
print "merge_cluster_ol: $merge_cluster_ol\n";
print "merge_overlap   : $merge_overlap\n";
print "min_cluster_size: $min_cluster_size\n";
print "bitscore cutoff : $cm_min_bitscore\n";
print "bitscore sig    : $cm_bitscore_sig\n";
print "evalue cutoff   : $cm_max_eval\n";
print "infernal version : 1_1\n" if($newInf);
print "new inf value  = $newInf \n";
################################################################################

my %cm_hitlists = ();
my $round_max   = 0;
my ( $class_Frags, $class_Size, $class_Names );

if ($evaluate) {
  ( $class_Frags, $class_Size, $class_Names ) = GraphClust::getEvalClassInfo("$rootDir/FASTA");
}

## remove soft partition to recompute all filterd cmsearch hits
system("rm -f $rootDir/RESULTS/partitions/final_partition.soft") if ($redo_all);

## reading existing soft partition, prevents to filter cmsearch hits again and again
## soft partition exactly contains filtered hits for fixed biscore cutoff/evalue
my $exist_part;
my %exist_part_used = ();
if ( -e "$rootDir/RESULTS/partitions/final_partition.soft" && -e "$rootDir/RESULTS/partitions/final_partition.used_cmsearch" ) {
  $exist_part = GraphClust::read_partition("$rootDir/RESULTS/partitions/final_partition.soft");

  my $used_cm = GraphClust::read_partition("$rootDir/RESULTS/partitions/final_partition.used_cmsearch");
  map { $exist_part_used{ $_->[0] } = 1 } @{$used_cm};
}

foreach my $file ( glob( "$rootDir" . "/CLUSTER/*.cluster/CMSEARCH/*.tabresult" ) ) {

  next if $file !~ /CLUSTER\/(\d+)\.(\d+)\.cluster\//;

  my $round = $1;
  my $index = $2;

  my $key = "$round.$index";
  next if ( $round_last && $round > $round_last );

  $round_max = $round if ( $round > $round_max );

  print "read $key... $file\n" if ($verbose);

  ## reuse already filtered tabresults from cmsearch, exactly what is in soft partition
  ## $key (eg. 1.10) becomes ->{NAME} of each hit
  my $cm_hits;
  if ( !exists $exist_part_used{$key} ) {
    if($newInf){
	 $cm_hits = GraphClust::read_CM_tabfile_ext_1_1( $file, $cm_min_bitscore, $cm_max_eval, $cm_bitscore_sig, $key );
    }
    else{
	 $cm_hits = GraphClust::read_CM_tabfile_ext_1_0( $file, $cm_min_bitscore, $cm_max_eval, $cm_bitscore_sig, $key );
     }
    $exist_part_used{$key} = 1;
  } else {
    my @exist_hits = grep { $_->[4] eq $key } @{$exist_part};
    $cm_hits = part2frags( \@exist_hits );
  }

  ## file is created when running, but empty
  ## keep all clusters with at least 2 hits
  next if ( !@{$cm_hits} || @{$cm_hits} < 2 );

  ## in greylist mode:
  ## check if we have non-greylist hits in cluster, otherwise skip cluster completely
  if ( -e "$rootDir/FASTA/data.greylist" ) {
    my %gl = ();
    open( GL, "$rootDir/FASTA/data.greylist" );
    map { my $key = $_; chomp($key); $gl{$key} = 1 } <GL>;
    close(GL);

    my $no_gl_hit = 0;
    map { $_->{SEQID} =~ /SEQ(\d+)/; my $id = $1; $no_gl_hit = 1 if ( !exists $gl{$id} ); } @{$cm_hits};

    if ( $no_gl_hit == 0 ) {
      print " ...skip cluster as we have only greylist elements left\n";
      next;
    }
  }

  $cm_hitlists{"$round.$index"}->{INDEX} = $index;
  $cm_hitlists{"$round.$index"}->{ROUND} = $round;
  $cm_hitlists{"$round.$index"}->{CLASS} = 0; ## used for do_overlap_matrix, is equal to CLASS_CLUSTER for tabfile hits only
  $cm_hitlists{"$round.$index"}->{HITS}    = [];
  $cm_hitlists{"$round.$index"}->{TABFILE} = $file;

  ## add infos from MODEL
  if ($evaluate) {

    my @tmp = readpipe("cat $rootDir/CLUSTER/$key.cluster/bestSubtrees");
    $tmp[0] =~ /class\s+(\S+).*names\s+(\S+)/;
    $cm_hitlists{$key}->{MODEL_CLASS} = $1;
    $cm_hitlists{$key}->{MODEL_NAMES} = $2;

  }

  ## add cm_hits to $cm_hitlists{$key}->{HITS}
  $cm_hitlists{$key}->{HITS} = $cm_hits;
}

#####################

my $ts = $merge_cluster_ol * 100;

system("mkdir -p $rootDir/EVAL/partitions");
system("mkdir -p $rootDir/EVAL/f_measure");
system("mkdir -p $rootDir/EVAL/overlap");
system("mkdir -p $rootDir/RESULTS/partitions");

my @allR = ();    ## store result line for each round
my @allF = ();    ## store result line for each round

system("touch $rootDir/RESULTS/partitions/final_partition.hard.merged");

$round_max = $round_last if ($round_last);

foreach my $round ( 1 .. $round_max ) {

  next if ( $round < $round_max && ( $do_summary_only || $round_last ) );

  print "\nround:$round " . "\n";
  my @keys = grep { $cm_hitlists{$_}->{ROUND} <= $round } sort { $cm_hitlists{$a}->{ROUND} <=> $cm_hitlists{$b}->{ROUND} || $cm_hitlists{$a}->{INDEX} <=> $cm_hitlists{$b}->{INDEX} } keys %cm_hitlists;

  ## map for old to new cluster names
  ## start: each cluster maps to itself, no merging etc
  my $old2new_map_ORIG = ();
  map { $old2new_map_ORIG->{$_} = $_ } @keys;

  my $all_hits = [];
  foreach my $clus (@keys) {
    push( @{$all_hits}, @{ $cm_hitlists{$clus}->{HITS} } );
  }
  @{$all_hits} = sort { $a->{SEQID} cmp $b->{SEQID} } @{$all_hits};

  my ( $sc0, $f0, $sc1, $f1, $sc2, $f2, $sc3, $f3 );

  ##############################
  ## soft partition
  my $old2new_map_SOFT = filterClusters( $all_hits, $old2new_map_ORIG );
  my $clus2idx = get_clus2idx( \%cm_hitlists, $old2new_map_SOFT );

  if ($evaluate) {
    print "\nSOFT PART until round $round:\n" if ($verbose);
    $f0 = get_f_measure( $all_hits, $clus2idx, $evaluate_min_overlap );

    write_res_summary( \%cm_hitlists, $f0, "$rootDir/EVAL/stage8.summary.soft" );

    if ( !$do_summary_only && $do_randi ) {
      $sc0 = get_rand_index( $all_hits, $clus2idx, $evaluate_min_overlap );
      print "RES: RI=" . $sc0->[0] . "\n" if ($verbose);
    }

   ## store clusClass for orginal hitlists, needed for merge map
   foreach my $key ( keys %cm_hitlists ) {
     $cm_hitlists{$key}->{CLASS} = $cm_hitlists{$key}->{HITS}->[0]->{CLASS_CLUSTER};
   }
  }
  system("rm -f $rootDir/RESULTS/partitions/final_partition.used_cmsearch");
  write_partition( $all_hits, $clus2idx, $old2new_map_SOFT, "$rootDir/EVAL/partitions/$round.soft" );

  next if ($do_summary_only);

  ###########################
  ## hard partition

  print "\nmake unique hits out of " . scalar( @{$all_hits} ) . "..." if ($verbose);
  my $hits = makeUniqueHits( $all_hits, $merge_overlap );
  @{$hits} = sort { $a->{SEQID} cmp $b->{SEQID} } @{$hits};
  print " result: " . scalar( @{$hits} ) . " hits\n" if ($verbose);

  ## hard BEST partition
  my $old2new_map_BEST = filterClusters( $hits, $old2new_map_ORIG );
  $clus2idx = get_clus2idx( \%cm_hitlists, $old2new_map_BEST );
  if ($evaluate) {
    print "\nHARD BEST PART until round $round:\n" if ($verbose);
    $f1 = get_f_measure( $hits, $clus2idx, $evaluate_min_overlap );
    if ($do_randi) {
      $sc1 = get_rand_index( $hits, $clus2idx, $evaluate_min_overlap );
      print "RES: RI=" . $sc1->[0] . "\n" if ($verbose);
    }
   ## store clusClass for orginal hitlists, needed for merge map
   foreach my $key ( keys %cm_hitlists ) {
     $cm_hitlists{$key}->{CLASS} = $cm_hitlists{$key}->{HITS}->[0]->{CLASS_CLUSTER};
   }
  }
  write_partition( $hits, $clus2idx, $old2new_map_BEST, "$rootDir/EVAL/partitions/$round.hard.best" );

  ## merged partitions
  ## overlap matrix is based on %cm_hitlists, i.e. on soft, not on hard partition of hits
  my $old2new_map_MERGED = do_overlap_matrix( $old2new_map_SOFT, "$rootDir/EVAL/overlap/$round.cm_overlap", $old2new_map_ORIG, $merge_overlap, $merge_cluster_ol );
  $old2new_map_MERGED = filterClusters( $hits, $old2new_map_MERGED );
  $clus2idx = get_clus2idx( \%cm_hitlists, $old2new_map_MERGED );

  if ($evaluate) {
    print "\nHARD MERGED PART until round $round:\n" if ($verbose);
    $f2 = get_f_measure( $hits, $clus2idx, $evaluate_min_overlap );
    if ($do_randi) {
      $sc2 = get_rand_index( $hits, $clus2idx, $evaluate_min_overlap );
      print "RES: RI=" . $sc2->[0] . "\n" if ($verbose);
    }

  }
  write_partition( $hits, $clus2idx, $old2new_map_MERGED, "$rootDir/EVAL/partitions/$round.hard.merged" );

  next if ( !$evaluate );

  ## ORACLE partition
  my $old2new_map_ORACLE = do_ORACLE_merge( \%cm_hitlists, $old2new_map_SOFT);
  $old2new_map_ORACLE = filterClusters( $hits, $old2new_map_ORACLE );
  $clus2idx = get_clus2idx( \%cm_hitlists, $old2new_map_ORACLE );
  print "\nHARD ORACLE PART until round $round:\n" if ($verbose);
  $f3 = get_f_measure( $hits, $clus2idx, $evaluate_min_overlap );
  if ($do_randi) {
    $sc3 = get_rand_index( $hits, $clus2idx, $evaluate_min_overlap );
    print "RES: RI=" . $sc3->[0] . "\n" if ($verbose);
  }
  write_partition( $hits, $clus2idx, $old2new_map_ORACLE, "$rootDir/EVAL/partitions/$round.hard.oracle" );

  write_f_measure( $f0, "$rootDir/EVAL/f_measure/round_$round.SOFT.fm" );
  write_f_measure( $f1, "$rootDir/EVAL/f_measure/round_$round.BEST.fm" );
  write_f_measure( $f2, "$rootDir/EVAL/f_measure/round_$round.MERGED.fm" );
  write_f_measure( $f3, "$rootDir/EVAL/f_measure/round_$round.ORACLE.fm" );

  write_res_summary( \%cm_hitlists, $f2, "$rootDir/EVAL/stage8.summary.merged" );
  write_res_summary( \%cm_hitlists, $f3, "$rootDir/EVAL/stage8.summary.oracle" );

  my $f0_a = join( " ", @{ ${$f0}[ $#{$f0} ] } );
  my $f1_a = join( " ", @{ ${$f1}[ $#{$f1} ] } );
  my $f2_a = join( " ", @{ ${$f2}[ $#{$f2} ] } );
  my $f3_a = join( " ", @{ ${$f3}[ $#{$f3} ] } );

  ## write file for all f measures
  my $all_F_str = "round $round HARD " . scalar( @{$hits} ) . " BEST $f1_a MERGED $f2_a ORACLE $f3_a SOFT $f0_a\n";
  ## write file for all rand scores
  if ($do_randi) {
    my $all_R_str = "round $round  HARD " . $sc1->[2] . " BEST "
      . $sc1->[0] . " " . $sc1->[1] . " MERGED " . $sc2->[0] . " "
      . $sc2->[1] . " ORACLE " . $sc3->[0] . " " . $sc3->[1] . " SOFT "
      . $sc0->[0] . " " . $sc0->[1] . " " . "\n";
    push( @allR, $all_R_str );
  }
  push( @allF, $all_F_str );
} ## foreach round

###################################
## write final data

## no cluster at all
exit if ($round_max == 0);

## soft part data
system("cp $rootDir/EVAL/partitions/$round_max.soft $rootDir/RESULTS/partitions/final_partition.soft");

open( OUT, ">$rootDir/RESULTS/partitions/final_partition.used_cmsearch" );
map { print OUT $_ . "\n" } keys %exist_part_used;
close(OUT);

if ( !$do_summary_only ) {
## hard best part data
  system("cp $rootDir/EVAL/partitions/$round_max.hard.best $rootDir/RESULTS/partitions/final_partition.hard.best");
## merged partition data
  system("cp $rootDir/EVAL/partitions/$round_max.hard.merged $rootDir/RESULTS/partitions/final_partition.hard.merged");
  system("cp $rootDir/EVAL/overlap/$round_max.cm_overlap $rootDir/RESULTS/partitions/final_overlap.matrix");
  system("cp $rootDir/EVAL/overlap/$round_max.cm_overlap.merge_map $rootDir/RESULTS/partitions/final_overlap.map");
}

if ( $evaluate && !$do_summary_only ) {

  my $r_all_file = "$rootDir/EVAL/rand_index.all.$ts";
  my $f_all_file = "$rootDir/EVAL/f_measure.all.$ts";

  open( ALLR, ">$r_all_file" );
  open( ALLF, ">$f_all_file" );

  map { print ALLF $_ } @allF;
  map { print ALLR $_ } @allR;

  close(ALLF);
  close(ALLR);

  #GraphClust::makeCol($r_all_file);
  #GraphClust::makeCol( $f_all_file );
}

GraphClust::SUBSECTION(" finished glob_results.pl");

exit;

################################################################################

## rand index all clusters, one instance in more than 1 cluster allowed
## rand index merged clusters, unique per merged cluster
##                             same instance in more than 1 cluster still allowed

sub do_overlap_matrix {
  my $map_soft               = $_[0];
  my $ol_file               = $_[1];
  my $merge_map             = $_[2];
  my $merge_ol              = $_[3];
  my $merge_cluster_overlap = $_[4];

  #my @keys = @{$keys_};
  my @keys = keys %{$map_soft};
  my %ol_map;

  @keys = sort { $cm_hitlists{$a}->{CLASS} <=> $cm_hitlists{$b}->{CLASS} } @keys;

  open( CMO, ">$ol_file" );
  print CMO "# MODEL " . join( " ", @keys ) . "\n";
  print CMO "MODEL CLASS " . join( " ", map { $cm_hitlists{$_}->{CLASS} } @keys ) . "\n";

  foreach my $i ( 0 .. $#keys ) {

    my @ols = ();
    my $best_ol = [ 0, 0, 0, 0, 0 ];    ## (ratio, #keys_i, #keys_j, #ol, #all)

    foreach my $j ( 0 .. $#keys ) {

      ## get overlaps, but do not ignore strand info
      my $listOLs = GraphClust::fragment_overlap( $cm_hitlists{ $keys[$i] }->{HITS}, $cm_hitlists{ $keys[$j] }->{HITS}, $merge_ol, 0 );

      ## get uniq overlap keys wrt. i, as hits could overlap each other
      my %uniq_ols_j = ();
      my %uniq_ols_i = ();
      foreach my $ol ( sort { $b->[2] <=> $a->[2] } @{$listOLs} ) {
        if ( !exists $uniq_ols_i{ $cm_hitlists{ $keys[$i] }->{HITS}->[ $ol->[0] ]->{KEY} } ) {

          $uniq_ols_j{ $cm_hitlists{ $keys[$j] }->{HITS}->[ $ol->[1] ]->{KEY} } = 1;
          $uniq_ols_i{ $cm_hitlists{ $keys[$i] }->{HITS}->[ $ol->[0] ]->{KEY} } = 1;
        }
      }

      my $ol  = keys %uniq_ols_j;
      my $num = scalar( @{ $cm_hitlists{ $keys[$i] }->{HITS} } );

      my $str;

      if ( $ol > 0 && $num > 0 ) {
        $str = sprintf( "%.1f", ( $ol / $num ) );
        $str .= "($ol/$num)";

        #print $keys[$i]." ".$keys[$j]." ol $ol/$num\n";
      } else {
        $str = "0";
      }

      # if ($keys[$j] eq "3.42" && $keys[$i] eq "3.5"){
      #  print "hallo\n";

# foreach my $ol (@{$listOLs}){
#  print "ol: ".join(" ",@{$ol})."\n";
#
#  print $cm_hitlists{ $keys[$i] }->{HITS}->[$ol->[0]]->{KEY}." ".$cm_hitlists{ $keys[$j] }->{HITS}->[$ol->[1]]->{KEY}."\n";
#}
#}

      $best_ol = [ ( $ol / $num ), $keys[$i], $keys[$j], $ol, $num ] if ( "$i" ne "$j" && $num > 0 && ( $ol / $num ) > $best_ol->[0] );

      push( @ols, $str );

    }
    $ol_map{ $best_ol->[1] . "#" . $best_ol->[2] } = $best_ol if ( $best_ol->[0] >= $merge_cluster_overlap );
    print CMO $keys[$i] . " " . $cm_hitlists{ $keys[$i] }->{CLASS} . " " . join( " ", @ols ) . "\n";

  }
  close(CMO);

  #GraphClust::makeCol($ol_file);

  ## sort all pairs by overlap ratio
  my @merge_pairs = sort { $ol_map{$b}->[0] <=> $ol_map{$a}->[0] } keys %ol_map;

  foreach my $pair (@merge_pairs) {

    my $parent = $merge_map->{ $ol_map{$pair}->[2] };
    my $child  = $merge_map->{ $ol_map{$pair}->[1] };

    next if ( $child eq $parent );

    #print "par $parent child $child ol ".$ol_map{$pair}->[0]."\n";
    ## continue only if bidirectional overlap > $merge_threshold
    #next if (!exists $ol_map{$parent."#".$child} );

    ## update previously into $child merged cluster with new parent
    foreach my $old_child ( keys %{$merge_map} ) {
      $merge_map->{$old_child} = $parent if ( $merge_map->{$old_child} eq $child );
    }

    $merge_map->{$child} = $parent;
  }

  open( MER, ">$ol_file.merge_map" );
  my @parents = values %{$merge_map};
  @parents = sort { $a cmp $b } unique(@parents);

  foreach my $par (@parents) {
    my @childs = sort { $a cmp $b } grep { $merge_map->{$_} eq $par } keys %{$merge_map};
    print MER $par . "\t" . join( " ", @childs ) . "\n";

  }

  close(MER);
  GraphClust::makeCol("$ol_file.merge_map");
  return $merge_map;
}

sub filterClusters {
  my $partition = $_[0];
  my $mergeMap  = $_[1];

  my @parentClusters = sort { $a <=> $b } unique( values %{$mergeMap} );

  my %mergeMap_filtered = ();

  foreach my $clus (@parentClusters) {

    my @hits;
    map { push(@hits,$_) if (exists $mergeMap->{ $_->{NAME} } && $mergeMap->{ $_->{NAME} } eq $clus)} @{$partition};
    #my @hits = grep { $mergeMap->{ $_->{NAME} } eq $clus } grep{ exists $mergeMap->{ $_->{NAME} } } @{$partition};

    if ( @hits >= $min_cluster_size ) {
      map { $mergeMap_filtered{$_} = $mergeMap->{$_} if ( $mergeMap->{$_} eq $clus ) } keys %{$mergeMap};
    }

  }

  return \%mergeMap_filtered;
}

## function that associates a running number to all final clusters according
## merge map
##
## clusters: hash with clusters and cmsearch hits
## merge_map is hash with keys: old cluster names; values: new clsuter name
## e.g. 1.1 -> 1.2
## return a hash with: key=final cluster name; value: idx
sub get_clus2idx {
  my $clusters  = $_[0];
  my $merge_map = $_[1];

  my %clus2idx = ();
  my @keys;

  ## @keys only contains the unique final/merged cluster names
  @keys = sort { $clusters->{$a}->{ROUND} <=> $clusters->{$b}->{ROUND} || $clusters->{$a}->{INDEX} <=> $clusters->{$b}->{INDEX} } unique( values %{$merge_map} );
  my %merg2idx = ();
  map { $merg2idx{ $keys[$_] } = ( $_ + 1 ) } 0 .. $#keys;

  ## map each cluster to the merged/final cluster id
  @keys = keys %{$merge_map};
  map { $clus2idx{$_} = $merg2idx{ $merge_map->{$_} } } @keys;

  return \%clus2idx;
}

sub do_ORACLE_merge {
  my $clusters  = $_[0];
  my $merge_map = $_[1];

  my %merge_map_NEW = ();
  my @keys = sort { $clusters->{$a}->{ROUND} <=> $clusters->{$b}->{ROUND} } keys %{$merge_map};
  my %classMap = ();
  ## first get one cluster-key for each class
  map { $classMap{ $clusters->{$_}->{CLASS} } = $_ } @keys;
  ## than map all clusters with the same class to that cluster-key
  map { if ( $clusters->{$_}->{CLASS} ne "0" ) { $merge_map_NEW{$_} = $classMap{ $clusters->{$_}->{CLASS} } } else { $merge_map_NEW{$_} = $_ } } @keys;

  return \%merge_map_NEW;
}

sub write_partition {
  my $partition = $_[0];    ## either all cmseach hits or unique cm search hits
  my $clus2idx  = $_[1];    ## new name -> idx
  my $merge_map = $_[2];    ## old -> new name
  my $outfile   = $_[3];

  print " write partition...";
  open( OUT, ">$outfile" );

  my @clus = sort { $a <=> $b } unique( map { $clus2idx->{$_} } keys %{$merge_map} );

  foreach my $idx (@clus) {

    my @hits = grep { exists $clus2idx->{ $_->{NAME} } && $clus2idx->{ $_->{NAME} } == $idx } @{$partition};

    foreach my $hit ( sort { $b->{BITSCORE} <=> $a->{BITSCORE} } @hits ) {
      my $idx   = $clus2idx->{ $hit->{NAME} };
      my $merge = $merge_map->{ $hit->{NAME} };

      ## seq-id-hit-key, score, parent-merge-cluster, child-merge-cluster, parent-merge-cluster-idx
      print OUT join( " ", $hit->{KEY}, $hit->{BITSCORE}, $hit->{EVALUE}, $merge, $hit->{NAME}, $idx );
      ## seq-true-class, seq-pred-class, class-idx

      if ($evaluate) {
        my $ckey = "-";
        $ckey = $hit->{CLASS_KEY} if ( $hit->{CLASS_KEY} ne "" );
        print OUT " " . join( " ", $hit->{CLASS}, $hit->{CLASS_CLUSTER}, $ckey, $hit->{CLASS_OL}, $hit->{CLASS_KEY_DUPL}, $class_Names->{ $hit->{CLASS} } );
      }
      print OUT "\n";
    }
  }
  close(OUT);

  GraphClust::makeCol("$outfile");
  print " finished!\n";
}

sub get_rand_index {
  my $partition    = $_[0];
  my $clus2idx     = $_[1];
  my $min_class_OL = $_[2];

  my %classHash_ = ();    ## used for missing class members
    # my ( $classFrags, $classSize ) = getEvalClassInfo("$rootDir/FASTA/class.hash.scan");

  map { $classHash_{ $_->{KEY} } = $_->{VALUE} } @{$class_Frags};
  my $classHash_size = keys %classHash_;

  my $max = 0; ## store max class index, used for correct class 0 idx for Rand index
  map { $max = $clus2idx->{$_} if ( $max < $clus2idx->{$_} ) } keys %{$clus2idx};

  my $outDir = "$rootDir/EVAL/_tmp_$$";
  system("mkdir -p $outDir");
  my $true_file = "$outDir/true.part";
  my $pred_file = "$outDir/pred.part";

  my @clus = sort { $a <=> $b } unique( values %{$clus2idx} );
  my $score = 0;

  ## check if we have clusters to analyze
  if (@clus) {

    open( TRUE, ">$true_file" );
    open( PRED, ">$pred_file" );

    my %classes = ();

    foreach my $idx (@clus) {

      my @hits = grep { exists $clus2idx->{ $_->{NAME} } && $clus2idx->{ $_->{NAME} } == $idx } @{$partition};

      ## assume CLASS is already assigned correctly, done by get_f_measure
#my $clus_class = GraphClust::evaluate_cm_hits( $classFrags, \@hits, $min_class_OL );
      my $clus_class = $hits[0]->{CLASS};

      foreach my $hit (@hits) {

        my $true = $hit->{CLASS};

        $true = $max + 1 if ( $true == 0 );

        print TRUE $true . "\n";
        print PRED $clus2idx->{ $hit->{NAME} } . "\n";

        delete $classHash_{ $hit->{CLASS_KEY} };
        $classes{ $hit->{CLASS_CLUSTER} } = 1;
      }
    }

    ## add missing elements from found classes
    my $midx = 0;
    foreach my $class ( keys %classes ) {

      my @tmp = grep { $classHash_{$_} eq $class && $class != 0 } keys %classHash_;

      #next;
      foreach my $hit_miss (@tmp) {
        my $true = $classHash_{$hit_miss};
        my $pred = $max + 1 + $true + $midx;

        print TRUE $true . "\n";
        print PRED $pred . "\n";

        $midx++;
      }

    }
    close(TRUE);
    close(PRED);

    my @call = readpipe( $CONFIG{PATH_OCTAVE} . "/octave -q $bin_dir/ari.m $true_file $pred_file 2> /dev/null" );
    chomp(@call);
    $score = $call[0];
  }
  $score = sprintf( "%.3f", $score );

  system("rm -R $outDir");
  my $num_clus = unique( values %{$clus2idx} );
  return ( [ $score, $num_clus, scalar( @{$partition} ), $classHash_size ] );
}

sub get_f_measure {
  my $partition    = $_[0];    ## hits
  my $clus2idx     = $_[1];
  my $min_class_OL = $_[2];

  print " calc f-measure...";
  my @res = ();
  my @clus = sort { $a <=> $b } unique( values %{$clus2idx} );

  foreach my $idx (@clus) {

    my @hits = grep { exists $clus2idx->{ $_->{NAME} } && $clus2idx->{ $_->{NAME} } == $idx } @{$partition};

    my $tp = 0;
    my $fp = 0;

    my $clus_class = GraphClust::evaluate_cm_hits( $class_Frags, \@hits, $min_class_OL );

    my @clusters = unique( map { $_->{NAME} } @hits );
    my @class    = unique( map { $_->{CLASS_CLUSTER} } @hits );

    print "\nWarning! Mergin of different classes! classes=" . join( ":", @class ) . " clus=" . join( ":", @clusters ) . "\n"
      if ( @class > 1 );

    foreach my $hit (@hits) {
      if ( $hit->{CLASS} ne "0" && $hit->{CLASS} eq $clus_class ) {
        $tp++ if ( !$hit->{CLASS_KEY_DUPL} );
      } else {
        $fp++ if ( ( $hit->{CLASS} eq "0" && $evaluate_class0_as_fp ) || $hit->{CLASS} ne "0" );
      }
    }

    #  ## True negatives depends on number of input seqs, not fragments
    #  ## this is correct as we mesaure F/R on scan over input set
    my $fn = 0;
    my $tn = 0;

    if ( $clus_class ne "0" ) {
      ## class size is based on class.size column 1 (all input class signals)
      $fn = $class_Size->{$clus_class} - $tp;
      ## use class 0 for true negatives
      $tn = scalar( $class_Size->{"0"} ); ## tn is only for accuracy, could be skipped
    } else {
      $tp = 0;
      $fp = 0;
    }

    my ( $f, $a, $p, $r, $s ) = f_measure( $tp, $tn, $fp, $fn );

    ## idx 0-6
    my @ent = ( $idx, "cluster", join( ":", @clusters[ 0 .. min( $#clusters, 10 ) ] ), "CLASS", $clus_class, "HITS", scalar(@hits) );
    ## idx 7-22
    push( @ent, ( "F", $f, "ACC", $a, "PREC", $p, "RECALL", $r, "TP", $tp, "FP", $fp, "FN", $fn, "TN", $tn ) );
    ## idx 23-24
    push( @ent, ( "ORIG_CLASS", join( ":", @class ) ) );

    push( @res, \@ent );
  }

  my @all = ( scalar(@res), "clusters", scalar(@res), "classes", 0, "hits", 0, "F", 0, "ACC", 0, "PREC", 0, "RECALL", 0 );
  my @cols         = ( 6, 8, 10, 12, 14 );    ## HITS F ACC PREC RECALL
  my %class_unique = ();
  foreach my $idx ( 0 .. $#res ) {
    map { $all[$_] += $res[$idx]->[$_] } @cols;
    $class_unique{ $res[$idx]->[4] } = 1 if ( $res[$idx]->[4] ne "0" );
  }

  ## check if we have clusters to analyze
  if (@res) {
    map { $all[$_] /= scalar(@res) } @cols[ 1 .. 4 ];   ## not for "hits"-column
    map { $all[$_] = sprintf( "%.3f", $all[$_] ) } @cols[ 1 .. 4 ]; ## not for "hits"-column
  }

  $all[4] = keys %class_unique;
  print "finished\n";
  print "RES: " . join( " ", @all ) . "\n" if ($verbose);
  push( @res, \@all );

  return \@res;
}

sub f_measure {
  my ( $tp, $tn, $fp, $fn ) = @_;

  if ( $tp != 0 ) {
    my $p = $tp / ( $tp + $fp );    ###precision
    my $r = $tp / ( $tp + $fn );    ###recall
    my $f = 2 * ( ( $p * $r ) / ( $p + $r ) );    ###f-measure
    my $s = 0;

    #my $s = $tn / ( $tn + $fp );  ###specifity
    my $a = ( $tp + $tn ) / ( $tp + $tn + $fp + $fn );    ###accuracy
    $f = sprintf( "%.3f", $f );
    $a = sprintf( "%.3f", $a );
    $p = sprintf( "%.3f", $p );
    $r = sprintf( "%.3f", $r );
    $s = sprintf( "%.3f", $s );
    return ( $f, $a, $p, $r, $s );

  } else {
    my $t = sprintf( "%.3f", 0 );
    return ( $t, $t, $t, $t, $t );
  }
}

sub write_f_measure {
  my $res  = $_[0];
  my $file = $_[1];

  open( OUT, ">$file" );

  map { print OUT join( " ", @{$_} ) . "\n" } @{$res};

  close(OUT);
  GraphClust::makeCol($file);
}

sub write_res_summary {
  my $hitlists = $_[0];
  my $fm       = $_[1];
  my $outfile  = $_[2];


  print " write summary... ";
  ## sort by iteration and then by f score, skip last entry (contains ALL info)
  my @fm_sort =
    sort { ( $a->[2] =~ /^(\d+)\.\d+/ )[0] <=> ( $b->[2] =~ /^(\d+)\.\d+/ )[0] || $b->[8] <=> $a->[8] } @{$fm}[ 0 .. $#{$fm} - 1 ];

  open( my $OUT_SUM, ">$outfile" );

  my $f_sum      = 0;
  my $curr_round = 0;

  my @sum_round = ();
  my $clus_count = 0; ## num clusters for average

  foreach my $idx ( 1 .. @fm_sort ) {

    my @res = @{ $fm_sort[ $idx - 1 ] };

    #print "IDX$idx ".join("#",@res)."\n";

    my $model_class;
    if ( $res[2] !~ /\:/ ) {
      $model_class = $hitlists->{ $res[2] }->{MODEL_CLASS}
    } else {
      $model_class = "MERGED:" . $res[24];
    }

    print $OUT_SUM " CLUSTER " . $res[0] . " CLASS " . $res[4] . " CLASS_NAME " . $class_Names->{ $res[4] } . " FMEASURE " . $res[8] . " MODEL_CLASS " . $model_class . " CENTER " . $res[2];
    print $OUT_SUM " HIT_NUM " . $res[6] . " " . join( " ", @res[ 15 .. 24 ] );
    print $OUT_SUM "\n";

    if ($evaluate_class0_as_fp || ($res[4] != "0" && !$evaluate_class0_as_fp) ){
      $f_sum += $res[8];
      $clus_count++;
    }


    ## do summary line for all results and for each interation

    if ( $idx - 1 == $#fm_sort || floor( ( $res[2] =~ /^(\d+)\.\d+/ )[0] ) ne floor( ( $fm_sort[$idx]->[2] =~ /^(\d+)\.\d+/ )[0] ) ) {

      $res[2] =~ /^(\d+)\.\d+/;
      $curr_round = $1;
      my $str = "ROUND $curr_round CENTERS " . sprintf( "%2d", $idx );

      my @centers_curr = @fm_sort[ 0 .. $idx - 1 ];
      @centers_curr = sort { $b->[8] <=> $a->[8] } @centers_curr;

      my %classes;
      map { $classes{ $_->[4] } = 1 if ( $_->[4] ne "0" ) } @fm_sort[ 0 .. ( $idx - 1 ) ];

      $str .= " HIT_CLASS_UNIQ " . ( keys %classes ) . " FMEASURE_AVG " . sprintf( "%.3f", $f_sum / $clus_count );

      my $cl_unique_ratio = "0.000";
      $cl_unique_ratio = sprintf( "%.3f", ( keys %classes ) / ( keys %{$class_Size} ) ) if ( keys %classes );

      $str .= " CL_FINAL " . @centers_curr . " CL_UNIQ_REL " . $cl_unique_ratio;

      my @f_steps = ( 0.9, 0.8, 0.7, 0.5 );
      my $fstep_str = " BETTER_FM";
      my $tf = 0;
      my @classes;

      foreach my $fstep (@f_steps) {

        @classes = ();
        my %cl2fm = ();
        foreach my $c (@centers_curr) {
          $cl2fm{$c->[4]} = $c->[8] if (!exists $cl2fm{$c->[4]} || $cl2fm{$c->[4]}<$c->[8] );
          push( @classes, $c->[4] ) if ( $c->[8] >= $fstep );
        }
        $fstep_str .= "  $fstep  " . @classes . " (" . unique(@classes) . ")";

        $tf = 0;
        map{ $tf += $cl2fm{$_} }unique(@classes);
      }

      $fstep_str .= " FM_AVG_BETTER_05 ".sprintf( "%.3f", $tf / unique(@classes) );

      $str .= "$fstep_str";
      push( @sum_round, $str );

    }    ## if new round

  }    ## foreach fmeasure

  print $OUT_SUM "\nSUMMARY\n\n";
  print $OUT_SUM join( "\n", @sum_round )."\n";

  close($OUT_SUM);
  GraphClust::makeCol($outfile);
  print "finished\n";
}

sub makeUniqueHits {
  my $hits    = $_[0];
  my $classOL = $_[1];

  my @hits_unique = ();

  ## sorting is required, but done already before
  ## @{$hits} = sort { $a->{SEQID} cmp $b->{SEQID} } @{$hits};

  my $hits_curr_start_idx = 0;
  my $hits_curr_end_idx   = 0;

  my $curr_frag;

  while ( $hits_curr_start_idx < @{$hits} ) {

    $curr_frag = $hits->[$hits_curr_start_idx];

    $hits_curr_end_idx = $hits_curr_start_idx;

    while ( $hits_curr_end_idx < @{$hits} - 1 && $curr_frag->{SEQID} eq $hits->[ $hits_curr_end_idx + 1 ]->{SEQID} ) { $hits_curr_end_idx++ }

    ## start..end is now for one seqid
    #    print "init: ".($hits_curr_end_idx - $hits_curr_start_idx)."\n";
    ## only one hit for current SEQID
    if ( $hits_curr_end_idx - $hits_curr_start_idx == 0 ) {

      push( @hits_unique, $curr_frag );

    } else {

      ## todo: check correctness!!!

      my @tmp       = @{$hits}[ $hits_curr_start_idx .. $hits_curr_end_idx ];
      my $curr_hits = \@tmp;

      foreach my $hit ( @{$curr_hits} ) {

        # print $hit->{KEY}."\n";
        my $tmp = [$hit];

        my $ols = GraphClust::fragment_overlap( $tmp, $curr_hits, $classOL );

        if ( !@{$ols} ) {
          push( @hits_unique, $hit );
          next;
        }

      #   print "OL:".$ols->[0]->[0]." ".$ols->[0]->[1]." ".$ols->[0]->[2]."\n";
        ## are there any overlapping hits with better score
        my @better = ();

        map { push( @better, $curr_hits->[ $_->[1] ] ) if ( $curr_hits->[ $_->[1] ]->{BITSCORE} > $hit->{BITSCORE} ) } @{$ols};

        ## than just continue and get them later
        next if ( @better > 0 );

        ## get hits with equal bitscore, but we have to select exactly one of them
        map { push( @better, $curr_hits->[ $_->[1] ] ) if ( $curr_hits->[ $_->[1] ]->{BITSCORE} == $hit->{BITSCORE} ) } @{$ols};
        ## add current hit
        push( @better, $hit );
        ## fix (arbitrary) order to check which of the equal hits to take
        @better = sort { $a->{NAME} cmp $b->{NAME} } @better;
        ## take hit if the first in the order is the current one
        push( @hits_unique, $hit ) if ( $better[0]->{NAME} eq $hit->{NAME} );
      }    ## foreach curr_hits
    }    ## else multiple hits on SEQID

    $hits_curr_start_idx = $hits_curr_end_idx + 1;
  }

  return \@hits_unique;
}

sub part2frags {
  my $part = $_[0];

  my @frags = ();

  foreach my $p ( @{$part} ) {

    my $f = GraphClust::str2frag( $p->[0] );
    $f->{BITSCORE} = $p->[1];
    $f->{EVALUE}   = $p->[2];
    $f->{NAME}     = $p->[4];

    push( @frags, $f );
  }

  @frags = sort { $a->{SEQID} cmp $b->{SEQID} } @frags;
  return \@frags;
}

#sub makeUniqueHitsNew {
#  my $hits    = $_[0];
#  my $classOL = $_[1];
#
#  my @hits_unique = ();
#
#  foreach my $hit ( @{$hits} ) {
#
#    # print $hit->{KEY}."\n";
#    my $tmp = [$hit];
#
#    my $ols = GraphClust::fragment_overlap( $tmp, $hits, $classOL );
#
#    if ( !@{$ols} ) {
#      push( @hits_unique, $hit );
#      next;
#    }
#
#    #   print "OL:".$ols->[0]->[0]." ".$ols->[0]->[1]." ".$ols->[0]->[2]."\n";
#    ## are there any overlapping hits with better score
#    my @better = ();
#
#    map { push( @better, $hits->[ $_->[1] ] ) if ( $hits->[ $_->[1] ]->{BITSCORE} > $hit->{BITSCORE} ) } @{$ols};
#
#    ## than just continue and get them later
#    next if ( @better > 0 );
#
#    ## get hits with equal bitscore, but we have to select exactly one of them
#    map { push( @better, $hits->[ $_->[1] ] ) if ( $hits->[ $_->[1] ]->{BITSCORE} == $hit->{BITSCORE} ) } @{$ols};
#    ## add current hit
#    push( @better, $hit );
#    ## fix (arbitrary) order to check which of the equal hits to take
#    @better = sort { $a->{NAME} cmp $b->{NAME} } @better;
#    ## take hit if the first in the order is the current one
#    push( @hits_unique, $hit ) if ( $better[0]->{NAME} eq $hit->{NAME} );
#  }
#
#  return \@hits_unique;
#}

## function that associates a running number to all final clusters according
## used mode
##
## clusters: hash with clusters and cmsearch hits
## mode: ALL MERGED ORACLE
## merge_map is hash with keys: old cluster names; values: new clsuter name
## e.g. 1.1 -> 1.2
## return a hash with: key=final cluster name; value: idx
#sub get_clus2idx {
#  my $clusters  = $_[0];
#  my $mode      = $_[1];
#  my $merge_map = $_[2];
#
#  my %clus2idx = ();
#  my @keys;
#
#  if ( $mode eq "ALL" ) {
#    @keys = sort { $clusters->{$a}->{ROUND} <=> $clusters->{$b}->{ROUND} || $clusters->{$a}->{INDEX} <=> $clusters->{$b}->{INDEX} } keys %{$merge_map};
#
#    map { $clus2idx{ $keys[$_] } = ( $_ + 1 ) } 0 .. $#keys;
#
#  } elsif ( $mode eq "MERGED" ) {
#
#    ## @keys only contains the unique final/merged cluster names
#    @keys = sort { $clusters->{$a}->{ROUND} <=> $clusters->{$b}->{ROUND} || $clusters->{$a}->{INDEX} <=> $clusters->{$b}->{INDEX} } unique( values %{$merge_map} );
#    my %merg2idx = ();
#    map { $merg2idx{ $keys[$_] } = ( $_ + 1 ) } 0 .. $#keys;
#
#    ## map each cluster to the merged/final cluster id
#    @keys = sort keys %{$merge_map};
#    map { $clus2idx{ $keys[$_] } = $merg2idx{ $merge_map->{ $keys[$_] } } } 0 .. $#keys;
#
#  } elsif ( $mode eq "ORACLE" ) {
#    @keys = sort { $clusters->{$a}->{ROUND} <=> $clusters->{$b}->{ROUND} } keys %{$merge_map};
#    my %classMap = ();
#    ## first get one idx for each class
#    map { $classMap{ $clusters->{ $keys[$_] }->{CLASS} } = ( $_ + 1 ) } 0 .. $#keys;
#    ## than map all clusters with the same class to that idx
#    map { $clus2idx{ $keys[$_] } = $classMap{ $clusters->{ $keys[$_] }->{CLASS} } } 0 .. $#keys;
#  }
#
#  return \%clus2idx;
#}
