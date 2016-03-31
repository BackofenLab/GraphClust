#!/usr/bin/perl

use strict;
use Getopt::Long;
use Pod::Usage;
use POSIX qw(ceil floor);
use Cwd qw(abs_path getcwd);
use File::Path;

use FindBin;
use lib "$FindBin::Bin";

use Array::Utils qw(:all);
use List::Util qw/ min max /;

use GraphClust;

use vars qw(%CONFIG);
*CONFIG = \%GraphClust::CONFIG;

my $in_SGE_TASK_ID;    ## used for results stage
my $in_eval_mode;
my $in_verbose = 0;
my $in_root_dir;

GetOptions(
  "task-id=i"  => \$in_SGE_TASK_ID,
  "evaluate"   => \$in_eval_mode,
  "verbose"    => \$in_verbose,
  "root-dir=s" => \$in_root_dir,
);

my $bin_dir = abs_path($FindBin::Bin);

################################################################################
## Load options from root-dir/config file
%GraphClust::CONFIG = readConfigFile( "$in_root_dir/config", 1 );
SUBSECTION("options loaded for alignCenter.pl");
printConfig( \%CONFIG ) if ($in_verbose);

my $infernal_path = $CONFIG{PATH_INFERNAL};
my $rnaz_path     = $CONFIG{PATH_RNAZ};
my $vrna_path     = $CONFIG{PATH_VRNA};
my $tmp_path      = $CONFIG{PATH_TMP};

my $evaluate = $CONFIG{evaluate};
$evaluate = $in_eval_mode if ($in_eval_mode);

my $results_top_num = $CONFIG{results_top_num};

$in_root_dir = abs_path($in_root_dir);

my $part_file = "$in_root_dir/RESULTS/partitions/final_partition.soft";
if ( !-e "$in_root_dir/FASTA/data.greylist" && uc( $CONFIG{results_partition_type} ) =~ /HARD/ ) {
  $part_file = "$in_root_dir/RESULTS/partitions/final_partition.hard.merged";
}

system("perl $bin_dir/glob_results.pl --root $in_root_dir --all") if ( !-e $part_file );

print $CONFIG{results_partition_type} . ": use partition: $part_file\n";

## summary contains evaluation info for used partition
my $summary;
if ($evaluate) {

  my $sum_file = "$in_root_dir/EVAL/stage8.summary.soft";
  if ( !-e "$in_root_dir/FASTA/data.greylist" && uc( $CONFIG{results_partition_type} ) =~ /HARD/ ) {
    $sum_file = "$in_root_dir/EVAL/stage8.summary.merged";
  }

  if ( -e $sum_file ) {
    ## use this function to read summary file as well
    $summary = GraphClust::read_partition($sum_file);
    print "here summary $sum_file $summary\n" . @{$summary} . "\n";
    $evaluate = 0 if ( !@{$summary} );
  } else {
    $evaluate = 0;
  }
}

## part: array of [SEQ700#26#65#+,42.40,-,1.2,1.1,1,63,63]
my $part = GraphClust::read_partition($part_file);

my @res_clus_sort = sort { $a <=> $b } unique( map { $_->[5] } @{$part} );

my @fa_scan = GraphClust::read_fasta_file("$in_root_dir/FASTA/data.fasta.scan");

my @res_todo = ();
push( @res_todo, $in_SGE_TASK_ID ) if ($in_SGE_TASK_ID);
push( @res_todo, ( 1 .. @res_clus_sort ) ) if ( !$in_SGE_TASK_ID );

foreach my $res_idx (@res_todo) {

  my %clus_hits = ();
  my $clus_idx  = $res_clus_sort[ $res_idx - 1 ];

  foreach my $p ( grep { $_->[5] eq $clus_idx } @{$part} ) {
    my $key = $p->[0];
    $clus_hits{$key}         = {};
    $clus_hits{$key}->{TYPE} = "CMSEARCH";
    $clus_hits{$key}->{PART} = $p;
    $clus_hits{$key}->{KEY}  = $key;
  }

  print "Warning! Used partition $part_file for cluster $clus_idx contains multiple fragments with same location!\n\n"
    if ( keys %clus_hits != grep { $_->[5] eq $clus_idx } @{$part} );

  my @clus_keys  = keys %clus_hits;
  my $clus_frags = GraphClust::list2frags( \@clus_keys );
  map { $clus_hits{ $_->{KEY} }->{FRAG} = $_ } @{$clus_frags};

  my @orig_clus = unique( map { $clus_hits{$_}->{PART}->[4] } keys %clus_hits );
  my $clus_dir = "$in_root_dir/RESULTS/$clus_idx";

  print "cluster idx: " . $clus_idx . "\n";
  print "cluster hits: " . ( keys %clus_hits ) . "\n";
  print "custer orig: " . ( join " ", @orig_clus ) . "\n";
  print "cluster dir:" . $clus_dir . "\n";
  mkdir($clus_dir);

  ## read in model ids of a final(merged) cluster, could be >1 orig clusters in case of merging
  my %model_ids = ();
  foreach my $orig_cl (@orig_clus) {
    my @model_fa = GraphClust::read_fasta_file("$in_root_dir/CLUSTER/$orig_cl.cluster/MODEL/model.tree.fa");
    map { $model_ids{$_} = 1 } @{ $model_fa[1] };
  }

  ## annotate %clus_hits with TYPE=MODEL or TYPE=BLASTCLUST
  my ( $model_map, $bc_map ) = getHitMap( $in_root_dir, \%clus_hits, \%model_ids );

  ## write out current cluster as partition
  open( PART, ">$clus_dir/cluster.part" );
  map { print PART join( " ", @{ $clus_hits{$_}->{PART} } ) . "\n"; } keys %clus_hits;
  close(PART);
  GraphClust::makeCol("$clus_dir/cluster.part");

  ############################################################################
  ## write cluster.all file with detailed cluster infos
  open( OUT, ">$clus_dir/cluster.all" );

  open( EVAL, ">$clus_dir/cluster.all.eval" ) if ($evaluate);

  ## write model ids to cluster file
  foreach my $key ( sort { $clus_hits{$b}->{PART}->[1] <=> $clus_hits{$a}->{PART}->[1] } grep { $clus_hits{$_}->{TYPE} eq "MODEL" } keys %clus_hits ) {
    my $score = $clus_hits{$key}->{PART}->[1];
    my $clus  = $clus_hits{$key}->{PART}->[4];
    print OUT $key . " RESULT $clus_idx CM_SCORE $score MODEL $clus ";
    print OUT $fa_scan[2]->{ $clus_hits{$key}->{FRAG}->{SEQID} } . "\n";

    ## add evaluation info if avaliable
    if ($evaluate) {
      my @eval = @{ $clus_hits{$key}->{PART} };
      print EVAL $key . " RESULT $clus_idx CM_SCORE $score MODEL $clus ";
      print EVAL "CLUS_CLASS " . $eval[7] . " CLASS " . $eval[6] . " CLASS_NAME " . $eval[11] . " CLASS_KEY " . $eval[8] . " CLASS_OL " . $eval[9] . " ";
      print EVAL $fa_scan[2]->{ $clus_hits{$key}->{FRAG}->{SEQID} } . "\n";
    }

  }

  ## add model seqs with no cm hit as well -> score = 0, MODEL_NOHIT tag
  foreach my $mod ( keys %{$model_map} ) {
    next if ( exists $model_map->{$mod}->{HIT} );
    print OUT $model_map->{$mod}->{FRAG}->{KEY} . " RESULT $clus_idx CM_SCORE 0 MODEL_NOHIT 0 " . $fa_scan[2]->{ $model_map->{$mod}->{FRAG}->{SEQID} } . "\n";
  }

  ## write all other hit seqs
  foreach my $key ( sort { $clus_hits{$b}->{PART}->[1] <=> $clus_hits{$a}->{PART}->[1] } grep { $clus_hits{$_}->{TYPE} ne "MODEL" } keys %clus_hits ) {
    my $score = $clus_hits{$key}->{PART}->[1];
    my $clus  = $clus_hits{$key}->{PART}->[4];
    print OUT $key . " RESULT $clus_idx CM_SCORE $score " . $clus_hits{$key}->{TYPE} . " $clus ";
    print OUT $fa_scan[2]->{ $clus_hits{$key}->{FRAG}->{SEQID} } . "\n";

    ## add evaluation info if avaliable
    if ($evaluate) {
      my @eval = @{ $clus_hits{$key}->{PART} };
      print EVAL $key . " RESULT $clus_idx CM_SCORE $score " . $clus_hits{$key}->{TYPE} . " $clus ";
      print EVAL "CLUS_CLASS " . $eval[7] . " CLASS " . $eval[6] . " CLASS_NAME " . $eval[11] . " CLASS_KEY " . $eval[8] . " CLASS_OL " . $eval[9] . " ";
      print EVAL $fa_scan[2]->{ $clus_hits{$key}->{FRAG}->{SEQID} } . "\n";
    }
  }

  close(OUT);
  GraphClust::makeCol("$clus_dir/cluster.all");
  if ($evaluate) {
    close(EVAL);
    GraphClust::makeCol("$clus_dir/cluster.all.eval");
  }

  ######################################################
  ## make alignments for top5 and selected top seqs (config: $results_top_num)

  ## create top5 alignment
  my $fa_top = alignTopResults( \%clus_hits, \@fa_scan, 5, $clus_dir, "top5" );

  ## create selected topXX, XX > 5 required
  if ( $results_top_num > 5 && @clus_keys > 5 ) {
    $fa_top = alignTopResults( \%clus_hits, \@fa_scan, $results_top_num, $clus_dir, "top" );
  }

  ## write BED file for larger top set
  fasta2BED( $fa_top, "$clus_dir/cluster.bed", "CLUSTER_$clus_idx" );

  ############################################################################
  ## write fasta file for all hits in cluster
  my $fa_file_all = "$clus_dir/" . "cluster.all.fa";
  writeClusterFasta( \@clus_keys, \%clus_hits, \@fa_scan, $fa_file_all );

  ############################################################################
  ## write stats file for global results stats
  my @orig_seqs = ();
  foreach my $key ( grep { $clus_hits{$_}->{PART}->[1] > 0 } sort { $clus_hits{$b}->{PART}->[1] <=> $clus_hits{$a}->{PART}->[1] } keys %clus_hits ) {
    my $header = $fa_scan[2]->{ $clus_hits{$key}->{FRAG}->{SEQID} };
    $header =~ /ORIGID\s+(\S+)/;
    push( @orig_seqs, $1 );
  }

  my @ids_unique = unique(@orig_seqs);

  ## mpi,she,zsc,sci,svm,...
  my $top5_rnaz = [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ];

  if ( -e "$clus_dir/locarna.top5/results/result.aln" ) {
    $top5_rnaz = GraphClust::getRNAzScores( "$clus_dir/locarna.top5/results/result.aln", "$rnaz_path/RNAz", "-n" );
    print "RNAz: " . join( ":", @{$top5_rnaz} ) . "\n";
  }

#    my $seqs_bc = grep { $hits_all{$_}->[0] > 0 && $hits_all{$_}->[1] eq "BLASTCL" } sort { $hits_all{$b}->[0] <=> $hits_all{$a}->[0] } keys %hits_all;
  open( STAT, ">$clus_dir/cluster.stats" );

  print STAT "CLUSTER $clus_idx SEQS " . scalar(@clus_keys) . " ";
  print STAT "IDS_UNIQUE " . scalar(@ids_unique) . " MODELS " . scalar(@orig_clus) . " MPI_TOP5 $top5_rnaz->[0] SCI_TOP5 $top5_rnaz->[9] ";

  if ($evaluate) {
    my @res = grep { $_->[1] == $clus_idx && $_->[0] eq "CLUSTER" } @{$summary}[ 0 .. $#res_clus_sort ];
    print STAT "CLUSTER_CLASS " . $res[0]->[3] . " CLUSTER_NAME " . $res[0]->[5] . " F_SCORE " . $res[0]->[7];
  }

  # print STAT "IDS_UNIQUE_LIST ".join(";",@ids_unique); ## too long for column
  print STAT "\n";

  close(STAT);

}    ## for all @res_todo

exit;

################################################################################

sub alignTopResults {
  my $clusHits = $_[0];
  my $faScan   = $_[1];
  my $top_num  = $_[2];
  my $resDir   = $_[3];
  my $name     = $_[4];

  my @sorted_hits = grep { $clusHits->{$_}->{PART}->[1] > 0 } sort { $clusHits->{$b}->{PART}->[1] <=> $clusHits->{$a}->{PART}->[1] } keys %{$clusHits};

  my $max = @sorted_hits;
  $max = $top_num if ( $max > $top_num );
  my @hitsTop = @sorted_hits[ 0 .. ( $max - 1 ) ];

  ## write out fasta with selected seqs, check direction of cm hit (+1,-1)
  my $fa_file_top = "$resDir/" . "cluster.$name.fa";

  writeClusterFasta( \@hitsTop, $clusHits, $faScan, $fa_file_top );

  my @hit_set = GraphClust::read_fasta_file($fa_file_top);
  my $useLocP = 1;

  if ( !-e "$resDir/locarna.$name/results/result.aln" && @{ $hit_set[1] } > 1 ) {
    GraphClust::mlocarna_center( $fa_file_top, "$resDir/locarna.$name", "", $useLocP );
    GraphClust::aln2alifold( "$resDir/locarna.$name/results/result.aln", $tmp_path, $vrna_path );
    system("cp $resDir/locarna.$name/results/result.aln.ps $resDir/cluster.$name.aln.ps");
    system("cp $resDir/locarna.$name/results/result.aln.alirna.ps $resDir/cluster.$name.alirna.ps");

    system("perl $bin_dir/mloc2stockholm.pl --split_input yes --con_struct $resDir/locarna.$name/results/result.aln.alifold -file $resDir/locarna.$name/results/result.aln");
    system("$infernal_path/cmbuild -F $resDir/locarna.$name/results/result.aln.cm $resDir/locarna.$name/results/result.aln.sth");
    system("cp  $resDir/locarna.$name/results/result.aln.cm $resDir/cluster.$name.cm");
  }

  return $fa_file_top;
}

sub getTrueLocation {
  my $loc = $_[0];    ## data.locations loc identified by graphFasta.pl
  my $seq = $_[1];    ## seq of data.scan.fasta to get available length
  my $hit = $_[2];    ## fragment datastructure of cm hit

  return $loc if ( $loc =~ /MISS/ );

  ##hg19.chr1:949858-949920:+
  my @ent = split( ":", $loc );
  my @pos = split( "-", $ent[1] );

  ## GraphClust locations always have smaller start than end, but check anyway
  my $loc_len = abs( $pos[1] - $pos[0] );

  return $loc if ( $loc_len != length($seq) );

  my $new_strand = "";
  ## + + or - -
  $new_strand = "+" if ( $hit->{STRAND} eq $ent[2] );
  ## + - or - +
  $new_strand = "-" if ( $hit->{STRAND} ne $ent[2] );

  ## example for all 4 cases:
  ##
  ## chr1:100-250:-
  ## hit : 33- 87:+
  ## true:163-217:-

  ## chr1:100-250:+
  ## hit : 33- 87:-
  ## true:133-187:-

  ## chr1:100-250:+
  ## hit : 33- 87:+
  ## true:133-187:+

  ## chr1:100-250:-
  ## hit : 33- 87:-
  ## true:163-217:+

  ## overwrite old start and end
  if ( $ent[2] eq "+" ) {
    $pos[0] = $pos[0] + $hit->{START} - 1;
    $pos[1] = $pos[0] + $hit->{STOP} - $hit->{START};
  } else {
    $pos[0] = $pos[1] - $hit->{END} - 1;
    $pos[1] = $pos[0] + $hit->{STOP} - $hit->{START};
  }

  $loc = $ent[0] . ":" . $pos[0] . "-" . $pos[1] . ":" . $new_strand;

  return $loc;
}

sub writeClusterFasta {
  my $hits         = $_[0];
  my $clusHits     = $_[1];
  my $faScan       = $_[2];
  my $out_filename = $_[3];

  open( FA, ">$out_filename" );

  my %locations = ();
  open( LOCS, "$in_root_dir/FASTA/data.locations" );
  while ( my $line = <LOCS> ) {
    chomp($line);
    my @ent = split( " ", $line );
    $locations{ $ent[0] } = $ent[1];
  }
  close(LOCS);

  foreach my $frag ( @{$hits} ) {

    my $seq;
    my $start = $clusHits->{$frag}->{FRAG}->{START} - 1;
    my $len = $clusHits->{$frag}->{FRAG}->{STOP} - $clusHits->{$frag}->{FRAG}->{START} + 1;

    my $seqid = $clusHits->{$frag}->{FRAG}->{SEQID};
    if ( $clusHits->{$frag}->{FRAG}->{STRAND} eq "+" ) {
      $seq = substr( $faScan->[0]->{$seqid}, $start, $len );
    } elsif ( $clusHits->{$frag}->{FRAG}->{STRAND} eq "-" ) {
      $seq = substr( $faScan->[0]->{$seqid}, $start, $len );
      $seq =~ tr/AUGC/UACG/;
      $seq = reverse($seq);
    } else {
      die "Strand Error for frag $frag! Exit...\n\n";
    }

    my $part = $clusHits->{$frag}->{PART};
    my $fr   = $frag;
    my $loc  = $locations{ $clusHits->{$frag}->{FRAG}->{SEQID} };

    $loc = getTrueLocation( $loc, $faScan->[0]->{$seqid}, $clusHits->{$frag}->{FRAG} );

    $fr =~ s/#/_/g;
    print FA ">$fr" . " RESULT " . $part->[5] . " SCORE " . $part->[1] . " EVALUE " . $part->[2];
    print FA " CLUSTER " . $part->[4] . " LOC $loc ";

    ## add evaluation info if avaliable
    if ( @{$part} >= 12 ) {
      print FA "CLUS_CLASS " . $part->[7] . " CLASS " . $part->[6] . " CLASS_NAME " . $part->[11] . " CLASS_KEY " . $part->[8] . " CLASS_OL " . $part->[9] . " ";
    }

    print FA "$faScan->[2]->{$seqid}\n";
    print FA "$seq\n";
  }
  close(FA);
}

sub getHitMap {
  my $rootDir    = $_[0];
  my $clus_href  = $_[1];
  my $model_href = $_[2];

  opendir( BLASTCL, "$rootDir/FASTA/" );
  my @blast_cl = readdir(BLASTCL);
  closedir(BLASTCL);

  my %blast_cl = ();

  foreach my $file (@blast_cl) {
    next if ( $file !~ /intermediate.blast_clusters.(\d+)/ );
    open( IN, "$rootDir/FASTA/$file" );
    while ( my $line = <IN> ) {
      chomp($line);
      my @ent = split( " ", $line );

      foreach my $idx ( 0 .. $#ent ) {
        my @push = grep { $_ =~ /SEQ\d+#\d+#\d+#/ } map { $ent[$_] if ( $_ != $idx ) } 0 .. $#ent;

        #print "push($ent[$idx])".join(":",@push)."\n";
        if ( exists( $blast_cl{ $ent[$idx] } ) ) {
          push( @{ $blast_cl{ $ent[$idx] } }, @push );
        } else {
          $blast_cl{ $ent[$idx] } = \@push;
        }
      }

    }

    close(IN);

  }

  my @bc_keys  = keys %blast_cl;
  my $bc_frags = GraphClust::list2frags( \@bc_keys );

  my $hit_frags = [];
  map { push( @{$hit_frags}, $clus_href->{$_}->{FRAG} ) } keys %{$clus_href};
  @{$hit_frags} = sort { $a->{SEQID} cmp $b->{SEQID} } @{$hit_frags};
  my $bc_ols = GraphClust::fragment_overlap( $hit_frags, $bc_frags, 0.51 );

  my %blast_hits = ();

  foreach my $ol ( @{$bc_ols} ) {
    $blast_hits{ $hit_frags->[ $ol->[0] ]->{KEY} } = 1;
    $clus_href->{ $hit_frags->[ $ol->[0] ]->{KEY} }->{TYPE} = "BLASTCLUST";
  }

  ###################################
  my $all_frags = GraphClust::read_fragments("$rootDir/FASTA/data.map");

  my %model_map   = ();
  my @model_frags = ();

  map { $model_map{ $_->{VALUE} }->{FRAG} = $_ if ( exists $model_href->{ $_->{VALUE} } ) } @{$all_frags};
  map { push( @model_frags, $model_map{$_}->{FRAG} ) } keys %model_map;

  @model_frags = sort { $a->{SEQID} cmp $b->{SEQID} } @model_frags;

  my $model_ols = GraphClust::fragment_overlap( \@model_frags, $hit_frags, 0.51 );

  foreach my $ol ( @{$model_ols} ) {
    $clus_href->{ $hit_frags->[ $ol->[1] ]->{KEY} }->{TYPE} = "MODEL";
    $model_map{ $model_frags[ $ol->[0] ]->{VALUE} }->{HIT} = $hit_frags->[ $ol->[1] ]->{KEY};
  }

  return ( \%model_map, \%blast_hits );
}

sub fasta2BED {

  my $fa_file  = $_[0];
  my $bed_file = $_[1];
  my $info     = $_[2];

  my @fa = GraphClust::read_fasta_file($fa_file);

  open( BED, ">$bed_file" );
  foreach my $id ( @{ $fa[1] } ) {
    my $header = $fa[2]->{$id};

    print $header. "\n";

    #$header =~ /range=(\S+)\:(\d+)\-(\d+)/;
    #my $chr   = 0;
    #my $start = 0;
    #my $end   = 0;
    #$chr   = $1;
    #$start = $2 - 1;
    #$end   = $3;

    ##hg19.chr1:949858-949920:+
    $header =~ /LOC\s+(\S+)/;

    my $loc = $1;

    if ( $loc =~ /MISS/ ) {

      #print BED "chr1\t1\t2\t$info\_$id\t0\t+\n";
      next;
    }

    print "loc:$loc\n";
    my @ent = split( ":",  $loc );
    my @chr = split( /\./, $ent[0] );
    print "ent:" . join( "#", @ent ) . "#\n";
    print "chr:" . join( "#", @chr ) . "#\n";
    my @pos = split( "-", $ent[1] );
    my $strand = $ent[2];

    print BED $chr[1] . "\t" . $pos[0] . "\t" . $pos[1] . "\t$info\_$id\t0\t$strand\n";

    #print "$chr\t$start\t$end\t$info\_$id\t0\t+\n";

  }
  close(BED);

  GraphClust::makeCol($bed_file);

}
