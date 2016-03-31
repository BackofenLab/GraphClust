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
use lib $GraphClust::CONFIG{PATH_LOCARNA} . "/../lib/perl";
use MLocarna;    ## for reliabilities, experimental tree alignemnts with locP

my $in_cluster_dir;
my $in_eval_mode;
my $in_verbose = 0;
my $in_root_dir;
my $debug = 0;
## not used but we need it for generic call script for sge/threads
my $in_taskid = 0;

GetOptions(
  "cluster-dir=s" => \$in_cluster_dir,
  "evaluate"      => \$in_eval_mode,
  "verbose"       => \$in_verbose,
  "root-dir=s"    => \$in_root_dir,
  "task-id=i"     => \$in_taskid,
  "debug"         => \$debug,
);

################################################################################
## Load options from root-dir/config file
%GraphClust::CONFIG = readConfigFile( "$in_root_dir/config", 1 );
SUBSECTION("options loaded for alignCenter.pl");
printConfig( \%CONFIG ) if ($in_verbose);

my $locarna_path  = $CONFIG{PATH_LOCARNA};
my $rnaz_path     = $CONFIG{PATH_RNAZ};
my $vrna_path     = $CONFIG{PATH_VRNA};
my $cmfinder_path = $CONFIG{PATH_CMFINDER};
my $tmp_path      = $CONFIG{PATH_TMP};

my $OPTS_locarna_paligs  = $CONFIG{OPTS_locarna_paligs};
my $OPTS_locarna_maligs  = $CONFIG{OPTS_locarna_maligs};
my $OPTS_locarna_p_model = $CONFIG{OPTS_locarna_p_model};

my $center_subtree_min = $CONFIG{center_subtree_min};
my $center_subtree_max = $CONFIG{center_subtree_max};
my $center_tree_type   = $CONFIG{center_tree_type};
my $center_model_type  = $CONFIG{center_model_type};
my $center_tree_aligs  = $CONFIG{center_tree_aligs};
my $center_skip_unstable = $CONFIG{center_skip_unstable};

my $tree_aligs_local = 0;
$tree_aligs_local = 1 if ( $OPTS_locarna_paligs =~ /local/ || $OPTS_locarna_paligs =~ /normalized/ );
$center_model_type = 1 if ((!-e "$cmfinder_path"."cmfinder" || !-e "$cmfinder_path/cmfinder") || ($cmfinder_path eq "false" && $center_model_type == 5) );

my $bin_dir = abs_path($FindBin::Bin);

################################################################################

my $center_fa_file = "$in_cluster_dir/center.fa";
my $dp_dir         = "$in_cluster_dir/dp";

################################################################################
## model msa

## get real knn, can be less than $nspdk_knn_center!
## @{$ids_knn} is not sorted, is in kernel sim order!
my $ids_knn = GraphClust::readSubset( "$in_cluster_dir/center.ids", 1 );

## write names file, 1 id per line, sorted by id!!!
open( OUT, ">$in_cluster_dir/names" );
map { print OUT $_ . "\n"; } sort { $a <=> $b } @{$ids_knn};
close(OUT);

my $tree_dir = "$in_cluster_dir/TREE";
system("mkdir -p $tree_dir");


## make matrix file from kernel sim file (we always have the kernel sim in new GraphClust)
if ( -e "$in_cluster_dir/center.sim" ) {
  sim2matrix( "$in_cluster_dir/center.sim", $ids_knn, "$tree_dir/matrix.kernel" );
  matrix_blacklist_overlap_frags( "$tree_dir/matrix.kernel", "$in_cluster_dir/center.ids", "$in_root_dir/FASTA/data.map", "$tree_dir/matrix.kernel.ol", 0.1, 6 );
} elsif ( $center_tree_type == 3 ) {
  $center_tree_type = 1;
  print "\nWarning! Change center_tree_type to 1 as I can find no kernel sim matrix!\n\n";
}

## extract locarna score matrix as well
if ( $center_tree_type == 1 || $center_tree_type == 2 ) {

  locarnaAligs2matrix( "$in_cluster_dir/paligs", scalar( @{$ids_knn} ), $tree_aligs_local, $tree_dir );
  matrix_blacklist_overlap_frags( "$tree_dir/matrix.locarna", "$in_cluster_dir/center.ids", "$in_root_dir/FASTA/data.map", "$tree_dir/matrix.locarna.ol", 0.1, 0 );
  system("cp $tree_dir/matrix.locarna.ol $tree_dir/matrix.tree");

} elsif ( $center_tree_type == 3 ) {

  system("cp $tree_dir/matrix.kernel.ol $tree_dir/matrix.tree")

}

## create tree according to sim matrix (either from kernel or from locarana aligs)
print "\n...compute guide tree!\n";
my $tree_file = "$tree_dir/tree";
matrix2tree( "$tree_dir/matrix.tree", "$in_cluster_dir/names", $tree_dir, $tree_file );

die "No tree file found ($tree_file)! Exit...\n\n" if ( !-e $tree_file );

## evaluate trees and matrices
if ($in_eval_mode) {
  evalTreeMatrix( "$in_cluster_dir/EVAL_NEW", $tree_dir, $tree_file, "$in_cluster_dir/names" );
}

my $use_sim_stree = 0;
my $use_prob_alig = 0;
$use_prob_alig = 1 if ( $center_tree_aligs == 2 );

my $mloc_opts = $OPTS_locarna_maligs;
$mloc_opts = $OPTS_locarna_p_model if $use_prob_alig;

## prefolding seqs with config params, does not work for locarna-p
if ( !$use_prob_alig ) {
  my $rnafold_opts   = $CONFIG{OPTS_RNAfold};
  my $rnaplfold_opts = $CONFIG{OPTS_RNAfold};
  my $plfold_minlen  = $CONFIG{GLOBAL_plfold_minlen};

  my $fold_opts = "--vrna-path $vrna_path --fasta $center_fa_file --tgtdir $dp_dir --new-only ";
  $fold_opts .= "--switch-len $plfold_minlen ";
  $fold_opts .= "--rnafold-opts \"$rnafold_opts\" ";
  $fold_opts .= "--rnaplfold-opts \"$rnaplfold_opts\" ";

  system("perl $bin_dir/foldFasta.pl $fold_opts");
}

computeTreeAligsLocP( "$in_cluster_dir/maligs.loc_new", $center_fa_file, $tree_file, $center_subtree_min, $center_subtree_max, $use_sim_stree, $use_prob_alig, $mloc_opts, $dp_dir, "$in_cluster_dir/SUBTREES" );

GraphClust::SUBSECTION("calc scores for subtrees and rank them...");
my $cluster_id = "MISS";
if ( $in_cluster_dir =~ /\/(\d+)\.(\d+)\.cluster/ ) {
  $cluster_id = "$1.$2";
}

my $subtrees = rankSubtrees( "$in_cluster_dir/SUBTREES", "$tree_dir", $in_eval_mode );

my $best_subtrees_file = "$in_cluster_dir/bestSubtrees";
open( BEST, ">$best_subtrees_file" );

my $rank           = 0;
my @subtrees_final = ();

foreach my $t ( sort { $b->{SCORE_SIMPLE} <=> $a->{SCORE_SIMPLE} } @{$subtrees} ) {

  $rank++;
  print BEST getNodeInfo( $t, $cluster_id, $rank ) . "\n";

  ## no model if subtree contains overlapping sequences (frags)
  next if ( $t->{OVERLAP} );

  ## no model if we are in greylist mode and node does not contain greylist/non-gl pair
  next if ( !$t->{GL_HIT} && -e "$in_root_dir/FASTA/data.greylist" );

  ## skip unstructured subtrees
  next if ( $t->{MFE} >= 0 && $center_skip_unstable);
  #next if ( $t->{RNAZ}->[4] >= 0 && $center_skip_unstructured );

  ## otherwise potential model alignment
  push( @subtrees_final, $t );
}

close(BEST);

GraphClust::makeCol("$best_subtrees_file");

if ( !@subtrees_final ) {
  print "\nNo subtree/cluster with model criteria found in tree $in_cluster_dir/tree! No model will be build. Exit...\n\n";
  exit;
}

################################################################################
## build cmsearch model from best subtree

GraphClust::SUBSECTION( "Found " . scalar(@subtrees_final) . " models. Build model alignment now for best subtree..." );
my $model_dir = "$in_cluster_dir/MODEL";
system("mkdir -p $model_dir");

my @fa_center = GraphClust::read_fasta_file($center_fa_file);
my $subset_fa = GraphClust::writeSubsetFasta( \@fa_center, $subtrees_final[0]->{NAMES}, "$model_dir/model.tree.fa", 1 );

## default file name = $center_model_type == 1
my $tree_aln_file = $subtrees_final[0]->{FILE};
system("cp $tree_aln_file $model_dir/best_subtree.aln");
system("cp $tree_aln_file.ps $model_dir/best_subtree.aln.ps");

## select corresponding files for final model alig
## realign with locarna ($center_model_type=2) or locarna-p ($center_model_type=3|4|5)
if ( $center_model_type == 2 ) {
  ## realign with locarna and diff opts 'OPTS_locarna_model'
  print "\n...realign best subtree with locarna...\n";
  GraphClust::mlocarna_center( $center_fa_file, $model_dir, $dp_dir, 0 );
  $tree_aln_file = "$model_dir/results/result.aln";

} elsif ( $center_model_type == 3 || $center_model_type == 4 || $center_model_type == 5 ) {

  if ( $center_tree_aligs != 2 && !-e "$model_dir/results/result.aln" ) {
    ## realign with LocarnaP as subtrees were aligned with locarna (with $center_tree_aligs==1 we already have locarnaP alignments)
    print "\n...realign best subtree with locarnaP...\n";
    GraphClust::mlocarna_center( $subset_fa, $model_dir, $dp_dir, 1 );
    $tree_aln_file = "$model_dir/results/result.aln";
  }

}

## in case we have realigned best subtree we create the new alifold files (annotated colored postscript)
if ( !( -e "$tree_aln_file.ps" && -e "$tree_aln_file.alifold" ) ) {
  GraphClust::aln2alifold( $tree_aln_file, $tmp_path, $vrna_path );
}

## default is best subtree alignment
my $final_tree_prefix = "$model_dir/model.tree.aln";
system("cp $tree_aln_file $final_tree_prefix");
system("cp $tree_aln_file.ps $final_tree_prefix.ps");
system("cp $tree_aln_file.alirna.ps $final_tree_prefix.alirna.ps");
system("cp $tree_aln_file.alifold $final_tree_prefix.alifold");

## reliabillity signal
if ( -e "$tree_aln_file.rel_signal" ) {
  system("cp $tree_aln_file.rel_signal $final_tree_prefix.rel_signal");
  system("cp $tree_aln_file.rel_plot.pdf $final_tree_prefix.rel_plot.pdf")
   if ( -e "$tree_aln_file.rel_plot.pdf" ); ## in case R is not working correctly, plot is not there
}

## prepare command to make stk file from aln file
my $cmd_stk = "$bin_dir/mloc2stockholm.pl -file $final_tree_prefix";
$cmd_stk .= " -split_input yes -con_struct $final_tree_prefix.alifold";

## use predicted signal region if we have aligned with locarna-p and we want it
if ( -e "$final_tree_prefix.rel_signal" && ( $center_model_type == 4 || $center_model_type == 5 ) ) {
  $cmd_stk .= " -interval_only yes -rel_signal $final_tree_prefix.rel_signal";
}

## create stk file from selected (and realigned) subtree
## name is $model_dir/model.sth from script, rename to stk
system($cmd_stk);
system("mv $final_tree_prefix.sth $model_dir/model.tree.stk");

## CMFINDER
if ( $center_model_type == 5 ) {

  ## number of additional frags is set in MASTER, currently 3*$nspdk_knn_center
  my $ids_ext = GraphClust::readSubset( "$in_cluster_dir/center.ids.ext", 1 );
  my @fa_ext_all = GraphClust::read_fasta_file("$in_cluster_dir/center.fa.ext");
  my $fa_ext_merged = GraphClust::mergeFrags( \@fa_ext_all );
  GraphClust::writeSubsetFasta( $fa_ext_merged, $fa_ext_merged->[1], "$model_dir/cmfinder.fa", 1 );
  system("rm -f $model_dir/model.cmfinder.stk");
  system_call( $cmfinder_path . "cmfinder --g 1.0 -a $model_dir/model.tree.stk $model_dir/cmfinder.fa $model_dir/model.cmfinder.cm > $model_dir/model.cmfinder.stk", 1 );
  system("rm -f $model_dir/model.cmfinder.cm");

  ## can be empty in case cmfinder does not find anything
  if ( -e "$model_dir/model.cmfinder.stk" && !-z "$model_dir/model.cmfinder.stk" ) {
    open( my $STK_IN, "$model_dir/model.cmfinder.stk" );
    my $model_stk = GraphClust::read_stockholm($STK_IN);
    close($STK_IN);
    open( my $ALN_OUT, ">$model_dir/model.cmfinder.aln" );
    GraphClust::write_clustal( $ALN_OUT, $model_stk );
    close($ALN_OUT);
    GraphClust::aln2alifold( "$model_dir/model.cmfinder.aln", $tmp_path, $vrna_path );

    open( my $ALI, "$model_dir/model.cmfinder.aln.alifold" );
    my @cmfinder_energy = <$ALI>;
    chomp(@cmfinder_energy);
    close($ALI);

    my $cmf_energy = 0;
    if ( $cmfinder_energy[2] =~ /\(\s*(\S+)\s+=\s+(\S+)\s+\+\s+(\S+)\)/ ) {
      $cmf_energy = $1;
    }

    ## check if cmfinder model has still any structure, otherwise choose initial model from tree
    if ( $cmf_energy < 0 || !$center_skip_unstable) {
      $final_tree_prefix = "$model_dir/model.cmfinder.aln";
      print "CMfinder model energy: $cmf_energy\n";
    }
  }
}

## copy files for final model stk from tree files
system("cp $final_tree_prefix $model_dir/model.aln");
system("cp $final_tree_prefix.ps $model_dir/model.aln.ps");
system("cp $final_tree_prefix.alirna.ps $model_dir/model.aln.alirna.ps");

$final_tree_prefix =~ s/\.aln//;
system("cp $final_tree_prefix.stk $model_dir/model.stk");


## just as info we keep this file
system("cp $in_cluster_dir/SGE_log_MSA/cmd.opts $in_cluster_dir/stage7.cmd");
## cleanup files
if ( !$debug ) {

  system("rm $tree_dir/*.ol");
  system("rm $tree_dir/*-list");
  system("rm -r $in_cluster_dir/SUBTREES");
  system("rm -r $in_cluster_dir/maligs.loc_new");
  system("rm -r -f $in_cluster_dir/dp");
  system("rm -r -f $in_cluster_dir/pp");
  system("rm -r -f $model_dir/results");
}

exit;

################################################################################
################################################################################
## subs

sub rankSubtrees {
  my $subtreeDir = $_[0];
  my $treeDir    = $_[1];
  my $evaluate   = $_[2];

  my @subtrees = ();

  open( TREE, "$treeDir/tree" );
  my $tree_str = <TREE>;
  close(TREE);
  chomp($tree_str);

  ## get sim-tree structure, $tree is hash-ref, keys are nodeIDs in postorder from 0
  my ( $nodes, $tree ) = GraphClust::newick_tree_to_postorder2($tree_str);
  my ( $childs2, $nodes2 ) = GraphClust::getNodeLeafs( $tree, 0 );

  foreach my $file ( glob("$subtreeDir/*.aln") ) {

    my %node = ();
    my $aln  = GraphClust::read_clustalw_alnloh($file);

    $node{FILE}      = $file;
    $node{BRANCHSUM} = 0;
    $node{BRANCH}    = 0;
    $node{DENSITY}   = 0;
    $node{NAMES}     = [ map { $_->{name} } @{$aln} ];
    $node{NODEID}    = 0;
    $node{NLEAFS}    = @{$aln};

    ## fill node with scores
    NodeScore( \%node, $aln, $childs2, "$treeDir/matrix.tree" );

    if ($evaluate) {
      ( $node{CLASS}, $node{QUAL}, $node{QUAL_ABS} ) = evalNode( \%node );
    } else {
      $node{CLASS}    = $node{NAMES};
      $node{QUAL}     = "0.0";
      $node{QUAL_ABS} = "0.0";
    }

    push( @subtrees, \%node );

  }
  return \@subtrees;
}

sub NodeScore {
  my $node             = $_[0];
  my $node_aln         = $_[1];
  my $tree_ids         = $_[2];
  my $tree_matrix_file = $_[3];

  $node->{SCI} = 0;
  $node->{COV_EI} = 1.0;

  ## alfold consensus MFE
  $node->{MFE} = 0;
  if ( -e $node->{FILE} . ".alifold" ) {
    my @alifold = readpipe( "cat " . $node->{FILE} . ".alifold" );
    $alifold[2] =~ /\(\s*(\S+)\s*=/;
    $node->{MFE} = $1;
  }

  ## MPI
  my $mpi = 1;

  my @getmpi = readpipe("perl $bin_dir/scoreAln.pl -i $node->{FILE} -f CLUSTALW -s mpi");
  if (@getmpi) {
    $mpi = $getmpi[0];
    chomp($mpi) if ($mpi);
  }
  $node->{MPI} = $mpi;

  ## REL
  $node->{REL} = [ 1, 1, 1, 1 ];
  if ( -e $node->{FILE} . ".rel_signal" ) {

    open( IN, $node->{FILE} . ".rel_signal" );
    my @rel_dat = <IN>;
    close(IN);

    ## line $rel_dat[3] is: 'REL_SCORES 48.92:38.57:53.42:42.12' -> split 2nd filed into array-ref
    $node->{REL} = [ split( ":", ( split( " ", $rel_dat[3] ) )[1] ) ];
  }

  ## RNAz
  $node->{RNAZ} = [ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 ];
  my $use_rnaz = 1;
  if ($use_rnaz) {
    $node->{RNAZ} = GraphClust::getRNAzScores( $node->{FILE}, "$rnaz_path/RNAz", "-n" );
#    if ( $node->{MFE} >= 0 && $node->{RNAZ}->[9] > 0 ) {
#      $node->{MFE} = $node->{RNAZ}->[4];
#      $node->{SCI} = $node->{RNAZ}->[9];
#    } elsif ( $node->{RNAZ}->[3] < 0 && $node->{RNAZ}->[4] < 0 ) {
#      $node->{SCI} = abs( $node->{RNAZ}->[4] / $node->{RNAZ}->[3] );
#    } else {
#      $node->{SCI} = 0.1;
#    }

    ## overwrite alifold mfe with RNAz consensus mfe if we have strange values
    if ( $node->{MFE} >= 0 && $node->{RNAZ}->[4] < 0 ){
      $node->{MFE} =  $node->{RNAZ}->[4];
    }

    if ( $node->{RNAZ}->[3] < 0 && $node->{MFE} < 0 && $node->{RNAZ}->[6] < 0 && $node->{RNAZ}->[5] < 0 ) {
      $node->{SCI} = abs( $node->{MFE} / $node->{RNAZ}->[3] );
    } else {
      $node->{SCI} = 0.01;
      $node->{MFE} = 0.0;
    }

    if ($node->{RNAZ}->[6] < 0){
      ## new 0.4
      $node->{COV_EI} = (abs($node->{RNAZ}->[6])) ** 0.5;
      ## new 0.7
      $node->{COV_EI} = 1 if ($node->{COV_EI} < 1);
    }

  }

  ## score
  $node->{SCORE}     = 0;
  $node->{SCORE_OLD} = 0;
  $node->{CENTERC}   = 0;

  $node->{LEN_MEAN} = 0;
  $node->{LEN_SD}   = 0;
  $node->{LEN_ALN}  = length( $node_aln->[0]->{seq} );
  $node->{SVM}      = 0;

  ## OVERLAP
  $node->{OVERLAP} = check_node_overlap( $node->{NAMES}, "$in_root_dir/FASTA/data.map", 0.1 );

  ## Greylist
  $node->{GL_HIT} = 0;
  if ( -e "$in_root_dir/FASTA/data.greylist" ) {
    my $node_gl_hit = check_node_greylist( $node->{NAMES}, "$in_root_dir/FASTA/data.greylist" );
    $node->{GL_HIT}  = $node_gl_hit;
    $node->{OVERLAP} = 0;
  }

  ## score simple is just average of score submatrix for current subtree
  my %leaf2idx = ();    ## stores score-matrix row for each leaf-id

  @{$tree_ids} = sort { $a <=> $b } @{$tree_ids};

  map { $leaf2idx{ $tree_ids->[$_] } = $_ } 0 .. $#{$tree_ids};
  my @leaf_idx = map { $leaf2idx{$_} if ( exists $leaf2idx{$_} ); } @{ $node->{NAMES} };

  my $sc = GraphClust::getMatrixSum( $tree_matrix_file, \@leaf_idx );
  $node->{SCORE} = $sc;
  print "\nmatrix avg = $sc " . ( $sc / ( $node->{NLEAFS} - 1 ) ) . "\n";
  print $node->{NAMES} . " " . $node->{FILE} . "\n";

  ## get locarna subtree score
  $node->{SCORE_LOC} = 1;
  if ( -e $node->{FILE} . ".matrix" ) {
    %leaf2idx = ();
    my $sub_ids = GraphClust::readSubset( $node->{FILE} . ".subtree_ids", 1 );
    @{$sub_ids} = sort { $a <=> $b } @{$sub_ids};
    map { $leaf2idx{ $sub_ids->[$_] } = $_ } 0 .. $#{$sub_ids};
    @leaf_idx = ();
    @leaf_idx = map { $leaf2idx{$_} if ( exists $leaf2idx{$_} ); } @{ $node->{NAMES} };
    my $sc_loc = GraphClust::getMatrixSum( $node->{FILE} . ".matrix", \@leaf_idx );

    #$node->{SCORE_LOC} =  ($sc_loc)/( $node->{NLEAFS}-1);
    #$node->{SCORE_LOC} =  $sc_loc/($center_subtree_max-@{ $node->{NAMES} }+1);

    if ( -e $node->{FILE} . ".rel_signal" ) {
     #$sc_loc = $sc_loc / ( $node->{NLEAFS} - 1 );
     $sc_loc = $sc_loc ;
    }
    $node->{SCORE_LOC} = $sc_loc;
    print "score locarna matrix = $sc_loc \n";
  }

  #print "idx:".join(":",@leaf_idx)."\n";
  #print "matrix_score:".$sc."\n";

  ## FINAL SCORE
  if ( $center_tree_type == 3 ) {
    ## tree based on kernel sim matrix
    $sc = $sc / ( $node->{NLEAFS} - 1 );
  } else {
    $sc = $sc / ( $node->{NLEAFS} - 1 );
  }
  $node->{SIZE_BAL} = 1 / ( ( $node->{MPI} * 100 )**( 1 / $node->{NLEAFS} ) );
  ## full_05_new_04
  $node->{SCORE_SIMPLE} = $node->{SCORE_LOC} * $node->{SIZE_BAL} * $node->{SCI}**0.7 * $node->{COV_EI};
  ##new
  #$node->{SCORE_SIMPLE} =  $node->{SCORE_LOC} * $node->{SIZE_BAL} * $node->{SCI} * $node->{SIZE_BAL} * $node->{COV_EI};

}

sub evalNode {
  my $node = $_[0];

  my ( $class_Frags, $class_Size, $class_Names, $class_data );
  ( $class_Frags, $class_Size, $class_Names, $class_data ) = GraphClust::getEvalClassInfo("$in_root_dir/FASTA");

  ## classes for evaluate
  my @leafs = @{ $node->{NAMES} };

  my @cl = ();
  map {
    if ( exists $class_data->{$_} ) {
      push( @cl, $class_data->{$_} )
    } else {
      push( @cl, 0 );
    }
  } @leafs;

  my ( $quals_n, $quals_a ) = GraphClust::subtreeQual( \@cl, $center_subtree_min, $center_subtree_max );

  return ( \@cl, sprintf( "%.2f", $quals_n->[0] ), sprintf( "%.2f", $quals_a->[0] ) );
}

sub check_node_overlap {
  my $node_ids    = $_[0];   ## leaf ids in one node
  my $id_map_file = $_[1];   ## usually data.map file line: "122 SEQ23#20#100#+"
  my $black_overlap = $_[2]; ## minimal required overlap of two frags to become blacklisted

  ## map id to frag
  my $frags   = GraphClust::read_fragments($id_map_file);
  my %id2frag = ();
  map { $id2frag{ $_->{VALUE} } = $_ } @{$frags};

  ## matrix frags
  my $node_frags = [];
  map { push( @{$node_frags}, $id2frag{$_} ) } @{$node_ids};

  ## matrix frag ols
  @{$node_frags} = sort { $a->{SEQID} cmp $b->{SEQID} } @{$node_frags};
  my $matrix_ols = GraphClust::fragment_overlap( $node_frags, $node_frags, $black_overlap, 1 );

  my $node_overlap = 0;
  $node_overlap = 1 if ( @{$matrix_ols} );

  return $node_overlap;
}

sub check_node_greylist {
  my $node_ids = $_[0];    ## leaf ids in one node
  my $gl_file  = $_[1];

  my %gl = ();
  open( GL, $gl_file );
  map { my $key = $_; chomp($key); $gl{$key} = 1 } <GL>;
  close(GL);

  my $hits    = 0;
  my $gl_hits = 0;
  map { $hits++ if ( !exists $gl{$_} ); $gl_hits++ if ( exists $gl{$_} ) } @{$node_ids};

  if ( $hits > 0 && $gl_hits > 0 ) {
    return 1;
  } else {
    return 0;
  }
}

sub getNodeInfo {
  my $node = $_[0];
  my $nidx = $_[1];
  my $rank = $_[2];

  my $info_str = "";

  $info_str .= "cent $nidx rank $rank";
  $info_str .= " score " . sprintf( "%.3f", $node->{SCORE} );
  $info_str .= " cent " . $node->{CENTERC};

  $info_str .= " SVM " . $node->{SVM};

  $info_str .= " class " . join( ":", @{ $node->{CLASS} } );
  $info_str .= " names " . join( ":", @{ $node->{NAMES} } );

  $info_str .= " sci " . sprintf( "%.2f", $node->{SCI} );
  $info_str .= " MPI " . $node->{MPI};
  $info_str .= " REL " . join( ":", @{ $node->{REL} } );
  $info_str .= " mfe " . $node->{MFE};
  $info_str .= " dens " . $node->{DENSITY};
  $info_str .= " BSum " . $node->{BRANCHSUM};
  $info_str .= " nodeID " . $node->{NODEID};
  $info_str .= " nleafs " . $node->{NLEAFS};

  my $tmp = $1 if ( $node->{FILE} =~ /intermediate_(\d+\.aln)/ || $node->{FILE} =~ /intermediate_(\d+_\d+\.aln)/ );
  $info_str .= " file " . $1;
  $info_str .= " RNAZ " . join( ":", @{ $node->{RNAZ} } );
  $info_str .= " LEN_MEAN " . $node->{LEN_MEAN};
  $info_str .= " LEN_SD " . $node->{LEN_SD};
  $info_str .= " LEN_ALN " . $node->{LEN_ALN};
  $info_str .= " QUAL " . $node->{QUAL};
  $info_str .= " QUAL_ABS " . $node->{QUAL_ABS};
  $info_str .= " SC_SIMPLE " . sprintf( "%.2f", $node->{SCORE_SIMPLE} );
  $info_str .= " OL " . $node->{OVERLAP};
  $info_str .= " GL " . $node->{GL_HIT};
  $info_str .= " LOC " . sprintf( "%.2f", $node->{SCORE_LOC} );
  $info_str .= " COV_EI " . sprintf( "%.2f", $node->{COV_EI} );
  $info_str .= " SIZE_BAL " . sprintf( "%.2f", $node->{SIZE_BAL} );

  return $info_str;

}

sub evalTreeMatrix {
  my $evalDir   = $_[0];
  my $treeDir   = $_[1];
  my $treeFile  = $_[2];
  my $namesFile = $_[3];

  mkdir($evalDir);
  my $tmpDir = "$tmp_path/tmp_tree_$$";
  mkdir($tmpDir);

  ## locarna eval
  if ( -e "$treeDir/matrix.locarna.ol" ) {
    system("perl $bin_dir/matrixClass.pl $in_root_dir/FASTA/class.hash $namesFile $treeDir/matrix.locarna.ol 0 | column -t > $evalDir/matrix.locarna.ol.class");
    matrix2tree( "$treeDir/matrix.locarna.ol", $namesFile, $tmpDir, "$evalDir/matrix.locarna.ol.tree" );
    system("perl $bin_dir/treeClass.pl $in_root_dir/FASTA/class.hash $evalDir/matrix.locarna.ol.tree > $evalDir/matrix.locarna.ol.tree.class");
  }

  ## kernal eval
  system("perl $bin_dir/matrixClass.pl $in_root_dir/FASTA/class.hash $namesFile $treeDir/matrix.kernel.ol 6 | column -t > $evalDir/matrix.kernel.ol.class");
  matrix2tree( "$treeDir/matrix.kernel.ol", $namesFile, $tmpDir, "$evalDir/matrix.kernel.ol.tree" );
  system("perl $bin_dir/treeClass.pl $in_root_dir/FASTA/class.hash $evalDir/matrix.kernel.ol.tree > $evalDir/matrix.kernel.ol.tree.class");

  ## final tree, similar to either locarna or kernel tree
  system("perl $bin_dir/treeClass.pl $in_root_dir/FASTA/class.hash $treeFile > $evalDir/tree.class");
  system("rm -r -f $tmpDir");
}

sub matrix2tree {
  my $matrix_file  = $_[0];
  my $names_file   = $_[1];
  my $tree_dir     = $_[2];
  my $tree_outfile = $_[3];

  mkdir($tree_dir);
  system("cat $matrix_file | awk '{for(i=1;i<NR;i++){print NR,i,\$(i)}}' > $tree_dir/tree.score-list");
  system("perl $bin_dir/rnaclustScores2Dist.pl --quantile 1.0 < $tree_dir/tree.score-list > $tree_dir/tree.dist-list");
  system("$bin_dir/pgma $names_file $tree_dir/tree.dist-list > $tree_outfile");
}

sub locarnaAligs2matrix {
  my $paligsDir  = $_[0];
  my $num_seqs   = $_[1];
  my $localAligs = $_[2];
  my $outDir     = $_[3];

  if ( !-e "$outDir/matrix.locarna" ) {

    if ( !$localAligs ) {

      system_call( "perl $bin_dir/getAdditionalScores.pl -paligs $paligsDir -n $num_seqs -outc $outDir/matrix.centers -outw $outDir/matrix.width --outl $outDir/matrix.locarna",
        $in_verbose );

    } else {

      system_call(
"perl $bin_dir/getAdditionalScoresLocal.pl -paligs $paligsDir -n $num_seqs -s $outDir/center.fa -outc $outDir/matrix.centers -outw $outDir/matrix.width --outl $outDir/matrix.locarna",
        $in_verbose
      );

    }
  } else {
    print "\nUse existing matrix file $outDir/matrix.locarna!\n\n";
  }

  die "\n\nError! Cannot find $outDir/matrix.locarna! Check run of getAdditionalScores(Local).pl! Exit...\n\n" if ( !-e "$outDir/matrix.locarna" );
}

sub sim2matrix {
  my $sim_file = $_[0];
  my $ids      = $_[1];    ## array ref
  my $out_file = $_[2];

  open( IN, "$sim_file" );
  my $t = <IN>;
  close(IN);
  chomp($t);
  my @ent = split( " ", $t );

  my $knn = @{$ids}; ## knn x knn matrix to build, sim file contains kernel similarities
  die "Sim matrix line has incorrect number of entries! Expected " . ( $knn * ( $knn - 1 ) / 2 ) . " found " . @ent . " Exit...\n\n"
    if ( @ent != $knn * ( $knn - 1 ) / 2 );

  ## fix order of frag ids
  my @new_order = sort { $a <=> $b } @{$ids};

  ## map frag-id to matrix col, start with col 0
  my %id2col = ();
  map { $id2col{ $new_order[ $_ - 1 ] } = $_ - 1 } 1 .. @new_order;

  ## final matrix
  my @matrix = ();
  foreach my $e (@ent) {
    my @vals = split( ":", $e );
    $matrix[ $id2col{ $vals[0] } ][ $id2col{ $vals[1] } ] = $vals[2];
    $matrix[ $id2col{ $vals[1] } ][ $id2col{ $vals[0] } ] = $vals[2];
  }

  ## diagonal
  map { $matrix[ $_ - 1 ][ $_ - 1 ] = 0 } 1 .. $knn;

  ## write out
  open( OUT, ">$out_file" );
  map { print OUT join( " ", @{$_} ) . "\n" } @matrix;
  close(OUT);
}

## min_seqs : further evaluate subtrees of sim-tree with >$min_seqs leafs
## max_seqs: further evaluate subtrees of sim-tree with <$max_seqs leafs
## use_sim_tree: uses tree structure from sim tree or not for realigning subtree
sub computeTreeAligsLocP {
  my $tgtDir            = $_[0];
  my $fa_tree_file      = $_[1];
  my $tree_file         = $_[2];
  my $min_seqs          = $_[3];
  my $max_seqs          = $_[4];
  my $use_sim_tree      = $_[5];
  my $use_probabilistic = $_[6];
  my $mloc_Opts         = $_[7];
  my $dpDir             = $_[8];
  my $resDir            = $_[9];

  my $new_only = 1;    ## compute only new aligs

  mkdir($tgtDir);

  open( TREE, "$tree_file" );
  my $tree_str = <TREE>;
  close(TREE);
  chomp($tree_str);

  ## get sim-tree structure, $tree is hash-ref, keys are nodeIDs in postorder from 0
  ## nodes: postordered list of node-labels (leaf-names & '$$nodesym')
  my ( $nodes, $tree ) = GraphClust::newick_tree_to_postorder2($tree_str);
  my ( $childs2, $nodes2 ) = GraphClust::getNodeLeafs( $tree, 0 );
  my $subtrees_simtree_href = getSubtrees( $tree, $min_seqs, $max_seqs ); ## key == value for max subtrees
  print "main tree nodes:" . join( ":", @{$nodes} ) . "\n";
  print "main tree nodes:" . join( ":", @{$nodes2} ) . "\n";

  my %imMap = ();
  my $imIdx = 1;
  foreach my $idx ( 1 .. @{$nodes2} ) {
    next if $tree->{ $nodes2->[ $idx - 1 ] }->{LEAF};
    $imMap{ $nodes2->[ $idx - 1 ] } = $imIdx;
    $imIdx++;
  }

  print "\nmap node-id (preorder) -> intermediate-id\n";
  foreach my $key ( sort keys %imMap ) {
    print "imMap $key -> " . $imMap{$key} . "\n";
  }
  mkdir("$tgtDir/log");
  mkdir($resDir);
  my @fa_tree = GraphClust::read_fasta_file($fa_tree_file);

  foreach my $id ( sort { $a <=> $b } keys %{$subtrees_simtree_href} ) {

    ## continue only for maximal trees
    next if ( $subtrees_simtree_href->{$id} != $id );

    print "\ntree $id used " . $subtrees_simtree_href->{$id} . "\n";

    my $subtree_dir = "$tgtDir/$id";
    mkdir("$tgtDir/$id");

    my ( $childs, $nodes ) = GraphClust::getNodeLeafs( $tree, $id );

    my $subset_fa = GraphClust::writeSubsetFasta( \@fa_tree, $childs, "$tgtDir/$id/$id.fa", 1 );

    if ($use_sim_tree) {
      print "use sim-tree structure for subtree node $id\n";
      my $subtree_str = subtree2Newick( $tree, $id, 1 );
      print "sim-tree for subtree node $id: $subtree_str\n";
      open( TREE, ">$tgtDir/$id/$id.subtree" );
      print TREE $subtree_str . ";";
      close(TREE);
    }

    ## reuse dot plots, does not work locarnaP currently!
    if ( !$use_probabilistic ) {
      my $loc_pp_dir = "$subtree_dir/input";
      mkdir($loc_pp_dir);
      foreach my $key ( @{$childs} ) {
        system("ln -f -s $dpDir/$key $loc_pp_dir/$key") if ( -e "$dpDir/$key" );
      }
    }

    my $call = "$locarna_path/mlocarna --tgtdir $subtree_dir ";
    $call .= "--treefile $tgtDir/$id/$id.subtree " if ($use_sim_tree);
    $call .= "$mloc_Opts --verbose ";

    if ($use_probabilistic) {
      $call .= "--probabilistic " if ( $call !~ /--probabilistic/ );
      system_call( $call . " $tgtDir/$id/$id.fa 1>$tgtDir/log/$id.locarna_call 2>$tgtDir/log/$id.locarna_call_err", $in_verbose )
        if ( !-e "$tgtDir/$id/results/result.aln" || !$new_only );
    } else {
      system_call( $call . " $tgtDir/$id/$id.fa 1>$tgtDir/log/$id.locarna_call 2>$tgtDir/log/$id.locarna_call_err", $in_verbose )
        if ( !-e "$tgtDir/$id/results/result.aln" || !$new_only );
    }

    print "eval subtree $id from simtree\n";
    open( TREE, "$tgtDir/$id/results/result.tree" );
    my $eval_tree_str = <TREE>;
    chomp($eval_tree_str);
    my ( $nodeList, $evalTree ) = GraphClust::newick_tree_to_postorder2($eval_tree_str);
    my ( $leafs_subtree, $nodes_subtree ) = getNodeLeafs( $evalTree, 0 );

    print "here eval subtree $id from simtree nodes " . join( ":", @{$nodes} ) . "\n";

    my $iidx = 0;

    ## go through nodes of subtree in postorder
    foreach my $node ( @{$nodes_subtree} ) {

      #print "node $node\n";
      ## ignore leafs
      next if ( $evalTree->{$node}->{LEAF} );
      $iidx++ if ( !$evalTree->{$node}->{LEAF} );

      ## ignore subtrees out of range, basically only too small ones
      next if ( $evalTree->{$node}->{NLEAFS} < $min_seqs || $evalTree->{$node}->{NLEAFS} > $max_seqs );

      ## analyse subtree
      #      my ( $leafsH, $nodesH ) = getNodeLeafs( $evalTree, $node );
      my ( $childs3, $nodes3 ) = GraphClust::getNodeLeafs( $tree, $node + $id );
      print "eval node $node (" . $id . ") intermediate $iidx orig-node:" . ( $node + $id ) . " leafs " . join( ":", @{$childs3} ) . "\n";

      ## intermediate name: either original intermediate idx if same subtree structure is used
      ## or if new structure is used: "intermediate_rootIntermID_idx.aln"
      my $im_filename;
      if ($use_sim_tree) {
        $im_filename = "intermediate_" . $imMap{ $node + $id } . ".aln";
      } else {
        $im_filename = "intermediate_$id\_$iidx.aln";
      }
      print "node=$node id=$id iidx:$iidx  im name: $im_filename\n";
      system("cp $tgtDir/$id/intermediates/intermediate$iidx.aln $resDir/$im_filename");

      ## write out parent leafs for matrix file
      open( PL, ">$resDir/$im_filename.subtree_ids" );
      print PL join( " ", sort @{$childs} ) . "\n";
      close(PL);

      ## probabilistic analysis
      if ($use_probabilistic) {
        MLocarna::forget_normalized_names();
        my @names_loh = ();
        map { my %t = (); $t{name} = $_; push( @names_loh, \%t ); } @$childs;
        MLocarna::register_normalized_seqnames( \@names_loh );
        locarnaP_aln2bmrels( "$tgtDir/$id/intermediates/intermediate$iidx.aln", "$tgtDir/$id/probs", 0, "$tgtDir/$id/intermediates" );
        my @rel_scores;
        my $call_eval = readpipe("$locarna_path/mlocarna --evaluate $resDir/$im_filename --tgtdir $tgtDir/$id $tgtDir/$id.fa");
        $call_eval =~ /RELIABILITY 1\/COL\s+(\S+)\%/;
        push( @rel_scores, $1 );
        $call_eval =~ /RELIABILITY 2\/COL\s+(\S+)\%/;
        push( @rel_scores, $1 );
        $call_eval =~ /RELIABILITY 1\/CCOL\s+(\S+)\%/;
        push( @rel_scores, $1 );
        $call_eval =~ /RELIABILITY 2\/CCOL\s+(\S+)\%/;
        push( @rel_scores, $1 );
        $call_eval =~ /MAX REL\. STRUCT\.\s+(\S+)/;
        my $rel_str = $1;
        open( LOG, ">$tgtDir/log/$im_filename.eval_out" );
        print LOG $call_eval;
        close(LOG);

        my @rel = readpipe("perl $locarna_path/reliability-profile.pl -v -fit-once-on --structure-weight=1  --fit-penalty 0.01 --beta 200 --show-sw --out=$resDir/$im_filename.rel_plot.pdf $tgtDir/$id/intermediates/intermediate$iidx.aln $tgtDir/$id/intermediates/intermediate$iidx.aln.bmreliability");
        open( REL, ">$resDir/$im_filename.rel_signal" );
        print REL @rel[ 0 .. 2 ];
        print REL "REL_SCORES " . join( ":", @rel_scores ) . "\n";
        print REL "REL_STRUCT " . $rel_str . "\n";
        close(REL);
        system("cp $tgtDir/$id/intermediates/intermediate$iidx.aln.bmreliability $resDir/$im_filename.bm_rel");
      }

      if ( -e "$tgtDir/$id/results/result.matrix" ) {
        system("cp $tgtDir/$id/results/result.matrix $resDir/$im_filename.matrix");
      }

      ##ali fold
      my $tmp_dir = "$tgtDir/alifold_$$";
      my $currDir = getcwd;
      mkdir($tmp_dir);
      chdir($tmp_dir);
      my @call_alifold = readpipe("$vrna_path/RNAalifold -r --noLP --color --aln $tgtDir/$id/intermediates/intermediate$iidx.aln 2>/dev/null");
      my $aln_cons_struct = $call_alifold[1];    ## store cons-struct
      chomp($aln_cons_struct);
      $aln_cons_struct =~ /([\(\).]*)(.*)/;      ## extract cons-struct
      $aln_cons_struct = $1;
      open( CONS, ">$resDir/$im_filename.alifold" );
      print CONS $call_alifold[0];
      print CONS "$aln_cons_struct\n";

      print CONS $2;
      close(CONS);
      system("mv alirna.ps $resDir/$im_filename.alirna.ps");
      system("mv aln.ps $resDir/$im_filename.ps");
      chdir($currDir);
      system("rm -R $tmp_dir");
    }

  }
}

sub getSubtrees {
  my $tree     = $_[0];
  my $min_seqs = $_[1];
  my $max_seqs = $_[2];

  my %used_IDs = ();

  foreach my $nodeID ( sort { $a <=> $b } keys %{$tree} ) {

    my ( $childs1, $nodes1 ) = GraphClust::getNodeLeafs( $tree, $nodeID );
    print "check node $nodeID leaf:" . $tree->{$nodeID}->{LEAF} . " nleafs:" . $tree->{$nodeID}->{NLEAFS} . " childs:" . join( ":", @{$childs1} ) . "\n";
    next if ( exists $used_IDs{$nodeID} );

    next if ( $tree->{$nodeID}->{LEAF}
      || $tree->{$nodeID}->{NLEAFS} < $min_seqs
      || $tree->{$nodeID}->{NLEAFS} > $max_seqs );

    print "use Node $nodeID\n";
    my ( $childs, $nodes ) = GraphClust::getNodeLeafs( $tree, $nodeID );

    $used_IDs{$nodeID} = $nodeID;

    foreach my $node ( @{$nodes} ) {
      $used_IDs{$node} = $nodeID if ( !$tree->{$node}->{LEAF} && !exists $used_IDs{$node} && $tree->{$node}->{NLEAFS} >= $min_seqs && $tree->{$node}->{NLEAFS} <= $max_seqs );
    }

    print "node $nodeID leafs " . join( ":", @{$childs} ) . " subnodes " . join( ":", @{$nodes} ) . "\n";

  }

  return \%used_IDs;
}

sub subtree2Newick {
  my $tree   = $_[0];
  my $nodeID = $_[1];
  my $root   = $_[2];

  my $str = "";

  if ( !$tree->{$nodeID}->{LEAF} ) {
    my @childs = @{ $tree->{$nodeID}->{CHILDS} };

    if ( $tree->{ $childs[0] }->{LEAF} ) {
      $str .= $tree->{ $childs[0] }->{NAME};
    } else {
      $str .= "(" . subtree2Newick( $tree, $childs[0] ) . ")";
    }

    if ( $tree->{ $childs[1] }->{LEAF} ) {
      $str .= "," . $tree->{ $childs[1] }->{NAME};
    } else {
      $str .= ",(" . subtree2Newick( $tree, $childs[1] ) . ")";
    }
  }
  $str = "(" . $str . ")" if ($root);

  return $str;
}

sub matrix_blacklist_overlap_frags {
  my $matrix_file    = $_[0];   ## one row per line, cols space seperated
  my $matrix_id_file = $_[1];   ## space-seperated file with all ids in one row!
  my $id_map_file = $_[2];   ## usually data.map file line: "122 SEQ23#20#100#+"
  my $matrix_out_file = $_[3];    ## filename for output matrix
  my $black_overlap = $_[4]; ## minimal required overlap of two frags to become blacklisted
  my $prec = $_[5];          ## precision for all values in matrix

  ## map id to frag
  my $frags   = GraphClust::read_fragments($id_map_file);
  my %id2frag = ();
  map { $id2frag{ $_->{VALUE} } = $_ } @{$frags};

  ## matrix ids
  open( IDS, $matrix_id_file );
  my @ids_matrix = <IDS>;
  close(IDS);
  chomp( $ids_matrix[0] );
  @ids_matrix = split( " ", $ids_matrix[0] );
  @ids_matrix = sort { $a <=> $b } @ids_matrix;

  ## matrix frags
  my $matrix_frags = [];
  map { push( @{$matrix_frags}, $id2frag{$_} ) } @ids_matrix;

  ## matrix frag ols
  @{$matrix_frags} = sort { $a->{SEQID} cmp $b->{SEQID} } @{$matrix_frags};
  my $matrix_ols = GraphClust::fragment_overlap( $matrix_frags, $matrix_frags, $black_overlap, 1 );

  ## store all overlaps in hash
  my %idOls = ();
  foreach my $ol ( @{$matrix_ols} ) {
    $idOls{ $matrix_frags->[ $ol->[0] ]->{VALUE} . "#" . $matrix_frags->[ $ol->[1] ]->{VALUE} } = $ol->[2];
  }

  ## get min of matrix
  open( MAT, $matrix_file );
  my $matrix_min = 10**10;
  while ( my $line = <MAT> ) {
    chomp $line;
    my @ent = split( " ", $line );

    map { $matrix_min = $_ if ( $_ < $matrix_min && $_ != 0 ) } @ent;
  }
  close(MAT);
  $matrix_min = 0 if ( $matrix_min > 0 && $matrix_min < 1 );
  print "matrix min $matrix_min\n";

  my $matrix_blacklist_value = $matrix_min - abs( $matrix_min / 2 );

  ## create new matrix
  my @newMAT = ();
  open( MAT, $matrix_file );

  my $row = 0;
  while ( my $line = <MAT> ) {

    chomp $line;
    my @ent = split( " ", $line );

    foreach my $col ( 0 .. $#ent ) {

      ## check if for current (row,col) exists overlap
      if ( exists $idOls{ $ids_matrix[$row] . "#" . $ids_matrix[$col] } || exists $idOls{ $ids_matrix[$col] . "#" . $ids_matrix[$row] } ) {
        $ent[$col] = $matrix_blacklist_value;
      } else {
        $ent[$col] = sprintf( "%1." . $prec . "f", $ent[$col] );
      }

    }
    push( @newMAT, join( " ", @ent ) );
    $row++;
  }
  close(MAT);

  ## write new matrix
  open( OUT, ">$matrix_out_file" );
  map { print OUT $_ . "\n" } @newMAT;
  close(OUT);
  GraphClust::makeCol($matrix_out_file);

}

sub locarnaP_aln2bmrels {
  my ( $alnfile, $probsdir, $consistency_transformation, $out_dir ) = @_;

  my $bmprobs;
  my $amprobs;

  my $aln = MLocarna::read_aln($alnfile);

  if ($consistency_transformation) {
    $bmprobs = MLocarna::read_bm_probs("$probsdir/bmprobs-cbt");
    $amprobs = MLocarna::read_am_probs("$probsdir/amprobs-cbt");
  } else {
    $bmprobs = MLocarna::read_bm_probs("$probsdir/bmprobs");
    $amprobs = MLocarna::read_am_probs("$probsdir/amprobs");
  }

  my ( $bmrels_seq_ref, $bmrels_str_ref, $amrels_ref );

  ( $bmrels_seq_ref, $bmrels_str_ref, $amrels_ref ) = MLocarna::compute_reliability( $aln, $bmprobs, $amprobs );

  ## nur ende des names nehmen
  $alnfile =~ /\/([^\/]+)$/;
  my $name = $1;
  MLocarna::write_bm_reliabilities( "$out_dir/$name.bmreliability", $bmrels_seq_ref, $bmrels_str_ref );
}

