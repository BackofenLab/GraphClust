#!/usr/bin/perl

## GraphClust MASTER script
## Author: Steffen Heyne heyne@informatik.uni-freiburg.de
## POD: see end of document

use strict;

########################################
## threading support
##
use Config;
my $THREADS_ENABLED = $Config{useithreads};
use threads;    # qw( async );
use threads::shared;
use Thread::Queue qw( );

use File::Temp qw(tempdir);
use IO::Handle;
use Pod::Usage;
use Getopt::Long;
use POSIX qw(ceil floor);
use Cwd qw(abs_path getcwd);
use FindBin;

use lib "$FindBin::Bin";
use GraphClust;
use vars qw(%CONFIG);
*CONFIG = \%GraphClust::CONFIG;
use Array::Utils qw(:all);
use List::Util qw/ min max /;

STDOUT->autoflush(1);

my $in_ROOTDIR;
my $in_FASTA_FILE;
my $in_USE_SGE = 0;
my $in_SGE_KILL;
my $in_config;
my $in_stage_clean;    ## not working, needs rewriting!!!
my $in_stage_end = -1;
my $in_verbose;
my $in_collect_results;
my $in_NO_SGE_SSH = 0;
my $in_grey_list_file;
my $in_debug = 0; ## 1: turns on storing of approx_knn,knn and ascii feature files
my $in_man;
my $in_help;
my $in_num_threads;
my $in_SGE_PE_name = "*";
my $newInf;

GetOptions(
  "root-dir=s"  => \$in_ROOTDIR,
  "fasta=s"     => \$in_FASTA_FILE,
  "sge"         => \$in_USE_SGE,
  "sge-kill"    => \$in_SGE_KILL,
  "sge-pe=s"    => \$in_SGE_PE_name,
  "no-ssh"      => \$in_NO_SGE_SSH,
  "config=s"    => \$in_config,
  "threads=i"   => \$in_num_threads,
  "stage-end=i" => \$in_stage_end,
  "results"     => \$in_collect_results,
  "v|verbose"   => \$in_verbose,
  "gl=s"        => \$in_grey_list_file,
  "debug"       => \$in_debug,
  "help"        => \$in_help,
  "man"         => \$in_man,
  "infernal_1_1" => \$newInf,
) || pod2usage(1);

pod2usage( -exitstatus => 0, -verbose => 1 ) if ($in_help);
pod2usage( -exitstatus => 0, -verbose => 2 ) if ($in_man);

if ( !$in_ROOTDIR ) {
  print "\nPlease provide root dir with --root-dir <DIR> option!\n\n";
  pod2usage( -exitstatus => -1, -verbose => 0 );
}

if ( !-e "$in_ROOTDIR/FASTA/data.fasta" && !$in_FASTA_FILE ) {
  print "\nPlease provide fasta file with --fasta <FILE> option!\n\n";
  pod2usage( -exitstatus => -1, -verbose => 0 );
}

################################################################################

SECTION("MASTER Graph Clustering");
SECTION("stage 0: initialize directories and data files");

$in_stage_end = 9 if ( $in_stage_end == -1 );

my $BIN_DIR = abs_path($FindBin::Bin);

$in_ROOTDIR = abs_path($in_ROOTDIR);
mkdir($in_ROOTDIR);

$in_FASTA_FILE = abs_path($in_FASTA_FILE);
$in_config = abs_path($in_config) if ($in_config);

my $CURRDIR = getcwd;
$in_ROOTDIR =~ /.*\/([^\/]+)$/;
my $ROOT_NAME = $1;

my $FASTA_DIR   = "$in_ROOTDIR/FASTA";
my $GSPAN_DIR   = "$in_ROOTDIR/GSPAN";
my $SVECTOR_DIR = "$in_ROOTDIR/SVECTOR";
my $CLUSTER_DIR = "$in_ROOTDIR/CLUSTER";
my $SGE_ERR_DIR = "$in_ROOTDIR/SGE_LOG";
my $EVAL_DIR    = "$in_ROOTDIR/EVAL";
my $RESULTS_DIR = "$in_ROOTDIR/RESULTS";

mkdir($FASTA_DIR);
mkdir($GSPAN_DIR);
mkdir($SVECTOR_DIR);
mkdir($CLUSTER_DIR);
mkdir($SGE_ERR_DIR);
mkdir($EVAL_DIR);
mkdir($RESULTS_DIR);
mkdir("$EVAL_DIR/times");
mkdir("$EVAL_DIR/svector");
mkdir("$EVAL_DIR/cluster");

print "used ROOT dir: " . $in_ROOTDIR . "\n";
print "used BINARY dir: " . $BIN_DIR . "\n";
print "used FASTA dir: " . $FASTA_DIR . "\n";
print "used GSPAN dir: " . $GSPAN_DIR . "\n";
print "used SVECTOR dir: " . $SVECTOR_DIR . "\n";
print "used CLUSTER dir: " . $CLUSTER_DIR . "\n";
print "used EVAL dir: " . $EVAL_DIR . "\n";
print "used RESULTS dir: " . $RESULTS_DIR . "\n";

################################################################################
## cleanup if provided

if ($in_stage_clean) {
  if ( $in_stage_clean !~ /1|2|3|4|5|6|7|8/ ) {
    print "\nClean stage $in_stage_clean not recognized! Please use only 1-8! Exit...\n\n";
  }
  SECTION("CLEANUP");
  cleanStage($in_stage_clean);
  SECTION("END MASTER SCRIPT Graph Clustering");
  exit;
}

################################################################################

my $OPTS_nspdk;
my $OPTS_nspdk_centers;
my $OPTS_fasta2shrep_gspan;
my $OPTS_RNAfold;
my $OPTS_RNAplfold;

my $GLOBAL_group_size;
my $GLOBAL_plfold_minlen;
my $GLOBAL_iterations;
my $GLOBAL_hit_blacklist_overlap;
my $GLOBAL_num_clusters;

my $evaluate;

my $nspdk_knn_center;
my $nspdk_nhf;
my $nspdk_nhf_max;
my $nspdk_nhf_step;
my $nspdk_fcs;

my $input_add_revcompl;
my $input_blastclust_id;
my $input_blastclust_len;
my $input_seq_min_length;
my $input_win_shift;
my $input_win_size;

################################################################################
## set config

SUBSECTION("USED Options for MASTER");
if ( -e "$in_ROOTDIR/config" && !$in_config ) {

  print "Use config from $in_ROOTDIR/config!\n";
  %CONFIG = readConfigFile("$in_ROOTDIR/config");

} elsif ( !-e "$in_ROOTDIR/config" && $in_config ) {

  print "Use config from $in_config!\n";
  %CONFIG = readConfigFile($in_config);
  writeConfig("$in_ROOTDIR/config");

} elsif ( !-e "$in_ROOTDIR/config" && !$in_config ) {
  print "Use default config from MASTER script!\n";
  writeConfig("$in_ROOTDIR/config");
} elsif ( -e "$in_ROOTDIR/config" && $in_config ) {
  print "\nRoot Dir already exists with user config file! Used config: $in_ROOTDIR/config\n";
  %CONFIG = readConfigFile("$in_ROOTDIR/config");
}

writeConfig("$in_ROOTDIR/config.start") if ( !-e "$in_ROOTDIR/config.start" );

setConfig();
printConfig( \%CONFIG );

###############################################################################
# end handler for temporary directory
# adds an error handler that deletes the directory in case of error
# SIGUSR{1/2} are sent by the sge prior to the uncatchable SIGKILL if the
# option -notify was set
###############################################################################
$SIG{'INT'}     = 'end_handler';
$SIG{'TERM'}    = 'end_handler';
$SIG{'ABRT'}    = 'end_handler';
$SIG{'USR1'}    = 'end_handler';
$SIG{'USR2'}    = 'end_handler';
$SIG{'__DIE__'} = 'end_handler';
$SIG{'KILL'}    = 'end_handler';

################################################################################
## tmp
my $tmp;
$tmp = $CONFIG{PATH_TMP} or $tmp = '/var/tmp/';
my $tmp_template = 'MASTER_GRAPHCLUST-XXXXXX';

# CLEANUP => 1 : automatically delete at exit
$tmp = tempdir( $tmp_template, DIR => $tmp, CLEANUP => 1 );

################################################################################
## threads/sge

my $SGE_PE_THREADS = 1;
my $NUM_THREADS = 1;

if ($in_USE_SGE && $in_num_threads>1){
  $SGE_PE_THREADS = $in_num_threads;
} elsif ($in_num_threads>1){
  $NUM_THREADS = $in_num_threads;
}

## set openMP environment variable to requested number of threads
$ENV{OMP_NUM_THREADS} = $NUM_THREADS;

$NUM_THREADS = 1 if ( !$THREADS_ENABLED || $in_SGE_KILL );

my %jobs_active;    ## active jobs/threads/sge tasks
my @workers;        ## worker threads if we use threads
my $q = Thread::Queue->new();    ## thread queue

if ( $NUM_THREADS > 1 ) {
  for ( 1 .. $NUM_THREADS ) {
    push @workers, async {
      while ( defined( my $call = $q->dequeue() ) ) {
        my $ret = call_thread($call);

        if ($ret != 0){
          print "Error during thread call:\n$call\n";
        }

      }
    };
  }
} elsif ( $in_USE_SGE || -e "$in_ROOTDIR/joblist.sge" ) {
  my $sge_last_jobs = [];
  $sge_last_jobs = GraphClust::read_partition("$in_ROOTDIR/joblist.sge")
    if ( -e "$in_ROOTDIR/joblist.sge" );
  print "\n";
  foreach my $job ( @{$sge_last_jobs} ) {
    print " read sge job " . join( ":", @{$job} ) . "\n" if ($in_verbose);
    $jobs_active{ $job->[0] } = $job;
  }

}

################################################################################
## sge kill

if ($in_SGE_KILL) {

  SECTION("Kill all SGE jobs");

  print "\n No SGE jobs found!\n" if ( !keys %jobs_active );

  foreach my $job ( keys %jobs_active ) {
    print "delete sge job id=" . $jobs_active{$job}->[3] . "...\n";
    my $ret = sge_call( "qdel " . $jobs_active{$job}->[3] );
    print "SGE reported: " . join( "", @{$ret} );
    delete $jobs_active{$job};
  }

  end_handler(0);
}

################################################################################
## do final results on request

if ($in_collect_results) {
  collect_results();
  SECTION( "END MASTER SCRIPT GraphClust " . $CONFIG{VERSION_INFO} . " ($ROOT_NAME)" );
  end_handler(0);
}


################################################################################
## stage 1 - prepare fasta file
SECTION("stage 1: preprocessing input sequences");

my $DATA_prefix = "data";

if ( !-e "$FASTA_DIR/data.fasta" ) {
  die "Please provide fasta file with --fasta <FILE>! Exit...\n\n" if ( !$in_FASTA_FILE );
  my $fasta_opts = "--prefix $DATA_prefix --tgtdir $FASTA_DIR ";
  $fasta_opts .= "--winsize $input_win_size ";
  $fasta_opts .= "--winshift $input_win_shift ";
  $fasta_opts .= "--bc-len $input_blastclust_len ";
  $fasta_opts .= "--bc-id $input_blastclust_id ";
  $fasta_opts .= "--min-len $input_seq_min_length ";
  $fasta_opts .= "--revcompl " if ($input_add_revcompl);
  $fasta_opts .= "--add-rc-sig " if ( !$CONFIG{cm_top_only} );
  $fasta_opts .= "--evaluate " if ($evaluate);
  $fasta_opts .= "--no-bc " if ( !$CONFIG{input_blastclust} );
  $fasta_opts .= "--gl $in_grey_list_file " if ($in_grey_list_file);
  $fasta_opts .= "--gl-num " . $CONFIG{nspdk_greylist_num} . " " if ( $CONFIG{nspdk_greylist_num} );

  system_call( "perl $BIN_DIR/graphFasta.pl $fasta_opts $in_FASTA_FILE ", $in_verbose, "$EVAL_DIR/times/time.stage.0.master" );
}

## used also in stage 6 for center.fa
my @fa       = GraphClust::read_fasta_file("$FASTA_DIR/data.fasta");
my $num_seqs = @{ $fa[1] };
print "Number of sequences in FASTA/data.fasta: " . $num_seqs . "\n";

die "Fasta file $FASTA_DIR/data.fasta contains only $num_seqs sequences! Exit...\n\n" if ( $num_seqs <= 2 );

end_handler(0) if ( $in_stage_end <= 1 );

if ( $num_seqs < $nspdk_knn_center ) {
  print "\nPlease select less 'nspdk_knn_center' as you only have $num_seqs sequences in you data set! Exit...\n\n";
  end_handler(0);
}

################################################################################
## 2) generate single gspan files from $in_FASTA_FILE

my $num_groups = ceil( $num_seqs / $GLOBAL_group_size );

SECTION("stage 2: Graph Generation (parallel)");

if ( !-e "$GSPAN_DIR/gspan.DONE" ) {
  my $CMD_fasta2shrep = [];
  $CMD_fasta2shrep->[0] = "perl $BIN_DIR/fasta2shrep_gspan.pl";
  $CMD_fasta2shrep->[1] = " --fasta $FASTA_DIR/data.fasta -o $GSPAN_DIR $OPTS_fasta2shrep_gspan --group $GLOBAL_group_size --tmp " . $CONFIG{PATH_TMP};
  my $job_name = "stage 2 (fasta2gspan)";
  my $job_uuid = "stage2";
  my $sge_status = job_call( $job_name, "$BIN_DIR/fasta2shrep_gspan.sge", $CMD_fasta2shrep, $num_groups, $SGE_ERR_DIR, $in_USE_SGE, "$GSPAN_DIR/sge_log", "$EVAL_DIR/times/time.stage.1", 1, "", $NUM_THREADS, $job_uuid );

  system("touch $GSPAN_DIR/tmp.no_match");
  system("cat $GSPAN_DIR/*.no_match > $FASTA_DIR/$DATA_prefix.no_match");
  system("rm $GSPAN_DIR/tmp.no_match");

  if ( $sge_status->[0] == 1 && $sge_status->[1] == 0 ) {
    system("touch $GSPAN_DIR/gspan.DONE");
  }

} else {
  print "graph generation for $num_seqs sequences done\n";
}

end_handler(0) if ( $in_stage_end <= 2 || !-e "$GSPAN_DIR/gspan.DONE" );

################################################################################
## 3) NSPDK

SECTION("stage 3: make sparse vector for all groups (parallel)");

if ( !-e "$SVECTOR_DIR/svector.groups.DONE" ) {
  my $CMD_gspanGroups = [];    ## additional arg for sge script
  $CMD_gspanGroups->[0] = "$GSPAN_DIR";       ## additional arg for sge script
  $CMD_gspanGroups->[1] = "$SVECTOR_DIR";     ## additional arg for sge script
  #$CMD_gspanGroups->[2] = "$in_debug";       ## additional arg for sge script
  $CMD_gspanGroups->[2] = "0";                ## never turn on debug mode, only useful in rare cases
  #$CMD_gspanGroups->[3] = "$BIN_DIR/EDeN";
  $CMD_gspanGroups->[3] = "$BIN_DIR/NSPDK";
  $CMD_gspanGroups->[4] = "$OPTS_nspdk";
  #--action FEATURE --binary_file_type --file_type GRAPH -r 3 -d 3 --graph_type DIRECTED --kernel_type NSPDK -i 1.group.gspan

  my $job_name = "stage 3 (gspanGroups.NSPDK)";
  my $job_uuid = "stage3";
  my $sge_status = job_call( $job_name, "$BIN_DIR/gspanGroups.NSPDK.sge", $CMD_gspanGroups, $num_groups, $SGE_ERR_DIR, $in_USE_SGE, "$SVECTOR_DIR/sge_log", "$EVAL_DIR/times/time.stage.3", 1, "", $NUM_THREADS, $job_uuid );

  if ( $sge_status->[0] == 1 && $sge_status->[1] == 0 ) {
    system("touch $SVECTOR_DIR/svector.groups.DONE");
  }

} else {
  print "svector for $num_groups groups already done\n";
}

end_handler(0) if ( $in_stage_end <= 3 || !-e "$SVECTOR_DIR/svector.groups.DONE" );

################################################################################
## 4) make single sparse vector (SVECTOR)
SECTION("Stage 4: make single sparse vector (serial)");

if ( !-e "$SVECTOR_DIR/data.svector.DONE" && !-e "$SVECTOR_DIR/data.svector.1" ) {
  system("rm $SVECTOR_DIR/svector.groups.DONE");
  print "\n\nError! data.svector cannot created without feature data!\n";
  print "Please re-run MASTER script in order to generate gspan features again!\n";
  end_handler(1);
}

if ( !-e "$SVECTOR_DIR/data.svector.DONE" ) {
  my $time_file = "$EVAL_DIR/times/time.stage.4.master";
  system( "echo  \`date\` > $time_file; echo " . time . " >> $time_file" );
  system_call( "\\rm -f $SVECTOR_DIR/data.svector; for i in \$(seq 1 $num_groups); do cat $SVECTOR_DIR/data.svector.\$i >> $SVECTOR_DIR/data.svector; done",
    "$in_verbose,$EVAL_DIR/times/time.master.4" );

  #if ($in_debug) {
    ##system_call( "\\rm -f $SVECTOR_DIR/feat.data.svector; for i in \$(seq 1 $num_groups); do cat $SVECTOR_DIR/feat.data.svector.\$i >> $SVECTOR_DIR/feat.data.svector; done");
    ##system("rm $SVECTOR_DIR/feat.data.svector.*");
  #}

  ## remove binary features of each group
  system_call("for i in \$(seq 1 $num_groups); do rm -f $SVECTOR_DIR/data.svector.\$i; done");
  ## remove gspan files
##  system_call("\\rm -f $GSPAN_DIR/*.group.gspan.* ");

  system("touch $SVECTOR_DIR/data.svector.DONE");
  system( "echo  \`date\` >> $time_file; echo " . time . " >> $time_file " );
} else {
  print "svector for $num_groups groups already done\n";
}

end_handler(0) if ( $in_stage_end <= 4 || !-e "$SVECTOR_DIR/data.svector.DONE" || !-e "$SVECTOR_DIR/data.svector" );

################################################################################
## ITERATIVE GraphClust starts here
## one round iterates through stages 5-8
## once all infernal scans are finshed, all hits are removed from the input set
## and everything is done again


## set automatic dense center overalpping values
## overlap is increased if we find <=50% of $GLOBAL_num_clusters
## max overlap between dense centers is 50% of dense center size
my $nspdk_mi_max  = int($nspdk_knn_center/2);
my $nspdk_mi_step = max(1, int($nspdk_mi_max/5)) ;
my $nspdk_mi      = 0;


foreach my $CI ( 1 .. $GLOBAL_iterations ) {

  SECTION("Iterative GraphClust: start round $CI/$GLOBAL_iterations");
  #sleep(1);

  if ( !-e "$EVAL_DIR/times/time.round.$CI" ) {    ## only for time measurment
    system("echo \`date\` > $EVAL_DIR/times/time.round.$CI");
    system( "echo " . time . " >> $EVAL_DIR/times/time.round.$CI" );
  }

################################################################################
## set params for this round
##  - increase number of hash functions for nspdk clustering in each round
##  - check how many models were found during last round
##    and increase dense region overlap (nspdk_mi) if no models were found

  my $clusters_last_round = [];

  if ( $CI > 1 ) {

    ## new: only do if we haven't found any cluster
    #$nspdk_nhf += $nspdk_nhf_step if ( $nspdk_nhf + $nspdk_nhf_step <= $nspdk_nhf_max );

    ## get list of center-idx from last round (soft partition) which lead to cluster (> results_min_cluster_size)
    $clusters_last_round = foundClusters( $CI - 1 );
    my $min_diff = 0;

    SUBSECTION("Parameters for round $CI");

    my $center_last_round = 0;
    ## fast_cluster file not exists if $in_stage_end = 5 (also used with special blacklist)
    if ( -e "$SVECTOR_DIR/data.svector." . ( $CI - 1 ) . ".fast_cluster" ) {
      my @t = readpipe( "cat $SVECTOR_DIR/data.svector." . ( $CI - 1 ) . ".fast_cluster" );
      chomp(@t);
      $center_last_round = @t;
    }

    print "\nModels found in last round: " . @{$clusters_last_round} . "\n";
    if ( @{$clusters_last_round} ) {
      print "\n" . join( "\n", map { ( $CI - 1 ) . ".$_.cluster/MODEL" } @{$clusters_last_round} ) . "\n";
    }

    my $tr1 = 0;

    if ( $center_last_round <= $GLOBAL_num_clusters / 2 && $nspdk_mi + $nspdk_mi_step <= $nspdk_mi_max && $nspdk_mi_step > 0 ) {
      $nspdk_mi += $nspdk_mi_step;
      print "\nToo few dense centers found in last round!\n";
      print "Set new overlap of dense regions: $nspdk_mi\n";
    } elsif ( !@{$clusters_last_round} ) { $tr1 += 1 }

    if ( @{$clusters_last_round} <= $center_last_round * (3/5) && $nspdk_nhf + $nspdk_nhf_step <= $nspdk_nhf_max && $nspdk_nhf_step > 0 ) {
      $nspdk_nhf += $nspdk_nhf_step;
      print "\nOnly few clusters found in last round!\n";
      print "Set new number of hash functions: $nspdk_nhf\n";
    } elsif (!@{$clusters_last_round}) { $tr1 += 1 }

    if ( $tr1 >= 2 ) {
      print "No cluster found in last round!\n";
      print "\nDense region overlap and number of hash functions is at maximum! (mi=$nspdk_mi, nhf=$nspdk_nhf )\n";
      print "No new iteration is started!\n\n";
      SECTION("Trigger set to stop terative GraphClust now!");
      system("echo  \`date\` >> $EVAL_DIR/times/time.round.$CI"); ## only for time measurment
      system( "echo " . time . " >> $EVAL_DIR/times/time.round.$CI" ); ## only for time measurment
      last;
    }

  }

  print "\n Parameters for round $CI: nspdk_nhf = $nspdk_nhf - nspdk_mi = $nspdk_mi\n";

  if ( -e "$in_ROOTDIR/$CI.round.DONE" ) {
    print "\n Round $CI already finished! Switch to next round...\n";
    next;
  }

################################################################################
## make white/blacklist for current round from CM hits from last round

  SUBSECTION("Create blacklist from all clusters/ cmsearch hits.");
  ## we always have blacklist paramter for NSPDK stage 5 -> need always a file
  system("touch $SVECTOR_DIR/data.svector.blacklist.$CI") if ( !-e "$SVECTOR_DIR/data.svector.blacklist.$CI" );
  ## todo: !!! tempBlacklist is not working, needs refactor but currently we can live without it
  system("touch $SVECTOR_DIR/blacklist.no_model");

  if ( $CI > 1 && !-e "$SVECTOR_DIR/data.svector.blacklist.$CI.special" ) {
    ## make blacklist for this round from all Infernal hits so far
    ## collect all cmsearch hits from last round as blacklist
    ## hits are based on hard merged partition file, i.e. on last glob_results.pl call
    makeBlacklist( $CI, "$SVECTOR_DIR/round." . ( $CI - 1 ) . ".hits", $GLOBAL_hit_blacklist_overlap );
  }

  ## combine temporary blacklist and hit list from last round to new blacklist for this round
  ## iteration 1 diff from other iterations
  if ( $CI > 1 ) {
    system( "cat $SVECTOR_DIR/blacklist.no_model $SVECTOR_DIR/round." . ( $CI - 1 ) . ".hits $FASTA_DIR/$DATA_prefix.no_match > $SVECTOR_DIR/data.svector.blacklist.$CI" );
  } else {
    system("cat $SVECTOR_DIR/blacklist.no_model $FASTA_DIR/$DATA_prefix.no_match > $SVECTOR_DIR/data.svector.blacklist.$CI");
  }

  ## blacklist size
  my $blacklist_curr_size = 0;
  my @bl_num = readpipe("cat $SVECTOR_DIR/data.svector.blacklist.$CI");
  chomp(@bl_num);
  $blacklist_curr_size = unique(@bl_num);
  print "\nFINAL blacklist size for round $CI: $blacklist_curr_size\n";

  ## $in_stage_end=10 -> special mode: only use NSPDK and predict candidate cluster,
  ## but do nothing else, stop each iteration after stage 5
  if ( -e "$SVECTOR_DIR/data.svector.blacklist.$CI.special" ) {
    system("cp $SVECTOR_DIR/data.svector.blacklist.$CI.special $SVECTOR_DIR/data.svector.blacklist.$CI");
    $in_stage_end = 10;
    $blacklist_curr_size = `cat $SVECTOR_DIR/data.svector.blacklist.$CI.special | wc -l`;
  }

################################################################################
## do nothing for this round if we have to few sequences left

  my $num_seqs_left = $num_seqs - $blacklist_curr_size;
  if ( $num_seqs_left < 2 ) {
    ## do not call NSPDK again
    print "\n...only $num_seqs_left sequences left for clustering.\n";
    print " Skip all following iterations!\n\n";
    SECTION("Trigger set to stop terative GraphClust now!");
    system("echo  \`date\` >> $EVAL_DIR/times/time.round.$CI"); ## only for time measurment
    system( "echo " . time . " >> $EVAL_DIR/times/time.round.$CI" ); ## only for time measurment
    last;
  }

################################################################################
## 5) call NSPDK on sparse vectors to get candiate clusters via MIN-HASH

  SECTION("Round $CI Stage 5: calculate candidate cluster with NSPDK (min-hash)");
  if ( !-e "$SVECTOR_DIR/data.svector.fast_cluster.$CI.DONE" ) {

    my $job_name = "stage 5 (fastCluster.NSPDK)";
    my $job_uuid = "stage5-$CI";

    ## fix nspdk binsize to 1000 instances to be independent of dataset size
    $OPTS_nspdk_centers =~ s/-msb\s+\S+//;
    my $nspdk_max_binsize = 1000;
    if ( $num_seqs_left > $nspdk_max_binsize ) {
      my $msb = sprintf( "%.4f", ( $nspdk_max_binsize / $num_seqs_left ) );
      $msb = 0.0001 if ( $msb < 0.0001 );
      $OPTS_nspdk_centers .= " -msb $msb";
      print "use msb $msb (num seqs left:$num_seqs_left binsize:$nspdk_max_binsize)\n";
    } else {
      $OPTS_nspdk_centers .= " -msb 1";
    }

    ## nspdk uses vector filename for output (append output type like .fast_cluster .knn .approx_knn)
    ## we need this here to use different names in case of stage_end=5,
    ## i.e. there could be multiple running nspdk instances on the same feature vector
    system("ln -f -s $SVECTOR_DIR/data.svector $SVECTOR_DIR/data.svector.$CI");

    my $CMD_fastClusterNSPDK = [];
    $CMD_fastClusterNSPDK->[0] = "$BIN_DIR/NSPDK";
    $CMD_fastClusterNSPDK->[1] =
      "-no-cache -rs $CI -fsb $SVECTOR_DIR/data.svector.$CI -bl $SVECTOR_DIR/data.svector.blacklist.$CI $OPTS_nspdk_centers -knn $nspdk_knn_center -ss " . $GLOBAL_num_clusters . " -nhf $nspdk_nhf -mi $nspdk_mi -fcs $nspdk_fcs ";

## NEW EDEN
# ~/workspace64/GraphClust/scripts/EDeN --action CLUSTER --binary_file_type --file_type SPARSE_VECTOR --cluster_type DENSE_CENTERS -R 1 --eccess_neighbour_size_factor 5 --num_nearest_neighbours 15 --sample_size 5 --num_hash_functions 300 --max_intersection_size 0 --fraction_center_scan 0.5 --max_size_bin 0.5 --shared_neighborhood --num_repeat_hash_functions 2 --force_approximate  -i 1.group.gspan.feature

# K quick shift clustering
#~/workspace64/GraphClust/scripts/EDeN --action CLUSTER --binary_file_type --file_type SPARSE_VECTOR --cluster_type K_QUICK_SHIFT -R 1 --eccess_neighbour_size_factor 5 --num_nearest_neighbours 15 --num_hash_functions 300 --max_size_bin 0.5 --shared_neighborhood --num_repeat_hash_functions 2 --force_approximate --cluster_threshold 1  -i 1.group.gspan.feature

# learning
#~/workspace64/GraphClust/scripts/EDeN --action TRAIN --binary_file_type --file_type SPARSE_VECTOR -t target.test -m model.name -i 1.group.gspan.feature
# ~/workspace64/GraphClust/scripts/EDeN --action FEATURE_SCALED --binary_file_type --file_type SPARSE_VECTOR  -m model.name -i 1.group.gspan.feature

    ## add special options for greylist clustering
    if ( -e "$FASTA_DIR/$DATA_prefix.greylist" ) {
      #my $greylist_size = `cat $FASTA_DIR/$DATA_prefix.greylist | wc -l`;
      #chomp($greylist_size);
      $CMD_fastClusterNSPDK->[2] = " -gl $FASTA_DIR/$DATA_prefix.greylist -fcs 1 -otknn "; # -knn ".($nspdk_knn_center+$greylist_size)." ";
    }

    ## add special debug options for NSPDK
    #if ($in_debug) {
    #  $CMD_fastClusterNSPDK->[3] = " -oaknn ";
    #}

    ## get size of feature vector and estimate memory requirement
    my $size = -s "$SVECTOR_DIR/data.svector";
    $size = ceil( $size / 1000000000 ) + ceil( $nspdk_nhf * 0.005 );
    $size = $size * ( $num_seqs_left / $num_seqs ) + 1;
    $size = sprintf( "%.1f", $size );
    print "used qsub h_vmem=$size" . "G\n";
    print "use threads $NUM_THREADS use SGE-PE threads $SGE_PE_THREADS\n";
    my $qsub_opts = " -l h_vmem=$size" . "G ";
    $qsub_opts .= " -v OMP_STACKSIZE=20M  -v OMP_NUM_THREADS=$SGE_PE_THREADS ";
    $qsub_opts .= " -pe \"$in_SGE_PE_name\" 1-$SGE_PE_THREADS " if ($SGE_PE_THREADS>1); ## use '-hard -l c_op2356=1' to have always same cpu model (FREIBURG only)

    ## we wait until NSPDK is finished, do not wait in special blacklist case == 10 (only for internal debug)
    my $sge_wait = 1;
    $sge_wait = 0 if ( $in_stage_end == 10 );

    my $sge_status =
      job_call( $job_name, "$BIN_DIR/fastCluster.NSPDK.sge", $CMD_fastClusterNSPDK, 1, $SGE_ERR_DIR, $in_USE_SGE, "$SVECTOR_DIR/sge_log_stage5.$CI", "$EVAL_DIR/times/time.stage.5.$CI", $sge_wait, $qsub_opts, 1, $job_uuid );

    if ( $sge_status->[0] == 1 && $sge_status->[1] == 0 ) {

      ## prefilter centers for interesting ones,
      ## i.e. centers which contain non greylist ids
      if ( -e "$FASTA_DIR/$DATA_prefix.greylist" ) {
        system("cp $SVECTOR_DIR/data.svector.$CI.fast_cluster $SVECTOR_DIR/data.svector.$CI.fast_cluster_orig");
        ## filterGreylistCenters("$FASTA_DIR/$DATA_prefix.greylist","$SVECTOR_DIR/data.svector.$CI.fast_cluster",$nspdk_knn_center);
        filterGreylistCenters( "$FASTA_DIR/$DATA_prefix.greylist", "$SVECTOR_DIR/data.svector.$CI.approx_knn", $nspdk_knn_center );
        system("cp $SVECTOR_DIR/data.svector.$CI.approx_knn.gl_filter $SVECTOR_DIR/data.svector.$CI.fast_cluster");
      }
      system("rm $SVECTOR_DIR/data.svector.$CI");
      system("touch $SVECTOR_DIR/data.svector.fast_cluster.$CI.DONE");

      if ($evaluate) {
        system_call(
          "$BIN_DIR/clusterQuality.sh $FASTA_DIR/class.hash $SVECTOR_DIR/data.svector.$CI.fast_cluster $SVECTOR_DIR/data.svector.$CI.fast_cluster $nspdk_knn_center $EVAL_DIR/svector $FASTA_DIR/class.size $FASTA_DIR/$DATA_prefix.map > $EVAL_DIR/svector/$CI.centers_qual",
          $in_verbose
        );
        system_call("mv $EVAL_DIR/svector/cluster.counts $EVAL_DIR/svector/$CI.cluster.counts");
        system_call("mv $EVAL_DIR/svector/cluster.class.all $EVAL_DIR/svector/$CI.cluster.class.all");
        system_call("mv $EVAL_DIR/svector/cluster.class $EVAL_DIR/svector/$CI.cluster.class");
        system_call("mv $EVAL_DIR/svector/cluster.class.frags $EVAL_DIR/svector/$CI.cluster.class.frags");
        system_call("mv $EVAL_DIR/svector/class.size $EVAL_DIR/svector/$CI.class.size");

        GraphClust::evalSVECTOR( "$EVAL_DIR/svector", $CI, "$EVAL_DIR/stage5.svector_qual.final" );
      }

    }

  } else {
    print "Round $CI fast_cluster already done\n";
  }

  next if ( $in_stage_end <= 5 || !-e "$SVECTOR_DIR/data.svector.fast_cluster.$CI.DONE" );

  ## special blacklist end, switch to next iteration
  if ( $in_stage_end == 10 ) {
    system("touch $in_ROOTDIR/$CI.round.DONE");
    system("echo  \`date\` >> $EVAL_DIR/times/time.round.$CI"); ## only for time measurment
    system( "echo " . time . " >> $EVAL_DIR/times/time.round.$CI" ); ## only for time measurment

    next;
  }

################################################################################
## 6-8: process each dense Center to build final cluster model

  ## get real number of dense centers/clusters found by nspdk ($CONFIG{GLOBAL_num_clusters} could be higher)
  my $num_clusters_curr = readpipe("wc -l $SVECTOR_DIR/data.svector.$CI.fast_cluster");
  $num_clusters_curr =~ /(\d+)\s.*$/;
  $num_clusters_curr = $1;

  SECTION("Round $CI Stage 6-8: cluster $num_clusters_curr centers with LocARNA (parallel)");
  my %toDo_models = ();
  map { $toDo_models{"$CI.$_"} = $_ } ( 1 .. $num_clusters_curr );

  my @time_file = readpipe("cat $EVAL_DIR/times/time.round.$CI");
  my $time      = $time_file[1];
  chomp($time);

  system("echo 'round $CI clusters $GLOBAL_num_clusters clusters_last ".scalar(@{$clusters_last_round})." clusters_curr $num_clusters_curr nhf $nspdk_nhf mi $nspdk_mi' >> $EVAL_DIR/params.round");
  my $cluster_error = 0;

  my %job_task_finished     = ();
  my $trigger_new_partition = 0;

  print "\nWait for " . ( keys %jobs_active ) . " SGE jobs to finish ...\n\n" if ( $in_USE_SGE && keys %jobs_active );

  while ( keys %toDo_models ) {

    print "\n>>> Round $CI models left: " . ( keys %toDo_models ) . " (after " . sprintf( "%.1f", ( time() - $time ) / 60 ) . " min) - center overlap=$nspdk_mi nhf=$nspdk_nhf\n"
      if ($in_verbose);

    foreach my $clus_idx ( sort { $toDo_models{$a} <=> $toDo_models{$b} } keys %toDo_models ) {

      ####################################################################################################
      ## model is completely done
      if ( -e "$CLUSTER_DIR/$clus_idx.cluster/cmsearch.DONE" ) {
        print "Round $CI cluster $clus_idx stages 6-8 already finished!\n";
        delete $toDo_models{$clus_idx};
        next;
      }

      ## alignCenter.pl params for locarna paligs
      $clus_idx =~ /\d+\.(\d+)/;
      my $clus_idx_ci      = $1;
      my $curr_cluster_dir = "$CLUSTER_DIR/$clus_idx.cluster";

      ## alignCenter.pl params for locarna model msa
      my $params_maligs = "--root-dir $in_ROOTDIR ";
      $params_maligs .= "--cluster-dir $curr_cluster_dir ";
      $params_maligs .= "--evaluate " if ($evaluate);
      $params_maligs .= "--verbose " if ($in_verbose);
      $params_maligs .= "--debug " if ($in_debug);

      ## cmsearcher opts
      my $curr_cmsearch_dir = "$curr_cluster_dir/CMSEARCH";
      my $params_cmsearcher = "--root $in_ROOTDIR -tgtdir $curr_cmsearch_dir ";
      $params_cmsearcher .= "-stk $curr_cluster_dir/MODEL/model.stk -db $FASTA_DIR/data.fasta.scan ";
      $params_cmsearcher .= "--verbose " if ($in_verbose);
      $params_cmsearcher .= "--infernal_1_1 " if ($newInf);

      ##########################################################################################
      ## pre 6 stage: fold used seqs

      if ( !-e "$CLUSTER_DIR/$clus_idx.cluster/pp.DONE" ) {

        system("mkdir -p $CLUSTER_DIR/$clus_idx.cluster");

        my $ids_aref = GraphClust::readSubset( "$SVECTOR_DIR/data.svector.$CI.fast_cluster", $clus_idx_ci, $nspdk_knn_center );
        GraphClust::writeSet( $ids_aref, "$CLUSTER_DIR/$clus_idx.cluster/center.ids" );
        GraphClust::writeSubsetFrags( \@fa, $ids_aref, "$CLUSTER_DIR/$clus_idx.cluster/center.frags", "SEQ" );
        GraphClust::writeSubsetFasta( \@fa, $ids_aref, "$CLUSTER_DIR/$clus_idx.cluster/center.fa", 1 );

        my $ids_all = GraphClust::readSubset( "$SVECTOR_DIR/data.svector.$CI.fast_cluster", $clus_idx_ci, $nspdk_knn_center * 3 );
        GraphClust::writeSet( $ids_all, "$CLUSTER_DIR/$clus_idx.cluster/center.ids.ext" );
        GraphClust::writeSubsetFasta( \@fa, $ids_all, "$CLUSTER_DIR/$clus_idx.cluster/center.fa.ext", 1 );

        my $knn_sim = GraphClust::readSubset( "$SVECTOR_DIR/data.svector.$CI.fast_cluster_sim", $clus_idx_ci );

        open( OUT, ">$CLUSTER_DIR/$clus_idx.cluster/center.sim" );
        print OUT join( " ", @{$knn_sim} ) . "\n";
        close(OUT);

        system("touch $CLUSTER_DIR/$clus_idx.cluster/pp.DONE");
      }

      ####################################################################################################
      ## stage 6: compute all vs. all pairwise alignments for knn frags of dense center

      my $pp_dir = "$CLUSTER_DIR/$clus_idx.cluster/pp";
      my $dp_dir = "$CLUSTER_DIR/$clus_idx.cluster/dp";

      ## do nothing in stage 6 if we use kernel sim matrix to build guide tree
      if ( $CONFIG{center_tree_type} == 3 && -e "$SVECTOR_DIR/$DATA_prefix.svector.$CI.fast_cluster_sim" ) {
        system("touch $CLUSTER_DIR/$clus_idx.cluster/cluster.DONE");
      }

      if ( !-e "$CLUSTER_DIR/$clus_idx.cluster/cluster.DONE" ) {
        my $in_rnafold_opts   = $CONFIG{OPTS_RNAfold};
        my $in_rnaplfold_opts = $CONFIG{OPTS_RNAfold};
        my $vrna_path         = $CONFIG{PATH_VRNA};

        my $fold_opts = "--vrna-path $vrna_path --fasta $CLUSTER_DIR/$clus_idx.cluster/center.fa --tgtdir $dp_dir --new-only ";
        $fold_opts .= "--switch-len $GLOBAL_plfold_minlen ";
        $fold_opts .= "--rnafold-opts \"$OPTS_RNAfold\" " if ($in_rnafold_opts);
        $fold_opts .= "--rnaplfold-opts \"$OPTS_RNAplfold\" " if ($in_rnaplfold_opts);

        ## SUBSECTION(" Round $CI cluster $clus_idx stage pre 6 : Fold fragments with RNAfold/RNAplfold... ");
        if ( !-e $dp_dir || !-e $pp_dir ) {
          system_call( "perl $BIN_DIR/foldFasta.pl $fold_opts", 0 );

          system("mkdir -p $pp_dir");
          my $idx = 0;
          my $ids_aref = GraphClust::readSubset( "$CLUSTER_DIR/$clus_idx.cluster/center.ids", 1, $nspdk_knn_center );
          foreach my $key ( sort { $a <=> $b } @{$ids_aref} ) {
            $idx++;
            system("rm $pp_dir/$idx") if ( -e "$pp_dir/$idx" );
            die "Dotplot dp/$key does not exists! Exit...\n\n" if ( !-e "$dp_dir/$key" );
            system("ln -s $dp_dir/$key $pp_dir/$idx");
          }
        }

        my $paligs_sge_log = "$CLUSTER_DIR/$clus_idx.cluster/SGE_log";
        my $job_name       = "Round $CI cluster $clus_idx stage 6";
        my $job_uuid       = "stage6-$clus_idx";

        ## just to get number of aligs
        my $curr_knn_aref = GraphClust::readSubset( "$CLUSTER_DIR/$clus_idx.cluster/center.ids", 1, -1 );
        my $curr_knn = @{$curr_knn_aref};

        my $aligs            = ( $curr_knn * ( $curr_knn - 1 ) ) / 2;
        my $aligs_per_task   = 10;
        my $sge_paligs_tasks = ceil( $aligs / $aligs_per_task );

        my $CMD_alignRange = [];
        $CMD_alignRange->[0] = "$aligs";
        $CMD_alignRange->[1] = "$sge_paligs_tasks";
        $CMD_alignRange->[2] = "perl $BIN_DIR/graphClust_AlignRange.pl";
        $CMD_alignRange->[3] = "--root-dir $in_ROOTDIR --tgtdir $CLUSTER_DIR/$clus_idx.cluster -dpdir $pp_dir --locarna-path " . $CONFIG{PATH_LOCARNA} . "/locarna ";

        my $sge_status =
          job_call( $job_name, "$BIN_DIR/graphClust_AlignRange.sge", $CMD_alignRange, $sge_paligs_tasks, $SGE_ERR_DIR, $in_USE_SGE, $paligs_sge_log, "$EVAL_DIR/times/time.stage.6.$clus_idx", 0, "", $NUM_THREADS, $job_uuid, $job_task_finished{$job_uuid} );

        $job_task_finished{$job_uuid} = $sge_status->[2];

        if ( $sge_status->[0] == 1 && $sge_status->[1] == 0 ) {
          system("touch $CLUSTER_DIR/$clus_idx.cluster/cluster.DONE");
          delete $job_task_finished{$job_uuid};
        } elsif ( $sge_status->[1] == 1 ) {
          ## stage 6 sge finished with error
          print "Round $CI cluster $clus_idx stage 6: SGE job generated some error! Skip cluster $clus_idx...\n";
          delete $toDo_models{$clus_idx};
          delete $job_task_finished{$job_uuid};
          $cluster_error++;
        }
      }
      next if ( !-e "$CLUSTER_DIR/$clus_idx.cluster/cluster.DONE" );

      if ( $in_stage_end == 6 ) {
        print "Round $CI stage 6 end requested for cluster $clus_idx.\n";
        delete $toDo_models{$clus_idx};
        next;
      }

      ####################################################################################################
      ## stage 7: submit model MSA to SGE
      ## SGE_log dir is stored in $CLUSTER_DIR/$clus_idx.cluster/SGE_log_MSA (Because $MODEL_DIR/$clus_idx.cluster.align is only present if MSA was found!)
      ## file $MODEL_DIR/$clus_idx.clutser.align/model_build.DONE is present if a MSA was found to build model
      ## old: file $MODEL_DIR/$clus_idx.model.DONE is present if everything was correct with SGE job
      ## new: file $CLUSTER_DIR/$clus_idx.cluster/model.DONE is present if job alignCenter.pl finished correctly

      if ( !-e "$CLUSTER_DIR/$clus_idx.cluster/model.DONE" ) {

        my $job_name = "Round $CI cluster $clus_idx stage 7";
        my $job_uuid = "stage7-$clus_idx";

        my $CMD_stage7 = [];
        $CMD_stage7->[0] = "perl $BIN_DIR/alignCenter.pl";
        $CMD_stage7->[1] = "$params_maligs";

        my $sge_status =
          job_call( $job_name, "$BIN_DIR/alignCenter.sge", $CMD_stage7, 1, $SGE_ERR_DIR, $in_USE_SGE, "$CLUSTER_DIR/$clus_idx.cluster/SGE_log_MSA", "$EVAL_DIR/times/time.stage.7.$clus_idx", 0, "", $NUM_THREADS, $job_uuid );

        if ( $sge_status->[0] == 1 && $sge_status->[1] == 0 ) {
          system("touch $CLUSTER_DIR/$clus_idx.cluster/model.DONE");
          system("touch $CLUSTER_DIR/$clus_idx.cluster/model_build.DONE") if ( -e "$curr_cluster_dir/MODEL/model.stk" );

          if ( -e "$CLUSTER_DIR/$clus_idx.cluster/model_build.DONE" ) {
            print "Round $CI cluster $clus_idx stage 7: Model found and scan all seqs in stage 8!\n";
          }

        } elsif ( $sge_status->[1] == 1 ) {
          ## stage 7 sge finished with error
          print "Round $CI cluster $clus_idx stage 7: SGE job generated some error! Skip cluster $clus_idx...\n";
          delete $toDo_models{$clus_idx};
          $cluster_error++;
        }
      }

      next if ( !-e "$CLUSTER_DIR/$clus_idx.cluster/model.DONE" );

      if ( !-e "$CLUSTER_DIR/$clus_idx.cluster/model_build.DONE" ) {
        print "Round $CI cluster $clus_idx stage 7: No Model found!\n";
        delete $toDo_models{$clus_idx};
        next;
      }

      if ( $in_stage_end == 7 ) {
        print "Round $CI stage 7 stop requested for cluster $clus_idx.\n";
        delete $toDo_models{$clus_idx};
        next;
      }

      ####################################################################################################
      ## stage 8: infernal stage

      if ( !-e "$CLUSTER_DIR/$clus_idx.cluster/cmsearch.DONE" ) {

        my $job_name = "Round $CI cluster $clus_idx stage 8";
        my $job_uuid = "stage8-$clus_idx";

        my $CMD_stage8 = [];
        $CMD_stage8->[0] = "perl $BIN_DIR/gc_cmsearch.pl ";
        $CMD_stage8->[1] = "$params_cmsearcher ";



        my $sge_status = job_call( $job_name, "$BIN_DIR/gc_cmsearch.sge", $CMD_stage8, 1, $SGE_ERR_DIR, $in_USE_SGE, "$curr_cmsearch_dir/SGE_log", "$EVAL_DIR/times/time.stage.8.$clus_idx", 0, "", $NUM_THREADS, $job_uuid );

        if ( $sge_status->[0] == 1 && $sge_status->[1] == 0 ) {
          ## stage 8 sge finished without error
          system_call("touch $CLUSTER_DIR/$clus_idx.cluster/cmsearch.DONE");
          delete $toDo_models{$clus_idx};
          $trigger_new_partition = 1;

        } elsif ( $sge_status->[1] == 1 ) {
          ## stage 8 sge finished with error
          print "Round $CI cluster $clus_idx stage 8: SGE job generated some error! Skip cluster $clus_idx...\n";
          delete $toDo_models{$clus_idx};
          $cluster_error++;
        }
      }    ## fi stage 8

    }    ## foreach keys %toDo_models

    ## write center and f-measure summary only

if ( $evaluate && $trigger_new_partition ) {
      my $CMD_glob_res= "perl $BIN_DIR/glob_results.pl --root $in_ROOTDIR --summary-only --last-ci $CI ";
      $CMD_glob_res .= " --infernal_1_1 " if($newInf);
      #system_call("perl $BIN_DIR/glob_results.pl --root $in_ROOTDIR --summary-only --last-ci $CI $newInf");
      system_call($CMD_glob_res);
    }


#    if ( $evaluate && $trigger_new_partition ) {
 #     system_call("perl $BIN_DIR/glob_results.pl --root $in_ROOTDIR --summary-only --last-ci $CI ". $newInf);
  #  }

    $trigger_new_partition = 0;
    sleep(5);

  }    ## while %toDo_models not empty

  #### some stuff to finish current iteration

  #GraphClust::evalCLUSTER( "$EVAL_DIR/cluster", $CI, "$EVAL_DIR/stage7.cluster_qual.final", $evaluate );

#  system_call( "perl $BIN_DIR/glob_results.pl --root $in_ROOTDIR --last-ci $CI --summary-only ". $in_verbose, $newInf );

  my $CMD_glob_res= "perl $BIN_DIR/glob_results.pl --root $in_ROOTDIR --last-ci $CI --summary-only  $in_verbose ";
  $CMD_glob_res .= " --infernal_1_1 " if($newInf);
  #system_call( "perl $BIN_DIR/glob_results.pl --root $in_ROOTDIR --last-ci $CI --summary-only ", $in_verbose, $newInf );
  system_call($CMD_glob_res);


  die "\n\nStage end $in_stage_end requested! Exit...\n\n" if ( $in_stage_end <= 8 );
  die "\n\nToo many errors in stages 6-8! Exit...\n\n" if ( $cluster_error >= 1 );

  system("touch $in_ROOTDIR/$CI.round.DONE");
  system("echo  \`date\` >> $EVAL_DIR/times/time.round.$CI"); ## only for time measurment
  system( "echo " . time . " >> $EVAL_DIR/times/time.round.$CI" ); ## only for time measurment

}    ## foreach iteration/round

## create final results
collect_results() if ( $in_stage_end == 9 );    ## normal case

SECTION( "END MASTER SCRIPT " . $CONFIG{VERSION_INFO} . " ($ROOT_NAME)" );

end_handler(0);

################################################################################
## subs

sub collect_results {


	  SECTION("stage 9: Final RESULTS - collect all found clusters");

  my $part_file = "$RESULTS_DIR/partitions/final_partition.soft";
  if ( !-e "$in_ROOTDIR/FASTA/data.greylist" && uc( $CONFIG{results_partition_type} ) =~ /HARD/ ) {
    $part_file = "$RESULTS_DIR/partitions/final_partition.hard.merged";
  }

  my $forInf = "";
  $forInf = " --infernal_1_1 " if($newInf);
  system("perl $BIN_DIR/glob_results.pl --root $in_ROOTDIR --all $forInf ") if ( !-e $part_file    || uc( $CONFIG{results_partition_type} ) =~ /HARD/ );

#  system("perl $BIN_DIR/glob_results.pl --root $in_ROOTDIR --all ". $newInf ) if ( !-e $part_file || uc( $CONFIG{results_partition_type} ) =~ /HARD/ );
  if (!-e $part_file){
    SUBSECTION("No partition file found! Maybe 0 clusters found so far!");
    return;
  }

  ## part: array of [528,42.40,1.1,1.2,1,63,63]
  my $part = GraphClust::read_partition($part_file);
  my @res_clus_sort = sort { $a <=> $b } unique( map { $_->[5] } @{$part} );

  return if ( !@res_clus_sort );    ## do we have results at all

  if ( !-e "$RESULTS_DIR/results.DONE" ) {

    ## delete all old results as a new partition could change everything
    if ( !-e "$RESULTS_DIR/SGE_log/task.submitted" ) {
      print "\nAll old results will be deleted!\n";
      sleep(1);
      my @res_dir = readpipe("ls $RESULTS_DIR");
      chomp(@res_dir);
      map { system("rm -r -f $RESULTS_DIR/$_ ") if ( $_ =~ /\d+/ ) } @res_dir;
    }

    my $res_jobs = @res_clus_sort;
    my $job_name = "stage 9: Collecting results and make final clusters";
    my $job_uuid = "stage9";

    my $CMD_stage9 = [];
    $CMD_stage9->[0] = "perl $BIN_DIR/gc_results_cluster.pl";
    $CMD_stage9->[1] = "--root-dir $in_ROOTDIR ";
    $CMD_stage9->[1] .= "--verbose "  if ($in_verbose);
    $CMD_stage9->[1] .= "--evaluate " if ($evaluate);
    $CMD_stage9->[1] .= "--infernal_1_1 " if ($newInf);

    ## (finished?, error?, jobs_finished, jobs_all)
    my $sge_status = [ 0, 1, 0, 0 ];
    while ( $sge_status->[0] == 0 ) {
      $sge_status = job_call( $job_name, "$BIN_DIR/gc_jobscript.sge", $CMD_stage9, $res_jobs, $SGE_ERR_DIR, $in_USE_SGE, "$RESULTS_DIR/SGE_log", "$EVAL_DIR/times/time.stage.9", 1, "", $NUM_THREADS, $job_uuid );
    }
    if ( $sge_status->[0] == 1 && $sge_status->[1] == 0 ) {
      system("touch $RESULTS_DIR/results.DONE");
    } elsif ( $sge_status->[1] == 1 ) {
      ## stage 7 sge finished with error
      print "Results stage 9: SGE job generated some error!\n";
    }
  }

  my $stats_file = "$RESULTS_DIR/cluster.final.stats";
  system("rm -f $stats_file");
  foreach my $clus (@res_clus_sort) {
    system("cat $RESULTS_DIR/$clus/cluster.stats >> $stats_file");
  }

  GraphClust::makeCol($stats_file);

  SUBSECTION("Results done!");

  print "\nAll clusters can be found in:\n $ROOT_NAME/CLUSTERS\n\n";
  print "To re-compute results please delete first file:\n $ROOT_NAME/RESULTS/results.DONE\n\n";
  my $CMD_res = "MASTER_GraphCluster.pl --root $ROOT_NAME --results ";
  $CMD_res .= "--sge " if ($in_USE_SGE);
  print "Then invoke pipeline again with:\n $CMD_res\n";
}

sub filterGreylistCenters {
  my $gl_file     = $_[0];
  my $center_file = $_[1];
  my $knn_final   = $_[2];

  my %gl = ();
  open( GL, $gl_file );
  map { my $key = $_; chomp($key); $gl{$key} = 1 } <GL>;
  close(GL);

  open( FC, $center_file );
  my @centers = <FC>;
  chomp(@centers);
  close(FC);

  foreach my $idx ( 0 .. $#centers ) {
    my @tmp = split( " ", $centers[$idx] );
    $centers[$idx] = \@tmp[ 0 .. ( $knn_final - 1 ) ];
  }

  my @newCenters = ();

  foreach my $idx ( 0 .. $#centers ) {

    my %hits = ();
    map { $hits{$_} = 1 if ( !exists $gl{$_} ) } @{ $centers[$idx] };

    next if ( ( keys %hits ) == 0 );

    my @newCent = @{ $centers[$idx] };

    #    my @newCent = ($centers[$idx]->[0]);

    #  push(@newCent,keys %hits);
    #
    #  if (@newCent >= $knn_final){
    #   @newCent = @newCent[0..($knn_final-1)]
    # } else{

#   my $ti = 1;
#  while (@newCent < $knn_final){
#    push(@newCent,$centers[$idx]->[$ti]) if (!exists $hits{$centers[$idx]->[$ti]});
#    $ti++;
#  }

    #}
    push( @newCenters, \@newCent );
  }

  open( OUT, ">$center_file.gl_filter" );

  foreach my $cent (@newCenters) {
    print OUT join( " ", @{$cent} ) . "\n";
  }

  close(OUT);

  #  system("cp $center_file $center_file.orig");
  #  system("cp $center_file.gl_filter $center_file");
}

## start or check status of a specific job
##
## the error behavior is as follows:
## if the full job is not finished yet (= number of .error + number of finished
##   files is less than number of tasks found in task.job_num) and job is not
##   in %jobs_active anymore, then we try to resubmit job.
##   This means for wait=1 we have to
##   start MASTER again, for wait=0 the main loop for stages 6-8 restarts job
## JOB.ERROR files is used for normal case that the full job produced some error
## and we want to avoid a re-submitting (both for wait=1 and wait = 0)
sub job_call {
  my ( $JOB_NAME, $JOB_script, $JOB_opts, $JOB_num_tasks, $JOB_sge_errdir, $JOB_sge_use, $JOB_sge_logdir, $JOB_times_file, $JOB_wait, $QSUB_opts, $THREADS_NUM, $JOB_UUID, $finished_last ) = @_;

  ## (finished?, error?, jobs_finished, jobs_all, $job_resubmit)
  my $sge_status = [ 0, 1, 0, 0, 0 ];

  if ( !-e "$JOB_sge_logdir/task.submitted" && !-e "$JOB_sge_logdir/JOB.ERROR" ) {
    SUBSECTION("submit $JOB_NAME - tasks = $JOB_num_tasks - use SGE = $JOB_sge_use");
    ## master start time
    system("\\rm -f $JOB_sge_logdir/*");
    system("echo  \`date\` > $JOB_times_file.master");
    system( "echo " . time . " >> $JOB_times_file.master" );
    my $err = job_submit( $JOB_script, $JOB_opts, $JOB_num_tasks, $JOB_sge_errdir, $JOB_sge_use, $JOB_sge_logdir, $QSUB_opts, $THREADS_NUM, $JOB_UUID );

    if ($err) {
      print "Error during job call/submission of $JOB_script!\n";
      end_handler(1);
    }

    system("touch $JOB_sge_logdir/task.submitted") if ( !$err );
  }

  if ( -e "$JOB_sge_logdir/task.submitted" ) {

    my @time_file    = readpipe("cat $JOB_times_file.master");
    my $time_started = $time_file[1];
    chomp($time_started);

    $sge_status = job_check( $JOB_sge_logdir, $JOB_NAME, $JOB_wait, $THREADS_NUM, $time_started );

    print "$JOB_NAME JOB status: finished=" . $sge_status->[0] . " error=" . $sge_status->[1] . " jobs " . $sge_status->[2] . "/" . $sge_status->[3] . "\n"
      if ( $in_verbose || $sge_status->[0] > 0 || $sge_status->[1] > 0 || $sge_status->[2] != $finished_last || $sge_status->[4] > 0 );

    if ( $sge_status->[0] == 1 && $sge_status->[1] == 0 ) {
      ## master end time
      system("echo  \`date\` >> $JOB_times_file.master");
      system( "echo " . time . " >> $JOB_times_file.master" );
      ## job times
      procSGEtimes( $JOB_sge_logdir, "$JOB_times_file.job" );
      system("rm $JOB_sge_logdir/*; rm -r $JOB_sge_logdir");
    }
  }

  ## resubmit job if MASTER is restarted
  if ( $JOB_wait && $sge_status->[4] ) {
    print "\nSGE/Thread/Local job ($JOB_NAME) did not finished correctly!\n\n";
    print "\nIf MASTER/JOB was just interrupted, then call MASTER again!\n\n";
    end_handler(1);
  }

  if ( $JOB_wait && ( -e "$JOB_sge_logdir/JOB.ERROR" || $sge_status->[1] ) ) {
    print "\nSGE/Loal job ($JOB_NAME) did not finished correctly!\n\n";
    print "To restart that phase please remove file:\n'$JOB_sge_logdir/JOB.ERROR'\n\n";
    print "call to remove: 'rm $JOB_sge_logdir/JOB.ERROR'\n\nGRAPHCLUST stopped.\n\n";
    end_handler(1);
  }

  return $sge_status;
}

################################################################################
## sub for time summary of one job, either sge or local job
## assumes a file "times" in $sge_log_dir
## each line contains "real:user" time, all lines are summed up and avg time is computed
## jobs, total time, avg time are written ti $out_times_file
sub procSGEtimes {
  my $sge_log_dir    = $_[0];
  my $out_times_file = $_[1];

  my @times_file_a = <$sge_log_dir/times>;
  my $times_file   = $times_file_a[0];

  my $sge_jobs = readpipe("cat $sge_log_dir/task.jobs");
  chomp($sge_jobs);

  open( TIMES, "$times_file" );
  my @times = <TIMES>;
  close(TIMES);
  chomp(@times);

  my @total = ();
  my @user  = ();

  foreach my $time (@times) {
    if ( $time =~ /^(\S+):(\S+)$/ ) {
      push( @total, $1 );
      push( @user,  $2 );
    }
  }

  my $avg_total  = 0;
  my $avg_user   = 0;
  my $total_time = 0;
  my $user_time  = 0;
  map { $total_time += $_ } @total;
  map { $user_time  += $_ } @user;
  $avg_total = $total_time / scalar(@total) if (@total);
  $avg_user  = $user_time / scalar(@user)   if (@user);

#  $avg_total = $avg_total + ( $sge_jobs - @total ) * $avg_total if ( @total != $sge_jobs );
#  $avg_user  = $avg_user +  ( $sge_jobs - @user ) * $avg_user   if ( @user != $sge_jobs );

  open( OUT, ">$out_times_file" );
  print OUT "jobs $sge_jobs jobs_finished " . @total . "  real $total_time avg_real $avg_total  user $user_time avg_user $avg_user\n";
  close(OUT);
}

sub job_check {
  my $sge_log_dir    = $_[0];
  my $section_prefix = $_[1];
  my $wait           = $_[2];
  my $use_threads    = $_[3];
  my $time_started   = $_[4];

  ## get number of tasks of job
  open( LOG, "$sge_log_dir/task.jobs" ) or die "Error in job_check! Cannot open $sge_log_dir/task.jobs. Exit...\n\n";
  my $sge_jobs = <LOG>;
  close(LOG);
  chomp $sge_jobs;
  die "Error! Unknown number of SGE jobs in $sge_log_dir/task.jobs. Exit...\n\n" if ( !defined($sge_jobs) || !$sge_jobs );

  ## file TASKID is created by qsubCall (SGE|LOCAL)
  open( ID, "$sge_log_dir/task.id" ) or die "Error! Cannot open $sge_log_dir/task.id! Exit...\n\n";
  my @t = split( " ", <ID> );
  chomp @t;
  close(ID);
  my $job_ID   = $t[0];
  my $job_UUID = $t[1];
  die "Error! Unknown TASKID in $sge_log_dir/task.id. Exit...\n\n" if ( !defined($job_ID) || !$job_ID );

  ## check if threads/sge-job were killed but we think they are active (e.g. ctrl-c last run)
  if ( !exists $jobs_active{$job_UUID} ) {
    print "\n$section_prefix: Job is not active anymore! Try to resubmit...\n";
    system("rm -f $sge_log_dir/task.submitted");
    system("rm -f $sge_log_dir/*.error");
    system("rm -f $sge_log_dir/*.started");
    system("rm -f $SGE_ERR_DIR/$job_UUID.*");
    ## set status to resubmit job
    return [ 0, 0, 0, $sge_jobs, 1 ];
  }

  if ( $wait && $use_threads > 1 ) {
    SUBSECTION("$job_ID:$section_prefix:  Wait for $use_threads ...");
  } elsif ( $wait && $job_ID ne "LOCAL" ) {
    SUBSECTION("$section_prefix:  Wait for SGE-JOB ($job_ID) ...");
    print "\nNote: You can stop script here if not all jobs are finished yet...\n\n";
  }

  my @jobs_fin = ();
  my @jobs_err = ();
  my $fin_last;

  do {

    #system("ls $sge_log_dir/*.started &>/dev/null");
    #system("ls $sge_log_dir/*.finished &>/dev/null");

    #@jobs_fin = <$sge_log_dir/*.finished>;
    @jobs_fin = ();
    map{push(@jobs_fin,1) if (system("stat $sge_log_dir/task-$_.finished 1>/dev/null 2>/dev/null") == 0)} 1..$sge_jobs;

    @jobs_err = <$sge_log_dir/*.error>;

    if ( $wait && !defined($fin_last) || ($wait && @jobs_fin > $fin_last) ) {
      $| = 1;
      print "Jobs finished: " . (@jobs_fin) . "/$sge_jobs  - time elapsed " . sprintf( "%.0f", ( time() - $time_started ) ) . " seconds.\r";
      $fin_last = @jobs_fin;
    }

    sleep(5) if ($wait);

  } while ( $wait && ( scalar(@jobs_fin) + scalar(@jobs_err) ) < $sge_jobs );

  my $sge_error = 0;

  if ( scalar(@jobs_fin) + scalar(@jobs_err) >= $sge_jobs ) {

    print "Jobs finished: " . ( @jobs_fin + @jobs_err ) . "/$sge_jobs  - time elapsed " . sprintf( "%.0f", ( time() - $time_started ) ) . " seconds.\r" if ( !$wait );

    ## job is now not active anymore
    delete $jobs_active{$job_UUID};

    ## combine stdout/stderr from array jobs
    system("rm -f $SGE_ERR_DIR/$job_UUID.out; rm -f $SGE_ERR_DIR/$job_UUID.err");
    system("cat $SGE_ERR_DIR/$job_UUID.out.* >> $SGE_ERR_DIR/$job_UUID.out 2> /dev/null");
    system("cat $SGE_ERR_DIR/$job_UUID.err.* >> $SGE_ERR_DIR/$job_UUID.err 2> /dev/null");

    system("rm -f $SGE_ERR_DIR/$job_UUID.out.* ");
    system("rm -f $SGE_ERR_DIR/$job_UUID.err.* ");

    open( ERR, "$SGE_ERR_DIR/$job_UUID.err" );
    my @err = <ERR>;
    close(ERR);

    system("rm -f $sge_log_dir/JOB.ERROR");
    if (@jobs_err) {
      print "ERROR! SGE job finished with exit code != 0!\n";
      system("cat $sge_log_dir/*.error >> $sge_log_dir/JOB.ERROR");
      system("rm -f $sge_log_dir/*.finished");
      $sge_error = 1;
    }

    if (@err) {
      SECTION("ERROR OUTPUT BEGIN $section_prefix");
      print @err;
      SECTION("ERROR OUTPUT END $section_prefix");
      print "ERROR! SGE job has generated error output\n";
      system("cat $SGE_ERR_DIR/$job_UUID.err >> $sge_log_dir/JOB.ERROR");
      system("rm -f $sge_log_dir/*.finished");
      $sge_error = 1;
    }

    if ( !$sge_error ) {
      print "$section_prefix -> SGE job ($job_ID) finished without error\n";
    }

    system("rm -f $sge_log_dir/task.submitted");
    system("cat $sge_log_dir/*.finished > $sge_log_dir/times 2> /dev/null");
    system("rm -f $sge_log_dir/*.finished");
  }

  my $sge_finished = 0;
  $sge_finished = 1 if ( @jobs_fin == $sge_jobs );

  my @ret = ( $sge_finished, $sge_error, scalar(@jobs_fin), $sge_jobs, 0 );

  return \@ret;
}

## set MASTER variables from CONFIG
sub setConfig {

  $OPTS_nspdk             = $CONFIG{OPTS_nspdk};
  $OPTS_nspdk_centers     = $CONFIG{OPTS_nspdk_centers};
  $OPTS_fasta2shrep_gspan = $CONFIG{OPTS_fasta2shrep_gspan};
  $OPTS_RNAfold           = $CONFIG{OPTS_RNAfold};
  $OPTS_RNAplfold         = $CONFIG{OPTS_RNAplfold};

  $nspdk_knn_center = $CONFIG{nspdk_knn_center};
  $nspdk_nhf        = $CONFIG{nspdk_nhf};
  $nspdk_nhf_max    = $CONFIG{nspdk_nhf_max};
  $nspdk_nhf_step   = $CONFIG{nspdk_nhf_step};
  $nspdk_fcs        = $CONFIG{nspdk_fcs};

  $GLOBAL_plfold_minlen = $CONFIG{GLOBAL_plfold_minlen};
  $GLOBAL_group_size    = $CONFIG{GLOBAL_group_size};
  $GLOBAL_iterations    = $CONFIG{GLOBAL_iterations};
  $GLOBAL_hit_blacklist_overlap = $CONFIG{GLOBAL_hit_blacklist_overlap};
  $GLOBAL_num_clusters  = $CONFIG{GLOBAL_num_clusters};

  $evaluate      = $CONFIG{evaluate};

  $input_blastclust_id  = $CONFIG{input_blastclust_id};
  $input_blastclust_len = $CONFIG{input_blastclust_len};
  $input_seq_min_length = $CONFIG{input_seq_min_length};
  $input_add_revcompl   = $CONFIG{input_add_revcompl};
  $input_win_size       = $CONFIG{input_win_size};
  $input_win_shift      = $CONFIG{input_win_shift};
}

sub cleanStage {
  my $stage = $_[0];

  if ( $stage <= 8 ) {
    system("rm -f -R $SGE_ERR_DIR/cmsearcher.*");
    system("rm -f -R $EVAL_DIR/times/time.stage.8.*");
    system("rm -f -R $EVAL_DIR/stage8.*");
  }

  if ( $stage <= 7 ) {
    system("rm -f -R $SGE_ERR_DIR/alignCenter.*");
    system("rm -f -R $EVAL_DIR/cluster/bestClusters.scores.*");
    system("rm -f -R $EVAL_DIR/cluster/bestClusters.qual.*");
    system("rm -f -R $EVAL_DIR/times/time.stage.7.*");
    system("rm -f -R $EVAL_DIR/stage7.*");
  }

  if ( $stage <= 6 ) {
    print "Cleanup stage 6 Cluster dir $CLUSTER_DIR...\n";
    system("rm -f -R $CLUSTER_DIR/*");
    system("rm -f -R $SGE_ERR_DIR/rnaclustAlignRange.*");
    system("rm -f -R $SGE_ERR_DIR/alignCenter.*");
    system("rm -f -R $EVAL_DIR/times/time.stage.6.*");

  }

  if ( $stage <= 5 ) {
    print "Cleanup stage 5 centers (data.svector.fast_cluster) in $SVECTOR_DIR...\n";
    system("rm -f $SVECTOR_DIR/centers.DONE");
    system("rm -f $SVECTOR_DIR/*.data.svector.fast_cluster*");
    system("rm -f $EVAL_DIR/svector/*");
  }

  if ( $stage <= 4 ) {
    print "Cleanup stage 4 sparse vector (data.svector) in $SVECTOR_DIR...\n";
    system("rm -f $SVECTOR_DIR/*.data.svector.DONE");
    system("rm -f $SVECTOR_DIR/*.data.svector.*");
    system("rm -f $SVECTOR_DIR/*.data.svector");
    system("rm -f $SVECTOR_DIR/*.blacklist");
    system("rm -f $in_ROOTDIR/*.DONE");
  }

  if ( $stage <= 3 ) {
    print "Cleanup stage 3 sparse vector groups in $SVECTOR_DIR...\n";
    system("rm -f -R $SVECTOR_DIR/svector.*");
    system("rm -f -R $SVECTOR_DIR/SGE_log");
    system("rm -f -R $SGE_ERR_DIR/gspanGroups.NSPDK.sge.*");
  }

  if ( $stage <= 2 ) {
    print "Cleanup stage 2 gspan groups in $GSPAN_DIR...\n";
    system("rm -f -R $GSPAN_DIR/groups.DONE");
    system("rm -f $GSPAN_DIR/*.group.gspan*");
  }

  if ( $stage <= 1 ) {
    print "Cleanup stage 1 gspan graphs in $GSPAN_DIR...\n";
    system("rm -f -R $GSPAN_DIR/*");
    system("rm -f -R $SGE_ERR_DIR/*");
  }

}


## returns sorted array of center-idx (1..$num_clusters) which finally produced a model
sub foundClusters {
  my $round = $_[0];

  my $part = GraphClust::read_partition("$EVAL_DIR/partitions/$round.soft");

  my %clusters = ();

  map { $clusters{$1} = 1 if ( $_->[3] =~ /^$round\.(\d+)$/ ) } @{$part};

  my @res = sort { $a <=> $b } keys %clusters;

  return \@res;
}

## write blacklist for all already clustered seqs, found by nspdk, alignment, cmsearch
sub makeBlacklist {
  my $curr_ci   = $_[0];
  my $bl_name   = $_[1];
  my $bl_min_ol = $_[2];

  my @blacklist = ();

  ## not necessary !?

  my $forInf = "";
  $forInf = " --infernal_1_1 " if ($newInf);
  system_call( "$BIN_DIR/glob_results.pl --root $in_ROOTDIR --all --summary-only --last-ci " . ( $curr_ci - 1 ), $in_verbose, $forInf );
#  system_call( "perl $BIN_DIR/glob_results.pl --root $in_ROOTDIR --all --summary-only --last-ci " . ( $curr_ci - 1 ), $in_verbose, $newInf );

  my $fragsDATA = GraphClust::read_fragments("$FASTA_DIR/$DATA_prefix.names");
  my $finalPart = GraphClust::read_partition("$RESULTS_DIR/partitions/final_partition.soft");

  my $finalHits = [];
  map { push( @{$finalHits}, $_->[0] ) } @{$finalPart};
  my $fragsCMSEARCH = GraphClust::list2frags($finalHits);

  ## we blacklist for all hits fragments on reverse strand of hit as well!
  ## this is probably the only situation where we can ignore strand in
  ## GraphClust::fragment_overlap(); all other overlap checks are strand aware!
  ## $bl_min_ol can be larger as we base all overlap on shorter fragment in ol check
  my $ignore_strand = 1;
  my $overlaps = GraphClust::fragment_overlap( $fragsCMSEARCH, $fragsDATA, $bl_min_ol, $ignore_strand );

  foreach my $ol ( @{$overlaps} ) {
    push( @blacklist, $fragsDATA->[ $ol->[1] ]->{VALUE} );
  }

  my @uniq_bl = unique(@blacklist);

  my %gl = ();
  if ( -e "$FASTA_DIR/$DATA_prefix.greylist" ) {
    open( GL, "$FASTA_DIR/$DATA_prefix.greylist" );
    map { my $key = $_; chomp($key); $gl{$key} = 1 } <GL>;
    close(GL);
  }

  open( OUT, ">$bl_name" );
  foreach my $key ( sort { $a <=> $b } @uniq_bl ) {
    ## write all out but do not blacklist greylist ids, can be used in non-greylist mode as well
    print OUT $key . "\n" if ( !exists $gl{$key} );
  }
  close(OUT);

  print "\nblacklist size (infernal hits) for round $curr_ci: " . @uniq_bl . "\n";
}

sub job_submit {
  my ( $CALL_script, $CALL_OPTS, $SGE_JOBS, $SGE_ERRDIR, $SGE_USE, $SGE_LOGDIR, $QSUB_OPTS, $THREADS, $JOB_UUID ) = @_;

  # !! check refactor !!!

  ## 0 is script path 1-end are options to script
  my $call_str = join( " ", @{$CALL_OPTS} );

  SUBSECTION("QSUBCALL $call_str") if ($in_verbose);

  system("\\mkdir -p $SGE_ERRDIR");
  system("\\mkdir -p $SGE_LOGDIR");

  system("rm -f $SGE_ERRDIR/$JOB_UUID.*");

  my $err = 0;

  ## number of tasks for this job
  system("echo $SGE_JOBS > $SGE_LOGDIR/task.jobs");

  if ($SGE_USE) {
    ## use SGE to submit jobs

    ## qsub path is added in sge_call
    my $qsub_call = "qsub -S /bin/bash $QSUB_OPTS -t 1-$SGE_JOBS -o $SGE_ERRDIR/$JOB_UUID.out.\\\$TASK_ID.\\\$JOB_ID -e $SGE_ERRDIR/$JOB_UUID.err.\\\$TASK_ID.\\\$JOB_ID ";
    $qsub_call .= "$CALL_script $CURRDIR $SGE_LOGDIR 0 $call_str";

    my $qsub_call_ret = sge_call($qsub_call);

    my $ret_str = join( "", @{$qsub_call_ret} );

    ## check if qsub was correctly done
    if ( $ret_str !~ /Your job-array (\d+\.\d+-\d+:\d+) \(/ || $ret_str =~ /usage: ssh/ ) {
      $err = 1;
      print $ret_str. "\n";
    } else {
      $ret_str =~ /Your job-array (\d+)/;
      system("echo \"$1 $JOB_UUID\" > $SGE_LOGDIR/task.id");
      print "Your SGE job $1 was successfully submitted!\n";
      $jobs_active{$JOB_UUID} = [ $JOB_UUID, "SGE", $SGE_LOGDIR, $1 ];
    }

  } else {
    ## Threaded or unthreaded job

    system("echo \"LOCAL $JOB_UUID\" > $SGE_LOGDIR/task.id");

    if ( $THREADS > 1 ) {

      foreach my $t ( 1 .. $SGE_JOBS ) {

        my $call = "cd $CURRDIR; $CALL_script $CURRDIR $SGE_LOGDIR $t $call_str 1>$SGE_ERRDIR/$JOB_UUID.out.$t 2>$SGE_ERRDIR/$JOB_UUID.err.$t;";

        $q->enqueue($call);
      }

    } else {

      ## print"cd $CURRDIR; $SGE_SCRIPT $CURRDIR $SGE_LOGDIR $call_str 1>$SGE_ERRDIR/$script_name.oLOCAL.0 2>$SGE_ERRDIR/$script_name.eLOCAL.0\n";
#system("cd $CURRDIR; for i in \$(seq 1 $SGE_JOBS); do echo -n TASK \$i started...; $SGE_SCRIPT $CURRDIR $SGE_LOGDIR \$i $call_str 1>$SGE_ERRDIR/$JOB_UUID.out.\$i 2>$SGE_ERRDIR/$JOB_UUID.err.\$i; echo TASK \$i/$SGE_JOBS finished ; done") == 0
#  or $err = 1;

      system("cd $CURRDIR");
      foreach my $t ( 1 .. $SGE_JOBS ) {
        print "TASK $t started...";
        system("$CALL_script $CURRDIR $SGE_LOGDIR $t $call_str 1>$SGE_ERRDIR/$JOB_UUID.out.$t 2>$SGE_ERRDIR/$JOB_UUID.err.$t") == 0
          or $err = 1;
        print "TASK $t/$SGE_JOBS finished err=$err\n";
        if ($err){
          print "\n";
          system("cat $SGE_ERRDIR/$JOB_UUID.err.$t");
          print "\n";
          die;
        }
      }
    }

    $jobs_active{"$JOB_UUID"} = [ $JOB_UUID, "THREAD", $SGE_LOGDIR, $SGE_JOBS ];

  }    ## else threaded/unthreaded

  return $err;
}

## do a real SGE call
## call_str needs to start with the sge tool WITHOUT any path!
## eg. "qsub ...", the sge_bin_path is added here
##
## from MASTER_GraphCLust.pl command line options we use:
##
##  --no-ssh       -> if set, then no ssh call is used but we call sge command
##                    directly on current host
##
## from global config file we use here:
##
## 'SGE_HOSTNAME'  -> is the host which allows to start/enqueue new sge jobs
## 'PATH_SGE'      -> is the path to the sge binaries like qstat, qsub etc
##
## from current environment we use following variables:
##
## 'SGE_ROOT'      -> needs to be set that sge can be used at all
##        Attention: we assume that this env-variable is set at current host
##        we do not check this variable at host 'SGE_HOSTNAME' due to missing env
##        after ssh login!
## 'whoami'        -> to get username for SSH login to SGE_HOSTNAME
##
## $CURRDIR        -> we use current directory (where MASTER was started)
##                    also on SGE_HOSTNAME if ssh is used,
##                    i.e. filesystem needs to be transparent to sge
##                    (usually via NFS)
sub sge_call {
  my $call_str = $_[0];

  my $CURRUSER = `whoami`;
  chomp $CURRUSER;

  my $SGE_HOSTNAME = $CONFIG{SGE_HOSTNAME};
  my $SGE_BINPATH  = $CONFIG{PATH_SGE};
  my $SGE_ROOT     = $ENV{"SGE_ROOT"};

  if ($in_verbose) {
    print "curr user    : $CURRUSER :\n";
    print "SGE_HOSTNAME : $SGE_HOSTNAME :\n";
    print "SGE_PATH     : $SGE_BINPATH :\n";
    print "use SSH      : " . ( !$in_NO_SGE_SSH ) . " :\n";
    print "SGE_ROOT     : $SGE_ROOT :\n";
    print "SGE CALL     : $call_str :\n";
  }

  die "Parameter 'SGE_HOSTNAME' or 'PATH_SGE' not set correctly! Please provide hostname which is able to
    call qsub <JOBSCRIPT>\n Set option either in config file or edit GraphClust_pathconf.pm\n"
    if ( $SGE_HOSTNAME eq "false" || $SGE_BINPATH eq "false" || $SGE_BINPATH eq "" || $SGE_HOSTNAME eq "" );

  my @sge_call_ret = ();

  my $sge_call = "$SGE_BINPATH/$call_str ";

  if ( !$in_NO_SGE_SSH ) {
    @sge_call_ret =
      readpipe("ssh $CURRUSER\@$SGE_HOSTNAME 'export SGE_ROOT=$SGE_ROOT; cd $CURRDIR; $sge_call ' 2>&1 ");
  } else {
    @sge_call_ret = readpipe("$sge_call 2>&1 ");
  }
  return \@sge_call_ret;
}

sub sge_qstat {

  my $CURRUSER = `whoami`;
  chomp $CURRUSER;

  my $qstat_call = "qstat ";
  $qstat_call .= " | grep $CURRUSER";
  my @qstat_call_ret;

  my $qstat_call_ret = sge_call($qstat_call);

  my %curr_user_jobs = ();
  foreach my $qs ( @{$qstat_call_ret} ) {
    my @t = split( " ", $qs );
    $curr_user_jobs{ $t[0] } = \@t;
  }
  return \%curr_user_jobs;
}

sub call_thread {
  my $call = $_[0];

  my $err = 0;
  system($call) == 0 or $err = 1;

  return $err;
}

sub end_handler {
  GraphClust::SUBSECTION(" GraphClust Exit handling");
  print STDERR "signal ", $_[0], " caught.\n" if ( $_[0] );

  if (@workers) {
    print "Kill threads to finish...\n";
    my $num = $q->pending();
    $q->dequeue_nb($num) if ( $num > 0 );

    ## Wait for workers to end
    $_->detach() for @workers;
  } elsif ( $in_USE_SGE || $in_SGE_KILL ) {
    open( SGE, ">$in_ROOTDIR/joblist.sge" );
    foreach my $job ( keys %jobs_active ) {
      print SGE join( " ", @{ $jobs_active{$job} } ) . "\n";
    }
    close(SGE);
  } elsif ( -e "$in_ROOTDIR/joblist.sge" ) {
    system("rm -f $in_ROOTDIR/joblist.sge");
  }

  chdir($CURRDIR) if ($CURRDIR);
  File::Temp::cleanup();
  GraphClust::SECTION( $CONFIG{VERSION_INFO} . " FINISHED" );

  exit;
}

## old makeBlacklist

#  foreach my $tab_file ( glob("$CMSEARCH_DIR/*.cmsearch/*.tabresult") ) {
#
#    $tab_file =~ /\/(\d+)\.(\d+)\.cmsearch\//;
#    next if ( $1 >= $curr_ci );    ## curr_ci=2 -> use only lists from round 1
#
#    print "$tab_file\n";
#    my $hits = GraphClust::read_CM_tabfile_ext( $tab_file, $cm_min_bitscore, $cm_max_eval, $cm_bitscore_sig, "BL" );
#
#    my $overlaps = GraphClust::fragment_overlap( $hits, $frags_clus, $bl_min_ol );    ##
#
#   # map { print "frags order: $_ " . $frags_clus->[$_]->{SEQID} . "\n"; } 0 .. $#{$frags_clus};
#
#    ## is it good to not blacklist too small final clusters?
#    ## next if (@{$hits} < $min_cluster_size);
#
#    foreach my $ol ( @{$overlaps} ) {
#
#      # my @ent = split("~",$ol_key); ## 0: hit-key 1:frag-key
#      # print $ol_key." ".$overlaps->{$ol_key}." ".$frags_clus->{$ent[1]}."\n";
#
#      push( @blacklist, $frags_clus->[ $ol->[1] ]->{VALUE} );
#    }
#  }

### write additional temp blacklist for one round if some centers do not lead to a model
### to avoid that exactly the same nspdk-centers appear again and again
#sub makeTempBlacklist {
#  my $curr_ci      = $_[0];
#  my $found_models = $_[1];
#  my $bl_name      = $_[2];
#
#  my @models = @{$found_models};
#
#  my %noModel = ();    ## blacklisted seq-keys
#  my $idx     = 0;
#  open( SVEC, "$SVECTOR_DIR/data.svector.fast_cluster." . ( $curr_ci - 1 ) ) or die "Cannot open last SVECTOR $SVECTOR_DIR/data.svector.fast_cluster." . ( $curr_ci - 1 ) . "! Exit...\n\n";
#  my $mod_c = 0;       ## model found count, similar to initial size of @models (?)
#  my $mod_b = 0;
#
#  ## foreach dense center
#  while ( my $line = <SVEC> ) {
#
#    $idx++;
#    chomp($line);
#
#    ## center round.idx lead to model
#    if ( defined( $models[0] ) && $models[0] == $idx ) {
#      shift(@models);
#      $mod_c++;
#      next;
#    }
#
#    ## no model for round.idx
#    $mod_b++;
#    next if ( $mod_b > $blacklist_centers_max );
#    next if ( $blacklist_center_frac <= 0 );
#
#    print "no model " . ( $curr_ci - 1 ) . ".$idx: $line\n";
#    my @ent = split( " ", $line );
#
#    $blacklist_center_frac = 1 if ( $blacklist_center_frac > 1 );
#    my $num_bl = int( @ent * $blacklist_center_frac );
#    print "num blacklist entries in center: " . $num_bl . "\n";
#    map { $noModel{$_} = ( $curr_ci - 1 ) . ".$idx" } @ent[ 0 .. ( $num_bl - 1 ) ];
#
#  }
#  close(SVEC);
#
#  ## check if exist no_model-blacklist from last round, yes, then read in old blacklist
#  if ( -e "$bl_name" && $mod_c > 0 ) {
#
#    print "temp blacklisted key this round: " . ( keys %noModel ) . "\n";
#
#    if ( !$blacklist_curr_only ) {
#      open( IN, "$bl_name" );
#      print "blacklist_curr_only = 0 -> read and add old temp blacklist...\n";
#      while ( my $line = <IN> ) {
#        chomp($line);
#        my @ent = split( " ", $line );
#        $noModel{ $ent[0] } = $ent[1];
#      }
#      close(IN);
#    }
#    print "temp blacklisted all: " . ( keys %noModel ) . "\n";
#  }
#
#  my $num_nokeys = keys %noModel;
#  if ( $mod_c > 0 && $mod_c <= $blacklist_max_models ) {
#    print "write new temp blacklist...\n";
#    open( OUT, ">$bl_name" );
#    foreach my $key ( sort { $a <=> $b } keys %noModel ) {
#      print OUT $key . "\n";    # . " " . $noModel{$key} . "\n";
#    }
#    close(OUT);
#  } else {
#    print "write EMPTY temp blacklist...\n";
#    system("rm -f $bl_name; touch $bl_name");
#  }
#}

## old eval svm  code from the beginning of iteration

## create special blacklist for evaluation
#if ( $CI > 1 && $evaluate_svm ) {
#  print "Create special SVM evaluation blacklist...\n";
#  ## calc this round$CI blacklist from subntree scores from last round
#  GraphClust::evalSVM_blacklist( "$FASTA_DIR", "$EVAL_DIR/stage7.cluster_qual.final.scores.round" . ( $CI - 1 ), "$SVECTOR_DIR/blacklist.svm.round$CI" );
#  ## merge final bl from last round and new one
#  system( "cat $SVECTOR_DIR/blacklist.svm.round$CI $SVECTOR_DIR/data.svector.blacklist." . ( $CI - 1 ) . " > $SVECTOR_DIR/data.svector.blacklist.$CI" );
#}

#use strict;
#use warnings;
#
#use threads       ;#qw( async );
#use Thread::Queue qw( );
#
#my $num_workers    = 5;
#my $num_work_units = 10;
#
#my $q = Thread::Queue->new();
#
## Create workers
#my @workers;
#for (1..$num_workers) {
#   push @workers, async {
#      while (defined(my $unit = $q->dequeue())) {
#         print("$unit\n");
#      }
#   };
#}
#
## Create work
#for (1..$num_work_units) {
#   $q->enqueue($_);
#}
#
## Tell workers they are no longer needed.
#$q->enqueue(undef) for @workers;
#
## Wait for workers to end
#$_->join() for @workers;

#  $script_name =~ /.*\/([^\/]+)$/;

################################################################################
################################################################################
## POD / manpage

=head1 NAME

GraphClust - structural clustering of RNA

=head1 SYNOPSIS

MASTER_GraphClust.pl --root <DIR> --fasta <FILE> [options]

See MASTER_GraphClust.pl --help or MASTER_GraphClust.pl --man for more
details!

=head1 DESCRIPTION

B<GraphCLust> performs a sequence-structure based clustering of RNA
sequences. Especially it is possible to apply GraphClust to verly large
datasets with thounsands of RNA sequences.

=head1 OPTIONS

=head2 Required Options for new data

=over 4

=item  B<--root <DIR>>

Target directory. All output files are written to this directory.
If GraphClust was interrupted at least a root dir should be given to
continue GraphClust.

=item  B<--fasta <FILE>>

File which contain your input sequences. There are several parameters
which influences the preprocessing of your data (see README).
The default case splits your input in 150nt fragments with 50%
overlap. BlastClust is used to filter out near identical fragments.

=back

=head2 Options

=over 4

=item  B<--config <FILE>>

Use the specified config file instead of default settings for all
paramters for a new root directory. See examples/config.default_global
for an example config file. NOTE: if the root directory already exists
and contains a config file than always the existing config file is
used! Please edit this file if you ant to change parameters.

=item  B<--sge>

Turn on parallel computing mode by using SGE (Sun Grid Engine). Please
see README to set all SGE paramters correctly. You need a working SGE
installation!

=item  B<--threads <NUM>>

Run GraphClust with the provided number of threads. Please use either
threads or SGE!

=item B<--results>

Invoke results stage for all found clusters so far. Results can be
found in ROOT/RESULTS/.

=back

=head2 Special Options

=over 4

=item  B<--gl <FILE>>

Use file as "greylist". Computes only clusters which inlcude given
sequences.

=item  B<--sge-kill>

Tries to kill all queued or running SGE jobs which were submitted for
the given root directory.

=item  B<--no-ssh>

Do not use ssh to submit jobs in SGE mode. Instead we call qsub
directly on current host.

=item  B<-v> B<--verbose>

Print more output about what is going on.

=back

=head2 Getting Help

=over 4

=item  B<--help>

Brief help message

=item  B<--man>

Detailed help message

=back

=head1 INFORMATION

Please see the README file in the source package for more informtion.

For updates and news please vist our webpage:
http://bioinf.uni.freiburg.de/Software/GraphClust


=head1 AUTHORS

Steffen Heyne, Fabrizio Costa, Dominic Rose and Rolf Backofen

=head1 REFERENCES

Steffen Heyne, Fabrizio Costa, Dominic Rose, and Rolf Backofen.
GraphClust: alignment-free structural clustering of local RNA secondary
structures. Bioinformatics, 28 no. 12 pp. i224-i232, 2012.

Please cite our paper if you use GraphClust!

=head1 CONTACT

heyne@informatik.uni-freiburg.de

=cut
