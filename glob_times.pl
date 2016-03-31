#!/usr/bin/perl -w
use warnings;
use strict;

use Getopt::Long;
use Cwd qw(abs_path getcwd);
use FindBin;
use lib "$FindBin::Bin";

use Array::Utils qw(:all);

use GraphClust;
use vars qw(%CONFIG);
*CONFIG = \%GraphClust::CONFIG;

my $rootDir;
my $results_file;
my $verbose;

my $merge_threshold = 0.51;

usage()
  unless GetOptions(
  "root=s"  => \$rootDir,
  "verbose" => \$verbose
  );

my $bin_dir = abs_path($FindBin::Bin);

################################################################################
## Load options from root-dir/config file
%GraphClust::CONFIG = readConfigFile( "$rootDir/config", 1 );
SUBSECTION("options loaded for glob_results.pl");
printConfig( \%CONFIG ) if ($verbose);

my $evaluate = $GraphClust::CONFIG{evaluate};
print "eval:$evaluate\n";
################################################################################

my @time_files = sort grep { $_ !~ /\~$/ } glob( "$rootDir" . "/EVAL/times/*" );

my $all_file = "$rootDir/EVAL/time_report";

my @master_times = ();
my @job_times    = ();

###############################################################################
## stage 0
foreach my $file ( grep { $_ =~ /time\.stage\.0/ } @time_files ) {

  my $s0t = read_time_syscall($file);

  my $str = "STAGE 0 0 0 preprocessing ";
  $str .= "MASTER " . join( " ", @{$s0t} ) . " - - $file";

  print $str. "\n";
  push( @master_times, $str );
}

#writeTimes($all_file,\@master_times,\@job_times);

###############################################################################
## stage 1

foreach my $file ( grep { $_ =~ /time\.stage\.1/ } @time_files ) {

  my $master = 0;
  $master = 1 if ( $file =~ /master/ );

  if ($master) {
    my $s1tm = read_time_master("$file");
    my $str  = "STAGE 1 0 0 fasta2gspan ";
    $str .= "MASTER " . ( @{$s1tm}[0] ) . " - - " . ( @{$s1tm}[1] ) . " " . ( @{$s1tm}[2] ) . " $file";
    push( @master_times, $str );
    print $str. "\n";
  } else {

    my $s1tj = read_time_job("$file");
    my $str  = "STAGE 1 0 0 fasta2gspan ";
    $str .= "JOB " . join( " ", @{$s1tj}[ 9, 11, 3 ] ) . " - - $file";
    push( @job_times, $str );
    print $str. "\n";
  }

}

writeTimes( $all_file, \@master_times, \@job_times );

###############################################################################
## stage 3

foreach my $file ( grep { $_ =~ /time\.stage\.3/ } @time_files ) {

  my $master = 0;
  $master = 1 if ( $file =~ /master/ );

  if ($master) {
    my $s1tm = read_time_master("$file");
    my $str  = "STAGE 3 0 0 NSPDK_feat ";
    $str .= "MASTER " . ( @{$s1tm}[0] ) . " - - " . ( @{$s1tm}[1] ) . " " . ( @{$s1tm}[2] ) . " $file";
    push( @master_times, $str );
    print $str. "\n";
  } else {

    my $s1tj = read_time_job("$file");
    my $str  = "STAGE 3 0 0 NSPDK_FEAT ";
    $str .= "JOB " . join( " ", @{$s1tj}[ 9, 11, 3 ] ) . " - - $file";
    push( @job_times, $str );
    print $str. "\n";
  }

}

writeTimes( $all_file, \@master_times, \@job_times );

###############################################################################
## stage 4

foreach my $file ( grep { $_ =~ /time\.stage\.4/ } @time_files ) {

  my $master = 0;
  $master = 1 if ( $file =~ /master/ );

  if ($master) {
    my $s1tm = read_time_master("$file");
    my $str  = "STAGE 4 0 0 SVECTOR_file ";
    $str .= "MASTER " . ( @{$s1tm}[0] ) . " - - " . ( @{$s1tm}[1] ) . " " . ( @{$s1tm}[2] ) . " $file";
    push( @master_times, $str );
    print $str. "\n";
  } else {

    my $s1tj = read_time_job("$file");
    my $str  = "STAGE 4 0 0 SVECTOR_file ";
    $str .= "JOB " . join( " ", @{$s1tj}[ 9, 11, 3 ] ) . " $file";
    push( @job_times, $str );
    print $str. "\n";
  }

}

writeTimes( $all_file, \@master_times, \@job_times );

###############################################################################
## stage 5

foreach my $file ( grep { $_ =~ /time\.stage\.5/ } @time_files ) {

  my $master = 0;
  $master = 1 if ( $file =~ /master/ );

  $file =~ /time\.stage\.5\.(\d+)/;
  my $round = $1;
  if ($master) {
    my $s1tm = read_time_master("$file");
    my $str  = "STAGE 5 $round 0  NSPDK_cluster ";
    $str .= "MASTER " . ( @{$s1tm}[0] ) . " - - " . ( @{$s1tm}[1] ) . " " . ( @{$s1tm}[2] ) . " $file";
    push( @master_times, $str );
    print $str. "\n";
  } else {

    my $s1tj = read_time_job("$file");
    my $str  = "STAGE 5 $round 0 NSPDK_cluster ";
    $str .= "JOB " . join( " ", @{$s1tj}[ 9, 11, 3 ] ) . " - - $file";
    push( @job_times, $str );
    print $str. "\n";
  }

}

writeTimes( $all_file, \@master_times, \@job_times );

###############################################################################
## stage 6
##           round idx  paligs maligs rnasoup procTree model_locarna cmbuild cmsearch
## stage 6-8     1   1

foreach my $file ( grep { $_ =~ /time\.stage\.6\./ } @time_files ) {

  my $master = 0;
  $master = 1 if ( $file =~ /master/ );
#  $master = 2 if ( $file =~ /\.0$/ );

  $file =~ /time\.stage\.6\.(\d+)\.(\d+)/;
  my $round = $1;
  my $idx   = $2;

  if ( $master == 1 ) {
    my $s1tm = read_time_master("$file");
    my $str = "STAGE 6 $round $idx loc_paligs ";
    $str .= "MASTER " . ( @{$s1tm}[0] ) . " - - " . ( @{$s1tm}[1] ) . " " . ( @{$s1tm}[2] ) . " $file";
    push( @master_times, $str );
    print $str. "\n";

  } elsif ( $master == 2 ) {
#    my $s1tm = read_time_syscall("$file");

    #my $s1tm = read_time_master("$file");
#    my $str = "STAGE 6 $round $idx loc_paligs ";
#    $str .= "JOB " . join( " ", @{$s1tm} ) . " - - $file";
#
    #$str .= "MASTER ".join(" ",@{$s1tm})." $file";
#    push( @job_times, $str );
#    print $str. "\n";
  } else {

    my $s1tj = read_time_job("$file");
    my $str  = "STAGE 6 $round $idx loc_paligs ";
    $str .= "JOB " . join( " ", @{$s1tj}[ 9, 11, 3 ] ) . " - - $file";
    push( @job_times, $str );
    print $str. "\n";
  }

}

writeTimes( $all_file, \@master_times, \@job_times );

###############################################################################
## stage 7
##           round idx  paligs maligs rnasoup procTree model_locarna cmbuild cmsearch
## stage 6-8     1   1

foreach my $file ( grep { $_ =~ /time\.stage\.7\./ } @time_files ) {

  my $master = 0;
  $master = 1 if ( $file =~ /master/ );
  #$master = 2 if ( $file =~ /\.0$/ );

  $file =~ /time\.stage\.7\.(\d+)\.(\d+)/;
  my $round = $1;
  my $idx   = $2;

  if ( $master == 1 ) {
    my $s1tm = read_time_master("$file");
    my $str  = "STAGE 7 $round $idx loc_paligs ";
    $str .= "MASTER " . ( @{$s1tm}[0] ) . " - - " . ( @{$s1tm}[1] ) . " " . ( @{$s1tm}[2] ) . " $file";
    push( @master_times, $str );
    print $str. "\n";
  } #elsif ( $master == 2 ) {
   # my $s1tm = read_time_syscall("$file");
   # my $str  = "STAGE 6 $round $idx loc_paligs ";
   # $str .= "JOB " . join( " ", @{$s1tm} ) . " - - $file";
   # push( @job_times, $str );
   # print $str. "\n";
 # }
   else {

    my $s1tj = read_time_job("$file");
    my $str  = "STAGE 7 $round $idx tree_model ";
    $str .= "JOB " . join( " ", @{$s1tj}[ 9, 11, 3 ] ) . " - - $file";
    push( @job_times, $str );
    print $str. "\n";
  }

}

writeTimes( $all_file, \@master_times, \@job_times );

###############################################################################
## stage 8
##           round idx  paligs maligs rnasoup procTree model_locarna cmbuild cmsearch
## stage 6-8     1   1

foreach my $file ( grep { $_ =~ /time\.stage\.8\./ } @time_files ) {

  my $master = 0;
  $master = 1 if ( $file =~ /master/ );

  $file =~ /time\.stage\.8\.(\d+)\.(\d+)/;
  my $round = $1;
  my $idx   = $2;

  if ($master) {
    my $s1tm = read_time_master("$file");
    my $str  = "STAGE 8 $round $idx loc_paligs ";
    $str .= "MASTER " . ( @{$s1tm}[0] ) . " - - " . ( @{$s1tm}[1] ) . " " . ( @{$s1tm}[2] ) . " $file";
    push( @master_times, $str );
    print $str. "\n";
  } else {

    my $s1tj = read_time_job("$file");
    my $str  = "STAGE 8 $round $idx cmsearch ";
    $str .= "JOB " . join( " ", @{$s1tj}[ 9, 11, 3 ] ) . " - - $file";
    push( @job_times, $str );
    print $str. "\n";
  }

}

writeTimes( $all_file, \@master_times, \@job_times );

GraphClust::makeCol($all_file);

makeSummary( "$rootDir/EVAL/time.summary", \@master_times, \@job_times );

################################################################################

sub makeSummary {
  my $file = $_[0];
  my $masT = $_[1];
  my $jobT = $_[2];

  my @t      = @{$masT};
  my @master = ();
  map { my @a = split( " ", $_ ); push( @master, \@a ) } @t;

  @t = @{$jobT};
  my @jobs = ();
  map { my @a = split( " ", $_ ); push( @jobs, \@a ) } @t;

  my $max_round = 0;

  map { $max_round = $_->[2] if ( $_->[2] > $max_round ) } @jobs;

  print "\n\nmax round $max_round\n";

  my %iterCenters = ();
  $iterCenters{0} = 0;

  map { $iterCenters{$_} = 0 } 1 .. $max_round;
  map { $iterCenters{ $_->[2] } = $_->[3] if ( $_->[2] > 0 && $_->[3] > $iterCenters{ $_->[2] } ) } @jobs;

  my @cols = ( "iter", "num_center", "sum_center", "gspan", "NSPDK_feat", "NSPDK_clus", "loc_paligs", "tree_model", "cmsearch", "sum_iter", "sum_all", "avg_per_center" );

  # my @row = (0) x @cols;
  # my @mat = (\@row) x ($max_round+1);
  # map{print join(":",@{$_})."\n"} @mat;

  my @mat = ();

  open( OUT, ">$file" );
  print OUT join( " ", @cols ) . "\n";
  my $sum_center = 0;

  foreach my $iter ( sort { $a <=> $b } keys %iterCenters ) {

    $sum_center += $iterCenters{$iter};

    my @row = ( $iter, $iterCenters{$iter}, $sum_center );

    ## gspan
    my @sel = grep { $_->[1] == 1 && $_->[2] == $iter } @jobs;
    if (@sel) {
      push( @row, $sel[0]->[6] );
    } else {
      push( @row, 0 );
    }

    ## nspdk_feat
    @sel = grep { $_->[1] == 3 && $_->[2] == $iter } @jobs;
    if (@sel) {
      push( @row, $sel[0]->[6] );
    } else {
      push( @row, 0 );
    }

    ## nspdk_clus
    @sel = grep { $_->[1] == 5 && $_->[2] == $iter } @jobs;
    if (@sel) {
      push( @row, $sel[0]->[6] );
    } else {
      push( @row, 0 );
    }

    ## loc_paligs
    @sel = grep { $_->[1] == 6 && $_->[2] == $iter } @jobs;
    if (@sel) {

      my $sum = 0;

      foreach my $s (@sel) {
        $sum += $s->[6];

        #    print $s->[6]."\n";
      }

      push( @row, $sum );
    } else {
      push( @row, 0 );
    }

    ## tree
    @sel = grep { $_->[1] == 7 && $_->[2] == $iter } @jobs;
    if (@sel) {

      my $sum = 0;

      foreach my $s (@sel) {
        $sum += $s->[6];

        #   print $s->[6]."\n";
      }

      push( @row, $sum );
    } else {
      push( @row, 0 );
    }

    ## cmsearch
    @sel = grep { $_->[1] == 8 && $_->[2] == $iter } @jobs;
    if (@sel) {

      my $sum = 0;

      foreach my $s (@sel) {
        $sum += $s->[6];

        #   print $s->[6]."\n";
      }

      push( @row, $sum );
    } else {
      push( @row, 0 );
    }

    my $sum = 0;
    map { $sum += $_ } @row[ 3 .. 8 ];
    push( @row, $sum );

    #$sum = 0;
    #map{$sum += $_} @row[3..8];
    #push(@row,$sum);

    push( @mat, \@row );
  }

  push( @{ $mat[0] }, $mat[0]->[9] );
  foreach my $idx ( 1 .. $#mat ) {
    my @row = @{ $mat[$idx] };

    my $sum = $mat[ $idx - 1 ][10] + $mat[$idx][9];
    push( @{ $mat[$idx] }, $sum );

    push( @{ $mat[$idx] }, sprintf( "%.1f", $mat[$idx][10] / $mat[$idx][2] ) );

  }

  foreach my $row_ (@mat) {
    print join( " ", @{$row_} ) . "\n";
    print OUT join( " ", @{$row_} ) . "\n";
  }

  close(OUT);
  GraphClust::makeCol($file);
}

sub writeTimes {
  my $file = $_[0];
  my $masT = $_[1];
  my $jobT = $_[2];

  open( ALL, ">$file" );

  print ALL "STAGE s i c NAME TYPE SERIAL_TIME AVG_JOB_TIME NUM_JOBS START END\n\n";

  my @t  = @{$masT};
  my @ts = ();
  map { my @a = split( " ", $_ ); push( @ts, \@a ) } @t;
  @ts = sort { $a->[1] <=> $b->[1] || $a->[2] <=> $b->[2] } @ts;
  foreach my $t (@ts) {
    print ALL join( " ", @{$t} ) . "\n";
    print join( " ", @{$t} ) . "\n";
  }

  @t  = @{$jobT};
  @ts = ();
  map { my @a = split( " ", $_ ); push( @ts, \@a ) } @t;
  @ts = sort { $a->[2] <=> $b->[2] || $a->[3] <=> $b->[3] || $a->[1] <=> $b->[1] } @ts;
  foreach my $t (@ts) {
    print ALL join( " ", @{$t} ) . "\n";
    print join( " ", @{$t} ) . "\n";
  }

  close(ALL);
}

sub read_time_master {
  my $tf = $_[0];

  my @cont = readpipe("cat $tf");
  chomp(@cont);

  my $time = ( $cont[3] - $cont[1] );
  $cont[0] =~ s/ /\_/g;
  $cont[2] =~ s/ /\_/g;
  my @ret = ( $time, $cont[0], $cont[2] );

  return ( \@ret );

}

sub read_time_master_old {
  my $tf = $_[0];

  my $cont = readpipe("cat $tf");
  chomp($cont);
  my @ent = split( " ", $cont );

  my $time = ( $ent[2] );
  my @ret  = ($time);

  return ( \@ret );
}

sub read_time_job {
  my $tf = $_[0];

  my $cont = readpipe("cat $tf");
  chomp($cont);
  my @ent = split( " ", $cont );

  return ( \@ent );
}

sub read_time_syscall {
  my $tf = $_[0];

  my $cont = readpipe("cat $tf");
  chomp($cont);
  my @ent = split( " ", $cont );

  #print $cont.":\n";
  my @ret = ( $ent[6], $ent[7], $ent[1] );
  return \@ret;
}
