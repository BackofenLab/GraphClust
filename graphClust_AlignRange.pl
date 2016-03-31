#!/usr/bin/perl -w
#
# process the pairs in the given range of indices/codes
# i.e.: align with locarna
# results are stores in files from-to.ar.bz2
#
use strict;
use threads;
use Thread::Semaphore;
use Getopt::Long;
use Pod::Usage;

use FindBin;
use lib "$FindBin::Bin";
use GraphClust;

use vars qw(%CONFIG);
*CONFIG = \%GraphClust::CONFIG;

sub align_range;
sub parse_ranges;
sub decode;
sub min;
sub create_arfile;

my $help = 0;
my $man;
my $max_ar_size     = 50000;       # maximal number of alignments stored in one archive
my $range           = '';
my $locarna_options = "";
my $locarna_path    = "locarna";
my $tgtdir          = "";
my $dpdir           = "";
my $numberThreads   = 1;
my $verbose         = 0;
my $in_root_dir;

GetOptions( "root-dir:s"     => \$in_root_dir,
			"verbose:i"      => \$verbose,
			"help"           => \$help,
			"man"            => \$man,
			"range:s"        => \$range,
			"locarna-opts=s" => \$locarna_options,
			"locarna-path=s" => \$locarna_path,
			"tgtdir:s"       => \$tgtdir,
			"dpdir:s"        => \$dpdir,
			"cpu:i"          => \$numberThreads
) || pod2usage(2);

if ($man) { pod2usage( -exitstatus => 0, -verbose => 2 ); }
if ( !$in_root_dir || $help || $dpdir eq "" || $tgtdir eq "" || $range !~ /^\d+\-\d+$/ ) { pod2usage(1); }

%GraphClust::CONFIG = readConfigFile( "$in_root_dir/config", 1 );
$locarna_options  = $CONFIG{OPTS_locarna_paligs};

 if ( $locarna_options =~ /local/ || $locarna_options =~ /normalized/ ){
  $locarna_options .= " --local-output --pos-output ";  
 }

#print "used loc_opts:$locarna_options\n";

my $outdir = $tgtdir . "/paligs";
align_range($range);

sub align_range {
  my @jobs                = parse_ranges(@_);
  my $file_from           = 0;
  my $file_to             = min( $#jobs, $max_ar_size - 1 );
  my $arfilename : shared = "";

  system_command("mkdir -p $outdir");

  create_arfile( \$arfilename, $file_from, $file_to, \@jobs );

  ## Call locarna - distribute calls on the available CPUs
  my @threads;
  my $i         = 0;
  my $semaphore = Thread::Semaphore->new();
  for ( my $idx = 0 ; $idx <= $#jobs ; $idx++ ) {
	my $decode_n;
	my $decode_m;
	if ( $idx - $file_from + 1 > $max_ar_size ) {    ## we need a new arfile
	  system_command("cd $outdir; bzip2 -f $arfilename");
	  $file_from = $idx;
	  $file_to = min( $#jobs, $file_from + $max_ar_size - 1 );
	  create_arfile( \$arfilename, $file_from, $file_to, \@jobs );
	}

	if ( $numberThreads == 1 ) {
	  call_locarna( $semaphore, $idx, $decode_n, $decode_m, $arfilename, @jobs );
	} else {

	  $threads[$i] = threads->create( \&call_locarna, $semaphore, $idx, $decode_n, $decode_m, $arfilename, @jobs );
	  $i++;

	  ## Wait for threads to exit
	  if ( $#threads + 1 == $numberThreads ) {
		for ( $i = 0 ; $i < $numberThreads ; $i++ ) {
		  my $r = $threads[$i]->join();
		}
		$i       = 0;
		@threads = ();    # Free for new threads
	  }
	}
  }

  ## Wait until remaining threads are ready
  for ( $i = 0 ; $i <= $#threads ; $i++ ) {
	my $r = $threads[$i]->join();
  }
  @threads = ();

  system_command("cd $outdir; bzip2 -f $arfilename");

}

sub call_locarna {
  my ( $semaphore, $idx, $decode_n, $decode_m, $arfilename, @jobs ) = @_;

  decode( $jobs[$idx], \$decode_n, \$decode_m );
  my $n           = $decode_n;
  my $m           = $decode_m;
  my $clustalfile = "/tmp/align_range_$n.$m.$$";

  #	if($locarna_options =~ /--kbest/){
  #	    system_command("$locarna_path $locarna_options $dpdir/$n $dpdir/$m > $clustalfile");
  #	}else{
  #print "$locarna_path $locarna_options $dpdir/$n $dpdir/$m > $clustalfile\n";
  system_command("$locarna_path $locarna_options $dpdir/$n $dpdir/$m > $clustalfile");

  #system_command("$locarna_path $locarna_options --clustal $clustalfile $dpdir/$n $dpdir/$m >/dev/null");
  #	}
  ## Printing to AR file must be synchronized!
  # only one thread at a time has access to AR
  $semaphore->down();
  open( AR, ">>$outdir/$arfilename" ) or die "ERROR ($0): Cannot open file '$arfilename' for appending!\n";
  print AR "\@$jobs[$idx]";
  open( CL, "$clustalfile" );
  my $line = <CL>;
  chomp $line;
  $line =~ s/CLUSTAL W --- LocARNA - Local Alignment of RNA --- //;
  print AR " ";
  print AR "$line\n";

  while ( $line = <CL> ) {
	chomp $line;
	print AR "$line\n";
  }
  close CL;
  unlink "$clustalfile";
  close AR;
  $semaphore->up();
  ##

  return 1;
}

sub create_arfile {
  my ( $ref_arfilename, $file_from, $file_to, $ref_jobs ) = @_;
  $$ref_arfilename = "$ref_jobs->[$file_from]-$ref_jobs->[$file_to].ar";
  if ( -e "$outdir/$$ref_arfilename" ) {
	rename "$outdir/$$ref_arfilename", "$outdir/$$ref_arfilename.bkp";
  }
  return;
}

sub min {
  my ( $x, $y ) = @_;
  return ( $x <= $y ) ? $x : $y;
}

sub decode {
  my ( $idx, $ref_decode_n, $ref_decode_m ) = @_;
  my $n;
  my $m;

  $n = 1;
  while ( $idx > $n ) {
	$idx -= $n;
	$n++;
  }
  $m = $idx;
  $n++;
  $$ref_decode_n = $n;
  $$ref_decode_m = $m;
  return;
}

sub parse_ranges {
  my @ranges = @_;
  my @list   = ();

  foreach my $range (@ranges) {
	if ( $range =~ /(\d+)-(\d+)/ ) {
	  my $fr = $1;
	  my $to = $2;
	  for ( my $i = $fr ; $i <= $to ; $i++ ) {
		push @list, $i;
	  }
	} else {
	  push @list, $range;
	}
  }
  return @list;
}

sub system_command {
  my ($sysCommand) = @_;

  if ($verbose) { print "CALL $sysCommand\n"; }
  system("$sysCommand") == 0
	or die "ERROR ($0): Could not call system command '$sysCommand'!\n\n";
  return 1;
}

__END__

=head1 NAME

rnaclustAlignRange.pl

=head1 SYNOPSIS

rnaclustAlignRange.pl [options] 

Options:

    --range      <i-j>				specifies pairwise alignments which will be calculated (e.g. --range 1-100)
						(required)

    --tgtdir     <dir>				target directory
						(required)

    --dpdir      <dir>                          Directory containing the dotplots
                                                (required)s

    --cpu        <n>                            Number of CPUs available on your machine. Calls to locarna will 
						be distributed equally between those CPUs. (default: 1)	
						(optional)

    --locarna-opts <"locarna options">          options passed directly to locarna (must be given as one string)
						(optional)

    --locarna-path <path>			path to locarna (default: )
						(optional)

    --help   					print this help message
						(optional)

    --man       	    	        	full documentation
						(optional)

=head1 DESCRIPTION

    Runs the pairwise locarna alignments of the dotplots in [--dpdir <dir>] specified by [--range <range>]. Output is written to [--tgtdir <tgtdir>]. 

=cut

