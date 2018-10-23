#!/usr/bin/env perl
use strict; 
use warnings;
use Getopt::Long;
use Pod::Usage;
use POSIX qw(ceil floor);
use Cwd qw(abs_path getcwd);
use File::Path;
use Data::Dumper;
use Array::Utils qw(:all);
use List::Util qw/ min max /;

if ($ARGV[0] eq "-help") {
    print "Usage: gc_res.pl\n";
    exit 0;
}
my $in_eval_mode;
my $in_verbose = 0;
my $in_root_dir;
$in_root_dir = "";

my ($myFasta, $myResDir,$clusName, $OPTS_locarna_model) = @ARGV;


## summary contains evaluation info for used partition

my @fa_scan = read_fasta_file($myFasta);

my $fa_top = alignTopResults( $myFasta,$myResDir,$clusName);


exit;

################################################################################
####subs






sub read_fasta_file {
  my ($file,$make_unique) = @_;
  my $FUNCTION = "read_fasta_file in Sequences.pm";

  my $id        = "";
  my $seqstring = "";
  my %fasta     = ();
  my %header    = ();
  my @order     = ();
  my $line      = "";
  my %meta      = ();
  my %seq_meta  = ();
  my $uniq_count = 0;

  open( IN_HANDLE, "<$file" ) || die "ERROR in $FUNCTION: " . "Couldn't open the following file in package Tool," . " sub read_fasta_file: $file\n";


  while ( $line = <IN_HANDLE> ) {
    chomp($line);

    # header (can contain one space after > symbol)
    if ( $line =~ /^\>\s?(\S+)\s*([\S*\s*]*)/ ) {
      if ($id) {
        if ( defined $fasta{$id} and ( $fasta{$id} ne $seqstring ) ) {
#          $uniq_count++;
#          $id .= "_$uniq_count";
#          print "Warning! Make Seq-id unique! now $id\n";
         die "ERROR in $FUNCTION: " . "multiple sequence id '$id', consider using function " . "read_fasta_with_nonunique_headers instead of read_fasta_file";
        }
        $seqstring =~ s/\s*//g;    ## do not allow spaces in sequence
        $fasta{$id} = $seqstring;
        $meta{$id}  = {%seq_meta};    ## anonymous hash reference

        $seqstring = "";
        undef(%seq_meta);
      }

      if ($make_unique){
        $uniq_count++;
        $id = $1."_$uniq_count";
      }else{
        $id = $1;
      }

      my $head = $2;
      $head = "" if ( !$head );
      $header{$id} = $head;
      push( @order, $id );
    } elsif ( $line =~ /(.+)\s+(#\S+)\s*$/ && $id ) {

      if ( exists $seq_meta{$2} ) {
        $seq_meta{$2} .= $1;
      } else {
        $seq_meta{$2} = $1;
      }

    } else {
      $seqstring .= $line if ($id);
    }
  }

  if ($id) {
    if ( defined $fasta{$id} and ( $fasta{$id} ne $seqstring ) ) {
      #$uniq_count++;
      #$id .= "_$uniq_count";
      #print "Warning! Make Seq-id unique! now $id\n";
      die "ERROR in $FUNCTION: " . "multiple sequence id '$id', consider using function " . "read_fasta_with_nonunique_headers instead of read_fasta_file";
    }
    $seqstring =~ s/\s*//g;    ## do not allow spaces in sequence
    $fasta{$id} = $seqstring;
    $meta{$id}  = +{%seq_meta};    ## anonymous hash reference
    $seqstring  = "";
    undef(%seq_meta);

  }

  return ( \%fasta, \@order, \%header, \%meta );
}




sub alignTopResults {
  my $fa_file   = $_[0];
  my $resDir   = $_[1];
  my $name     = $_[2];

  my $currPath =  getcwd ;
  

  my $useLocP = 1;

    mlocarna_center( $fa_file, "$resDir/locarna.$name", "", $useLocP );
    aln2alifold( "$resDir/locarna.$name/results/result.aln", $resDir, "" );
    system("cp $resDir/locarna.$name/results/result.aln cluster.aln");
    system("cp $resDir/locarna.$name/results/result.aln.ps cluster.aln.ps");
    system("cp $resDir/locarna.$name/results/result.aln.alirna.ps cluster.alirna.ps");


    system("gm convert cluster.aln.ps cluster.aln.png");
    system("gm convert cluster.alirna.ps cluster.alirna.png");

    system("mloc2stockholm.pl --split_input yes --con_struct $resDir/locarna.$name/results/result.aln.alifold -file $resDir/locarna.$name/results/result.aln");
    system("cmbuild -F $resDir/locarna.$name/results/result.aln.cm $resDir/locarna.$name/results/result.aln.sth");
    system("cp $resDir/locarna.$name/results/result.aln.sth result.aln.sth");
    # system("R-scape --outdir $resDir/locarna.$name/results/ $resDir/locarna.$name/results/result.aln.sth");
    # system("cp $resDir/locarna.$name/results/result.aln_1.R2R.sto.pdf  cluster.result.aln_1.R2R.sto.pdf");
    system("cp  $resDir/locarna.$name/results/result.aln.cm $resDir/cluster.$name.cm");

  return $fa_file;
}

sub mlocarna_center {
    my $fasta    = $_[0];
    my $dir      = $_[1];
    my $dpDir    = $_[2];
    my $use_locP = $_[3];

    my $loc_pp_dir = "$dir/input";
    system("mkdir -p $loc_pp_dir");

    my @fa = read_fasta_file($fasta);
    foreach my $key ( keys %{ $fa[0] } ) {
        system("ln -f -s $dpDir/$key $loc_pp_dir/$key")
          if ( -e "$dpDir/$key" && !$use_locP );

    }
    #my $OPTS_locarna_model = "-p 0.001 --max-diff-am 50 --tau 50  --max-diff 100 --alifold-consensus-dp --indel-open -400 --indel -200 --struct-weight 180";
    system(
"mlocarna $OPTS_locarna_model  --skip-pp --verbose --tgtdir $dir $fasta > $dir/locarna.out 2>&1"
    );
}


sub aln2alifold {
    my $aln_file  = $_[0];
    my $tmp_path  = $_[1];
    my $vrna_path = $_[2];

    ## alifold result to get consensus structure string for infernal and some nice pictures
    my $tmp_dir = "$tmp_path/alifold_$$";
    my $currDir = getcwd;

    system("mkdir -p $tmp_dir");

    chdir($tmp_dir);
    print("Change path to $tmp_dir\n");
    print("RNAalifold -r --noLP --color --aln $currDir/$aln_file 2>/dev/null\n");
    my @call_alifold = readpipe( "RNAalifold -r --noLP --color --aln $currDir/$aln_file 2>/dev/null" );
    print "call_alifold:",join("\n",@call_alifold),"\n";

    my $aln_cons_struct = $call_alifold[1];    ## store cons-struct
    print "aln_cons_struct = $aln_cons_struct \n";
    chomp($aln_cons_struct);
    $aln_cons_struct =~ /([\(\).]*)(.*)/;      ## extract cons-struct
    $aln_cons_struct = $1;
    open( CONS, ">$currDir/$aln_file.alifold" );
    print CONS $call_alifold[0];
    print CONS "$aln_cons_struct\n";
    print CONS $2;
    close(CONS);
  system("mv alirna.ps $currDir/$aln_file.alirna.ps");
  system("mv aln.ps $currDir/$aln_file.ps");
  chdir($currDir);

  #system("rm -R $tmp_dir");
}
