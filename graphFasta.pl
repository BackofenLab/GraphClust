#!/usr/bin/perl

## GraphClust input preprocessing script
##
## Author: Steffen Heyne (heyne@informatik.uni-freiburg.de)
##
## This script creates the final input for GraphCLust based on config parameters
##

use strict;
use Getopt::Long;
use Pod::Usage;
use POSIX qw(ceil floor);
use Cwd qw(abs_path getcwd);
use File::Path;
use List::Util 'shuffle';
use List::Util qw/ min max /;

use FindBin;
use lib "$FindBin::Bin";
use GraphClust;

use vars qw(%CONFIG);
*CONFIG = \%GraphClust::CONFIG;

my $tgtdir;
my $in_fasta;
my $in_prefix      = "data";
my $SEQPREFIX      = "SEQ";
my $add_rc         = 0;
my $max_length     = 0;        ## 0 means no splitting
my $min_seq_length = 5;
my $blastclust_len = 0.9;
my $blastclust_id  = 90;
my $max_N_stretch  = 15;       ## allowed length of stretch of N's
my $in_winShift    = 100;      ## window shift in percent
my $in_evaluate;
my $in_no_blastclust_filter;
my $in_grey_list;
my $in_nspdk_greylist_num  = 0;
my $in_add_rc_class_signal = 0;
my $in_class_files_only    = 0;
## grey list: add as file with orig seq-ids, fasta contains all seqs at once
## during file creation also create grey list with frag-idx which matching
## seqid from grey-list file

GetOptions(
  "tgtdir=s"   => \$tgtdir,
  "prefix=s"   => \$in_prefix,
  "revcompl"   => \$add_rc,
  "winsize=i"  => \$max_length,
  "bc-id=i"    => \$blastclust_id,
  "bc-len=f"   => \$blastclust_len,
  "min-len=i"  => \$min_seq_length,
  "winshift=i" => \$in_winShift,
  "evaluate"   => \$in_evaluate,
  "add-rc-sig" => \$in_add_rc_class_signal,
  "no-bc"      => \$in_no_blastclust_filter,
  "gl=s"       => \$in_grey_list,
  "gl-num=i"   => \$in_nspdk_greylist_num,
  "class-files-only" =>\$in_class_files_only,
);

$in_fasta = $ARGV[0];

################################################################################
## Load options from root-dir/config file
#%GraphClust::CONFIG = readConfigFile( "$rootDir/config", 1 );
#SUBSECTION("options loaded for glob_results.pl");

## this is a preconfigured option, so we can use default value
my $blastclust_path = $CONFIG{PATH_BLASTCLUST};
$blastclust_path = "" if ( !-e "$blastclust_path" || $blastclust_path eq "false" || $in_no_blastclust_filter);

## skip blastclust if it cannot be found
if ( !$in_no_blastclust_filter && (!$blastclust_path || !-e $blastclust_path."blastclust") ) {
  my $loc = `which blastclust`;
  chomp($loc);
  if ( !$loc || !-e $loc ){
    print "\nCannot find BlastClust! Skip BlastClust step!\n"; 
    $in_no_blastclust_filter = 1;
  } else{
    $blastclust_path = $loc;
  }
    
}

if ( !$tgtdir ) {
  $tgtdir = getcwd;
} else {
  $tgtdir = abs_path($tgtdir);
  File::Path::make_path($tgtdir);
  #system("rm -f -R $tgtdir/*");
}

die "Please provide fasta file!\n" if ( !$in_fasta );
die "Cannot find fasta file!\n "   if ( !-e $in_fasta );
die "Split length ($max_length) < 20! Please choos larger value! Exit...\n\n" if ( $max_length > 0 && $max_length < 20 );

#die "Split offset ($split_offset) < 5! Please choos larger value! Exit...\n\n" if ( $split_offset && $split_offset < 5 );
die "Shift = $in_winShift! Please choose value between 10 and 100! Exit...\n\n" if ( $in_winShift < 10 || $in_winShift > 100 );
$min_seq_length = 5 if ( $min_seq_length <= 5 );

## grey list
die "Please use either greylist file or a number of used grey list seqs! Exit...\n\n"
  if ( $in_nspdk_greylist_num && $in_grey_list );

$in_fasta = abs_path($in_fasta);
my @fa_path = split( "/", $in_fasta );
my $fa_name = $fa_path[$#fa_path];

$fa_name = $in_prefix if ($in_prefix);

$in_add_rc_class_signal = 1 if ($add_rc);

################################################################################
## read given fasta

## parse grey list

my $grey_href;

if ($in_grey_list) {

  my @tmp = GraphClust::read_fasta_with_nonunique_headers($in_grey_list);
  ## if fasta file then concat with in_fasta
  if ( @{ $tmp[0] } && @{ $tmp[1] } && length( $tmp[1]->[0] ) > 0 ) {

    system("cp $in_grey_list $tgtdir/orig.gl");
    system("cp $in_fasta $tgtdir/orig.in.fa");
    system("cat $tgtdir/orig.gl $tgtdir/orig.in.fa > $tgtdir/input.fa");
    $in_fasta              = "$tgtdir/input.fa";
    $in_nspdk_greylist_num = @{ $tmp[0] };

    print "\nInfo: Provided grey-list file $in_grey_list is fasta file! Treat " . $in_nspdk_greylist_num .
      " seqs as grey listed!\n\n";
    undef($in_grey_list);
  }
}

my @fa_in = GraphClust::read_fasta_file($in_fasta,1);

## check if header fasta file comes from GraphClust already, reuse ORIGID and ORIGHEAD
@fa_in = reuseFastaFormat( \@fa_in );

## normalize all input seqs
my %normalized_seqs = ();
foreach my $seq ( @{ $fa_in[1] } ) {

  $normalized_seqs{$seq} = uc( $fa_in[0]->{$seq} );
  $normalized_seqs{$seq} =~ s/\s*//g;
  $normalized_seqs{$seq} =~ tr/T/U/;
  $normalized_seqs{$seq} =~ s/[^AUGCN]/N/g;

}

#$fa_in[3] = $fa_in[0]; ## [3] is used for seq in scan fasta and seq length in classhash, but why???
$fa_in[0] = \%normalized_seqs;
my $num_seqs_in = @{ $fa_in[1] };

## todo: debug
#if ($in_evaluate) {
#  writeEvaluateFiles( \@fa_in, "$tgtdir" );
#}
#exit;

if ( $in_grey_list || $in_nspdk_greylist_num ) {
  $grey_href = getGreySeqIndices( \@fa_in, $in_grey_list, $in_nspdk_greylist_num );
}

## get array ref with locations
my $genome_locs = genomeLocations( \@fa_in );

## get array of array refs with all fragments, splitted on stretches with
## more than $max_N N's
## array-idx is same as in $fa_in
my $frags_splitN = splitMasked( \@fa_in, $max_N_stretch );

## todo: grey_href
writeFiles( \@fa_in, $frags_splitN, $genome_locs, $grey_href );

## split 50,100 with 50% shift
my @split_wins = ($max_length);

## we use only one window size currently, but we should evaluate
## to use multiple windows at the same time, but this complicated things
## with overlap checkings etc.
##
## use minwin, next smallerwin to check if fragment will be used
#foreach my $splitWin (@split_wins) {

my $splitWin = $split_wins[0];
my $frags_shift = splitLengthShift( $frags_splitN, $splitWin, $in_winShift, $max_N_stretch );

my $frags_keep = [ keys %{$frags_shift} ];
my $frags_reject;
if ( !$in_no_blastclust_filter ) {
  ( $frags_keep, $frags_reject ) = filterFasta_BlastClust( $frags_shift, $blastclust_id, $blastclust_len );
}

if ($grey_href) {
  addRejGreyFrags( $frags_keep, $frags_reject, $grey_href );
}

if ($add_rc) {
  addReverseCompl( $frags_shift, $frags_keep, $grey_href );
}

writeFrags( \@fa_in, $frags_shift, $frags_keep, $grey_href, "$tgtdir/$fa_name" );

if ($in_evaluate) {
  writeEvaluateFiles( \@fa_in, "$tgtdir" );
}

#}

print "\n#########################################################################\n";
print "graphFasta.pl STATS\n\n";
print "input seqs                                    : " . ( @{ $fa_in[1] } ) . "\n";
print "    fragments (after removing masked regions) : " . ( map { @{$_} } @{$frags_splitN} ) . "\n";
print "    fragments (after using windows/shift)     : " . ( keys %{$frags_shift} ) . "\n";
print "    fragments (after blastclust)              : " . ( @{$frags_keep} ) . "\n";
print "#########################################################################\n";

################################################################################
## subs

sub writeEvaluateFiles {
  my $fa  = $_[0];
  my $out = $_[1];

  ## count all, with overlapping frags/signals
  my %class_count = ();
  my %class_names = ();

  my @class_frags = ();
  $class_count{ "0" } = 0;
  foreach my $idx ( 1 .. @{ $fa->[1] } ) {

    my $id   = $fa->[1]->[ $idx - 1 ];
    my $head = $fa->[2]->{$id};

    ## multiple signal info in header
    if ( $head =~ /CLASS_MULTI\s+(\S+)/ ) {

      my @sigs = split( ";", $1 );

      foreach my $sig (@sigs) {

        my @ent = split( "#", $sig );
        next if ( @ent != 4 );

        $class_count{ $ent[0] }++;
        $class_names{ $ent[0] } = $ent[1];

        ## class 0 is noise/unknown class
        next if ( $ent[0] == 0 );

        my $sigFrags = classSignal2frag( $idx, $ent[0], $ent[2], $ent[3], $in_add_rc_class_signal );
        push( @class_frags, @{$sigFrags} );

        #print CH $sig_str;
      }

    } elsif ( $head =~ /CLASSIDX\s+(\d+)/ ) {
      ## single signal
      my $sig_idx  = $1;
      my $sig_name = "UNKNOWN";

      if ( $head =~ /CLASS\s+(\S+)/ ) {
        $sig_name = $1;
      }

      $class_count{$sig_idx}++;
      $class_names{$sig_idx} = $sig_name;

      ## class 0 is noise/unknown class
      next if ( $sig_idx == 0 );

      my $sig_start = 1;
      my $sig_stop  = length( $fa->[0]->{$id} );

      if ( $head =~ /CLASS_SIGNAL\s+(\d+)-(\d+)/ ) {
        $sig_start = $1;
        $sig_stop  = $2;
      }

      my $sigFrags = classSignal2frag( $idx, $sig_idx, $sig_start, $sig_stop, $in_add_rc_class_signal );
      push( @class_frags, @{$sigFrags} );

      #   print CH $sig_str;

    } else {
      $class_count{"0"}++;
    }
  }    ## foreach seq

  ## signals used for cmsearch during stage 8
  ## use orig signals as we scan on input seqs
  ## VALUE = CLASSIDX
  ## note: class.hash.scan can contain overlapping frags, depend on user enconding in input fasta!,
  ## i.e. also signals for evaluation can overlap but we care about it
  ## class 0 seqs are just estimated, number of class 0 FRAGS could be higher!
  open( CH, ">$out/class.hash.scan" );
  foreach my $fr (@class_frags) {
    print CH $fr->{VALUE} . " " . $fr->{KEY} . "\n";
  }
  close(CH);

  ## read in fragments produced from input fasta
  ## VALUE = seqidx
  my $frags = GraphClust::read_fragments("$out/$in_prefix.map");

  ## annotate frags with CLASS
  ##
  ## overlap 0.05 is a bit weak, but not really important here
  ## each fragment get class of signal with best overlap
  ## if there are multiple signals which overlap, they are storedin CLASS_MULTI!
  ##
  ## signals on minus strand filtered out if not present in .map,
  ## i.e. minus strand seqs not wantened via option input_add_revcompl!
  ## option 'ignore strand' can be given to fragment_overlap, but this needs to be implemented
  ## for evaluate_cm_hits!
  ##
  ## @class_frags does not contain class 0 signals!

  ## not yet sorted
  @class_frags = sort { $a->{SEQID} cmp $b->{SEQID} } @class_frags;
  my ($hitlist_class) = GraphClust::evaluate_cm_hits( \@class_frags, $frags, 0.05 );

  ##todo: check multiple signals in class.hash?? what to do with them?
  my %class_frag_count = ();
  my %class_uniq_count = ();
  open( CH, ">$tgtdir/class.hash" );
  foreach my $frag ( @{$frags} ) {

    #print $frag->{KEY}." ".$frag->{VALUE}." ".$frag->{CLASS}." ".$frag->{CLASS_KEY}." ".$frag->{CLASS_OL}."\n";
    print CH $frag->{VALUE} . " " . $frag->{CLASS} . "\n";
    $class_frag_count{ $frag->{CLASS} }++;

    foreach my $key (@{$frag->{CLASS_MULTI}}){
      $class_uniq_count{$key->[1]} = $key->[0];
      #$class_frag_count{ $key->[0] }++;
    }

  }
  close(CH);

  ## class size info based on original seqs/signals and after blastclust
  open( CS, ">$tgtdir/class.size" );

  $class_names{"0"} = "UNKNOWN";

  foreach my $idx ( sort { $a <=> $b } keys %class_count ) {
    my $frag_count = 0;
    $frag_count = $class_frag_count{$idx} if (exists $class_frag_count{$idx});

    my $uniq_count = 0;
    $uniq_count = grep { $class_uniq_count{$_} eq "$idx" } keys %class_uniq_count;
    $uniq_count /= 2 if ($in_add_rc_class_signal);
    print CS $class_count{$idx} . " $idx " .$frag_count ." ".$uniq_count." ".$class_names{$idx} . "\n";
  }
  close(CS);

  GraphClust::makeCol("$tgtdir/class.size");
  system("cat $tgtdir/class.size | sort -rgk1 > $tgtdir/class.size.sort");
  system("cat $tgtdir/class.size | sort -rgk3 > $tgtdir/class.size.sort_frags");
}

sub reuseFastaFormat {
  my $fa = $_[0];

  my %newSeqs  = ();
  my @newOrder = ();
  my %newHeads = ();
  my %newMeta  = ();

  foreach my $id ( @{ $fa->[1] } ) {

    my $head     = $fa->[2]->{$id};
    my $new_id   = $id;
    my $new_head = $head;

    if ( $head =~ /^$SEQPREFIX\d+#\d+#\d+#\S\s+ORIGID\s+(\S+)\s+ORIGHEAD\s+(.*)$/ ) {
      $new_id   = $1;
      $new_head = $2;
    }# elsif ( $head =~ /ORIGID\s+(\S+)\s+(.*)$/ ) {
     # $new_id   = $1;
     # $new_head = $2;
     # }

    $newSeqs{$new_id} = $fa->[0]->{$id};
    push( @newOrder, $new_id );
    $newHeads{$new_id} = $new_head;
    $newMeta{$new_id}  = $fa->[3]->{$id};
  }

  return ( \%newSeqs, \@newOrder, \%newHeads, \%newMeta );
}

sub addRejGreyFrags {
  my $frags_keep = $_[0];
  my $frags_rej  = $_[1];
  my $gl         = $_[2];

  foreach my $key ( @{$frags_rej} ) {

    my @ent = split( "#", $key );
    $ent[0] =~ /^$SEQPREFIX(\d+)/;
    my $sIdx = $1;

    if ( exists $gl->{$sIdx} && $gl->{$sIdx}==1){
      push( @{$frags_keep}, $key );
    }

  }
}

## create hash with 0/1 for each seq-idx according to input fasta
## 1 = TRUE grey list member
## 0 = NO grey list member
## returns hash-ref
sub getGreySeqIndices {
  my $fa      = $_[0];
  my $gl_file = $_[1];
  my $gl_num  = $_[2];

  my %gl = ();

  if ($gl_file) {
    ## gl_file provided
    open( IN, "$gl_file" );
    my @tmp = <IN>;
    close(IN);

    foreach my $ent (@tmp) {
      chomp $ent;
      $ent =~ />?(\S+)\s*/; ## use first non-space word, ignore '>' if we have eg full fasta header lines
      my $key = $1;
      $gl{$key} = 1 if ($key);
    }

  } else {
    ## gl_num provided
    $gl_num = 1 if ( $gl_num < 0 || $gl_num > @{ $fa->[1] } );
    foreach my $num ( 1 .. $gl_num ) {
      $gl{$num} = 1;
    }

  }

  ## check if we only have numbers as keys in greylist
  ## if yes, then treat hash keys as seq-indices according to input fasta
  ## else treat keys as pattern which is matched against seq-id & header
  my $match_seq_idx = 1;
  map { $match_seq_idx = 0 if ( $_ !~ /^\d+$/ ) } keys %gl;

  ## greylist hash: $key->1 = seq_idx is in grey list, $key->0 = NOT
  if ( !$match_seq_idx ) {
    ## %gl contains patterns which we check against all seqs in $fa
    ## each key is a non-space word, see above
    my $pattern = join( "|", keys %gl );

    my %gl_tmp = ();

    foreach my $seq_idx ( 1 .. @{ $fa->[1] } ) {

      my $tgt = $fa->[1]->[ $seq_idx - 1 ];    ## seq-id
      $tgt .= $fa->[2]->{$tgt} if ( $fa->[2]->{$tgt} ); ## append header if we have one

      if ( $tgt =~ /$pattern/ ) { ## if we have a match use seq idx for greylist
        $gl_tmp{$seq_idx} = 1;
      } else {
        $gl_tmp{$seq_idx} = 0;
      }
    }

    %gl = %gl_tmp;

  } else {
    ## %gl only contains seq-indices of of $fa
    foreach my $seq_idx ( 1 .. @{ $fa->[1] } ) {
      $gl{$seq_idx} = 0 if ( !exists $gl{$seq_idx} );
    }

  }

  return \%gl;
}

sub addReverseCompl {
  my $fragsAll  = $_[0];
  my $fragsKeep = $_[1];
  my $gl_href   = $_[2];

  my @add_keys = ();
  foreach my $key ( @{$fragsKeep} ) {

    ## do not add revcompl if greylist fragment
    my $fr = GraphClust::str2frag($key);
    $fr->{SEQID} =~ /$SEQPREFIX(\d+)$/;
    next if ( $gl_href && exists $gl_href->{$1} && $gl_href->{$1} > 0 );

    my $rc = $fragsAll->{$key};
    $rc =~ tr/AUGC/UACG/;
    $rc = reverse($rc);
    my $rcKey = $key;
    $rcKey =~ s/#\+/#\-/;
    $fragsAll->{$rcKey} = $rc;
    push( @add_keys, $rcKey );
  }

  push( @{$fragsKeep}, @add_keys );
}

sub writeFiles {
  my $fa       = $_[0];
  my $f_splitN = $_[1];
  my $locas    = $_[2];
  my $gl_href  = $_[3];

  open( FA, ">$tgtdir/orig.fasta" );
  open( SC, ">$tgtdir/$in_prefix.fasta.scan" );
  open( LO, ">$tgtdir/data.locations" );
  open( FR, ">$tgtdir/fragments.splitMasked.list" );

  foreach my $idx ( 1 .. @{ $fa->[1] } ) {

    my $id = $fa->[1]->[ $idx - 1 ];

    print FA ">$id " . $fa->[2]->{$id} . "\n" . $fa->[0]->{$id} . "\n";

    ## meta information from fasta
    map { print FA $fa->[3]->{$id}->{$_} . " $_\n" } keys %{ $fa->[3]->{$id} };

    print SC ">$SEQPREFIX$idx ORIGID $id ORIGHEAD " . $fa->[2]->{$id} . "\n" . $fa->[0]->{$id} . "\n";
    print LO "$SEQPREFIX$idx $locas->[$idx-1]\n";

    my $pos_abs_end   = 0;
    my $pos_abs_start = 0;

    foreach my $seq_frag ( @{ $f_splitN->[ $idx - 1 ] } ) {
      $pos_abs_end += length($seq_frag);
      $pos_abs_start = $pos_abs_end - length($seq_frag) + 1;

      if ( length($seq_frag) < $min_seq_length ) {
        print FR "$SEQPREFIX$idx#$pos_abs_start#$pos_abs_end SMALLER_THAN_$min_seq_length\_NT\n";
        next;
      }

      if ( $seq_frag =~ /N{$max_N_stretch,}/ ) {
        print FR "$SEQPREFIX$idx#$pos_abs_start#$pos_abs_end MASKED_FRAGMENT\n";
        next;
      }

      print FR "$SEQPREFIX$idx#$pos_abs_start#$pos_abs_end USED_FRAGMENT\n";

    }
  }
  close(FA);
  close(SC);
  close(LO);
  close(FR);

  GraphClust::makeCol("$tgtdir/fragments.splitMasked.list");
}

sub classSignal2frag {
  my $seq_idx    = $_[0];
  my $class_idx  = $_[1];
  my $sig_start  = $_[2];
  my $sig_stop   = $_[3];
  my $sig_add_rc = $_[4];

  my $frag_str = "";
  my @frags    = ();

  $frag_str .= "$SEQPREFIX$seq_idx";

  if ( $sig_start < $sig_stop ) {

    $frag_str .= "#$sig_start#$sig_stop";
    push( @frags, $frag_str . "#+" );

    if ($sig_add_rc) {
      push( @frags, $frag_str . "#-" );
    }

  } else {

    $frag_str .= "#$sig_stop#$sig_start";
    push( @frags, $frag_str . "#-" );

    if ($sig_add_rc) {
      push( @frags, $frag_str . "#+" );
    }

  }

  map { $frags[$_] = GraphClust::str2frag( $frags[$_] ); $frags[$_]->{VALUE} = $class_idx } 0 .. $#frags;
#  print "FRAG " . $frags[0]->{KEY} . " " . $frags[0]->{VALUE} . "\n";

  return \@frags;
}

sub filterFasta_BlastClust {
  my $frags  = $_[0];
  my $bc_id  = $_[1];
  my $bc_len = $_[2];

  my $num_seqs_curr;
  my @seqs_keep = sort { ( $a =~ /^$SEQPREFIX(\d+)#/ )[0] <=> ( $b =~ /^$SEQPREFIX(\d+)#/ )[0] || ( $a =~ /[^#]+#(\d+)/ )[0] <=> ( $b =~ /[^#]+#(\d+)/ )[0] || $a cmp $b } keys %{$frags};
  my $iter = 0;

  my @seqs_reject = ();

  do {

    $iter++;

    my $si = 1;
    open( FA, ">$tgtdir/intermediate_bc.fasta" );
    foreach my $key (@seqs_keep) {
      print FA ">$key\n";
      print FA $frags->{$key} . "\n";
      $si++;
    }
    close(FA);

# system("cp $tgtdir/intermediate_bc.fasta $tgtdir/intermediate.start.fa") if ( $iter == 1 );

    #print "### BLASTCLUST seqs ###\n";
    my $bc_opts = "-p F -b T -W 20 -S $bc_id -L $bc_len -o $tgtdir/intermediate.bclust.$iter -i $tgtdir/intermediate_bc.fasta";
    print "#################################################################################\n";
    print "blastclust $bc_opts 2>/dev/null 1>/dev/null";
    system( $blastclust_path. "blastclust $bc_opts 2>/dev/null" );

    open( BC, "$tgtdir/intermediate.bclust.$iter" );
    my @clusters = <BC>;
    close(BC);
    chomp(@clusters);

    $num_seqs_curr = @seqs_keep;
    @seqs_keep     = ();

    open( CL, ">$tgtdir/intermediate.blast_clusters.$iter" );
    foreach my $clus (@clusters) {
      my @ent = split( " ", $clus );
      push( @seqs_keep, $ent[0] );
      push( @seqs_reject, @ent[ 1 .. $#ent ] ) if ( @ent > 1 );
      print CL $clus . "\n" if ( scalar(@ent) > 1 );
    }
    close(CL);

    @seqs_keep = shuffle(@seqs_keep);

    print ">>> blastclust  iter = $iter " . scalar(@seqs_keep) . " / $num_seqs_curr / " . ( keys %{$frags} ) . "\n";

  } while ( scalar(@seqs_keep) + 5 <= $num_seqs_curr );

  system("rm $tgtdir/intermediate_bc.fasta");
  system("rm $tgtdir/intermediate.bclust.*");
  return ( \@seqs_keep, \@seqs_reject );
}

sub writeFrags {
  my $fa        = $_[0];
  my $frags     = $_[1];
  my $fragsKeep = $_[2];
  my $gl        = $_[3];
  my $prefix    = $_[4];

  ## files with kept fragments
  open( NA,   ">$prefix.names" );
  open( SEQS, ">$prefix.fasta" );
  open( LEN,  ">$prefix.lens" );
  open( MAP,  ">$prefix.map" );

  my $seq_idx = 1;
  my @keep_sorted = sort { ( $a =~ /^$SEQPREFIX(\d+)#/ )[0] <=> ( $b =~ /^$SEQPREFIX(\d+)#/ )[0] || ( $a =~ /[^#]+#(\d+)/ )[0] <=> ( $b =~ /[^#]+#(\d+)/ )[0] || $a cmp $b } @{$fragsKeep};

  my $write_gl = 0;
  if ( keys %{$gl} ) {
    ## sort to get grey list frags first
    @keep_sorted = sort {
      $gl->{ ( $b =~ /^$SEQPREFIX(\d+)#/ )[0] } <=> $gl->{ ( $a =~ /^$SEQPREFIX(\d+)#/ )[0] }
        || ( $a =~ /^$SEQPREFIX(\d+)#/ )[0] <=> ( $b =~ /^$SEQPREFIX(\d+)#/ )[0]
        || ( $a =~ /[^#]+#(\d+)/ )[0] <=>       ( $b =~ /[^#]+#(\d+)/ )[0]
        || $a cmp $b
    } @{$fragsKeep};
    $write_gl = 1;
    open( GL, ">$prefix.greylist" );
  }

  foreach my $key (@keep_sorted) {

    my @ent = split( "#", $key );
    $ent[0] =~ /^$SEQPREFIX(\d+)/;
    my $sIdx    = $1;
    my $orig_id = $fa->[1]->[ $sIdx - 1 ];

    print NA "$seq_idx $key ORIGID $orig_id ORIGHEAD " . $fa->[2]->{$orig_id} . "\n";
    print SEQS ">$seq_idx $key ORIGID $orig_id ORIGHEAD " . $fa->[2]->{$orig_id} . "\n";
    print SEQS $frags->{$key} . "\n";

    ## add meta information, but avoid for splitted seqs
    map { print SEQS $fa->[3]->{$orig_id}->{$_} . " $_\n" if ( length( $frags->{$key} ) == length( $fa->[3]->{$orig_id}->{$_} ) ) } keys %{ $fa->[3]->{$orig_id} };

    print LEN length( $frags->{$key} ) . "\n";
    print MAP $seq_idx . " " . $key . "\n";

    print GL $seq_idx . "\n" if ( $gl->{$sIdx} == 1 );

    $seq_idx++;
  }

  close(NA);
  close(SEQS);
  close(LEN);
  close(MAP);
  close(GL) if ($write_gl);

  ## file with all frags, not blastclust filtered
  my @all_sorted = sort { ( $a =~ /^$SEQPREFIX(\d+)#/ )[0] <=> ( $b =~ /^$SEQPREFIX(\d+)#/ )[0] || ( $a =~ /[^#]+#(\d+)/ )[0] <=> ( $b =~ /[^#]+#(\d+)/ )[0] || $a cmp $b } keys %{$frags};
  open( IM, ">$prefix.fasta.all_frags" );

  foreach my $key (@all_sorted) {
    $key =~ /^$SEQPREFIX(\d+)#/;
    my $idx = $1;
    my $id  = $fa->[1]->[ $idx - 1 ];
    print IM ">$key ORIGID $id ORIGHEAD " . $fa->[2]->{$id} . "\n";
    print IM $frags->{$key} . "\n";
  }
  close(IM);
}

sub splitLengthShift {
  my $seqs     = $_[0];
  my $winSize  = $_[1];
  my $winShift = $_[2];
  my $max_N    = $_[3];

  my %usedSeqs  = ();
  my $shift_abs = 0;

  if ( $winSize > 0 ) {
    $shift_abs = floor( ( $winSize / 100 ) * $winShift );
    $shift_abs = 1        if ( $shift_abs < 1 );
    $shift_abs = $winSize if ( $shift_abs > $winSize );
  }
  foreach my $seq_idx ( 1 .. @{$seqs} ) {

    my $frags = $seqs->[ $seq_idx - 1 ];

    my $pos_abs_end   = 0;
    my $pos_abs_start = 0;

    #my $seqKey = "$SEQPREFIX$seq_idx";

    foreach my $frag_idx ( 1 .. @{$frags} ) {

      my $seq_frag = $frags->[ $frag_idx - 1 ];

      $pos_abs_end += length($seq_frag);
      $pos_abs_start = $pos_abs_end - length($seq_frag) + 1;

      next if ( length($seq_frag) < $min_seq_length );
      next if ( $seq_frag =~ /N{$max_N,}/ );

      ## N's at start end end we can ignore always
      ## we need only adjust $pos_abs_start by $start_cue later
      ## end cueing can be ignored
      $seq_frag =~ s/^([N]+)//;
      my $start_cue = 0;
      $start_cue = length($1) if ( length($1) > 0 );
      $seq_frag =~ s/[N]+$//;

## old split/shift if something went wrong with the new code below :-)
#      my $pos = 1;
#      while ( length($seq_frag) >= $min_seq_length ) {
#
#        my $length_used = $winSize;
#        $length_used = length($seq_frag) if ( $winSize == 0 || length($seq_frag) <= $winSize + $min_seq_length );
#
#        my $curr_seq = substr( $seq_frag, 0, $length_used );
#
#        next if ( length($curr_seq) < $min_seq_length );
#
#        my $abs_pos = ( $pos_abs_start + $pos - 1 ) . "#" . ( $pos_abs_start + $pos + length($curr_seq) - 2 );
#        my $key = $SEQPREFIX . $seq_idx . "#" . $abs_pos;
#
#        # print "key: $key len_frag " . length($seq_frag) . " len_curr " . length($curr_seq) . " abs: " . $abs_pos . " abs_f: $pos_abs_start-$pos_abs_end \n";
#        $usedSeqs{ $key . "#+" } = $curr_seq;
#
#        last if ( $winSize == 0 );
#        last if ( length($seq_frag) <= $winSize + $min_seq_length );
#
#        substr( $seq_frag, 0, $shift_abs, "" );
#        $pos += $shift_abs;
#
#      }

      my $end_pos = length($seq_frag);
      $end_pos = $winSize if ( $winSize > 0 && length($seq_frag) > $winSize );
      my $start_idx = 0;
      my $used_len  = length($seq_frag);

      do {

        ## check if current frag is last, if then use winSize length from end
        ## if too small overhang then shrink frag a bit

        if ($winSize) {

          ## overlap of two fragments, defined by $winShift
          my $overlap = $winSize - $shift_abs;
          ## tolerance for certain fragmets (within 10% length diff)
          my $special_len = $winSize * 0.1;

          ## idea: take always frags of $winSize, only for last frags
          ## use special length and startidx
          ## debug
          my $last_len = length($seq_frag) - $start_idx;
          if ( length($seq_frag) <= $winSize ) {
            ## case: frag is between $min_seq_length and $winSize
            $used_len = length($seq_frag);

            # print "case 0 last $last_len used $used_len\n";
          } elsif ( length($seq_frag) - $start_idx <= $winSize * 0.5 ) {
            ## case: last is very short, take at least $min_seq_length, or 60% of $winSize
            $used_len = max( $min_seq_length, $winSize * 0.6 );
            $start_idx = length($seq_frag) - $used_len;

            #print "case 1 last $last_len used $used_len\n";
          } elsif ( length($seq_frag) - $start_idx - $overlap <= $shift_abs + $special_len ) {
            ## case: last frag is in tolerance range +- 10% $winSize
            $used_len = max( $winSize, length($seq_frag) - $start_idx - $overlap );
            $start_idx = length($seq_frag) - $used_len;

            #print "case 2 last $last_len used $used_len seq $seq_idx\n";
          } #elsif ( length($seq_frag) - $start_idx <= $winSize + $special_len ) {
            #$used_len = length($seq_frag) - $start_idx;
            #print "case 3 last $last_len used $used_len\n";
            #}
          else {
            ## case: if last frag is long enough then take just $winSize
            $used_len = $winSize;

            #print "case 5 last $last_len used $used_len\n";
          }
        }

        my $curr_seq = substr( $seq_frag, $start_idx, $used_len );

        my $abs_pos = ( $pos_abs_start + $start_cue + $start_idx ) . "#" . ( $pos_abs_start + $start_idx + length($curr_seq) - 1 );
        my $key = $SEQPREFIX . $seq_idx . "#" . $abs_pos;
        $usedSeqs{ $key . "#+" } = $curr_seq;

        $start_idx += $shift_abs;

#print "$key len frag ".length($seq_frag)." start_idx $start_idx ".($start_idx-$shift_abs)." end_pos $end_pos spec $used_len\n";
        } while ( $start_idx - $shift_abs + $used_len < length($seq_frag) ); }
  }

  return \%usedSeqs;
}

sub genomeLocations {
  my $fa = $_[0];

  my @locs = ();

  foreach my $idx ( 1 .. @{ $fa->[1] } ) {

    if ( $fa->[2]->{ $fa->[1]->[ $idx - 1 ] } =~ /GENOME_LOC (\S+\.chr\S+:\d+-\d+:\S+)/ ) {
      push( @locs, $1 );
      next;
    }

    ## str = id." ".header
    my $str = $fa->[1]->[ $idx - 1 ] . " " . $fa->[2]->{ $fa->[1]->[ $idx - 1 ] };
    $str = lc($str);

    my $genome = "undef";
    my $chr;
    my $start;
    my $end;
    my $strand;
    my $loc = "MISS";

    if ( $str =~ /chr(\S+)[:_](\d+)[-_](\d+)/ ) {
      $chr   = $1;
      $start = $2;
      $end   = $3;

      $strand = "+" if ( $start <= $end );
      if ( $start > $end ) {
        $strand = "-";
        ## we use always start before end, even for minus strand
        my $tmp = $start;
        $start = $end;
        $end   = $tmp;
      }

      if ( $str =~ /chr\S+[:_]\d+[-_]\d+[:_](\S+)/ ) {
        $strand = $1;
      }

      if ( $str =~ /strand=([+-])/ && $start <= $end ) {
        $strand = $1;
      }

      if ( $str =~ /(\S+)[\._]chr\S+[:_]\d+[-_]\d+/ ) {
        $genome = $1;
      }

      $loc = "$genome.chr$chr:$start-$end:$strand";

    }

    push( @locs, $loc );

  }

  return \@locs;

}

sub splitMasked {
  my $fa    = $_[0];
  my $max_N = $_[1];

  my @all_frag = ();

  foreach my $idx ( 0 .. $#{ $fa->[1] } ) {

    #print " next seq $idx\n";

    my @seq_a = split( "", $fa->[0]->{ $fa->[1]->[$idx] } );

    my $curr_frag   = "";
    my $curr_N_frag = "";
    my @seq_frags   = ();

    while ( @seq_a > 0 ) {

      my $curr_nt = shift(@seq_a);

      while ( $curr_nt eq "N" ) {
        $curr_N_frag .= $curr_nt;
        $curr_nt = shift(@seq_a);
      }

      while ( $curr_nt ne "N" && @seq_a ) {
        $curr_frag .= $curr_nt;
        $curr_nt = shift(@seq_a);
      }

      ## use always last nt, also "N"
      $curr_frag .= $curr_nt if ( !@seq_a );

      if ( length($curr_N_frag) <= $max_N ) {

        ## append at last fragment, get last idx first
        $seq_frags[ 0 < $#seq_frags ? $#seq_frags : 0 ] .= ( $curr_N_frag . $curr_frag );

      } elsif ( length($curr_N_frag) > $max_N ) {
        push( @seq_frags, $curr_N_frag );
        push( @seq_frags, $curr_frag ) if ( length($curr_frag) > 0 );
      }

      $curr_N_frag = "";
      $curr_N_frag = $curr_nt if ( $curr_nt eq "N" );

      $curr_frag = "";

    }    ## while curr seq

    my $len = 0;
    map { $len += length($_) } @seq_frags;

    die "split error! Different length after splitting! Exit..\n\n" if ( $len != length( $fa->[0]->{ $fa->[1]->[$idx] } ) );

    #    print join("\n",@seq_frags)."\n";
    push( @all_frag, \@seq_frags );

  }    ## foreach seq

  return \@all_frag;
}

########

## orig.fasta
##  seqA1 > 1#ORIG
##        > 2#ORIG

## Repeat Mask -> fragments
##        >FRAG=1:1-100
##        >FRAG=1:101-200 REPEAT
##        >FRAG=2:1-200
##        >FRAG=3:201-300 REPEAT

## Splitting/overlaps / discard repeats

##       >1 FRAG=1:1-50:F
##       >2 FRAG=1:25-75:F
##       >3 FRAG=1:50-100:F
## add reverse?
##       >4 FRAG=1:1-50:R ORIGID <ID> ORIGPOS hg19.chr3:100000-111000:+ ORIGHEAD
##       >5 FRAG=1:25-75:R
##       >6 FRAG=1:50-100:R

## index1:  1 FRAG=SEQ1:1-50:R ORIGID <ID> ORIGPOS hg19.chr3:100000-111000 ORIGHEAD

## index2:  SEQ1  <ORIG-ID>

#foreach my $frag ( sort { ($a =~ /#(\d+)/)[0] <=> ($b =~ /#(\d+)/)[0] } grep { $_ =~ /$idx#/ } keys %all_frag ) {
#  print "L:" . $idx . " " . $frag . " :" . $all_frag{$frag} . ":\n";
#}

# 1 1 1 1 100 +
# 2 1 2 101 200 +
# 3 1 3 201 300 +
# 4 1 4 1 50 +
# 5 1 5 51 10 +
# 6 1 6 101 151 +
# 7 2 1 1 100 +
# 8 2 2 101 200 +
# 9 2 3 201 300 +

#  if ($add_rc) {
#    my $rc = $fa_ids{$seq};
#    $rc =~ tr/AUGC/UACG/;
#    $rc = reverse($rc);
#    $fa_rc{$seq} = $rc;
#  }

#open( NA,   ">$tgtdir/intermediate.names" );
#open( SEQS, ">$tgtdir/intermediate.fasta" );
#open( LEN,  ">$tgtdir/intermediate.lens" );
#my $seq_idx = 1;
#
#my $split_count = 0;
#my $short_count = 0;
#
#$max_length = $min_seq_length if ( $max_length && $max_length < $min_seq_length );
#
#foreach my $id ( @{ $fa[1] } ) {
#
#  my $seq_full = $fa[0]->{$id};
#
#  #print "\n";
#
#  my @seqs_Nsplit = split( /N{$max_N,}/, $seq_full );
#
#  #foreach my $i (0..$#seqs_Nsplit){
#  #  print "$i SPLIT $id ".length($seqs_Nsplit[$i])."\n";
#  #}
#
#  my $n_split_idx = 0;
#
#  foreach my $n_split_seq (@seqs_Nsplit) {
#
#    my $seq = $n_split_seq;
#
#    my $split_idx = 0;
#    $n_split_idx++;
#
#    next if ( length($seq) == 0 );
#
#    $short_count++ if ( length($seq) < 5 );
#
#    while ( length($seq) >= $min_seq_length ) {
#
#      my $seq_out = $seq;
#      if ( $max_length && length($seq) > $max_length + $min_seq_length ) {
#
#        my $max_length_used = $max_length;
#
#        $max_length_used = $split_offset if ( $split_offset > 0 && $split_idx == 0 );
#
#        $seq_out = substr( $seq, 0, $max_length_used, "" );
#        $split_count++;
#        $split_idx++;
#      } else {
#        $seq_out = substr( $seq, 0, length($seq), "" );
#        $split_idx++ if ( $split_idx > 0 );
#      }
#
#      #    print "$id orig " . length( $fa_ids{$id} ) . " rest " . length($seq) . " len new " . length($seq_out) . "\n";
#
#      print NA "$seq_idx $id SPLIT $n_split_idx-$split_idx STR +1 " . $fa[2]->{$id} . "\n";
#      print SEQS ">$seq_idx ORIGID $id SPLIT $n_split_idx-$split_idx STR +1 " . $fa[2]->{$id} . "\n";
#      print SEQS $seq_out . "\n";
#      print LEN length($seq_out) . "\n";
#
#      $seq_idx++;
#
#      #     if ($add_rc) {
#      #       my $header = $fa_heads_ref->{$id};
#
#      #$header =~ s/CLASSIDX\s+(\d+)/CLASSIDX 0/;
#      #       print NA "$seq_idx $id SPLIT $split_idx STR -1 " . $header . "\n";
#      #       print SEQS ">$seq_idx ORIGID $id SPLIT $split_idx STR -1 " . $header . "\n";
#      #       print SEQS $fa_rc{$id} . "\n";
#      #       print LEN length( $fa_rc{$id} ) . "\n";
#      #       $seq_idx++;
#      #     }
#    }    ## while length($seq) >= $min_seq_length
#  }    ## foreach split
#}    ## foreach seq

#system("cp $tgtdir/intermediate.fasta $tgtdir/intermediate.start.fa");

#close(FA);
#close(NA);
#close(SEQS);
#close(LEN);

#$seq_idx = $seq_idx - 1;
#my $num_seqs    = $seq_idx;
#my $bclust_seqs = 0;
#my @seqs_keep   = ();
#
#my $iter = 0;
#
#@fa = GraphClust::read_fasta_file("$tgtdir/intermediate.fasta");
#
#print "seqs after split $seq_idx\n";
#
#while ( scalar(@seqs_keep) != $num_seqs ) {
#
#  $iter++;
#
#  $num_seqs = @seqs_keep;
#
#  #print "### BLASTCLUST seqs ###\n";
#  #print "blastclust -p F -b F -W 32 -S $blastclust_id -L $blastclust_len -o $tgtdir/intermediate.bclust.$iter -i $tgtdir/intermediate.fasta\n";
#  system( $blastclust_path. "blastclust -p F -b T -W 32 -S $blastclust_id -L $blastclust_len -o $tgtdir/intermediate.bclust.$iter -i $tgtdir/intermediate.fasta 2>&1 >/dev/null" );
#
#  open( BC, "$tgtdir/intermediate.bclust.$iter" );
#  my @clusters = <BC>;
#  chomp(@clusters);
#  close(BC);
#
#  @seqs_keep = ();
#  open( CL, ">$tgtdir/intermediate.blast_clusters.$iter" );
#  foreach my $clus (@clusters) {
#    my @ent = split( " ", $clus );
#    push( @seqs_keep, $ent[0] );
#    print CL $clus . "\n" if ( @ent > 1 );
#  }
#  close(CL);
#
#  @seqs_keep = shuffle(@seqs_keep);
#
#  open( FA, ">$tgtdir/intermediate.fasta" );
#
#  foreach my $id (@seqs_keep) {
#    print FA ">$id " . $fa[2]->{$id} . "\n";
#    print FA $fa[0]->{$id} . "\n";
#  }
#  close(FA);
#
#  print ">>> blastclust  iter = $iter " . scalar(@seqs_keep) . " / $num_seqs / $seq_idx\n";
#}

#open( LEN, ">$tgtdir/$fa_name.lens" );
#open( FA,  ">$tgtdir/$fa_name.fasta" );
#open( NA,  ">$tgtdir/$fa_name.names" );
#open( MAP, ">$tgtdir/$fa_name.map" );
#
#@seqs_keep = sort { $a <=> $b } (@seqs_keep);
#
#foreach my $idx ( 1 .. @seqs_keep ) {
#
#  my $id = $seqs_keep[ $idx - 1 ];
#
#  print FA ">$idx " . $fa[2]->{$id} . "\n";
#  print FA $fa[0]->{$id} . "\n";
#  my $name = $fa[2]->{$id};
#  $name =~ s/ORIGID //;
#  print NA "$idx $name\n";
#  print LEN length( $fa[0]->{$id} ) . "\n";
#  print MAP $idx . " " . $id . "\n";
#}
#close(FA);
#close(NA);
#close(LEN);
#close(MAP);
