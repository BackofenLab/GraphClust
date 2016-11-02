#!/usr/bin/perl
use warnings;
use strict;

## GraphClust stage 8 script for cmsearch

use List::Util qw/ min max /;
use FindBin;
use lib "$FindBin::Bin";
use GraphClust;

use vars qw(%CONFIG);
*CONFIG = \%GraphClust::CONFIG;

use Getopt::Long;

my $binDir = "$FindBin::Bin";

my $stk_file;
my $db;
my $tgtdir;
my $evaluate;
my $use_cmcalibrate;
my $in_root_dir;
my $verbose = 1;
my $newInf;
my $OPTS_cmsearch;
my $OPTS_cmbuild;
my $OPTS_cmcalibrate;




usage()
  unless GetOptions(
  "root=s"           => \$in_root_dir,
  "tgtdir=s"         => \$tgtdir,
  "stk=s"            => \$stk_file,
  "db=s"             => \$db,
  "calibrate"        => \$use_cmcalibrate,
  "verbose"          => \$verbose,
  "infernal_1_1"     => \$newInf
  );

################################################################################
## Load options from root-dir/config file
%GraphClust::CONFIG = readConfigFile( "$in_root_dir/config", 1 );
SUBSECTION("options loaded for gc_cmsearch.pl");
printConfig( \%CONFIG ) if ($verbose);



## infernal 1.1
if($newInf){
        $OPTS_cmsearch = " -g --nohmm --noali --cpu 1 ";
        $OPTS_cmbuild = "";
        $OPTS_cmcalibrate = " -L 0.01 --cpu 1 ";
}

## infernal 1.0.2
else{
        $OPTS_cmsearch = " -g --fil-no-hmm --noalign ";
        $OPTS_cmbuild = "";
        $OPTS_cmcalibrate = " -L 0.01 ";
}



my $infernal_path;
my $cm_top_only = 1;
$cm_top_only = $CONFIG{cm_top_only};
if($newInf){
	$infernal_path = $CONFIG{PATH_INFERNAL_1_1};
}
else{
	$infernal_path = $CONFIG{PATH_INFERNAL_1_0};
}
my $cm_bitscore_sig = $CONFIG{cm_bitscore_sig};
my $cm_min_bitscore = $CONFIG{cm_min_bitscore};

## check top_only option consistency
## if we have revcompl seqs than we should also scan on them
## remember: scan is done on data.fasta.scan which are the original input seqs
## and they dont contain the added revcompl seqs in general
if ( $CONFIG{input_add_revcompl} && $cm_top_only ) {
  print "\nATTENTION! Option 'cm_top_only' changed to 0 (false)!\n";
  print " This is necessary because option 'input_add_revcompl' is 1 (true)!\n\n";
  $cm_top_only = 0;
}

################################################################################

mkdir($tgtdir);

die "No input file specified....\n\n" if ( !$stk_file || !-e $stk_file);
die "No input database specified....\n" if ( !$db || !-e $db );

################################
### build the covariance model |
################################

my @file_line = split( /\//, $stk_file );
my $target_file = $file_line[ scalar(@file_line) - 1 ];

my $cmd_cmbuild = "$infernal_path/cmbuild $OPTS_cmbuild -F $tgtdir/$target_file.cm $stk_file";

system($cmd_cmbuild);

####################################
### calibrate the covariance model |
####################################

#### for infernal 1.1 calibration is mandatory
my $cmd_cal="";

if($newInf){

	$cmd_cal = "$infernal_path/cmcalibrate $OPTS_cmcalibrate $tgtdir/$target_file.cm";
	system($cmd_cal);
}
else{
	print "####################### elsi  mejja #######################\n";

	if ($use_cmcalibrate ) {

	  $cmd_cal = "$infernal_path/cmcalibrate $OPTS_cmcalibrate $tgtdir/$target_file.cm";
	  system($cmd_cal);
  	}
}

#####################################################################################################
### search the database with the calibrated covariance model and report the findings into a tabfile |
#####################################################################################################
### Each non-# prefixed line of this file corresponds to a hit, and each such line has 9 fields:    |
### <model name> the name of the CM used for the search, 					    |
### <target name> the name of the target sequence the hit was found in, 			    |
### <target coord - start> the start position of the hit in the target sequence, 		    |
### <target coord - stop> the end position of hit in the target sequence, 			    |
### <query coord - start> the start position of the hit in the query model, 			    |
### <query coord - stop> the end position of hit in the query sequence, 			    |
### <bit sc> the bit score of the hit, 								    |
### <E-value> the E-value of the hit (if available, ”-” if not),				    |
### <GC> the percentage of G and C residues in the hit within the target sequence.		    |
#####################################################################################################
## add '-o $e/$target_file.cm.alignments' if also alignment file should be created

my $tmp_tab_file = "$tgtdir/$target_file.cm.tab_tmp";
my $cmd_cms="";

if($newInf){
	$cmd_cms = "$infernal_path/cmsearch  $OPTS_cmsearch --tblout $tmp_tab_file ";
	##die "hastat mtnuma newInfi mej.... $cmd_cms \n\n";
}
else{
	$cmd_cms = "$infernal_path/cmsearch $OPTS_cmsearch --tabfile $tmp_tab_file ";
}
$cmd_cms .= "--toponly " if ($cm_top_only);
my $min_bitscore = min(10,$cm_min_bitscore);
$cmd_cms .= " -T $min_bitscore " if ($cm_bitscore_sig == 1);
$cmd_cms .= " $tgtdir/$target_file.cm $db ";

system("rm -f $tmp_tab_file");
system("rm -f $tgtdir/$target_file.cm.tabresult");

## run cmsearch
system($cmd_cms);

## postprocessing
system("mv $tmp_tab_file $tgtdir/$target_file.cm.tabresult ");
