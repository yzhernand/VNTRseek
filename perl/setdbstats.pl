#!/usr/bin/perl

# sets a number of db stat variables
#
# command line usage example:
#  ./setdbstats.pl reference_file reads_profiles_folder reference_folder reads_profile_folder_clean dbname dblogin dbpass dbhost
# where inputfile is the main cluster file
#

use strict;
use warnings;
use Cwd;

use FindBin;
use File::Basename;

use List::Util qw[min max];

use lib "$FindBin::RealBin/lib"; # must be same as install dir!

use vutil qw( get_config set_statistics );


my $argc = @ARGV;
if ($argc<6) { die "Usage: setdbstats.pl reference_file reads_profiles_folder reference_folder reads_profile_folder_clean dbsuffix msdir\n"; }

my $reffile = $ARGV[0];
my $readpf = $ARGV[1];
my $reffolder = $ARGV[2];
my $rpfc = $ARGV[3];
my $DBSUFFIX = $ARGV[4];
my $MSDIR = $ARGV[5];

# set these mysql credentials in vs.cnf (in installation directory)
my %run_conf = get_config($MSDIR . "vs.cnf");
my ( $LOGIN, $PASS, $HOST ) = @run_conf{qw(LOGIN PASS HOST)};

####################################

my $rc;
my $exstring;
my $input;

open($input, "-|", "wc -l $reffile | tail -1");
$rc = <$input>;
if ($rc =~ /(\d+)/) {
  set_statistics($DBSUFFIX, 'NUMBER_REF_TRS', $rc);
}
close($input);

open($input, "-|", "wc -l $readpf/*.indexhist | tail -1");
$rc = <$input>;
if ($rc =~ /(\d+)/) {
  set_statistics($DBSUFFIX, 'NUMBER_TRS_IN_READS', $rc);
}
close($input);

open($input, "-|", "wc -l $reffolder/reference.leb36.rotindex | tail -1");
$rc = <$input>;
if ($rc =~ /(\d+)/) {
  set_statistics($DBSUFFIX, 'NUMBER_REFS_TRS_AFTER_REDUND', $rc);
}
close($input);

open($input, "-|", "wc -l $rpfc/*.rotindex | tail -1");
$rc = <$input>;
if ($rc =~ /(\d+)/) {
  set_statistics($DBSUFFIX, 'NUMBER_TRS_IN_READS_AFTER_REDUND', $rc);
}
close($input);


my $readTRsWithPatternGE7 = 0;
my $totalReadsWithTRsPatternGE7 = 0;
my $totalReadsWithTRs = 0;
my $readTRsWPGE7AfterCyclicRedundancyElimination = 0;

open($input, "-|", "./ge7.pl $readpf/*.index");
$rc = <$input>;
if ($rc =~ /(\d+) (\d+) (\d+)/) {
   $readTRsWithPatternGE7 = $1;
   $totalReadsWithTRsPatternGE7 = $2;
   $totalReadsWithTRs = $3;
}
close($input);


open($input, "-|", "./ge7.pl $readpf/*.indexhist");
$rc = <$input>;
if ($rc =~ /(\d+) (\d+) (\d+)/) {
  $totalReadsWithTRs = $3;
}
close($input);


open($input, "-|", "cat $rpfc/*.rotindex | wc");
$rc = <$input>;
if ($rc =~ /(\d+) (\d+) (\d+)/) {
  $readTRsWPGE7AfterCyclicRedundancyElimination = $1;
}
close($input);


set_statistics($DBSUFFIX, "NUMBER_TRS_IN_READS_GE7", $readTRsWithPatternGE7);
set_statistics($DBSUFFIX, "NUMBER_READS_WITHTRS_GE7", $totalReadsWithTRsPatternGE7);
set_statistics($DBSUFFIX, "NUMBER_READS_WITHTRS", $totalReadsWithTRs);
set_statistics($DBSUFFIX, "NUMBER_READS_WITHTRS_GE7_AFTER_REDUND", $readTRsWPGE7AfterCyclicRedundancyElimination);


1;




