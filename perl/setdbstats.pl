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

my %stats;
my $exstring;
my $input;

open($input, "-|", "wc -l $reffile | tail -1");
my $rc = <$input>;
if ($rc =~ /(\d+)/) {
  $stats{NUMBER_REF_TRS} = $1;
}
close($input);

open($input, "-|", "wc -l $readpf/*.indexhist | tail -1");
$rc = <$input>;
if ($rc =~ /(\d+)/) {
  $stats{NUMBER_TRS_IN_READS} = $1;
}
close($input);

open($input, "-|", "wc -l $reffolder/reference.leb36.rotindex | tail -1");
$rc = <$input>;
if ($rc =~ /(\d+)/) {
  $stats{NUMBER_REFS_TRS_AFTER_REDUND} = $1;
}
close($input);

open($input, "-|", "wc -l $rpfc/*.rotindex | tail -1");
$rc = <$input>;
if ($rc =~ /(\d+)/) {
  $stats{NUMBER_TRS_IN_READS_AFTER_REDUND} = $1;
}
close($input);


$stats{NUMBER_TRS_IN_READS_GE7} = 0;
$stats{NUMBER_READS_WITHTRS_GE7} = 0;
$stats{NUMBER_READS_WITHTRS} = 0;
$stats{NUMBER_READS_WITHTRS_GE7_AFTER_REDUND} = 0;

$rc = qx(./ge7.pl $readpf/*.index 2>/dev/null);
if ($rc =~ /(\d+) (\d+) (\d+)/) {
   $stats{NUMBER_TRS_IN_READS_GE7} = $1;
   $stats{NUMBER_READS_WITHTRS_GE7} = $2;
   $stats{NUMBER_READS_WITHTRS} = $3;
}
close($input);


open($input, "-|", "./ge7.pl $readpf/*.indexhist");
$rc = <$input>;
if ($rc =~ /(\d+) (\d+) (\d+)/) {
  $stats{NUMBER_READS_WITHTRS} = $3;
}
close($input);


open($input, "-|", "cat $rpfc/*.rotindex | wc");
$rc = <$input>;
if ($rc =~ /(\d+) (\d+) (\d+)/) {
  $stats{NUMBER_READS_WITHTRS_GE7_AFTER_REDUND} = $1;
}
close($input);

set_statistics($DBSUFFIX, %stats);

1;