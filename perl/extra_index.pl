#!/usr/bin/perl

# prints sequences for BBB for pcr_dup step

use strict;
use warnings;
use Cwd;

# use strict;
use warnings;
use Cwd;
use DBI;

use FindBin;
use File::Basename;

use lib "$FindBin::Bin/vntr"; 
require "vutil.pm";

use vutil ('get_credentials');

my $sec; 
my $min;
my $hour;
my $mday;
my $mon;
my $year;
my $wday;
my $yday;
my $isdst;

my $FASTA = 1;

# Perl trim function to remove whitespace from the start and end of the string
sub trim($)
{
	my $string = shift;
	$string =~ s/^\s+//;
	$string =~ s/\s+$//;
	return $string;
}


($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
printf "\n\nstart: %4d-%02d-%02d %02d:%02d:%02d\n\n\n",$year+1900,$mon+1,$mday,$hour,$min,$sec;

my $argc = @ARGV;

if ($argc<3) { die "Usage: extra_index.pl expects 3 arguments!\n"; }

my $curdir =  getcwd;

my $folder = $ARGV[0];
my $DBNAME = $ARGV[1];
my $MSDIR = $ARGV[2];

# set these mysql credentials in vs.cnf (in installation directory)
my ($LOGIN,$PASS,$HOST) = get_credentials($MSDIR);


 # create folder
 my $exstring = "rm -rf $folder";
 system($exstring);
 mkdir($folder);


 # db setup
 my $dbh = DBI->connect("DBI:mysql:$DBNAME;host=$HOST", "$LOGIN", "$PASS"
                   ) || die "Could not connect to database: $DBI::errstr";

 my $sth;
 my $query;
 my $result;
 my $num;
 my $i;


 $sth = $dbh->prepare('SELECT map.refid, map.readid, replnk.sid, replnk.first, replnk.last, replnk.copynum, replnk.patsize, replnk.pattern,fasta_reads.dna from map INNER JOIN rank ON rank.refid=map.refid AND rank.readid=map.readid INNER JOIN rankflank ON rankflank.refid=map.refid AND rankflank.readid=map.readid INNER JOIN replnk ON replnk.rid=map.readid INNER JOIN fasta_reads on fasta_reads.sid=replnk.sid ORDER BY map.refid,map.readid;')
                or die "Couldn't prepare statement: " . $dbh->errstr;

 $sth->execute();
 $num = $sth->rows; 

 print "\n best best best records: $num\n";
 
 my $oldref = -1;
 my $oldread = -1;
 $i=0;
 my $nrefs=0;

 while ($i<$num) {
  my @data = $sth->fetchrow_array();
  if ($data[0] != $oldref) {
     if ($i!=0) { close(FILE); }
     $nrefs++;
     print "\n$nrefs";
     open FILE, ">$folder/$data[0].seq" or die $!;    
  }
  print FILE "$data[1] $data[2] $data[3] $data[4]";
  printf FILE " %.2lf",$data[5];
  $data[8] =~ s/\s+//g;
  print FILE " $data[6] $data[7] $data[8]\n";
  $oldref=$data[0];
  $oldread=$data[1];
  $i++;
 }
 
 close(FILE);

 $sth->finish;



print "\n\nProcessing complete (extra_index.pl), $nrefs files created.\n";

($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
printf "\n\nend: %4d-%02d-%02d %02d:%02d:%02d\n",$year+1900,$mon+1,$mday,$hour,$min,$sec;

 #
$dbh->disconnect();

1;




