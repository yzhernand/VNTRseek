#!/usr/bin/perl

my $RECORDS_PER_INFILE_INSERT = 100000;

use List::Util qw[min max];

use strict;
use warnings;
use Cwd;
use DBI;

use FindBin;
use File::Basename;

use lib "$FindBin::Bin/vntr"; 
require "vutil.pm";

use vutil ('get_credentials');
use vutil ('write_mysql');
use vutil ('stats_set');

my $sec; 
my $min;
my $hour;
my $mday;
my $mon;
my $year;
my $wday;
my $yday;
my $isdst;

sub nowhitespace($)
{
	my $string = shift;
	$string =~ s/\s+//g;
	return $string;
}

($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
print STDERR sprintf("\n\nstart: %4d-%02d-%02d %02d:%02d:%02d\n\n",$year+1900,$mon+1,$mday,$hour,$min,$sec);

my $curdir =  getcwd;

my $argc = @ARGV;
if ($argc<4) { die "Usage: run_flankcomp.pl inputfile dbname msdir tempdir\n"; }


my $inputfile = $ARGV[0];
my $DBNAME = $ARGV[1];
my $MSDIR = $ARGV[2];
my $TEMPDIR = $ARGV[3];

# set these mysql credentials in vs.cnf (in installation directory)
my ($LOGIN,$PASS,$HOST) = get_credentials($MSDIR);


my $clusters_processed = 0;
my $totalRefReps = 0;
my $totalReadReps = 0;


my $dbh = DBI->connect("DBI:mysql:$DBNAME;mysql_local_infile=1;host=$HOST", "$LOGIN", "$PASS"
	           ) || die "Could not connect to database: $DBI::errstr";

my $sth;
my $sth1;
my $sth2;
my $sth3;

my $mostReps = 0;
my $mostRefReps = 0;
my $maxRange = 0;


my $BREAK_SIZE = 4000;



#for my $i (keys %BELONG) {
#  print STDERR $i.":".$BELONG{$i}."\n";
#}
#exit(1);

# clear database cluster tables
$sth = $dbh->prepare('TRUNCATE TABLE clusters; ') or die "Couldn't prepare statement: " . $dbh->errstr;
     $sth->execute()             # Execute the query
            or die "Couldn't execute statement: " . $sth->errstr;

$sth = $dbh->prepare('TRUNCATE TABLE clusterlnk; ') or die "Couldn't prepare statement: " . $dbh->errstr;
     $sth->execute()             # Execute the query
            or die "Couldn't execute statement: " . $sth->errstr;


$sth = $dbh->prepare('ALTER TABLE clusterlnk DISABLE KEYS;')
                or die "Couldn't prepare statement: " . $dbh->errstr;
$sth->execute or die "Cannot execute: " . $sth->errstr();;
$sth->finish;

$sth = $dbh->prepare('SET AUTOCOMMIT = 0;')
                or die "Couldn't prepare statement: " . $dbh->errstr;
$sth->execute()             # Execute the query
            or die "Couldn't execute statement: " . $sth->errstr;
$sth->finish;

$sth = $dbh->prepare('SET FOREIGN_KEY_CHECKS = 0;')
                or die "Couldn't prepare statement: " . $dbh->errstr;
$sth->execute()             # Execute the query
            or die "Couldn't execute statement: " . $sth->errstr;
$sth->finish;

$sth = $dbh->prepare('SET UNIQUE_CHECKS = 0;')
                or die "Couldn't prepare statement: " . $dbh->errstr;
$sth->execute()             # Execute the query
            or die "Couldn't execute statement: " . $sth->errstr;
$sth->finish;



$sth2 = $dbh->prepare('INSERT INTO clusters(cid,minpat,maxpat,repeatcount,refcount) VALUES(?,?,?,?,?)')
                or die "Couldn't prepare statement: " . $dbh->errstr;


####################################
sub SetStatistics {

  my $argc = @_;
  if ($argc <2) { die "stats_set: expects 2 parameters, passed $argc !\n"; }

  my $NAME = $_[0];
  my $VALUE = $_[1];

  #print "$DBNAME,$LOGIN,$PASS,$NAME,$VALUE\n";
  return stats_set($DBNAME,$LOGIN,$PASS,$HOST,$NAME,$VALUE);
}
####################################



#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################


print STDERR "\nInserting into clusterlnk table...\n";

open FILE, "<$inputfile" or die $!;
#open QFILE, ">$inputfile.reads.fq" or die $!;


$sth3 = $dbh->prepare("LOAD DATA LOCAL INFILE '$TEMPDIR/clusterlnk_$DBNAME.txt' INTO TABLE clusterlnk FIELDS TERMINATED BY ',' LINES TERMINATED BY '\n';")
                or die "Couldn't prepare statement: " . $dbh->errstr;


my $TEMPFILE;
open ($TEMPFILE, ">$TEMPDIR/clusterlnk_$DBNAME.txt") or die $!;


# insert into clusterlnk
my $totalreps = 0;
while (<FILE>) { 
 $clusters_processed++;

 chomp;
 
 my @values = split(',', $_);

 foreach my $val (@values) {

     # insert clusterlnk entry

     $totalreps++;

     my $dir = '\'';
     if ($val =~ m/([\'\"])/) { $dir = $1; }

     print $TEMPFILE $clusters_processed,",",$val,",",$dir,",",0,",",0,"\n";

     if ($totalreps % $RECORDS_PER_INFILE_INSERT == 0) {
       close($TEMPFILE);
       $sth3->execute();
       open ($TEMPFILE, ">$TEMPDIR/clusterlnk_$DBNAME.txt") or die $!;
     }


 }


} # end of while loop


# finish clustelnk load infile
close($TEMPFILE);
$sth3->execute();
unlink("$TEMPDIR/clusterlnk_$DBNAME.txt");


# reenable indices
$sth = $dbh->prepare('ALTER TABLE clusterlnk ENABLE KEYS;')
                or die "Couldn't prepare statement: " . $dbh->errstr;
$sth->execute or die "Cannot execute: " . $sth->errstr();;
$sth->finish;


#cleanup
$sth2->finish;
$sth3->finish;



#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################

print STDERR "\nPrinting DNA and inserting into cluster table...\n";

# now print dna and quals (also insert into cluster table)
$sth = $dbh->prepare('SELECT rid,flankleft,sequence,flankright,pattern,copynum,direction FROM fasta_ref_reps INNER JOIN clusterlnk ON rid=-repeatid WHERE clusterid = ?')
                or die "Couldn't prepare statement: " . $dbh->errstr;
$sth1 = $dbh->prepare('SELECT rid, dna, first, last, pattern, copynum,direction FROM fasta_reads INNER JOIN replnk ON fasta_reads.sid=replnk.sid INNER JOIN clusterlnk ON rid=repeatid WHERE clusterid = ?')
                or die "Couldn't prepare statement: " . $dbh->errstr;


my $i;
seek(FILE,0,0);
$clusters_processed = 0;
while (<FILE>) {
 $clusters_processed++;

 chomp;
 
 my @values = split(',', $_);
   
 my $repeatcount = 0;
 my $refcount = 0;
 my $readcount = 0;
 my $minpat = 1000000;
 my $maxpat = 0;
 my $range = 0;


# process each line
{

  # for statistics
  foreach my $val (@values) {

     $val =~ s/[\'\"]//g;

     # insert clusterlnk entry
     if ($val<=0) {

         $refcount++;

         $totalRefReps++;

     } else {

         $readcount++;

         $totalReadReps++;

     }
 
     $repeatcount++;
  }


  # execute ref and read pulls
  $sth->execute($clusters_processed)
      or die "Couldn't execute statement: " . $sth->errstr;
  $sth1->execute($clusters_processed)
      or die "Couldn't execute statement: " . $sth1->errstr;


  my $numrefs = $sth->rows;
  my $numreads = $sth1->rows;

   
  # store refs for later use
  $i=0;
  open (RFILE, ">$TEMPDIR/refs_$DBNAME.txt") or die $!;
  while ($i<$numrefs) { 
     (my @data = $sth->fetchrow_array()) or die "Can't fetch reference row from database!";

     print RFILE "-".$data[0].$data[6].",".$data[1].",".$data[2].",".$data[3].",".$data[4]."\n";

     $minpat = min( $minpat, length($data[4]) );
     $maxpat = max( $maxpat, length($data[4]) );

     $i++;
  }
  close(RFILE);


  # print reads in blocks of BREAK_SIZE
  $i = 0;
  my $hcount = 0;
  while ($i<$numreads) { 

    if (($i % $BREAK_SIZE)==0) {
     $hcount++;
     
     print "@($clusters_processed\_$hcount):";
     print "\n**********************************************************************\n";

     # print refs each time 
     open (RFILE, "$TEMPDIR/refs_$DBNAME.txt") or die $!;
     while (<RFILE>) { print $_; }
     close(RFILE);
    }

    (my @data = $sth1->fetchrow_array()) or die "Can't fetch read row from database!";
    my $dna = nowhitespace($data[1]);
    print $data[0].$data[6].",".$data[2].",".$data[3].",".$dna.",".$data[4]."\n";

    $minpat = min( $minpat, length($data[4]) );
    $maxpat = max( $maxpat, length($data[4]) );

    $i++;
  }
  


} # end of process each line


 # do for 1st 10 clusters for now
 #if ($clusters_processed >= 20) { last; }

 # insert database records (cluster table)
 $sth2->execute($clusters_processed,$minpat,$maxpat,$repeatcount,$refcount)  # Execute the query
         or die "Couldn't execute statement: " . $sth2->errstr;

 # stats
 $range = int(($maxpat / $minpat - 1.0) * 100 + .5);
 $mostReps =  max( $mostReps, $repeatcount );
 $mostRefReps =  max( $mostRefReps, $refcount );
 $maxRange =  max( $maxRange, $range );

 
} # end of while loop


# delete temp ref file
unlink("$TEMPDIR/refs_$DBNAME.txt");


#cleanup
close(FILE);
#close(QFILE);

$sth->finish;
$sth1->finish;


# enable old settings
$sth = $dbh->prepare('SET AUTOCOMMIT = 1;')
                or die "Couldn't prepare statement: " . $dbh->errstr;
$sth->execute()             # Execute the query
            or die "Couldn't execute statement: " . $sth->errstr;
$sth->finish;

$sth = $dbh->prepare('SET FOREIGN_KEY_CHECKS = 1;')
                or die "Couldn't prepare statement: " . $dbh->errstr;
$sth->execute()             # Execute the query
            or die "Couldn't execute statement: " . $sth->errstr;
$sth->finish;

$sth = $dbh->prepare('SET UNIQUE_CHECKS = 1;')
                or die "Couldn't prepare statement: " . $dbh->errstr;
$sth->execute()             # Execute the query
            or die "Couldn't execute statement: " . $sth->errstr;
$sth->finish;



$dbh->disconnect();


# update the stats table
SetStatistics("CLUST_LARGEST_NUMBER_OF_TRS_IN_PROCLU_CLUSTER",$mostReps);
SetStatistics("CLUST_LARGEST_NUMBER_OF_REFS_IN_PROCLU_CLUSTER",$mostRefReps);
SetStatistics("CLUST_LARGEST_PATRANGE_IN_PROCLU_CLUSTER",$maxRange);

SetStatistics("CLUST_NUMBER_OF_PROCLU_CLUSTERS",$clusters_processed);
SetStatistics("CLUST_NUMBER_OF_REF_REPS_IN_CLUSTERS",$totalRefReps);
SetStatistics("CLUST_NUMBER_OF_READ_REPS_IN_CLUSTERS",$totalReadReps);



print STDERR "Processing complete -- processed $clusters_processed cluster(s).\n";

($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
print STDERR sprintf("\n\nend: %4d-%02d-%02d %02d:%02d:%02d\n\n",$year+1900,$mon+1,$mday,$hour,$min,$sec);



1;




