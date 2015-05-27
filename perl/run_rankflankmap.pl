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



my $updatedClustersCount = 0;
my $updatedRefsCount = 0;


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
if ($argc<5) { die "Usage: run_rankflankmap.pl inputfile  mapdir tmpdir dbname msdir\n"; }


my $inputfile = $ARGV[0];
my $mapdir = $ARGV[1];
my $tmp = $ARGV[2];
my $DBNAME = $ARGV[3];
my $MSDIR = $ARGV[4];


# set these mysql credentials in vs.cnf (in installation directory)
my ($LOGIN,$PASS,$HOST) = get_credentials($MSDIR);



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



my $clusters_processed = 0;


my $dbh = DBI->connect("DBI:mysql:$DBNAME;mysql_local_infile=1;host=$HOST;mysql_local_infile=1", "$LOGIN", "$PASS"
	           ) || die "Could not connect to database: $DBI::errstr";

my $sth;
my $sth1;
my $sth6;
my $sth8;

#goto AAA;

# clear map
$sth = $dbh->prepare('truncate table map;')
                or die "Couldn't prepare statement: " . $dbh->errstr;
$sth->execute()             # Execute the query
            or die "Couldn't execute statement: " . $sth->errstr;
$sth->finish;


# clear rankmap
$sth = $dbh->prepare('truncate table rankflank;')
                or die "Couldn't prepare statement: " . $dbh->errstr;
$sth->execute()             # Execute the query
            or die "Couldn't execute statement: " . $sth->errstr;
$sth->finish;


# disable indices
$sth = $dbh->prepare('ALTER TABLE map DISABLE KEYS;')
                or die "Couldn't prepare statement: " . $dbh->errstr;
$sth->execute()             # Execute the query
            or die "Couldn't execute statement: " . $sth->errstr;
$sth->finish;

$sth = $dbh->prepare('ALTER TABLE rankflank DISABLE KEYS;')
                or die "Couldn't prepare statement: " . $dbh->errstr;
$sth->execute()             # Execute the query
            or die "Couldn't execute statement: " . $sth->errstr;
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


# prepare statments
$sth6 = $dbh->prepare("LOAD DATA LOCAL INFILE '/$tmp/${DBNAME}_map.txt' INTO TABLE map FIELDS TERMINATED BY ',' LINES TERMINATED BY '\n';")
                or die "Couldn't prepare statement: " . $dbh->errstr;
$sth8 = $dbh->prepare("LOAD DATA LOCAL INFILE '/$tmp/${DBNAME}_rankflank.txt' INTO TABLE rankflank FIELDS TERMINATED BY ',' LINES TERMINATED BY '\n';")
                or die "Couldn't prepare statement: " . $dbh->errstr;




open FILE, "<$inputfile" or die "error openin for reading '$inputfile': $!";

opendir(DIR, $mapdir);
my @allfiles = readdir(DIR);
closedir(DIR);


 my $j=0;
 my $k=0;
 my $upload;
 my $uploadedrank=0;
 my $uploadedmap=0;

  # open out files
my $MAPFILE;
my $RFFILE;

  open $MAPFILE, ">/$tmp/${DBNAME}_map.txt" or die $!;
  open $RFFILE, ">/$tmp/${DBNAME}_rankflank.txt" or die $!;


 #foreach my $file (@files) {
 foreach my $file (@allfiles) {


if (index($file,"_map")) {

 $clusters_processed++;


 my $mfile = "$mapdir/$file";
 open MFILE, "<$mfile" or die $!; 


 my %READVECTOR = ();
 while (<MFILE>) {
  chomp;
  my @mfields = split('\t', $_);
  my $msize =  scalar @mfields;
  if ($msize >= 8) { 
   my $readid = $mfields[0];

   #print "\n\n".$_."\n\n";
   #exit(1);

   my @rfields = split(',',$mfields[7] );

   #my @temparray = ();

   my $bestscore = 0;
   my $bestref = "";
 
   foreach my $refstr (@rfields) {
    if ($refstr =~ /^-?(\d+):(\d+):(\d+)/) {

      #print "\n".$1." ".$2." ".$3."\n";
      #exit(1);

      # calculate score 
      my $score;

      if (( $mfields[2] + $mfields[3] ) == 0) # if no flanks, it will only be marked best if nothing else is available
         { $score = 0; }
      else
         #{ # asked to add by Dr. Benson to balance out small flanks

         #  my $A = ($2 + $3) / ( $mfields[2] + $mfields[3] );
         #  my $B = 1 / ( $mfields[2] + $mfields[3]);
         #  $score = 1 - max($A,$B);
         #}

         { $score = 1 - (($2 + $3) / ( $mfields[2] + $mfields[3] )); }

# filter to remove all flank scores below .9, added Nov 5, 2012
if ($score >= 0.90) {

      if ($score > $bestscore) {
	$bestref = $1;
	$bestscore = $score;
      }	elsif ($score == $bestscore) {
	if ($bestref eq "") { $bestref = $1; } else { $bestref .= (",".$1);}
      }
}
      # create a map entry in database (added 11/19/2010)
      #$sth6->execute($1,$readid)             # Execute the query
      #      or die "Couldn't execute statement: " . $sth6->errstr;
      $k++;
      print $MAPFILE "$1,$readid\n";

      if ($k % $RECORDS_PER_INFILE_INSERT == 0) {
        close($MAPFILE);
        $upload = $sth6->execute();
        $uploadedmap += $upload;

        open ($MAPFILE, ">/$tmp/${DBNAME}_map.txt") or die $!;
      }


    }      
   }

   # insert the rankflank
   if ($bestref ne "") {
     my @ranks = split(',', $bestref);
     my $ties = scalar(@ranks) - 1;
     foreach my $rstr (@ranks) {
      #$sth8->execute($readid,$rstr,$bestscore)             # Execute the query
      #      or die "Couldn't execute statement: " . $sth8->errstr;
      $j++;
      print $RFFILE "$rstr,$readid,$bestscore,$ties\n";

      if ($j % $RECORDS_PER_INFILE_INSERT == 0) {
        close($RFFILE);
        $upload = $sth8->execute();
        $uploadedrank += $upload;
        open ($RFFILE, ">/$tmp/${DBNAME}_rankflank.txt") or die $!;
      }

     }
   }

   #$READVECTOR{$readid} = @temparray;

   #print "\n";
  }
 }


 # load the files and remove the temp files


 print STDERR "\nprocessed: $clusters_processed";


 } # end of if (index($file,"_map")) {
 


} # end of foreach @files
close(FILE);



 SetStatistics("RANKFLANK_EDGES_INSERTED",$j);


# finish writing and loading out files
 close($RFFILE);
 close($MAPFILE);

 $upload = $sth6->execute();
 $uploadedmap += $upload;             

 $upload = $sth8->execute();
 $uploadedrank += $upload;            

 $sth6->finish;
 $sth8->finish;

 unlink("/$tmp/${DBNAME}_map.txt");
 unlink("/$tmp/${DBNAME}_rankflank.txt");



# enable indices
print STDERR "Enabling indices...\n";

$sth = $dbh->prepare('ALTER TABLE map ENABLE KEYS;')
                or die "Couldn't prepare statement: " . $dbh->errstr;
$sth->execute()             # Execute the query
            or die "Couldn't execute statement: " . $sth->errstr;
$sth->finish;

$sth = $dbh->prepare('ALTER TABLE rankflank ENABLE KEYS;')
                or die "Couldn't prepare statement: " . $dbh->errstr;
$sth->execute()             # Execute the query
            or die "Couldn't execute statement: " . $sth->errstr;
$sth->finish;



if ($uploadedmap != $k) { die "\nUploaded number of map entries($uploadedmap) not equal to the number of uploaded counter ($k), aborting!"; }
if ($uploadedrank != $j) { die "\nUploaded number of rankflank entries($uploadedrank) not equal to the number of uploaded counter ($j), aborting!"; }


AAA:

 # create temp table for deletions
 $sth = $dbh->prepare('CREATE TEMPORARY TABLE ranktemp (`refid` INT(11) NOT NULL, `readid` INT(11) NOT NULL, PRIMARY KEY (refid,readid)) ENGINE=INNODB;')
                or die "Couldn't prepare statement: " . $dbh->errstr;
 $sth->execute()             # Execute the query
            or die "Couldn't execute statement: " . $sth->errstr;
 $sth->finish;

 $sth = $dbh->prepare('ALTER TABLE ranktemp DISABLE KEYS;')
                or die "Couldn't prepare statement: " . $dbh->errstr;
 $sth->execute()             # Execute the query
            or die "Couldn't execute statement: " . $sth->errstr;
 $sth->finish;



 my $TEMPFILE;
 open ($TEMPFILE, ">$tmp/ranktemp_$DBNAME.txt") or die $!;



 print STDERR "Prunning (keep best ref for each read) from rankflank table...\n";
 my $query = "Select refid,readid,sid,score FROM rankflank INNER JOIN replnk ON rankflank.readid=replnk.rid ORDER BY readid, score;";
 $sth = $dbh->prepare($query);
 $sth->execute();
 my $num = $sth->rows;
 my $i=0;
 my $count=0;
 my $oldseq = -1;
 my $oldref = -1;
 my $oldread = -1;
 my $oldscore = -1.0;
 while ($i<$num) {
  my @data = $sth->fetchrow_array();
  if ($data[1] == $oldread && $data[3]!=$oldscore) {

    # delete old one
    print $TEMPFILE $oldref,",",$oldread,"\n";  

    $count++;
  }
  $oldref=$data[0];
  $oldread=$data[1];
  $oldseq=$data[2];
  $oldscore=$data[3];
  $i++;
 }

 $sth->finish;


 print STDERR "Prunning complete. Pruned $count rankflank records.\n";
 SetStatistics("RANKFLANK_REMOVED_SAMEREF",$count);


 print STDERR "Prunning all (one TR/same read) rankflank table...\n";
 $query = "Select refid,readid,sid,score FROM rankflank INNER JOIN replnk ON rankflank.readid=replnk.rid ORDER BY refid, sid, score, readid;"; # readid added for tie resolution to keep rank and rankflank entries more in sync
 $sth = $dbh->prepare($query);
 $sth->execute();
 $num = $sth->rows;
 $i=0;
 $count=0;
 $oldseq = -1;
 $oldref = -1;
 $oldread = -1;
 while ($i<$num) {
  my @data = $sth->fetchrow_array();
  if ($data[0] == $oldref && $data[2] == $oldseq) {

    # delete old one
    print $TEMPFILE $oldref,",",$oldread,"\n";  

    $count++;
  }
  $oldref=$data[0];
  $oldread=$data[1];
  $oldseq=$data[2];
  $i++;
 }

 $sth->finish;
 print STDERR "Prunning complete. Pruned $count rankflank records.\n";
 SetStatistics("RANKFLANK_REMOVED_SAMESEQ",$count);



 # load the file into tempfile
 close($TEMPFILE);
 $sth = $dbh->prepare("LOAD DATA LOCAL INFILE '$tmp/ranktemp_$DBNAME.txt' INTO TABLE ranktemp FIELDS TERMINATED BY ',' LINES TERMINATED BY '\n';")
                or die "Couldn't prepare statement: " . $dbh->errstr; 
 $sth->execute();
 $sth->finish;
 $sth = $dbh->prepare('ALTER TABLE ranktemp ENABLE KEYS;')
                or die "Couldn't prepare statement: " . $dbh->errstr;
 $sth->execute();
 $sth->finish;


 # delete from rankdflank based on temptable entries
 my $delfromtable = 0;
 $query = "DELETE FROM t1 USING rankflank t1 INNER JOIN ranktemp t2 ON ( t1.refid = t2.refid AND t1.readid = t2.readid );";
 $sth = $dbh->prepare($query);
 $delfromtable=$sth->execute();
 $sth->finish;


 # cleanup temp file
 unlink("$tmp/ranktemp_$DBNAME.txt");


 # set old settings
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



if ($delfromtable != $count) { die "Deleted number of entries($delfromtable) not equal to the number of deleted counter ($count), aborting!"; }


$dbh->disconnect();

print STDERR "\n\n";

print STDERR "Processing complete -- processed $clusters_processed cluster(s). Deleted from rankflank using temptable: $delfromtable\n";

($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
print STDERR sprintf("\n\nend: %4d-%02d-%02d %02d:%02d:%02d\n\n",$year+1900,$mon+1,$mday,$hour,$min,$sec);



1;




