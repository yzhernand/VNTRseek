#!/usr/bin/perl

my $RECORDS_PER_INFILE_INSERT = 100000;

use List::Util qw[min max];

#use strict;
use warnings;
use Cwd;
use DBI;

use FindBin;
use File::Basename;

use lib "$FindBin::RealBin/lib"; 
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
#my $PROCLU_PARAM = "$curdir/eucledian.dst $curdir/eucledian.dst 0 70 -se ";
my $PROCLU_PARAM = " $curdir/eucledian.dst 70 0 0 ";


my $files_to_process = 100;        # number of files to process in one batch
my $files_processed = 0;        # files processed
my %p; # associates forked pids with output pipe pids

my $MYLOCK=0;

my $argc = @ARGV;

if ($argc<8) { die "Usage: run_edges.pl reference_file edges_folder dbname msdir MINPROFSCORE NPROCESSORS PSEARCH TMPDIR\nn"; }


my $inputfile = $ARGV[0];
my $folder = $ARGV[1];
my $DBNAME = $ARGV[2];
my $MSDIR = $ARGV[3];
my $MINPROFSCORE = $ARGV[4];
my $cpucount = $ARGV[5];
my $PROCLU = "$curdir/$ARGV[6]";
my $tmp = $ARGV[7];

# set these mysql credentials in vs.cnf (in installation directory)
my ($LOGIN,$PASS,$HOST) = get_credentials($MSDIR);

my $clusters_processed = 0;
my $totalRefReps = 0;
my $totalReadReps = 0;

my $dbh = DBI->connect("DBI:mysql:$DBNAME;mysql_local_infile=1;host=$HOST", "$LOGIN", "$PASS"
	           ) || die "Could not connect to database: $DBI::errstr";

my $sth;
my $sth2;
my $query;
my $result;
my $num;
my $i;
my $clusterid;
my $clusteridold;
my $repid;
my $exstring;
my $leftname = "";
my $rightname = "";

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

#goto AAA;

 # create folder
 $exstring = "rm -f $folder -R";
 system($exstring);
 mkdir($folder);

 $query = "CREATE TEMPORARY TABLE tempmap (rid int  PRIMARY KEY);";
 $sth = $dbh->prepare($query);
 $sth->execute();
 $sth->finish;


$sth = $dbh->prepare('ALTER TABLE tempmap DISABLE KEYS;')
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



 $query = "INSERT INTO tempmap(rid) SELECT DISTINCT -refid FROM map;";
 $sth = $dbh->prepare($query);
 $sth->execute();
 $sth->finish;


$sth = $dbh->prepare('ALTER TABLE tempmap ENABLE KEYS;')
                or die "Couldn't prepare statement: " . $dbh->errstr;
$sth->execute()             # Execute the query
            or die "Couldn't execute statement: " . $sth->errstr;
$sth->finish;



 # print ref leb files
 my %LHASH = ();
 open(MYINPUTFILE, "<$inputfile") or die "\nCannot open file '$inputfile'!\n"; # open for inpu
 while (<MYINPUTFILE>) {  # read file into list
  if (/^(\d+)/) {
    $LHASH{$1} = $_;
    #print "aaa: ($inputfile): $1";
    #exit(0);
  }
 }
 close(MYINPUTFILE);

 $clusters_processed=0;
 $query = "SELECT repeatid,clusterid,profile,profilerc,patsize,copynum FROM clusterlnk LEFT OUTER JOIN replnk ON clusterlnk.repeatid=replnk.rid INNER JOIN tempmap on clusterlnk.repeatid=tempmap.rid ORDER BY clusterid;";
 $sth = $dbh->prepare($query);
 $sth->execute();
 $num = $sth->rows;
 $i=0;
 $clusteridold = -1;
 open(MYOUTFILE, ">$folder/dummy.txt") or die "\nCannot open file '$folder/dummy.txt'!\n"; # open for output
 while ($i<$num) {
  my @data = $sth->fetchrow_array();
  $clusterid=$data[1];
  $repid=$data[0];

  #print "repid: $repid"; exit(0);
  #print "$clusterid / $repid\n";
 
  if ($clusterid != $clusteridold) {
   $leftname = "refs.${clusterid}.leb36";
   close(MYOUTFILE);
   open(MYOUTFILE, ">$folder/$leftname") or die "\nCannot open file '$folder/$leftname'!\n"; # open for output
   
 
   $clusters_processed++;
   print "$clusters_processed\n";
  }

  {
   $repid = -1 * $repid;
   #$exstring = "echo -`grep $repid $inputfile`  >> $folder/$leftname";
   #system($exstring);

   if (exists $LHASH{$repid}) 
     { print MYOUTFILE "-".$LHASH{$repid}; } 
   else 
     {
       print STDERR "\nCould not lookup reference in hash ($repid)!\n";
       exit(1);
     }

  } 

  $clusteridold = $clusterid;
  $i++;
 }
 close MYOUTFILE;

$sth->finish;

print STDERR "Processing complete -- outputed $clusters_processed ref leb files. Creating input edges files ...\n\n";



 # create edges files
 $clusters_processed=0;
 $query = "SELECT clusterid,refid,readid FROM map,clusterlnk WHERE -refid=clusterlnk.repeatid ORDER BY clusterid, refid, readid;";
 $sth = $dbh->prepare($query);
 $sth->execute();
 $num = $sth->rows;
 $i=0;
 $clusteridold = -1;
 my $refid; 
 my $readid;
 while ($i<$num) {
  my @data = $sth->fetchrow_array();
  $clusterid=$data[0];
  $refid=$data[1];
  $readid=$data[2];

  if ($clusterid != $clusteridold) {
   $leftname = "reads.${clusterid}.edgein";
   close(EFILE); 
   open( EFILE,">$folder/$leftname") or die "Can't open $folder/$leftname!";
   $clusters_processed++;
   print "$clusters_processed\n";
  }
   print EFILE "-".$refid.",".$readid."\n";

   $clusteridold = $clusterid;
   $i++;
 }
 $sth->finish;

close(EFILE);


print STDERR "Processing complete -- outputed $clusters_processed edges input files.\n\n";

#AAA:


 $query = "DROP TABLE tempmap;";
 $sth = $dbh->prepare($query);
 $sth->execute();
 $sth->finish;

 $query = "CREATE TEMPORARY TABLE tempmap (rid int  PRIMARY KEY);";
 $sth = $dbh->prepare($query);
 $sth->execute();
 $sth->finish;


$sth = $dbh->prepare('ALTER TABLE tempmap DISABLE KEYS;')
                or die "Couldn't prepare statement: " . $dbh->errstr;
$sth->execute()             # Execute the query
            or die "Couldn't execute statement: " . $sth->errstr;
$sth->finish;


 $query = "INSERT INTO tempmap(rid) SELECT DISTINCT readid FROM map;";
 $sth = $dbh->prepare($query);
 $sth->execute();
 $sth->finish;


$sth = $dbh->prepare('ALTER TABLE tempmap ENABLE KEYS;')
                or die "Couldn't prepare statement: " . $dbh->errstr;
$sth->execute()             # Execute the query
            or die "Couldn't execute statement: " . $sth->errstr;
$sth->finish;



 # print read leb files
 $clusters_processed=0;

 $query = "select cid FROM clusters ORDER BY cid;";
 $sth2 = $dbh->prepare($query);
 $sth2->execute();
 $nclusters = $sth2->rows;
 $clusteridold = -1;

 while ($clusters_processed<$nclusters) {
  my @data2 = $sth2->fetchrow_array();

   $clusterid=$data2[0];

   print "$clusters_processed\n";


 #$query = "SELECT repeatid,clusterid,profile,profilerc,patsize,copynum FROM clusterlnk LEFT OUTER JOIN replnk ON clusterlnk.repeatid=replnk.rid INNER JOIN tempmap on clusterlnk.repeatid=tempmap.rid ORDER BY clusterid;";	
 $query = "SELECT repeatid,clusterid,profile,profilerc,patsize,copynum FROM clusterlnk LEFT OUTER JOIN replnk ON clusterlnk.repeatid=replnk.rid  WHERE clusterid=$clusterid AND clusterlnk.repeatid IN (SELECT rid FROM tempmap);";
 $sth = $dbh->prepare($query);
 $sth->execute();
 $num = $sth->rows;

 if ($num > 0) {
   $rightname = "reads.${clusterid}.leb36";
   close(MYFILE); 
   open (MYFILE, ">$folder/$rightname") or die "\nCannot open file '$folder/$rightname'!\n";;
 }

 $i=0;
 while ($i<$num) {
  my @data = $sth->fetchrow_array();
  $repid=$data[0];
  print "$clusterid / $repid\n";
 
  {
   my $copies = sprintf("%.2lf",$data[5]);
   my $pline = $repid ." ". $data[4] ." ". $copies  ." ". (length($data[2]) / 2) ." ". (length($data[3]) / 2) ." ". $data[2] ." ". $data[3] ." 0 0 0 0 |\n";
   print MYFILE $pline;
  }

  $i++;
 }

 $sth->finish;

 $clusteridold = $clusterid;
 $clusters_processed++;
}
$sth2->finish;



close(MYFILE); 



print STDERR "Processing complete -- outputed $clusters_processed read leb files.\n\n";

$dbh->disconnect();


#AAA:


# get a list of input files
opendir(DIR, $folder);
# the only extensions are .leb36
my @tarballs = grep(/reads\.(\d+)\.(?:leb36)$/, readdir(DIR));
closedir(DIR);
my $tarball_count = @tarballs;
print STDERR "$tarball_count supported files found in $folder\n";
#die "Exiting\n" if $tarball_count == 0;
$files_to_process = $tarball_count if $files_to_process > $tarball_count;


# enter dir
chdir($folder);


# fork as many new processes as there are CPUs
for (my $i = 0; $i < $cpucount; $i++) { $p{fork_proclu()} = 1 }

# wait for processes to finish and then fork new ones
while ((my $pid = wait) != -1) {


      # check return value
      my ($rc, $sig, $core) = ($? >> 8, $? & 127, $? & 128);
      if ($core){
         print STDERR "proclu process $pid dumped core\n";
         exit (1000);
      }elsif($sig == 9){
         print STDERR "proclu process $pid was murdered!\n";
         exit (1001);
      }elsif ($rc != 0){
         print STDERR  "proclu process $pid has returned $rc!\n";
         exit ($rc);
      }

        if ($p{$pid}) {
                # one instance has finished processing -- start a new one
                delete $p{$pid};
                $p{fork_proclu()} = 1;
        } else {
                die "************ Process $pid finished (not in hash)\n";
        }
}


print STDERR "Processing complete -- processed $files_processed cluster(s).\n";

($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
print STDERR sprintf("\n\nend: %4d-%02d-%02d %02d:%02d:%02d\n\n",$year+1900,$mon+1,$mday,$hour,$min,$sec);

# update database
#AAA:

 my %RHASH = ();
 my %SHASH = ();


 $dbh = DBI->connect("DBI:mysql:$DBNAME;mysql_local_infile=1;host=$HOST", "$LOGIN", "$PASS"
                   ) || die "Could not connect to database: $DBI::errstr";


 #$query = "UPDATE clusters SET profdensity=? WHERE cid=?;";
 #$sth = $dbh->prepare($query);


 # update cluster entries via temp file
 print STDERR "Updating cluster table...\n"; 


 $sth = $dbh->prepare('CREATE TEMPORARY TABLE clustemp (`cid` INT(11) NOT NULL PRIMARY KEY, `pd` float NOT NULL) ENGINE=INNODB;')
                or die "Couldn't prepare statement: " . $dbh->errstr;
 $sth->execute()             # Execute the query
            or die "Couldn't execute statement: " . $sth->errstr;
 $sth->finish;

 $sth = $dbh->prepare('ALTER TABLE clustemp DISABLE KEYS;')
                or die "Couldn't prepare statement: " . $dbh->errstr;
 $sth->execute();
 $sth->finish;

 $sth = $dbh->prepare("LOAD DATA LOCAL INFILE '$tmp/pcd_$DBNAME.txt' INTO TABLE clustemp FIELDS TERMINATED BY ',' LINES TERMINATED BY '\n';")
                or die "Couldn't prepare statement: " . $dbh->errstr;


 my $TEMPFILE;
 open ($TEMPFILE, ">$tmp/pcd_$DBNAME.txt") or die $!;


 opendir(DIR, $folder);
 @tarballs = grep(/(?:edges)$/, readdir(DIR));
 closedir(DIR);
 $tarball_count = @tarballs;

 $i=0;
 my $id="0";
 my $ds=0.0;
 my $pcd=0;
 my $pcdupd=0;
 foreach (@tarballs) {
   $i++;
   #print $i.". ". $_ ."\n";
   if (/(\d+)\.leb36\.edges/i) {
     $id = $1;
   }

   open (MYFILE, "$folder/$_") or die;
   while (<MYFILE>) { 


        if ($_ =~ /-(\d+)(['"]),(\d+),(\d+\.\d+)/i) {

	   # do a check here if the read-ref pair is in the map file
	   #$sth2->execute($1,$3);
 	   #$num = $sth2->rows;
	   #if ($num>=1) 
           { 
		#print $1.$2." ".$3." ".$4."\n";
		$ds = $4;	
                if ($ds >= $MINPROFSCORE) {

  		    if (exists $SHASH{$3}) {
			if ($ds > $SHASH{$3}) {
	                        $RHASH{$3} = $1;
                        	$SHASH{$3} = $ds;
			} elsif ($ds == $SHASH{$3}) {
	                        $RHASH{$3} .= (",".$1);
                        	$SHASH{$3} = $ds;
			}			
  		    } else {
			$RHASH{$3} = $1;
			$SHASH{$3} = $ds;
   		    }
                }
	   }		

	} elsif ($_ =~ /ave: (\d+\.\d+)/i) {
		#print $id. ": ". $1."\n";
		#$sth->execute($1,$id);

		$pcd++;
		print $TEMPFILE "$id,$1\n";

                if ($pcd % $RECORDS_PER_INFILE_INSERT == 0) {
                  close($TEMPFILE);
                  $pcdupd += $sth->execute();
                  open ($TEMPFILE, ">$tmp/pcd_$DBNAME.txt") or die $!;
                }


	}
    }
   close MYFILE;;

 }
 
 # finish insert and update based on temp table
 close($TEMPFILE);
 $pcdupd += $sth->execute();
 $sth->finish;
 unlink("$tmp/pcd_$DBNAME.txt");

 $sth = $dbh->prepare('ALTER TABLE clustemp ENABLE KEYS;')
                or die "Couldn't prepare statement: " . $dbh->errstr;
 $sth->execute();
 $sth->finish;

 my $updfromtable = 0;
 $query = "UPDATE clustemp p, clusters pp SET pp.profdensity=p.pd WHERE pp.cid = p.cid;";
 $sth = $dbh->prepare($query);
 $updfromtable=$sth->execute();
 $sth->finish;



 # populate rank table 
 print STDERR "Populating rank table...\n"; 
 $i=0; $j=0; $rankins=0;


 # prepare rank table
 $query = "TRUNCATE TABLE rank;";
 $sth = $dbh->prepare($query);
 $sth->execute();
 $sth->finish;
 

 $sth = $dbh->prepare('ALTER TABLE rank DISABLE KEYS;')
                or die "Couldn't prepare statement: " . $dbh->errstr;
 $sth->execute()             # Execute the query
            or die "Couldn't execute statement: " . $sth->errstr;
 $sth->finish;


 # insert rank entries via temp file
 open ($TEMPFILE, ">$tmp/ranktemp_$DBNAME.txt") or die $!;

 $query = "LOAD DATA LOCAL INFILE '$tmp/ranktemp_$DBNAME.txt' INTO TABLE rank FIELDS TERMINATED BY ',' LINES TERMINATED BY '\n';";
 $sth = $dbh->prepare($query)
                or die "Couldn't prepare statement: " . $dbh->errstr;

 while (($key,$value) =  each(%RHASH)) {
	$i++;
	#print "$key => $value (". $SHASH{$key}. ")\n";
        my @pieces = split(/,/, $value);
        my $ties = scalar(@pieces) - 1;
	foreach my $ps (@pieces) {

     	  $j++; 
	  print $TEMPFILE $ps,",",$key,",",$SHASH{$key},",",$ties,"\n";

          if ($j % $RECORDS_PER_INFILE_INSERT == 0) {
             close($TEMPFILE);
             $rankins += $sth->execute();
	     open ($TEMPFILE, ">$tmp/ranktemp_$DBNAME.txt") or die $!;
          }
	}
 }

 # finish insert
 close($TEMPFILE);
 $rankins += $sth->execute();
 $sth->finish;


 # enable keys in rank, set stats, print msg
 $sth = $dbh->prepare('ALTER TABLE rank ENABLE KEYS;')
                or die "Couldn't prepare statement: " . $dbh->errstr;
 $sth->execute()             # Execute the query
            or die "Couldn't execute statement: " . $sth->errstr;
 $sth->finish;


 SetStatistics("RANK_EDGES_OVERCUTOFF",$j);
 print STDERR "Inserted ($j) rank records for $i reads.\n"; 




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


 open ($TEMPFILE, ">$tmp/ranktemp_$DBNAME.txt") or die $!;


 print STDERR "Prunning (keep best ref TR for each read TR) from rank table...\n"; 
 $query = "Select refid,readid,sid,score FROM rank INNER JOIN replnk ON rank.readid=replnk.rid ORDER BY readid, score;";
 $sth = $dbh->prepare($query);
 $sth->execute();
 $num = $sth->rows;
 $i=0;
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
 print STDERR "Prunning complete. Pruned $count rank records.\n"; 
 SetStatistics("RANK_REMOVED_SAMEREF",$count);


 print STDERR "Prunning (one TR/same read) from rank table...\n"; 
 $query = "Select refid,readid,sid,score FROM rank INNER JOIN replnk ON rank.readid=replnk.rid ORDER BY refid, sid, score, readid;"; # readid added for tie resolution to keep rank and rankflank entries more in sync
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
 print STDERR "Prunning complete. Pruned $count rank records.\n"; 
 SetStatistics("RANK_REMOVED_SAMESEQ",$count);


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


 # delete from rankd based on temptable entries
 my $delfromtable = 0;
 $query = "DELETE FROM t1 USING rank t1 INNER JOIN ranktemp t2 ON ( t1.refid = t2.refid AND t1.readid = t2.readid );";
 $sth = $dbh->prepare($query);
 $delfromtable=$sth->execute();
 $sth->finish;


 # cleanup temp file
 unlink("$tmp/ranktemp_$DBNAME.txt");


# set old db settings
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


if ($updfromtable != $pcdupd) { die "Updated number of cluster entries($updfromtable) not equal to the number of inserted counter ($pcdupd), aborting! You might need to rerun from step 12."; }
if ($rankins != $j) { die "Inserted number of rank entries($rankins) not equal to the number of inserted counter ($j), aborting! You might need to rerun from step 12."; }
if ($delfromtable != $count) { die "Deleted number of entries($delfromtable) not equal to the number of deleted counter ($count), aborting! You might need to rerun from step 12."; }


$dbh->disconnect();


print STDERR "Finished. Deleted from rank using temptable: $delfromtable\n";

($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
print STDERR sprintf("\n\nend: %4d-%02d-%02d %02d:%02d:%02d\n\n",$year+1900,$mon+1,$mday,$hour,$min,$sec);


1;


############################ Procedures ###############################################################

sub fork_proclu {
        if ($files_processed >= $tarball_count) {
                return 0;
        }

        # wait for shared variables to unlock
        while ($MYLOCK) { }

        # lock shared vars
        $MYLOCK = 1;

        # use a predefined number of files
        my $until = $files_processed + $files_to_process - 1;
        $until = $tarball_count - 1 if $until > ($tarball_count - 1);
        print STDERR 'Processing files '. ($files_processed + 1) . ' to '. ($until + 1) . "\n";
        #my $output_prefix = "$root/$files_processed-$until";
        my @file_slice = @tarballs[($files_processed)..($until)];
        my $file_slice_count = @file_slice;
        $files_processed += $file_slice_count;
        my $proclu_string;

        # unlock shared vars
        $MYLOCK = 0;

        defined(my $pid = fork)
                or die "Unable to fork: $!\n";

        # child
        if ($pid == 0) {

                foreach (@file_slice) {
                  #print STDERR "\t" . $_ . "\n";
                 
                  my $edgesin = $_;
		  $edgesin =~ s/leb36/edgein/; 

		  my $reffile = $_;
		  $reffile =~ s/read/ref/; 

                  #>/dev/null";
                  #>$_.log";
                  $proclu_string = $PROCLU . " " .  $_ . " "  . $reffile . $PROCLU_PARAM  . "$edgesin > /dev/null";
                  #$proclu_string = $PROCLU . " " .  $_ . " "  . $reffile . $PROCLU_PARAM  . "$edgesin > ${edgesin}.proclu_log";
                  #print STDERR $proclu_string."\n";
                  #exit(1);
                  system($proclu_string);
                  if ( $? == -1 ) { die "command failed: $!\n"; }
                  else {
                    my $rc = ($? >> 8);
		    # empty clusters cause nonzero return code, cause some files missing
                    # if ( 0 != $rc ) { print "proclu returned $rc ( $proclu_string  )!"; exit($rc); }
                  }


                }
                #print STDERR "\n";

                # child must never return
                exit 0;

        # parent
        } else {
                return $pid;
        }

 return 0;
}




