#!/usr/bin/perl

my $RECORDS_PER_INFILE_INSERT = 100000;

use List::Util qw[min max];


use strict;
use warnings;
use Cwd;
use DBI;

sub nowhitespace($)
{
        my $string = shift;
        $string =~ s/\s+//g;
        return $string;
}

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



($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
print STDERR sprintf("\n\nstart: %4d-%02d-%02d %02d:%02d:%02d\n\n",$year+1900,$mon+1,$mday,$hour,$min,$sec);

my $curdir =  getcwd;

my $argc = @ARGV;
if ($argc<6) { die "Usage: run_variability.pl inputfile  mapdir dbname msdir minflank tempdir\n"; }


my $inputfile = $ARGV[0];
my $mapdir = $ARGV[1];
my $DBNAME = $ARGV[2];
my $MSDIR = $ARGV[3];
my $MIN_FLANK_REQUIRED = $ARGV[4];
my $TEMPDIR = $ARGV[5];

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


my $dbh = DBI->connect("DBI:mysql:$DBNAME;mysql_local_infile=1;host=$HOST", "$LOGIN", "$PASS"
	           ) || die "Could not connect to database: $DBI::errstr";

my $sth;
my $sth1;
my $sth7;
my $sth8;


# change settings to speedup updates and inserts
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



# update reserved field on entire table
$sth = $dbh->prepare('UPDATE clusterlnk SET reserved=0,reserved2=0;')
                or die "Couldn't prepare statement: " . $dbh->errstr;
$sth->execute()             # Execute the query
            or die "Couldn't execute statement: " . $sth->errstr;
$sth->finish;
$sth = $dbh->prepare('UPDATE map SET reserved=0,reserved2=0;')
                or die "Couldn't prepare statement: " . $dbh->errstr;
$sth->execute()             # Execute the query
            or die "Couldn't execute statement: " . $sth->errstr;
$sth->finish;
$sth = $dbh->prepare('truncate table vntr_support;')
                or die "Couldn't prepare statement: " . $dbh->errstr;
$sth->execute()             # Execute the query
            or die "Couldn't execute statement: " . $sth->errstr;
$sth->finish;


 # create temp table for updates
 $sth = $dbh->prepare('CREATE TEMPORARY TABLE mapr (`refid` INT(11) NOT NULL, `readid` INT(11) NOT NULL, PRIMARY KEY (refid,readid)) ENGINE=INNODB;')
                or die "Couldn't prepare statement: " . $dbh->errstr;
 $sth->execute()             # Execute the query
            or die "Couldn't execute statement: " . $sth->errstr;
 $sth->finish;

 $sth = $dbh->prepare('ALTER TABLE mapr DISABLE KEYS;')
                or die "Couldn't prepare statement: " . $dbh->errstr;
 $sth->execute()             # Execute the query
            or die "Couldn't execute statement: " . $sth->errstr;
 $sth->finish;


# prepare statments
$sth = $dbh->prepare('SELECT rid,flankleft,sequence,flankright,pattern,copynum,(lastindex-firstindex+1) as arlen FROM fasta_ref_reps WHERE rid = ?')
                or die "Couldn't prepare statement: " . $dbh->errstr;
$sth1 = $dbh->prepare('SELECT rid, dna, first, last, pattern, copynum FROM fasta_reads INNER JOIN replnk ON fasta_reads.sid=replnk.sid WHERE rid = ?')
                or die "Couldn't prepare statement: " . $dbh->errstr;
$sth7 = $dbh->prepare("LOAD DATA LOCAL INFILE '$TEMPDIR/mapr_$DBNAME.txt' INTO TABLE mapr FIELDS TERMINATED BY ',' LINES TERMINATED BY '\n';")
                or die "Couldn't prepare statement: " . $dbh->errstr;

$sth8 = $dbh->prepare('SELECT map.refid, map.readid from clusterlnk INNER JOIN map ON map.refid=-clusterlnk.repeatid WHERE map.bbb=1 AND clusterlnk.clusterid=? ORDER BY map.readid,map.refid;')
                or die "Couldn't prepare statement: " . $dbh->errstr;


#$sth9 = $dbh->prepare('UPDATE vntr_support SET support=support+1 WHERE refid=? AND copies=?;')
#                or die "Couldn't prepare statement: " . $dbh->errstr;
#$sth10 = $dbh->prepare('INSERT IGNORE INTO vntr_support(refid,copies,sameasref,support,copiesfloat,representative) VALUES (?,?,?,?,?,?);')
#                or die "Couldn't prepare statement: " . $dbh->errstr;



my(%VNTR_REF) = ();
my(%VNTR_COPIES) = ();
my(%VNTR_COPIESFLOAT) = ();
my(%VNTR_SUPPORT) = ();
my(%VNTR_SAMEASREF) = ();
my(%VNTR_REPRESENTATIVE) = ();


my $SUPPORT_INCREMENTED = 0;

my $processed=0;
my $TEMPFILE;
open ($TEMPFILE, ">$TEMPDIR/mapr_$DBNAME.txt") or die $!;

my $cl_processed=0;
my $TEMP_CLNK;
open ($TEMP_CLNK, ">$TEMPDIR/clnk_$DBNAME.txt") or die $!;


# this function adds support to reference copynumber count 
sub add_support {

 my $refid = $_[0];
 my $sameasref = $_[1];
 my $copies = $_[2];
 my $copiesfloat = $_[3];
 my $representative = $_[4];

 my $key = $refid."_".$copies;

 if (exists $VNTR_SUPPORT{$key}) {

   $VNTR_SUPPORT{$key}++;

 } else {

   $VNTR_REF{$key} = $refid;
   $VNTR_SAMEASREF{$key} = $sameasref;
   $VNTR_COPIES{$key} = $copies;
   $VNTR_COPIESFLOAT{$key} = $copiesfloat;
   $VNTR_REPRESENTATIVE{$key} = $representative;
   $VNTR_SUPPORT{$key} = 1;
 } 

 $SUPPORT_INCREMENTED++;

 return 0;
}

# shows that there is 0 support, ignores if already an entry for this reference
sub add_zero_support {

 my $refid = $_[0];
 my $sameasref = $_[1];
 my $copies = $_[2];
 my $copiesfloat = $_[3];

 my $key = $refid."_".$copies;

 if (!exists $VNTR_SUPPORT{$key}) {

   $VNTR_REF{$key} = $refid;
   $VNTR_SAMEASREF{$key} = $sameasref;
   $VNTR_COPIES{$key} = $copies;
   $VNTR_COPIESFLOAT{$key} = $copiesfloat;
   $VNTR_SUPPORT{$key} = 0;
 }

 return 0;
}


# this is a function that detects if reads have (somewhat) different number of copy numbers then reference(s)
sub VNTR_YES_NO {

 
 my(%refs) = %{$_[0]};
 my(%reads) = %{$_[1]};
 my(%readhash) = %{$_[2]};
 my(%newrefs) = %{$_[3]}; # this is to make sure we have an etnry in vntr_support for all BBB refs
 my $clusterid = $_[4];

 my $varyes = 0;
 my $change = 0;

 my(%REF_UPDATED) = ();

 #print "\nVNTR_YES_NO:";

 # for each read in map file
 while ( my ($key, @temp) = each(%readhash) ) {

	my $valread;
	my $valref;


   if (exists $reads{$key}) {
	
  	$valread = $reads{$key};

        #print STDERR "\n$key($valread) =>";


	# for each references listed in map file as being associated with the read
        #foreach my $val (@{$temp[0]}) {
        foreach my $val (@{ $readhash{$key}  }) {



	  if (exists $refs{$val}) {
		$valref = $refs{$val};
	  } else {
		print "Undefined reference `$val`. Aborting!\n";
            	exit(1);
  	  }



  	  #print STDERR " ".$val."($valref)";

	  # print "$valRead - $valRef\n";
    	  #if (($valread >= ($valref-.5)) && ($valread <= ($valref+.5))) {
    	  if (($valread > ($valref-.8)) && ($valread < ($valref+.8))) {

		# add support
		#add_support($val, 1, int(abs($valref) + 0.5), $valref );
		add_support($val, 1, 0, $valref, $key );
 

	  } else {

		$varyes = 1;

		# add support
		#add_support($val, 0, int(abs($valread) + 0.5) , 0.0);

                my $fchange;

                if ($valread > $valref)
                   { $fchange =  1 + int($valread - ($valref+.8)); } # truncation
                else 
                   { $fchange = -1 - int(($valref-.8) - $valread); }   # truncation

		add_support($val, 0, $fchange , 0.0, $key);

                # this ignored if there is already an entry for this ref
                # this would be incremented to 0 if there is a read found later on
                # add_zero_support($val, 1, int(abs($valref) + 0.5) );

		#print "\t$clusterid: $refid\n";		

		$change = int(abs($valread-$valref) + 0.5);

		if (!exists $REF_UPDATED{$val} || $change > $REF_UPDATED{$val}) {

  		  #$sth3->execute($change,$clusterid,$val)             # set the reserved field on reference
	          #   or die "Couldn't execute statement: " . $sth3->errstr;

	 	  $REF_UPDATED{$val}=$change;
		}

		#$sth7->execute(-$val,$key)             # set the reserved field on map read, (added 11/19/2010)
	        #    or die "Couldn't execute statement: " . $sth7->errstr;

                $processed++;

                print $TEMPFILE -$val,",",$key,"\n";

                if ($processed % $RECORDS_PER_INFILE_INSERT == 0) {
                       close($TEMPFILE);
                       $sth7->execute();
                       open ($TEMPFILE, ">$TEMPDIR/mapr_$DBNAME.txt") or die $!;
                }



	  }
	
	}
	#print "\n";


    } # end of exists

 } # end of loop

 # create self support for refs
 foreach my $val ( keys %newrefs ) {
   my $valref = $refs{$val};
   #add_zero_support($val, 1, int(abs($valref) + 0.5), $valref );
   add_zero_support($val, 1, 0, $valref );
 }


 # write to clustmp file (to update clusterlnk entries)
 foreach my $val ( keys %REF_UPDATED ) {
   $cl_processed++;
   $change = $REF_UPDATED{$val};
   print $TEMP_CLNK "$clusterid,$val,$change\n";
 }


 # print "\nNOT VARIABLE\n\n";
 return $varyes;

} # end of func





open FILE, "<$inputfile" or die $!;



while (<FILE>) { 
 $clusters_processed++;

 my @values = split(',', $_);

 my %REFCOPIES = ();
 my %ASKLENGTH = ();
 my %READCOPIES = ();
 my $repeatcount = 0;
 my $refcount = 0;
 my $readcount = 0;

 my $readlen;
 my $first;
 my $last;

 foreach my $val (@values) {

   my $dir = '\'';
   if ($val =~ m/([\'\"])/) { $dir = $1; } 

   chomp($val);
 
   $val =~ s/[\'\"]//g;

   $repeatcount++;

   # go though all refs for each read
   if ($val<=0) {

     $refcount++;	

     #$val = ($val<0)? -$val : $val; # ids are positive in database

     $sth->execute(-$val)             # Execute the query
            or die "Couldn't execute statement: " . $sth->errstr;

     my @data = $sth->fetchrow_array();

     if ($sth->rows == 0) {
            print STDERR "No record in database for entry `$val'. Aborting!\n\n";
 	    exit(1);
     }
   

     $ASKLENGTH{$val} = int($data[6]); # remember length of the ref array sequence
     $REFCOPIES{$val} = $data[5];
  
   } else {

     $readcount++;

     $sth1->execute($val)             # Execute the query
            or die "Couldn't execute statement: " . $sth1->errstr;

     my @data = $sth1->fetchrow_array();

     if ($sth1->rows == 0) {
            print STDERR "No record in database for entry `$val'. Aborting!\n";
 	    exit(1);
     }

     $readlen = length(nowhitespace($data[1]));

     $ASKLENGTH{$val} = $readlen; # remember length of the read

     $first = $data[2];
     $last = $data[3];
     if (($first-1) >= $MIN_FLANK_REQUIRED && ($readlen - $last) >= $MIN_FLANK_REQUIRED) {

        $READCOPIES{$val} = $data[5];
     }



   }

 }

 # get map data
 my %READVECTOR = ();
 my %REFHASH = ();
 my $readidold=-1;
 $sth8->execute($clusters_processed)             # Execute the query
            or die "Couldn't execute statement: " . $sth8->errstr;
 my $nmapped = $sth8->rows;
 for (my $i=0; $i<$nmapped; $i++) {
     my @data = $sth8->fetchrow_array();

     my $readid;
     my $refid;
     my $readlen;

     $refid = -1 * int($data[0]);
     $readid = int($data[1]);
    
     $REFHASH{$refid}=1;

print STDERR "\n\n$refid | $readid \n";

     if ($readid != $readidold) {
       $READVECTOR{$readid} = ();
     }     

     $readlen = $ASKLENGTH{$readid};

 
     if (!exists $ASKLENGTH{$readid})  {
       print STDERR "Entry for $readid does not exist!";
       exit(1);
     }
     if (!exists $ASKLENGTH{$refid})  {
       print STDERR "Entry for $refid does not exist!";
       exit(1);
     }



     #if ($ASKLENGTH{$refid} <= $readlen) {
       push( @{ $READVECTOR{$readid}  }, $refid  );
     #}

     $readidold=$readid;
 }
$sth8->finish;
 print STDERR "\nprocessed: $clusters_processed\n\n";
 


 # do for 1st 10 clusters for now
 # if ($clusters_processed >= 20) { last; }

 # insert database records (cluster table)
 #print "VYN:  ";
 #print  VNTR_YES_NO(\%REFCOPIES,\%READCOPIES,\%READVECTOR);

 my $vYes=VNTR_YES_NO(\%REFCOPIES,\%READCOPIES,\%READVECTOR,\%REFHASH,$clusters_processed);

 $updatedClustersCount += $vYes;


}
close(FILE);
$sth->finish;
$sth1->finish;





 # finish loading the map file into tempfile and switch indexes back on
 close($TEMPFILE);
 $sth7->execute();
 $sth7->finish;

 $sth7 = $dbh->prepare('ALTER TABLE mapr ENABLE KEYS;')
                or die "Couldn't prepare statement: " . $dbh->errstr;
 $sth7->execute();
 $sth7->finish;


 # update map based on temp table
 my $updfromtable = 0;
 my $query = "UPDATE mapr p, map pp SET reserved=1 WHERE pp.refid = p.refid AND pp.readid = p.readid;";
 $sth7 = $dbh->prepare($query);
 $updfromtable=$sth7->execute();
 $sth7->finish;


 # cleanup temp file
 unlink("$TEMPDIR/mapr_$DBNAME.txt");




 # write SUPPORT info to temp files to be loaded into vntr_support
 my $supcounter = 0;
 open ($TEMPFILE, ">$TEMPDIR/support_$DBNAME.txt") or die $!;
 foreach my $key ( keys %VNTR_REF ) {
    $supcounter++;
    print $TEMPFILE $VNTR_REF{$key},",",$VNTR_COPIES{$key},",",$VNTR_SAMEASREF{$key},",",$VNTR_SUPPORT{$key},",",$VNTR_COPIESFLOAT{$key};
    if (exists $VNTR_REPRESENTATIVE{$key}) { print $TEMPFILE ",",$VNTR_REPRESENTATIVE{$key}; }
    print $TEMPFILE "\n";
 }
 close($TEMPFILE);

 $sth = $dbh->prepare('ALTER TABLE vntr_support DISABLE KEYS;')
                or die "Couldn't prepare statement: " . $dbh->errstr;
 $sth->execute();
 $sth->finish;

 $sth = $dbh->prepare("LOAD DATA LOCAL INFILE '$TEMPDIR/support_$DBNAME.txt' INTO TABLE vntr_support FIELDS TERMINATED BY ',' LINES TERMINATED BY '\n';")
                or die "Couldn't prepare statement: " . $dbh->errstr;
 my $supInsert = $sth->execute()             # Execute the query
            or die "Couldn't execute statement: " . $sth->errstr;
 $sth->finish;

 $sth = $dbh->prepare('ALTER TABLE vntr_support ENABLE KEYS;')
                or die "Couldn't prepare statement: " . $dbh->errstr;
 $sth->execute();
 $sth->finish;

 # cleanup temp file
 unlink("$TEMPDIR/support_$DBNAME.txt");





 # create temp table for clusterlnk table updates and update clusterlnk table
 $sth = $dbh->prepare('CREATE TEMPORARY TABLE ctrlnk (`clusterid` INT(11) NOT NULL, `repeatid` INT(11) NOT NULL, `change` INT(11) NOT NULL,PRIMARY KEY (clusterid,repeatid)) ENGINE=INNODB;')
                or die "Couldn't prepare statement: " . $dbh->errstr;
 $sth->execute()             # Execute the query
            or die "Couldn't execute statement: " . $sth->errstr;
 $sth->finish;

 $sth = $dbh->prepare('ALTER TABLE ctrlnk DISABLE KEYS;')
                or die "Couldn't prepare statement: " . $dbh->errstr;
 $sth->execute()             # Execute the query
            or die "Couldn't execute statement: " . $sth->errstr;
 $sth->finish;

 close($TEMP_CLNK);

 $sth = $dbh->prepare("LOAD DATA LOCAL INFILE '$TEMPDIR/clnk_$DBNAME.txt' INTO TABLE ctrlnk FIELDS TERMINATED BY ',' LINES TERMINATED BY '\n';")
                or die "Couldn't prepare statement: " . $dbh->errstr;
 $sth->execute()             # Execute the query
            or die "Couldn't execute statement: " . $sth->errstr;
 $sth->finish;

 unlink("$TEMPDIR/clnk_$DBNAME.txt");

 $sth = $dbh->prepare('ALTER TABLE ctrlnk ENABLE KEYS;')
                or die "Couldn't prepare statement: " . $dbh->errstr;
 $sth->execute()             # Execute the query
            or die "Couldn't execute statement: " . $sth->errstr;
 $sth->finish;

 $sth = $dbh->prepare('UPDATE ctrlnk p, clusterlnk pp SET reserved=p.change WHERE pp.clusterid = p.clusterid AND pp.repeatid=p.repeatid;')
                or die "Couldn't prepare statement: " . $dbh->errstr;
 my $updCLNKfromfile = $sth->execute()             # Execute the query
            or die "Couldn't execute statement: " . $sth->errstr;
 $sth->finish;


SetStatistics("CLUST_NUMBER_OF_REFS_WITH_PREDICTED_VNTR", $updCLNKfromfile );
SetStatistics("CLUST_NUMBER_OF_CLUSTERS_WITH_PREDICTED_VNTR",$updatedClustersCount);


 # create temp table for cluster table updates and update cluster table
 $sth = $dbh->prepare('CREATE TEMPORARY TABLE ctr (`clusterid` INT(11) NOT NULL PRIMARY KEY, `varbl` INT(11) NOT NULL DEFAULT 0) ENGINE=INNODB;')
                or die "Couldn't prepare statement: " . $dbh->errstr;
 $sth->execute()             # Execute the query
            or die "Couldn't execute statement: " . $sth->errstr;
 $sth->finish;

 $sth = $dbh->prepare('ALTER TABLE ctr DISABLE KEYS;')
                or die "Couldn't prepare statement: " . $dbh->errstr;
 $sth->execute()             # Execute the query
            or die "Couldn't execute statement: " . $sth->errstr;
 $sth->finish;

 $sth = $dbh->prepare('INSERT INTO ctr SELECT clusterid, count(*) as vrefs FROM clusterlnk WHERE reserved>0 GROUP by clusterid ORDER BY vrefs DESC')
                or die "Couldn't prepare statement: " . $dbh->errstr;
 my $InsClusToFile = $sth->execute()             # Execute the query
            or die "Couldn't execute statement: " . $sth->errstr;
 $sth->finish;

 $sth = $dbh->prepare('ALTER TABLE ctr ENABLE KEYS;')
                or die "Couldn't prepare statement: " . $dbh->errstr;
 $sth->execute()             # Execute the query
            or die "Couldn't execute statement: " . $sth->errstr;
 $sth->finish;

 $sth = $dbh->prepare('UPDATE ctr p, clusters pp SET variability=varbl WHERE pp.cid = p.clusterid;')
                or die "Couldn't prepare statement: " . $dbh->errstr;
 my $UpdClusFromFile = $sth->execute()             # Execute the query
            or die "Couldn't execute statement: " . $sth->errstr;
 $sth->finish;



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




$dbh->disconnect();


if ($supInsert != $supcounter) { die "Inserted number of vntr_support entries($supInsert) not equal to the number of inserted counter ($supcounter), aborting!"; }
if ($updfromtable != $processed) { die "Updated number of map entries($updfromtable) not equal to the number of inserted counter ($processed), aborting!"; }
if ($updCLNKfromfile != $cl_processed) { die "Updated number of cluster entries($updCLNKfromfile) not equal to the number of inserted counter ($cl_processed), aborting!"; }
if ($UpdClusFromFile != $InsClusToFile) { die "Updated number of clusterlnk entries($UpdClusFromFile) not equal to the number of inserted counter ($InsClusToFile), aborting!"; }


print STDERR "\n\n";

print STDERR "Processing complete -- processed $clusters_processed cluster(s), support entries created = $supInsert.\n";

($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
print STDERR sprintf("\n\nend: %4d-%02d-%02d %02d:%02d:%02d\n\n",$year+1900,$mon+1,$mday,$hour,$min,$sec);



1;




