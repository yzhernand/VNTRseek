#!/usr/bin/perl

# removes *reads* that are mapped to multiple references (due to multiple TRs). Picks the best TR by profile, when same, picks best by flank. When same, removes both.
# !!!  makes an exception if references are on same chromosome and close together

use strict;
use warnings;
use Cwd;
use DBI;
use List::Util qw[min max];
use POSIX qw(strftime);

use FindBin;
use File::Basename;

use lib "$FindBin::Bin/vntr";
require "vutil.pm";

use vutil qw(get_credentials write_mysql stats_set);

print strftime( "\n\nstart: %F %T\n\n\n", localtime );

my $argc = @ARGV;

if ( $argc < 3 ) { die "Usage: map_dup.pl dbname msdir tempdir\n"; }

my $curdir            = getcwd;
my $maxRepeatsPerRead = -1;

my $DBNAME  = $ARGV[0];
my $MSDIR   = $ARGV[1];
my $TEMPDIR = $ARGV[2];

# set these mysql credentials in vs.cnf (in installation directory)
my ( $LOGIN, $PASS, $HOST ) = get_credentials($MSDIR);

my $dbh = DBI->connect( "DBI:mysql:$DBNAME;mysql_local_infile=1;host=$HOST",
    "$LOGIN", "$PASS" )
    or die "Could not connect to database: $DBI::errstr";

#my $sth3 = $dbh->prepare("UPDATE map SET bbb=0 WHERE refid=? AND readid=?;")
#                or die "Couldn't prepare statement: " . $dbh->errstr;

# create temp table for updates
my $sth
    = $dbh->prepare(
    'CREATE TEMPORARY TABLE mduptemp (`refid` INT(11) NOT NULL, `readid` INT(11) NOT NULL, PRIMARY KEY (refid,readid)) ENGINE=INNODB;'
    ) or die "Couldn't prepare statement: " . $dbh->errstr;
$sth->execute()    # Execute the query
    or die "Couldn't execute statement: " . $sth->errstr;
$sth->finish;

$sth = $dbh->prepare('ALTER TABLE mduptemp DISABLE KEYS;')
    or die "Couldn't prepare statement: " . $dbh->errstr;
$sth->execute()    # Execute the query
    or die "Couldn't execute statement: " . $sth->errstr;
$sth->finish;

$sth = $dbh->prepare('SET AUTOCOMMIT = 0;')
    or die "Couldn't prepare statement: " . $dbh->errstr;
$sth->execute()    # Execute the query
    or die "Couldn't execute statement: " . $sth->errstr;
$sth->finish;

$sth = $dbh->prepare('SET FOREIGN_KEY_CHECKS = 0;')
    or die "Couldn't prepare statement: " . $dbh->errstr;
$sth->execute()    # Execute the query
    or die "Couldn't execute statement: " . $sth->errstr;
$sth->finish;

$sth = $dbh->prepare('SET UNIQUE_CHECKS = 0;')
    or die "Couldn't prepare statement: " . $dbh->errstr;
$sth->execute()    # Execute the query
    or die "Couldn't execute statement: " . $sth->errstr;
$sth->finish;

my $TEMPFILE;
open( $TEMPFILE, ">", "$TEMPDIR/mduptemp_$DBNAME.txt" ) or die $!;

my $trsInRead_sth
    = $dbh->prepare(
    "SELECT map.readid,map.refid,(SELECT head FROM fasta_reads WHERE fasta_reads.sid=replnk.sid),rank.score,rankflank.score,fasta_ref_reps.head as refhead,fasta_ref_reps.firstindex,fasta_ref_reps.lastindex,(SELECT length(DNA) FROM fasta_reads WHERE fasta_reads.sid=replnk.sid) FROM map INNER JOIN replnk on replnk.rid=map.readid INNER JOIN rank ON rank.refid=map.refid AND rank.readid=map.readid INNER JOIN rankflank ON rankflank.refid=map.refid AND rankflank.readid=map.readid INNER JOIN fasta_ref_reps ON fasta_ref_reps.rid=map.refid WHERE replnk.sid=? AND bbb=1 ORDER BY rank.score ASC, rankflank.score ASC, map.refid ASC;"
    ) or die "Couldn't prepare statement: " . $dbh->errstr;

# first get list of reads with multiple TRs that are mapped to more than one reference (sorted by read id)
my $deleted      = 0;
my $ReadsDeleted = 0;
my $readsWithMultTRsMappedMultRefs_sth
    = $dbh->prepare(
    "SELECT replnk.sid,count(map.refid) FROM map INNER JOIN replnk on replnk.rid=map.readid WHERE bbb=1 GROUP BY replnk.sid HAVING count(map.refid)>1;"
    ) or die "Couldn't prepare statement: " . $dbh->errstr;
$readsWithMultTRsMappedMultRefs_sth->execute()
    or die "Cannot execute: " . $readsWithMultTRsMappedMultRefs_sth->errstr();
my $numReadsWithMultTRsMappedMultRefs = $readsWithMultTRsMappedMultRefs_sth->rows;
my $i                            = 0;
while ( $i < $numReadsWithMultTRsMappedMultRefs ) {
    my @data = $readsWithMultTRsMappedMultRefs_sth->fetchrow_array();
    $i++;

    $trsInRead_sth->execute( $data[0] );
    my $numTRsInRead = $trsInRead_sth->rows;

    my $j          = 0;
    my $oldp       = -1;
    my $oldf       = -1;
    my $oldreadid  = -1;
    my $oldrefid   = -1;
    my $oldrefhead = "";
    my $oldRfirst  = -1;
    my $oldRlast   = -1;
    my $isDeleted  = 0;
    while ( $j < $numTRsInRead ) {
        $j++;
        my @data2   = $trsInRead_sth->fetchrow_array();
        my $readlen = $data2[8];
        my $RefDiff
            = ( $oldRfirst == -1 || $oldRlast == -1 )
            ? 1000000
            : (
            max( $oldRlast, $data2[7] ) - min( $oldRfirst, $data2[6] ) + 1 );

        if ( $j == 1 ) {
            print "\n"
                . $i
                . ". read='"
                . $data2[2]
                . "' refs=("
                . $data[1] . ")";
        }

        # 1st entry can only be deleted due to NOT being on same chromosome and close together as 2nd entry (or more than $maxRepeatsPerRead entries exist)
        if ($j == $maxRepeatsPerRead
            && (  !( $data2[5] eq $oldrefhead && $RefDiff <= $readlen )
                || ( $numTRsInRead > $maxRepeatsPerRead ) )
            )
        {
            #$sth3->execute($oldrefid,$oldreadid);
            print $TEMPFILE $oldrefid, ",", $oldreadid, "\n";
            print "X";
            $deleted++;
            $isDeleted = 1;
        }

        print "\n\t"
            . $data2[0] . "->"
            . $data2[1]
            . " P=$data2[3] F=$data2[4] ";

        # delete every entry (except first which is deleted in another block) if more than $maxRepeatsPerRead exist
        if ( $j > 1 && $numTRsInRead > $maxRepeatsPerRead ) {

            #$sth3->execute($data2[1],$data2[0]);
            print $TEMPFILE $data2[1], ",", $data2[0], "\n";
            print "X";
            $deleted++;
            $isDeleted = 1;
        }

        # else we are on 2nd etnry, delete 2nd entry if previous entry score is equal to this entry
        #elsif ($j==2 && $data2[3]==$oldp && $data2[4]==$oldf) {
        elsif ( $j == $maxRepeatsPerRead ) {

            # make an exception (DO NOT DELETE) if references are on same chromosome and close together as previous entry
            if ( $data2[5] eq $oldrefhead && $RefDiff <= $readlen ) {

            }
            else {

                #$sth3->execute($data2[1],$data2[0]);
                print $TEMPFILE $data2[1], ",", $data2[0], "\n";
                print "X";
                $deleted++;
                $isDeleted = 1;
            }
        }

        $oldrefhead = $data2[5];
        $oldRfirst  = $data2[6];
        $oldRlast   = $data2[7];
        $oldp       = $data2[3];
        $oldf       = $data2[4];
        $oldreadid  = $data2[0];
        $oldrefid   = $data2[1];
    }

    if ($isDeleted) { $ReadsDeleted++; }

}

# $sth3->finish();
$trsInRead_sth->finish();
$readsWithMultTRsMappedMultRefs_sth->finish();

# load the file into tempfile
close($TEMPFILE);
$sth
    = $dbh->prepare(
    "LOAD DATA LOCAL INFILE '$TEMPDIR/mduptemp_$DBNAME.txt' INTO TABLE mduptemp FIELDS TERMINATED BY ',' LINES TERMINATED BY '\n';"
    ) or die "Couldn't prepare statement: " . $dbh->errstr;
$sth->execute();
$sth->finish;
$sth = $dbh->prepare('ALTER TABLE mduptemp ENABLE KEYS;')
    or die "Couldn't prepare statement: " . $dbh->errstr;
$sth->execute();
$sth->finish;

# update based on temp table
my $updfromtable = 0;
my $query
    = "UPDATE mduptemp p, map pp SET bbb=0 WHERE pp.refid = p.refid AND pp.readid = p.readid;";
$sth          = $dbh->prepare($query);
$updfromtable = $sth->execute();
$sth->finish;

# set old db settings
$sth = $dbh->prepare('SET AUTOCOMMIT = 1;')
    or die "Couldn't prepare statement: " . $dbh->errstr;
$sth->execute()    # Execute the query
    or die "Couldn't execute statement: " . $sth->errstr;
$sth->finish;

$sth = $dbh->prepare('SET FOREIGN_KEY_CHECKS = 1;')
    or die "Couldn't prepare statement: " . $dbh->errstr;
$sth->execute()    # Execute the query
    or die "Couldn't execute statement: " . $sth->errstr;
$sth->finish;

$sth = $dbh->prepare('SET UNIQUE_CHECKS = 1;')
    or die "Couldn't prepare statement: " . $dbh->errstr;
$sth->execute()    # Execute the query
    or die "Couldn't execute statement: " . $sth->errstr;
$sth->finish;

# cleanup temp file
unlink("$TEMPDIR/mduptemp_$DBNAME.txt");

if ( $updfromtable != $deleted ) {
    die
        "Updated bbb=0 number of entries($updfromtable) not equal to the number of deleted counter ($deleted), aborting! You might need to rerun from step 12.";
}

# update BBB on stats
$sth = $dbh->prepare(
    "UPDATE stats SET BBB=(SELECT count(*) FROM map WHERE bbb=1);");
$sth->execute();
$sth->finish();

$dbh->disconnect();

printf "\n\n%d entries deleted!\n", $deleted;
printf "%d reads deleted!\n",       $ReadsDeleted;

print strftime( "\n\nend: %F %T\n\n\n", localtime );

1;

