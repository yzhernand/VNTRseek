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

use lib "$FindBin::RealBin/lib";
use vutil qw(get_config get_dbh set_statistics get_trunc_query);

print strftime( "\n\nstart: %F %T\n\n\n", localtime );

my $argc = @ARGV;

if ( $argc < 3 ) { die "Usage: map_dup.pl dbname msdir tempdir\n"; }

my $curdir = getcwd;

# TODO Better default or calculate in advance
my $maxRepeatsPerRead = 2;

my $DBSUFFIX = $ARGV[0];
my $MSDIR    = $ARGV[1];
my $TEMPDIR  = $ARGV[2];

# set these mysql credentials in vs.cnf (in installation directory)
my %run_conf = get_config( $MSDIR . "vs.cnf" );
my ( $LOGIN, $PASS, $HOST ) = @run_conf{qw(LOGIN PASS HOST)};
my $dbh = get_dbh( $DBSUFFIX, $MSDIR . "vs.cnf" )
    or die "Could not connect to database: $DBI::errstr";

#my $sth3 = $dbh->prepare("UPDATE map SET bbb=0 WHERE refid=? AND readid=?;")
#                or die "Couldn't prepare statement: " . $dbh->errstr;

# create temp table for updates
my $query = q{CREATE TEMPORARY TABLE mduptemp (
    `refid` INT(11) NOT NULL,
    `readid` INT(11) NOT NULL,
    PRIMARY KEY (refid,readid)
    )};

$dbh->do($query)
    or die "Couldn't do statement: " . $dbh->errstr;

my $sth;
my $TEMPFILE;
if ( $run_conf{BACKEND} eq "mysql" ) {
    $dbh->do('ALTER TABLE mduptemp DISABLE KEYS;')
        or die "Couldn't do statement: " . $dbh->errstr;

    $dbh->do('SET AUTOCOMMIT = 0;')
        or die "Couldn't do statement: " . $dbh->errstr;

    $dbh->do('SET FOREIGN_KEY_CHECKS = 0;')
        or die "Couldn't do statement: " . $dbh->errstr;

    $dbh->do('SET UNIQUE_CHECKS = 0;')
        or die "Couldn't do statement: " . $dbh->errstr;
    open( $TEMPFILE, ">", "$TEMPDIR/mduptemp_$DBSUFFIX.txt" ) or die $!;
}
elsif ( $run_conf{BACKEND} eq "sqlite" ) {
    $dbh->do("PRAGMA foreign_keys = OFF");
    # warn "\nTurning off AutoCommit\n";
    $dbh->{AutoCommit} = 0;
}

my $trsInRead_sth
    = $dbh->prepare(q{SELECT COUNT(*)
    FROM map
    INNER JOIN replnk on replnk.rid=map.readid
    INNER JOIN rank ON rank.refid=map.refid AND rank.readid=map.readid
    INNER JOIN rankflank ON rankflank.refid=map.refid AND rankflank.readid=map.readid
    INNER JOIN fasta_ref_reps ON fasta_ref_reps.rid=map.refid
    WHERE replnk.sid=? AND bbb=1
    ORDER BY rank.score ASC, rankflank.score ASC, map.refid ASC})
    or die "Couldn't prepare statement: " . $dbh->errstr;
$trsInRead_sth->execute();
my $numTRsInRead = $trsInRead_sth->rows;

$trsInRead_sth
    = $dbh->prepare(q{SELECT map.readid,map.refid,(
        SELECT head
        FROM fasta_reads
        WHERE fasta_reads.sid=replnk.sid
    ),rank.score,rankflank.score,fasta_ref_reps.head AS refhead,fasta_ref_reps.firstindex,fasta_ref_reps.lastindex,(
        SELECT length(DNA)
        FROM fasta_reads
        WHERE fasta_reads.sid=replnk.sid
    )
    FROM map
    INNER JOIN replnk on replnk.rid=map.readid
    INNER JOIN rank ON rank.refid=map.refid AND rank.readid=map.readid
    INNER JOIN rankflank ON rankflank.refid=map.refid AND rankflank.readid=map.readid
    INNER JOIN fasta_ref_reps ON fasta_ref_reps.rid=map.refid
    WHERE replnk.sid=? AND bbb=1
    ORDER BY rank.score ASC, rankflank.score ASC, map.refid ASC})
    or die "Couldn't prepare statement: " . $dbh->errstr;

# first get list of reads with multiple TRs that are mapped to more than one reference (sorted by read id)
my $deleted      = 0;
my $ReadsDeleted = 0;
my $readsWithMultTRsMappedMultRefs_sth
    = $dbh->prepare(q{SELECT replnk.sid,count(map.refid)
    FROM map
    INNER JOIN replnk on replnk.rid=map.readid
    WHERE bbb=1
    GROUP BY replnk.sid
    HAVING count(map.refid)>1})
    or die "Couldn't prepare statement: " . $dbh->errstr;
$readsWithMultTRsMappedMultRefs_sth->execute()
    or die "Cannot execute: " . $readsWithMultTRsMappedMultRefs_sth->errstr();
my $insert_mduptemp_sth = $dbh->prepare(q{INSERT INTO mduptemp VALUES(?, ?)});
my $i = 0;
while ( my @data = $readsWithMultTRsMappedMultRefs_sth->fetchrow_array() ) {

    # TODO Use read length to calculate $maxRepeatsPerRead
    $i++;

    $trsInRead_sth->execute( $data[0] );
    my $j          = 0;
    my $oldp       = -1;
    my $oldf       = -1;
    my $oldreadid  = -1;
    my $oldrefid   = -1;
    my $oldrefhead = "";
    my $oldRfirst  = -1;
    my $oldRlast   = -1;
    my $isDeleted  = 0;
    while ( my @data2   = $trsInRead_sth->fetchrow_array() ) {
        $j++;
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
# TODO Make this compare all adjacent pairs (not just TRs 1 and 2 in the read) for longer reads.
        if ($j == 2
            && (  !( $data2[5] eq $oldrefhead && $RefDiff <= $readlen )
                || ( $numTRsInRead > $maxRepeatsPerRead ) )
            )
        {
            #$sth3->execute($oldrefid,$oldreadid);
            if ( $run_conf{BACKEND} eq "mysql" ) {
                print $TEMPFILE $oldrefid, ",", $oldreadid, "\n";
            }
            elsif ( $run_conf{BACKEND} eq "sqlite" ) {
                $insert_mduptemp_sth->execute($oldrefid, $oldreadid);
            }
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
            if ( $run_conf{BACKEND} eq "mysql" ) {
                print $TEMPFILE $data2[1], ",", $data2[0], "\n";
            }
            elsif ( $run_conf{BACKEND} eq "sqlite" ) {
                $insert_mduptemp_sth->execute($data2[1], $data2[0]);
            }
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
                if ( $run_conf{BACKEND} eq "mysql" ) {
                    print $TEMPFILE $data2[1], ",", $data2[0], "\n";
                }
                elsif ( $run_conf{BACKEND} eq "sqlite" ) {
                    $insert_mduptemp_sth->execute($data2[1], $data2[0]);
                }
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

my $numReadsWithMultTRsMappedMultRefs
    = $readsWithMultTRsMappedMultRefs_sth->rows;
# $sth3->finish();
$trsInRead_sth->finish();
$readsWithMultTRsMappedMultRefs_sth->finish();
if ( $run_conf{BACKEND} eq "mysql" ) {
    # load the file into tempfile
    close($TEMPFILE);
    $dbh->do(q{LOAD DATA LOCAL INFILE '$TEMPDIR/mduptemp_$DBSUFFIX.txt'
        INTO TABLE mduptemp FIELDS TERMINATED BY ',' LINES TERMINATED BY '\n'})
    or die "Couldn't do statement: " . $dbh->errstr;
    $dbh->do('ALTER TABLE mduptemp ENABLE KEYS;')
        or die "Couldn't do statement: " . $dbh->errstr;
}
elsif ( $run_conf{BACKEND} eq "sqlite" ) {
    $dbh->commit;
}


# update based on temp table
my $updfromtable = 0;
$query = q{UPDATE map SET bbb=0
    WHERE EXISTS (
    SELECT refid FROM mduptemp t2
    WHERE map.refid = t2.refid AND map.readid=t2.readid
)};
$updfromtable = $dbh->do($query)
    or die "Couldn't do statement: " . $dbh->errstr;

# set old db settings
if ( $run_conf{BACKEND} eq "mysql" ) {
    $sth = $dbh->do('SET AUTOCOMMIT = 1;')
        or die "Couldn't do statement: " . $dbh->errstr;

    $sth = $dbh->do('SET FOREIGN_KEY_CHECKS = 1;')
        or die "Couldn't do statement: " . $dbh->errstr;

    $sth = $dbh->do('SET UNIQUE_CHECKS = 1;')
        or die "Couldn't do statement: " . $dbh->errstr;
    # cleanup temp file
    unlink("$TEMPDIR/mduptemp_$DBSUFFIX.txt");
}
elsif ($run_conf{BACKEND} eq "sqlite") {
    $dbh->do("PRAGMA foreign_keys = ON");
    $dbh->{AutoCommit} = 1;
}


if ( $updfromtable != $deleted ) {
    die
        "Updated bbb=0 number of entries($updfromtable) not equal to the number of deleted counter ($deleted), aborting! You might need to rerun from step 12.";
}

# update BBB on stats
$dbh->do(q{UPDATE stats SET BBB=(SELECT count(*) FROM map WHERE bbb=1)})
    or die "Couldn't do statement: " . $dbh->errstr;

$dbh->disconnect();

printf "\n\n%d entries deleted!\n", $deleted;
printf "%d reads deleted!\n",       $ReadsDeleted;

print strftime( "\n\nend: %F %T\n\n\n", localtime );

1;

