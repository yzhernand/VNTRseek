#!/usr/bin/perl

# takes .index.seq files and calls pcr_dup program to remove duplicates

use strict;
use warnings;
use Cwd;
use DBI;
use List::Util qw[min max];

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

# Perl trim function to remove whitespace from the start and end of the string
sub trim($) {
    my $string = shift;
    $string =~ s/^\s+//;
    $string =~ s/\s+$//;
    return $string;
}

( $sec, $min, $hour, $mday, $mon, $year, $wday, $yday, $isdst )
    = localtime(time);
printf "\n\nstart: %4d-%02d-%02d %02d:%02d:%02d\n\n\n", $year + 1900,
    $mon + 1, $mday, $hour, $min, $sec;

my $argc = @ARGV;

if ( $argc < 6 ) {
    die
        "Usage: pcr_dup.pl indexfolder profcleanfolder dbname msdir cpucount tempdir ignorepcrdups\n";
}

my $curdir = getcwd;

my $indexfolder  = $ARGV[0];
my $pcleanfolder = $ARGV[1];
my $DBNAME       = $ARGV[2];
my $MSDIR        = $ARGV[3];
my $cpucount     = $ARGV[4];
my $TEMPDIR      = $ARGV[5];
my $IGNOREPCRDUP = $ARGV[6];

# set these mysql credentials in vs.cnf (in installation directory)
my ( $LOGIN, $PASS, $HOST ) = get_credentials($MSDIR);

####################################
sub SetStatistics {

    my $argc = @_;
    if ( $argc < 2 ) {
        die "stats_set: expects 2 parameters, passed $argc !\n";
    }

    my $NAME  = $_[0];
    my $VALUE = $_[1];

    #print "$DBNAME,$LOGIN,$PASS,$NAME,$VALUE\n";
    return stats_set( $DBNAME, $LOGIN, $PASS, $HOST, $NAME, $VALUE );
}

#goto AAA;

# process
print STDERR "reading: $indexfolder";

opendir( DIR, $indexfolder );
my @indexfiles = grep( /\.(?:seq)$/, readdir(DIR) );
closedir(DIR);

#my $i=0;
#foreach my $ifile (@indexfiles) {
# $i++;
# print "\n $i. ". $ifile."...";

# my $exstring = "./pcr_dup.exe $indexfolder/$ifile $indexfolder/$ifile.pcr_dup 2 3 > /dev/null";
# #my $exstring = "./pcr_dup.exe $indexfolder/$ifile $indexfolder/$ifile.pcr_dup 2 3 > $indexfolder/$ifile.pcr_log";
# system($exstring);
#}

my $files_to_process = 100;    # number of files to process in one batch
my $files_processed  = 0;      # files processed
my %p;                         # associates forked pids with output pipe pids

my $MYLOCK = 0;

my $tarball_count = @indexfiles;
print STDERR "$tarball_count supported files found in $indexfolder\n";

#die "Exiting\n" if $tarball_count == 0;
$files_to_process = $tarball_count if $files_to_process > $tarball_count;

# fork as many new processes as there are CPUs
for ( my $i = 0; $i < $cpucount; $i++ ) { $p{ fork_pcrdup() } = 1 }

# wait for processes to finish and then fork new ones
while ( ( my $pid = wait ) != -1 ) {
    if ( $p{$pid} ) {

        # one instance has finished processing -- start a new one
        delete $p{$pid};
        $p{ fork_pcrdup() } = 1;
    }
    else {
        die "************ Process $pid finished (not in hash)\n";
    }
}

print STDERR
    "Processing complete -- processed $files_processed cluster(s).\n";

AAA:

# load results
print STDERR "reading: $indexfolder";

my $dbh = DBI->connect( "DBI:mysql:$DBNAME;mysql_local_infile=1;host=$HOST",
    "$LOGIN", "$PASS" )
    || die "Could not connect to database: $DBI::errstr";

# first count the intersect before pcr dup
my $rrintersect = 0;
my $sth
    = $dbh->prepare(
    "SELECT count(*) FROM rank INNER JOIN rankflank ON rank.refid=rankflank.refid AND rank.readid=rankflank.readid;"
    ) or die "Couldn't prepare statement: " . $dbh->errstr;
$sth->execute() or die "Cannot execute: " . $sth->errstr();
my $num = $sth->rows;
if ($num) {
    my @data = $sth->fetchrow_array();
    $rrintersect = $data[0];
}
$sth->finish();
SetStatistics( "INTERSECT_RANK_AND_RANKFLANK_BEFORE_PCR", $rrintersect );

# deleteing PCR DUPS in database
my %PENTRIES = ();

# create temp table for updates
$sth
    = $dbh->prepare(
    'CREATE TEMPORARY TABLE pduptemp (`refid` INT(11) NOT NULL, `readid` INT(11) NOT NULL, PRIMARY KEY (refid,readid)) ENGINE=INNODB;'
    ) or die "Couldn't prepare statement: " . $dbh->errstr;
$sth->execute()    # Execute the query
    or die "Couldn't execute statement: " . $sth->errstr;
$sth->finish;

$sth = $dbh->prepare('ALTER TABLE pduptemp DISABLE KEYS;')
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
open( $TEMPFILE, ">$TEMPDIR/pduptemp_$DBNAME.txt" ) or die $!;

#my $query = "DELETE FROM rank WHERE refid=? AND readid=?;";
#$sth = $dbh->prepare($query);

#$query = "DELETE FROM rankflank WHERE refid=? AND readid=?;";
#my $sth1 = $dbh->prepare($query);

opendir( DIR, $indexfolder );
@indexfiles = grep( /\.(?:pcr_dup)$/, readdir(DIR) );
closedir(DIR);

my $i       = 0;
my $deleted = 0;
foreach my $ifile (@indexfiles) {

    $i++;

    if ( $ifile =~ /(\d+)\.seq\.pcr_dup/ ) {

        my $ref        = $1;
        my $filedelted = 0;
        print "\n $i. " . $ifile . "...";

        if ( open( FILE, "$indexfolder/$ifile" ) ) {

            my %RHASH = ();

            # added at 1.02, to eliminate most connected nodes preferentiably
            my %RCOUNTS = ();
            my %NEWIDS  = ();
            while (<FILE>) {
                if (/^compare: (\d+) (\d+) (\d+) (\d+)\|(\d+)/) {
                    $RCOUNTS{$1}++;
                    $RCOUNTS{$5}++;
                }
            }
            my @keys = sort { $RCOUNTS{$a} <=> $RCOUNTS{$b} or $a <=> $b }
                keys %RCOUNTS;
            my $iditer = 1;
            foreach my $key (@keys) { $NEWIDS{$key} = $iditer; $iditer++; }

            seek FILE, 0, 0;
            while (<FILE>) {
                if (/^compare: (\d+) (\d+) (\d+) (\d+)\|(\d+)/) {

                    #my $read = max($1,$5);

                    my $read;
                    if ( $NEWIDS{$1} > $NEWIDS{$5} ) {
                        $read = $1;
                    }
                    elsif ( $NEWIDS{$1} < $NEWIDS{$5} ) {
                        $read = $5;
                    }
                    else {
                        $read = max( $1, $5 );
                    }

                    if ( !exists $RHASH{$read} ) {

                        #$sth->execute($ref,$read);
                        #$sth1->execute($ref,$read);
                        print $TEMPFILE $ref, ",", $read, "\n";

                        #print "deleting: -$ref -> $read\n";
                        $deleted++;
                        $filedelted++;

                  #print "$1 (newid:$NEWIDS{$1}) - $5 (newid: $NEWIDS{$5})\n";
                  #exit(1);

                        $PENTRIES{ $ref . "_" . $read } = 1;

                    }

                    $RHASH{$read} = 1;

                }
            }

            close(FILE);
            print "($filedelted deleted)";

        }

    }

}

#$sth->finish;
#$sth1->finish;

# load the file into tempfile
close($TEMPFILE);
$sth
    = $dbh->prepare(
    "LOAD DATA LOCAL INFILE '$TEMPDIR/pduptemp_$DBNAME.txt' INTO TABLE pduptemp FIELDS TERMINATED BY ',' LINES TERMINATED BY '\n';"
    ) or die "Couldn't prepare statement: " . $dbh->errstr;
$sth->execute();
$sth->finish;
$sth = $dbh->prepare('ALTER TABLE pduptemp ENABLE KEYS;')
    or die "Couldn't prepare statement: " . $dbh->errstr;
$sth->execute();
$sth->finish;

# delete from rankdflank based on temptable entries
my $delfromtable = 0;
my $query
    = "DELETE FROM t1 USING rank t1 INNER JOIN pduptemp t2 ON ( t1.refid = t2.refid AND t1.readid = t2.readid );";
$sth          = $dbh->prepare($query);
$delfromtable = $sth->execute();
$sth->finish;
$query
    = "DELETE FROM t1 USING rankflank t1 INNER JOIN pduptemp t2 ON ( t1.refid = t2.refid AND t1.readid = t2.readid );";
$sth = $dbh->prepare($query);
$sth->execute();
$sth->finish;

# cleanup temp file
unlink("$TEMPDIR/pduptemp_$DBNAME.txt");

print "\n\nProcessing complete (pcr_dup.pl), deleted $deleted duplicates.\n";
SetStatistics( "RANK_REMOVED_PCRDUP",      $deleted );
SetStatistics( "RANKFLANK_REMOVED_PCRDUP", $deleted );

# for accounting of pcr dups
print STDERR "\nMaking a list of pcr_dup removed...\n";
if ( open( FILE, ">$pcleanfolder/result/$DBNAME.pcr_dup.txt" ) ) {
    $i = 0;
    for my $key ( sort keys %PENTRIES ) {
        $i++;
        print FILE "$i\t-" . $key . "\n";
    }
    print STDERR "PCR_DUP list complete with $i removed entries.\n";
    close(FILE);
}

# first count the intersect
$rrintersect = 0;
$sth
    = $dbh->prepare(
    "SELECT count(*) FROM rank INNER JOIN rankflank ON rank.refid=rankflank.refid AND rank.readid=rankflank.readid;"
    ) or die "Couldn't prepare statement: " . $dbh->errstr;
$sth->execute() or die "Cannot execute: " . $sth->errstr();
$num = $sth->rows;
if ($num) {
    my @data = $sth->fetchrow_array();
    $rrintersect = $data[0];
}
$sth->finish();
SetStatistics( "INTERSECT_RANK_AND_RANKFLANK", $rrintersect );

# now exclude ties, mark in map table and record the number
print STDERR "Updating BEST BEST BEST map entries...\n";
$sth = $dbh->prepare("UPDATE map SET bbb=0;");    # clear all bbb entries
$sth->execute();
$sth->finish;

$sth = $dbh->prepare("truncate table pduptemp;"); # clear all pduptemp entries
$sth->execute();
$sth->finish;

$sth = $dbh->prepare('ALTER TABLE pduptemp DISABLE KEYS;')
    or die "Couldn't prepare statement: " . $dbh->errstr;
$sth->execute();
$sth->finish;

$sth
    = $dbh->prepare(
    "INSERT INTO pduptemp SELECT map.refid, map.readid FROM map INNER JOIN rank ON rank.refid=map.refid AND rank.readid=map.readid INNER JOIN rankflank ON rankflank.refid=map.refid AND rankflank.readid=map.readid WHERE rank.ties=0 OR rankflank.ties=0;"
    );                                            # clear all pduptemp entries
$sth->execute();
$sth->finish;

$sth = $dbh->prepare('ALTER TABLE pduptemp ENABLE KEYS;')
    or die "Couldn't prepare statement: " . $dbh->errstr;
$sth->execute();
$sth->finish;

$query
    = "UPDATE pduptemp p, map pp SET bbb=1 WHERE pp.refid = p.refid AND pp.readid=p.readid;";
$sth = $dbh->prepare($query);
$i   = $sth->execute();
$sth->finish;

SetStatistics( "BBB_WITH_MAP_DUPS", $i );

# make a list of ties
print STDERR "Making a list of ties (references)...\n";
if ( open( FILE, ">$pcleanfolder/result/$DBNAME.ties.txt" ) ) {

    $query
        = "SELECT map.refid, max(bbb) as mbb, (select head from fasta_ref_reps where rid=map.refid) as chr,(select firstindex from fasta_ref_reps where rid=map.refid) as tind  FROM map INNER JOIN rank ON rank.refid=map.refid AND rank.readid=map.readid INNER JOIN rankflank ON rankflank.refid=map.refid AND rankflank.readid=map.readid GROUP BY map.refid HAVING mbb=0 ORDER BY chr, tind;";
    $sth = $dbh->prepare($query);
    $sth->execute();
    $num = $sth->rows;
    $i   = 0;
    while ( $i < $num ) {
        my @data = $sth->fetchrow_array();
        $i++;
        print FILE "$i\t-"
            . $data[0] . "\t"
            . $data[2] . "\t"
            . $data[3] . "\n";
    }
    $sth->finish;
    close(FILE);
}
print STDERR "Ties list complete with $i removed references.\n";

# make a list of ties
print STDERR "\nMaking a list of ties (entries)...\n";
if ( open( FILE, ">$pcleanfolder/result/$DBNAME.ties_entries.txt" ) ) {

    $query
        = "SELECT map.refid, map.readid,rank.ties,rankflank.ties  FROM map INNER JOIN rank ON rank.refid=map.refid AND rank.readid=map.readid INNER JOIN rankflank ON rankflank.refid=map.refid AND rankflank.readid=map.readid WHERE bbb=0 ORDER BY map.refid,map.readid;";
    $sth = $dbh->prepare($query);
    $sth->execute();
    $num = $sth->rows;
    $i   = 0;
    while ( $i < $num ) {
        my @data = $sth->fetchrow_array();
        $i++;
        print FILE "$i\t-"
            . $data[0] . "\t"
            . $data[1] . "\t"
            . $data[2] . "\t"
            . $data[3] . "\n";
    }
    $sth->finish;
    close(FILE);
}
print STDERR "Ties list complete with $i removed entries.\n";

# set old settings
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

if ( $delfromtable != $deleted ) {
    die
        "Deleted number of entries($delfromtable) not equal to the number of deleted counter ($deleted), aborting! You might need to rerun from step 12.";
}

$dbh->disconnect();

( $sec, $min, $hour, $mday, $mon, $year, $wday, $yday, $isdst )
    = localtime(time);
printf "\n\nend: %4d-%02d-%02d %02d:%02d:%02d\n", $year + 1900, $mon + 1,
    $mday, $hour, $min, $sec;

1;

############################ Procedures ###############################################################

sub fork_pcrdup {
    if ( $files_processed >= $tarball_count ) {
        return 0;
    }

    # wait for shared variables to unlock
    while ($MYLOCK) { }

    # lock shared vars
    $MYLOCK = 1;

    # use a predefined number of files
    my $until = $files_processed + $files_to_process - 1;
    $until = $tarball_count - 1 if $until > ( $tarball_count - 1 );
    print STDERR 'Processing files '
        . ( $files_processed + 1 ) . ' to '
        . ( $until + 1 ) . "\n";

    #my $output_prefix = "$root/$files_processed-$until";
    my @file_slice = @indexfiles[ ($files_processed) .. ($until) ];
    my $file_slice_count = @file_slice;
    $files_processed += $file_slice_count;
    my $exstring;

    # unlock shared vars
    $MYLOCK = 0;

    defined( my $pid = fork )
        or die "Unable to fork: $!\n";

    # child
    if ( $pid == 0 ) {

        foreach (@file_slice) {

            #print STDERR "\t" . $_ . "\n";

            my $exstring
                = "./pcr_dup.exe $indexfolder/${_} $indexfolder/${_}.pcr_dup 0 2 $IGNOREPCRDUP > /dev/null";
            system($exstring);

        }

        #print STDERR "\n";

        # child must never return
        exit 0;

        # parent
    }
    else {
        return $pid;
    }

    return 0;
}

