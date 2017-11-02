#!/usr/bin/perl

my $RECORDS_PER_INFILE_INSERT = 100000;

use List::Util qw[min max];

use strict;
use warnings;
use Cwd;
use DBI;
use POSIX qw(strftime);
use FindBin;
use File::Basename;

use lib "$FindBin::RealBin/lib";

use vutil qw(get_config get_dbh set_statistics get_trunc_query);

my $updatedClustersCount = 0;
my $updatedRefsCount     = 0;

sub nowhitespace($) {
    my $string = shift;
    $string =~ s/\s+//g;
    return $string;
}

warn strftime( "\n\nstart: %F %T\n\n", localtime );

my $curdir = getcwd;

my $argc = @ARGV;
if ( $argc < 5 ) {
    die
        "Usage: run_rankflankmap.pl inputfile  mapdir tmpdir dbsuffix msdir\n";
}

my $inputfile = $ARGV[0];
my $mapdir    = $ARGV[1];
my $tmp       = $ARGV[2];
my $DBSUFFIX  = $ARGV[3];
my $MSDIR     = $ARGV[4];

# set these mysql credentials in vs.cnf (in installation directory)
my %run_conf = get_config( $MSDIR . "vs.cnf" );
my ( $LOGIN, $PASS, $HOST ) = @run_conf{qw(LOGIN PASS HOST)};

my $clusters_processed = 0;

my $dbh = get_dbh( $DBSUFFIX, $MSDIR . "vs.cnf" )
    or die "Could not connect to database: $DBI::errstr";

my $sth;
my $sth1;
my $map_insert_sth;
my $rankflank_insert_sth;

#goto AAA;

# clear map
$dbh->do( get_trunc_query( $run_conf{BACKEND}, "map" ) )
    or die "Couldn't do statement: " . $dbh->errstr;

# clear rankflank
$dbh->do( get_trunc_query( $run_conf{BACKEND}, "rankflank" ) )
    or die "Couldn't do statement: " . $dbh->errstr;

# disable indices
if ( $run_conf{BACKEND} eq "mysql" ) {
    $dbh->do('ALTER TABLE map DISABLE KEYS;')
        or die "Couldn't do statement: " . $dbh->errstr;

    $dbh->do('ALTER TABLE rankflank DISABLE KEYS;')
        or die "Couldn't do statement: " . $dbh->errstr;

    $dbh->do('SET AUTOCOMMIT = 0;')
        or die "Couldn't do statement: " . $dbh->errstr;

    $dbh->do('SET FOREIGN_KEY_CHECKS = 0;')
        or die "Couldn't do statement: " . $dbh->errstr;

    $dbh->do('SET UNIQUE_CHECKS = 0;')
        or die "Couldn't do statement: " . $dbh->errstr;

    # prepare statments
    $map_insert_sth
        = $dbh->prepare(
        "LOAD DATA LOCAL INFILE '/$tmp/${DBSUFFIX}_map.txt' INTO TABLE map FIELDS TERMINATED BY ',' LINES TERMINATED BY '\n';"
        ) or die "Couldn't prepare statement: " . $dbh->errstr;
    $rankflank_insert_sth
        = $dbh->prepare(
        "LOAD DATA LOCAL INFILE '/$tmp/${DBSUFFIX}_rankflank.txt' INTO TABLE rankflank FIELDS TERMINATED BY ',' LINES TERMINATED BY '\n';"
        ) or die "Couldn't prepare statement: " . $dbh->errstr;
}
elsif ( $run_conf{BACKEND} eq "sqlite" ) {
    $dbh->do("PRAGMA foreign_keys = OFF");
    $dbh->{AutoCommit} = 0;
    $map_insert_sth = $dbh->prepare(
        qq{INSERT INTO map (refid, readid, reserved, reserved2)
        VALUES (?, ?, 0, 0)}
    ) or die "Couldn't prepare statement: " . $dbh->errstr;
    $rankflank_insert_sth
        = $dbh->prepare(qq{INSERT INTO rankflank VALUES (?, ?, ?, ?)})
        or die "Couldn't prepare statement: " . $dbh->errstr;
}

open FILE, "<$inputfile" or die "error opening for reading '$inputfile': $!";

opendir( DIR, $mapdir );
my @allfiles = readdir(DIR);
closedir(DIR);

my $j = 0;
my $k = 0;
my $upload;
my $uploadedrank = 0;
my $uploadedmap  = 0;

# open out files
my $MAPFILE;
my $RFFILE;

if ( $run_conf{BACKEND} eq "mysql" ) {
    open $MAPFILE, ">/$tmp/${DBSUFFIX}_map.txt"       or die $!;
    open $RFFILE,  ">/$tmp/${DBSUFFIX}_rankflank.txt" or die $!;
}

#foreach my $file (@files) {
foreach my $file (@allfiles) {

    if ( index( $file, "_map" ) ) {

        $clusters_processed++;

        my $mfile = "$mapdir/$file";
        open MFILE, "<$mfile" or die $!;

        my %READVECTOR = ();
        while (<MFILE>) {
            chomp;
            my @mfields = split( '\t', $_ );
            my $msize = scalar @mfields;
            if ( $msize >= 8 ) {
                my $readid = $mfields[0];

                #print "\n\n".$_."\n\n";
                #exit(1);

                my @rfields = split( ',', $mfields[7] );

                #my @temparray = ();

                my $bestscore = 0;
                my $bestref   = "";

                foreach my $refstr (@rfields) {
                    if ( $refstr =~ /^-?(\d+):(\d+):(\d+)/ ) {

                        #print "\n".$1." ".$2." ".$3."\n";
                        #exit(1);

                        # calculate score
                        my $score;

                        if ( ( $mfields[2] + $mfields[3] )
                            == 0
                            ) # if no flanks, it will only be marked best if nothing else is available
                        {
                            $score = 0;
                        }
                        else

                   #{ # asked to add by Dr. Benson to balance out small flanks

                         #  my $A = ($2 + $3) / ( $mfields[2] + $mfields[3] );
                         #  my $B = 1 / ( $mfields[2] + $mfields[3]);
                         #  $score = 1 - max($A,$B);
                         #}

                        {
                            $score = 1 - (
                                ( $2 + $3 ) / ( $mfields[2] + $mfields[3] ) );
                        }

               # filter to remove all flank scores below .9, added Nov 5, 2012
                        if ( $score >= 0.90 ) {

                            if ( $score > $bestscore ) {
                                $bestref   = $1;
                                $bestscore = $score;
                            }
                            elsif ( $score == $bestscore ) {
                                if ( $bestref eq "" ) { $bestref = $1; }
                                else { $bestref .= ( "," . $1 ); }
                            }
                        }

       # create a map entry in database (added 11/19/2010)
       #$map_insert_sth->execute($1,$readid)             # Execute the query
       #      or die "Couldn't execute statement: " . $map_insert_sth->errstr;
                        $k++;
                        if ( $run_conf{BACKEND} eq "mysql" ) {
                            print $MAPFILE "$1,$readid\n";

                            if ( $k % $RECORDS_PER_INFILE_INSERT == 0 ) {
                                close($MAPFILE);
                                $upload = $map_insert_sth->execute();
                                $uploadedmap += $upload;

                                open( $MAPFILE, ">/$tmp/${DBSUFFIX}_map.txt" )
                                    or die $!;
                            }
                        }
                        elsif ( $run_conf{BACKEND} eq "sqlite" ) {
                            $map_insert_sth->execute( $1, $readid );
                            $uploadedmap++;
                        }
                    }
                }

                # insert the rankflank
                if ( $bestref ne "" ) {
                    my @ranks = split( ',', $bestref );
                    my $ties = scalar(@ranks) - 1;
                    foreach my $rstr (@ranks) {

#$rankflank_insert_sth->execute($readid,$rstr,$bestscore)             # Execute the query
#      or die "Couldn't execute statement: " . $rankflank_insert_sth->errstr;
                        $j++;
                        if ( $run_conf{BACKEND} eq "mysql" ) {
                            print $RFFILE "$rstr,$readid,$bestscore,$ties\n";

                            if ( $j % $RECORDS_PER_INFILE_INSERT == 0 ) {
                                close($RFFILE);
                                $upload = $rankflank_insert_sth->execute();
                                $uploadedrank += $upload;
                                open( $RFFILE,
                                    ">/$tmp/${DBSUFFIX}_rankflank.txt" )
                                    or die $!;
                            }
                        }
                        elsif ( $run_conf{BACKEND} eq "sqlite" ) {
                            $rankflank_insert_sth->execute( $rstr, $readid,
                                $bestscore, $ties );
                            $uploadedrank++;
                        }

                    }
                }

                #$READVECTOR{$readid} = @temparray;

                #print "\n";
            }
        }

        # load the files and remove the temp files

        print STDERR "\nprocessed: $clusters_processed";

    }    # end of if (index($file,"_map")) {

}    # end of foreach @files
close(FILE);

if ( $run_conf{BACKEND} eq "mysql" ) {

    # finish writing and loading out files
    close($RFFILE);
    close($MAPFILE);

    $upload = $map_insert_sth->execute();
    $uploadedmap += $upload;

    $upload = $rankflank_insert_sth->execute();
    $uploadedrank += $upload;
    unlink("/$tmp/${DBSUFFIX}_map.txt");
    unlink("/$tmp/${DBSUFFIX}_rankflank.txt");
    print STDERR "Enabling indices...\n";

    # enable indices
    $dbh->do('ALTER TABLE map ENABLE KEYS;')
        or die "Couldn't do statement: " . $dbh->errstr;

    $dbh->do('ALTER TABLE rankflank ENABLE KEYS;')
        or die "Couldn't do statement: " . $dbh->errstr;
}

if ( $uploadedmap != $k ) {
    die
        "\nUploaded number of map entries($uploadedmap) not equal to the number of uploaded counter ($k), aborting!";
}
if ( $uploadedrank != $j ) {
    die
        "\nUploaded number of rankflank entries($uploadedrank) not equal to the number of uploaded counter ($j), aborting!";
}

AAA:

# create temp table for deletions
if ( $run_conf{BACKEND} eq "mysql" ) {
    $dbh->do(
        'CREATE TEMPORARY TABLE ranktemp (
        `refid` INT(11) NOT NULL,
        `readid` INT(11) NOT NULL,
        PRIMARY KEY (refid,readid)
        ) ENGINE=INNODB;'
    ) or die "Couldn't do statement: " . $dbh->errstr;

    $dbh->do('ALTER TABLE ranktemp DISABLE KEYS;')
        or die "Couldn't do statement: " . $dbh->errstr;
}
elsif ( $run_conf{BACKEND} eq "sqlite" ) {
    $dbh->do(
        'CREATE TEMPORARY TABLE ranktemp (
        `refid` integer NOT NULL,
        `readid` integer NOT NULL,
        PRIMARY KEY (`refid`, `readid`))'
    ) or die "Couldn't do statement: " . $dbh->errstr;
}

my $TEMPFILE;

if ( $run_conf{BACKEND} eq "mysql" ) {
    open( $TEMPFILE, ">$tmp/ranktemp_$DBSUFFIX.txt" ) or die $!;
}
elsif ( $run_conf{BACKEND} eq "sqlite" ) {
    $sth1 = $dbh->prepare(
        qq{
        INSERT INTO ranktemp VALUES (?, ?);
    }
    );
}
print STDERR
    "Prunning (keep best ref for each read) from rankflank table...\n";
my $query
    = "Select refid,readid,sid,score FROM rankflank INNER JOIN replnk ON rankflank.readid=replnk.rid ORDER BY readid, score;";
$sth = $dbh->prepare($query);
$sth->execute();
my $i        = 0;
my $count    = 0;
my $oldseq   = -1;
my $oldref   = -1;
my $oldread  = -1;
my $oldscore = -1.0;

while ( my @data = $sth->fetchrow_array() ) {
    if ( $data[1] == $oldread && $data[3] != $oldscore ) {

        # delete old one
        if ( $run_conf{BACKEND} eq "mysql" ) {
            print $TEMPFILE $oldref, ",", $oldread, "\n";
        }
        elsif ( $run_conf{BACKEND} eq "sqlite" ) {
            $sth1->execute( $oldref, $oldread );
        }

        $count++;
    }
    $oldref   = $data[0];
    $oldread  = $data[1];
    $oldseq   = $data[2];
    $oldscore = $data[3];
    $i++;
}

my $num = $sth->rows;
$sth->finish;

print STDERR "Prunning complete. Pruned $count rankflank records.\n";

print STDERR "Prunning all (one TR/same read) rankflank table...\n";
$query
    = "Select refid,readid,sid,score FROM rankflank INNER JOIN replnk ON rankflank.readid=replnk.rid ORDER BY refid, sid, score, readid;"
    ; # readid added for tie resolution to keep rank and rankflank entries more in sync
$sth = $dbh->prepare($query);
$sth->execute();
$i       = 0;
$count   = 0;
$oldseq  = -1;
$oldref  = -1;
$oldread = -1;
while ( my @data = $sth->fetchrow_array() ) {

    if ( $data[0] == $oldref && $data[2] == $oldseq ) {

        # delete old one
        if ( $run_conf{BACKEND} eq "mysql" ) {
            print $TEMPFILE $oldref, ",", $oldread, "\n";
        }
        elsif ( $run_conf{BACKEND} eq "sqlite" ) {
            $sth1->execute( $oldref, $oldread );
        }

        $count++;
    }
    $oldref  = $data[0];
    $oldread = $data[1];
    $oldseq  = $data[2];
    $i++;
}

$num = $sth->rows;
$sth->finish;

# load the file into tempfile
if ( $run_conf{BACKEND} eq "mysql" ) {
    close($TEMPFILE);
    $dbh->do(
        "LOAD DATA LOCAL INFILE '$tmp/ranktemp_$DBSUFFIX.txt' INTO TABLE ranktemp FIELDS TERMINATED BY ',' LINES TERMINATED BY '\n';"
    ) or die "Couldn't do statement: " . $dbh->errstr;
    $dbh->do('ALTER TABLE ranktemp ENABLE KEYS;')
        or die "Couldn't do statement: " . $dbh->errstr;

    # cleanup temp file
    unlink("$tmp/ranktemp_$DBSUFFIX.txt");
}

print STDERR "Prunning complete. Pruned $count rankflank records.\n";

# delete from rankflank based on temptable entries
my $delfromtable = 0;
if ( $run_conf{BACKEND} eq "mysql" ) {
    $query = qq{
    DELETE FROM t1 USING rankflank t1
    INNER JOIN ranktemp t2
        ON ( t1.refid = t2.refid AND t1.readid = t2.readid )
    };
}
elsif ( $run_conf{BACKEND} eq "sqlite" ) {
    $query = qq{
        DELETE FROM rankflank
        WHERE EXISTS (
            SELECT * FROM ranktemp t2
            WHERE rankflank.refid = t2.refid
                AND rankflank.readid = t2.readid
        )
    };
}
$sth          = $dbh->prepare($query);
$delfromtable = $sth->execute();

# $sth->finish;

# set old settings
if ( $run_conf{BACKEND} eq "mysql" ) {
    $sth = $dbh->do('SET AUTOCOMMIT = 1;')
        or die "Couldn't do statement: " . $dbh->errstr;
    $sth = $dbh->do('SET FOREIGN_KEY_CHECKS = 1;')
        or die "Couldn't do statement: " . $dbh->errstr;

    $sth = $dbh->do('SET UNIQUE_CHECKS = 1;')
        or die "Couldn't do statement: " . $dbh->errstr;
}
elsif ( $run_conf{BACKEND} eq "sqlite" ) {
    $dbh->do("PRAGMA foreign_keys = ON");
    $dbh->{AutoCommit} = 1;
}

if ( $delfromtable != $count ) {
    die
        "Deleted number of entries($delfromtable) not equal to the number of deleted counter ($count), aborting!";
}

$dbh->disconnect();
set_statistics(
    $DBSUFFIX,
    (   "RANKFLANK_EDGES_INSERTED"  => $j,
        "RANKFLANK_REMOVED_SAMEREF" => $count,
        "RANKFLANK_REMOVED_SAMESEQ" => $count,
    )
);

print STDERR "\n\n";

print STDERR
    "Processing complete -- processed $clusters_processed cluster(s). Deleted from rankflank using temptable: $delfromtable\n";

warn strftime( "\n\nend: %F %T\n\n", localtime );

1;

