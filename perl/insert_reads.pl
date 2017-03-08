#!/usr/bin/perl

my $RECORDS_PER_INFILE_INSERT = 100000;

sub read_file_line {
    my $fh = shift;

    if ( $fh and my $line = <$fh> ) {
        chomp $line;
        return $line;
    }
    return;
}

# set to 1 if using proclu 1.86 and lower
my $ADD_REDUNDANT_BACK = 0;

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

my $timestart;

my $count;

# Perl function to remove whitespace from the start and end of the string
sub trim($) {
    my $string = shift;
    $string =~ s/^\s+//;
    $string =~ s/\s+$//;
    return $string;
}

# Perl function to remove all whitespace
sub trimall($) {
    my $string = shift;
    $string =~ s/\s+//g;
    return $string;
}

# Perl flipc function to reverse direction of repeats
sub flipc($) {
    my $string = shift;

    $string =~ s/\'/`/g;
    $string =~ s/"/\'/g;
    $string =~ s/`/"/g;

    return $string;
}

sub dummyquals($) {
    my $dna = shift;
    my @arr = split( //, $dna );
    my $len = scalar @arr;
    for ( my $i = 0; $i < $len; $i++ ) {
        $arr[$i] = 'a';
    }
    return join( "", @arr );
}

( $sec, $min, $hour, $mday, $mon, $year, $wday, $yday, $isdst )
    = localtime(time);
printf STDERR "\n\nstart: %4d-%02d-%02d %02d:%02d:%02d\n\n\n", $year + 1900,
    $mon + 1, $mday, $hour, $min, $sec;

my $argc = @ARGV;

if ( $argc < 10 ) {
    die
        "Usage: insert_reads.pl  clusterfile indexfolder fastafolder rotatedfolder rotatedreffile strip_454_keytags dbname msdir tempdir ispairedreads\n";
}

my $curdir          = getcwd;
my $clusterfile     = $ARGV[0];
my $indexfolder     = $ARGV[1];
my $fastafolder     = $ARGV[2];
my $rotatedfolder   = $ARGV[3];
my $rotatedreffile  = $ARGV[4];
my $strip454        = $ARGV[5];
my $DBNAME          = $ARGV[6];
my $MSDIR           = $ARGV[7];
my $TEMPDIR         = $ARGV[8];
my $IS_PAIRED_READS = $ARGV[9];

# set these mysql credentials in vs.cnf (in installation directory)
my ( $LOGIN, $PASS, $HOST ) = get_credentials($MSDIR);

my $totalReads = 0;

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
####################################

my $dbh = DBI->connect( "DBI:mysql:$DBNAME;mysql_local_infile=1;host=$HOST",
    "$LOGIN", "$PASS" )
    || die "Could not connect to database: $DBI::errstr";
my %RHASH    = ();
my %SHASH    = ();
my %HEADHASH = ();
my $negcount = 0;

$timestart = time();

if ( 0 == $ADD_REDUNDANT_BACK ) {

    system("cp $clusterfile $rotatedfolder/allwithdups.clusters");

    # THIS IS NO LONGER NEEDED AS PROCLU PUTS THEM BACK TOGETHER SINCE 1.87
}
else {

# read clusters to see what values we will store (so we don't have to store all)
    print STDERR "\nreading clusterfile to hash clustered ids..." . "";
    open FILE, "<$clusterfile" or die $!;
    while (<FILE>) {

        my @values = split( ',', $_ );

        foreach my $val (@values) {

            my $dir = '\'';
            if ( $val =~ /\"/ ) { $dir = '"'; }
            $val =~ s/[\'\"]//g;
            $val = trim($val);

            $SHASH{"$val"} = $dir;

        }
    }
    close(FILE);

    print STDERR "..." . keys(%SHASH) . " entries inserted into hash. \n\n";

    # Reading cluster file to insert redundant repeat ids back into clusters
    print STDERR
        "\nCreating new cluster file (allwithdups.clusters) to insert redundant repeat ids back into clusters...\n\n";
    opendir( DIR, $rotatedfolder );
    my @redfiles = sort grep( /\.(?:rotindex)$/, readdir(DIR) );
    closedir(DIR);
    push( @redfiles, "ref" );

    foreach my $ifile (@redfiles) {
        print STDERR "\n" . $ifile . "...";
        if ( $ifile ne "ref" ) {
            open FILE, "<$rotatedfolder/$ifile" or die $!;
        }
        else {
            open FILE, "<$rotatedreffile" or die $!;
        }

        my $inscount = 0;

        while (<FILE>) {
            my @values = split( ' ', $_ );
            $count = @values;
            if ( $count > 0 ) {
                my $ind  = $values[0];
                my $rdir = '\'';
                if ( $ind =~ /\"/ ) { $rdir = '"'; }

                $ind =~ s/[\'\"]//g;

                if ( $count >= 2 && ( exists $SHASH{"$ind"} ) ) {
                    $inscount++;
                    my $dir = $SHASH{"$ind"};
                    my $flip = ( $dir eq $rdir ) ? 0 : 1;

                    if ($flip) {
                        $SHASH{"$ind"} = trim( flipc($_) );
                    }
                    else {
                        $SHASH{"$ind"} = trim($_);
                    }

                }
            }    # end of count if cond

        }
        close(FILE);
        print STDERR "(inserts: $inscount)";
    }

    # create new file with extra repeats
    open FILEO, ">$rotatedfolder/allwithdups.clusters" or die $!;
    open FILE, "<$clusterfile" or die $!;
    while (<FILE>) {

        my @values = split( ',', $_ );
        my $i = 0;
        foreach my $val (@values) {

            $i++;

            $val = trim($val);

            my $valor = $val;

            $val =~ s/[\'\"]//g;

            if ( exists $SHASH{"$val"} && length( $SHASH{"$val"} ) > 2 ) {
                $valor = $SHASH{"$val"};
                $valor =~ s/\s/,/g;
                $valor = $valor;
            }

            if ( $i > 1 ) { print FILEO ","; }
            print FILEO $valor;

        }

        print FILEO "\n";

    }
    close(FILE);
    close(FILEO);

}    # end of ADD_REDUNDANT_BACK

# load the RHASH now with new values added
# read clusters to see what values we will store (so we don't have to store all)
print STDERR
    "\n\nreading clusterfile allwithdups.clusters to hash clustered ids (with added rotated repeats)...";
open FILE, "<$rotatedfolder/allwithdups.clusters" or die $!;
$negcount = 0;
while (<FILE>) {

    #print STDERR $_;

    my @values = split( ',', $_ );

    foreach my $val (@values) {

        my $dir = '\'';
        if ( $val =~ /\"/ ) { $dir = '"'; }
        $val =~ s/[\'\"]//g;
        $val = trim($val);

        if ( $val > 0 ) {
            $RHASH{"$val"} = $dir;

            #print "setting `$val`\n";
        }
        else {
            $negcount++;
        }
    }

    #exit(1);
}
close(FILE);

print STDERR
    keys(%RHASH)
    . " positive entries inserted into hash. (plus $negcount neg reference ones not in hash) ("
    . ( time() - $timestart )
    . ") secs\n\n";

# clear  tables
print STDERR "\ntruncating database tables\n\n";
my $sth = $dbh->prepare('TRUNCATE TABLE replnk;')
    or die "Couldn't prepare statement: " . $dbh->errstr;
$sth->execute or die "Cannot execute: " . $sth->errstr();
$sth->finish;

$sth = $dbh->prepare('ALTER TABLE replnk DISABLE KEYS;')
    or die "Couldn't prepare statement: " . $dbh->errstr;
$sth->execute or die "Cannot execute: " . $sth->errstr();
$sth->finish;

$sth = $dbh->prepare('TRUNCATE TABLE fasta_reads;')
    or die "Couldn't prepare statement: " . $dbh->errstr;
$sth->execute or die "Cannot execute: " . $sth->errstr();
$sth->finish;

$sth = $dbh->prepare('ALTER TABLE fasta_reads DISABLE KEYS;')
    or die "Couldn't prepare statement: " . $dbh->errstr;
$sth->execute or die "Cannot execute: " . $sth->errstr();
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

# now only insert the sequences that we saw in clusters
#my $sth0 = $dbh->prepare('INSERT INTO fasta_reads (head) VALUES (?)')
#                or die "Couldn't prepare statement: " . $dbh->errstr;

#$sth = $dbh->prepare('INSERT INTO replnk(rid,sid,first,last,copynum,patsize,pattern,profile,profilerc,profsize) VALUES (?,?,?,?,?,?,?,?,?,?)')
#                or die "Couldn't prepare statement: " . $dbh->errstr;
$sth
    = $dbh->prepare(
    "LOAD DATA LOCAL INFILE '$TEMPDIR/replnk_$DBNAME.txt' INTO TABLE replnk FIELDS TERMINATED BY ',' LINES TERMINATED BY '\n';"
    ) or die "Couldn't prepare statement: " . $dbh->errstr;

my $TEMPFILE;
open( $TEMPFILE, ">$TEMPDIR/replnk_$DBNAME.txt" ) or die $!;

$timestart = time();
print STDERR
    "\nreading index file and storing relevant entries in database..."
    . "\n\n";
opendir( DIR, $indexfolder );
my @indexfiles = sort grep( /\.(?:index\.renumbered)$/, readdir(DIR) );
closedir(DIR);

my $indexcount = @indexfiles;
my $i          = 0;
my $inserted   = 0;
my $processed  = 0;
my $id;
my $head;
my $first;
my $last;
my $copy;
my $pat;
my $pattern;
my $fh1;
my $fh2;

foreach my $ifile (@indexfiles) {

    open( $fh1, "<$indexfolder/$ifile" ) or die $!;
    $i = 0;

    my $lfile = $ifile;
    $lfile =~ s/index/leb36/;

    open( $fh2, "<$indexfolder/$lfile" ) or die $!;

    print STDERR "\n" . $ifile . "-" . $lfile . "...";

    my $line1 = read_file_line($fh1);
    my $line2 = read_file_line($fh2);

    while ( $line1 && $line2 ) {
        if ( $line1
            =~ /^(\d+)\t(.+)\t(\d+)\t(\d+)\t(\d+\.\d)\t(\d+)\t([A-Z]+)$/ )
        {

            if ( exists $RHASH{"$1"} ) {
                $processed++;

                $id      = $1;
                $head    = $2;
                $first   = 0;
                $last    = 0;
                $copy    = 0.0;
                $pat     = 0;
                $pattern = "";

               #   if ($3 =~  /\s(\d+)\s(\d+)\s(\d+\.\d)\s(\d+)\s([A-Z]+)/ ) {
                {
                    $first   = $3;
                    $last    = $4;
                    $copy    = $5;
                    $pat     = $6;
                    $pattern = $7;
                    $i++;

                    $head = trim($head);

                    if ( exists $HEADHASH{"$head"} ) {
                    }
                    else {
          #$sth0->execute("$head") or die "Cannot execute: " . $sth->errstr();
                        $HEADHASH{"$head"} = $processed;
                    }

       #print $id." ".$head." ".$first." ".$last." ".$copy." ".$pat." "." \n";

                    my @values = split( ' ', $line2 );
                    if ( $values[0] != $id ) {
                        die
                            "id from index file ($id) does not match id from leb36 file ($values[0])";
                    }

                    my $profile   = $values[5];
                    my $profilerc = $values[6];
                    my $proflen   = length($profile) / 2;

#$sth->execute($id,$HEADHASH{"$head"},$first,$last,$copy,$pat,$pattern,$profile,$profilerc,$proflen) or die "Cannot execute: " . $sth->errstr();

                    print $TEMPFILE $id, ",", $HEADHASH{"$head"}, ",",
                        $first, ",", $last, ",", $pat, ",", $copy, ",",
                        $pattern, ",", $profile, ",", $profilerc, ",",
                        $proflen, "\n";

                    if ( $processed % $RECORDS_PER_INFILE_INSERT == 0 ) {
                        close($TEMPFILE);
                        $count = $sth->execute();
                        $inserted += $count;
                        open( $TEMPFILE, ">$TEMPDIR/replnk_$DBNAME.txt" )
                            or die $!;
                    }

                }

            }
            else {
                #print "skip..";
            }

            #if ($i>10) { last; }
        }

        #print $_;
        $line1 = read_file_line($fh1);
        $line2 = read_file_line($fh2);

    }
    close($fh1);
    close($fh2);
    print STDERR $i;

}

close($TEMPFILE);
$count = $sth->execute();
$inserted += $count;
unlink("$TEMPDIR/replnk_$DBNAME.txt");
$sth->finish;

if ( $inserted == keys(%RHASH) ) {
    print STDERR "\n\n..."
        . $inserted
        . " read repeat entries inserted into database ("
        . ( time() - $timestart )
        . ") secs.\n\n";
}
else {
    die "\n\nERROR: hash contains " .
        keys(%RHASH)
        . " entries, while index files only have $processed matching entries and only $inserted were inserted into the database. Aborting!\n\n";
}

#goto AAA;

# get the fasta/fastq files
$totalReads = 0;
$inserted   = 0;
$processed  = 0;
$timestart  = time();
print STDERR
    "\nreading gzipped fasta files and storing relevant entries in database..."
    . "\n\n";
opendir( DIR, $fastafolder );

# the only extensions are .tgz, .tar.gz, and .gz
#my @tarballs = grep(/\.(?:tgz|gz)$/, readdir(DIR));

my @tarballs;
@tarballs = sort grep( /^fasta.*\.(?:tgz|gz)$/, readdir(DIR) );
if ( @tarballs <= 0 ) {
    rewinddir(DIR);
    @tarballs = sort grep( /^fastq.*\.(?:tgz|gz)$/, readdir(DIR) );
}

if ( @tarballs <= 0 ) {
    die "\n\nNo fasta or fastq files found. Aborting!";
}

closedir(DIR);
my $tarball_count = @tarballs;

$sth
    = $dbh->prepare(
    "LOAD DATA LOCAL INFILE '$TEMPDIR/fastareads_$DBNAME.txt' INTO TABLE fasta_reads FIELDS TERMINATED BY '\t' LINES TERMINATED BY '\n';"
    ) or die "Couldn't prepare statement: " . $dbh->errstr;

my $headstr       = "";
my $dnastr        = "";
my $qualstr       = "";
my $HEADER_SUFFIX = "";

open( $TEMPFILE, ">$TEMPDIR/fastareads_$DBNAME.txt" ) or die $!;

FILE:
foreach (@tarballs) {
    print STDERR $_ . "";

    if ( $_ =~ /\.(?:tgz|tar\.gz)$/ ) {
        open( FASTA_IN, "tar xzfmoO '$fastafolder/$_' |" );
    }
    elsif ( $_ =~ /\.gz$/ ) {
        open( FASTA_IN, "gunzip -c '$fastafolder/$_' |" );
    }
    else {
        warn "File $_ has wrong extension. Skipping this file\n";
        next FILE;
    }

    # if this is a paired file, add .1 and .2 to all headers
    if    ( $IS_PAIRED_READS && $_ =~ /s_\d+_1/ ) { $HEADER_SUFFIX = ".1"; }
    elsif ( $IS_PAIRED_READS && $_ =~ /s_\d+_2/ ) { $HEADER_SUFFIX = ".2"; }
    else                                          { $HEADER_SUFFIX = ""; }

    my $first_line = <FASTA_IN>;

    # fasta
    if ( $first_line =~ /^>(.*)(?:\n|\r)/ ) {

        #print $first_line;
        $headstr = $1;
        $dnastr  = "";
        while (<FASTA_IN>) {

            #print $_
            #end of sequence, insert and reset vars
            if ( $_ =~ /^>(.*)(?:\n|\r)/ ) {
                $headstr = trim($headstr) . $HEADER_SUFFIX;
                chomp($dnastr);
                $dnastr = trimall($dnastr);
                my $dnabak = $dnastr;
                $totalReads++;
                if ( exists $HEADHASH{"$headstr"} ) {

                    $processed++;

                    if ( $strip454 eq "1" && $dnastr !~ s/^TCAG//i ) {
                        warn "Read does not start with keyseq TCAG : "
                            . $dnabak . " ("
                            . $headstr . ")\n";
                    }

                    print $TEMPFILE $HEADHASH{"$headstr"}, "\t", "$headstr",
                        "\t", "$dnastr", "\n";

                    if ( $processed % $RECORDS_PER_INFILE_INSERT == 0 ) {
                        close($TEMPFILE);
                        $count = $sth->execute();
                        $inserted += $count;
                        open( $TEMPFILE, ">$TEMPDIR/fastareads_$DBNAME.txt" )
                            or die $!;
                    }

                }

                #print "\nhead: "."`$headstr`";
                #print "\ndna: ".$dnastr;
                #last FILE;
                $headstr = $1;
                $dnastr  = "";
            }
            else {
                $dnastr .= $_;
            }
        }
        close FASTA_IN;

        # fastq
    }
    elsif ( $first_line =~ /^\@(.*)(?:\n|\r)/ ) {

        my $qmode = 0;

        #print $first_line;
        $headstr = $1;
        $dnastr  = "";
        $qualstr = "";
        while (<FASTA_IN>) {

            #print $_
            #end of sequence, insert and reset vars
            if ( $qmode == 2 && $_ =~ /^\@(.*)(?:\n|\r)/ ) {
                $qmode   = 0;
                $headstr = trim($headstr) . $HEADER_SUFFIX;
                $dnastr  = trimall($dnastr);
                my $dnabak = $dnastr;

                #$qualstr = trimall($qualstr);
                $qualstr = ""
                    ; # need to escape this field properly, currently causing problems
                $totalReads++;
                if ( exists $HEADHASH{"$headstr"} ) {

                    my $stripped = 1;

                    $processed++;

                    if ( $strip454 eq "1" && $dnastr !~ s/^TCAG//i ) {
                        $stripped = 0;
                        warn "Read does not start with keyseq TCAG : "
                            . $dnabak . " ("
                            . $headstr . ")\n";
                    }

                    if ( $strip454 eq "1" && $stripped ) {

                        print $TEMPFILE $HEADHASH{"$headstr"}, "\t",
                            "$headstr", "\t", "$dnastr", "\t",
                            substr( $qualstr, 4 ), "\n";

                    }
                    else {

                        print $TEMPFILE $HEADHASH{"$headstr"}, "\t",
                            "$headstr", "\t", "$dnastr", "\t", $qualstr, "\n";

                    }

                    if ( $processed % $RECORDS_PER_INFILE_INSERT == 0 ) {
                        close($TEMPFILE);
                        $count = $sth->execute();
                        $inserted += $count;
                        open( $TEMPFILE, ">$TEMPDIR/fastareads_$DBNAME.txt" )
                            or die $!;
                    }

                    $qmode = 0;
                }

                $headstr = $1;
                $dnastr  = "";
                $qualstr = "";

                # start reading quals
            }
            elsif ( $qmode == 0 && $_ =~ /^\+/ ) {
                $qmode = 1;
            }
            elsif ( $qmode != 2 ) {
                if ( $qmode == 0 ) {
                    $dnastr .= $_;
                }
                else {
                    $qualstr .= $_;
                    if ( length($qualstr) >= length($dnastr) ) { $qmode = 2; }
                }
            }
        }
        close FASTA_IN;

    }
    else {

        warn "File $_ is not a FASTA file. Skipping this file\n";
        next FILE;

    }

    # one more sequence is in the buffer, finish it off
    $headstr = trim($headstr) . $HEADER_SUFFIX;
    chomp($dnastr);
    $dnastr = trimall($dnastr);
    my $dnabak = $dnastr;
    $totalReads++;
    if ( exists $HEADHASH{"$headstr"} ) {

        my $stripped = 1;

        if ( $strip454 eq "1" && $dnastr !~ s/^TCAG//i ) {
            $stripped = 0;
            warn "Read does not start with keyseq TCAG : "
                . $dnabak . " ("
                . $headstr . ")\n";
        }

        $processed++;

        if ( $first_line =~ /^>(.*)(?:\n|\r)/ ) {
            print $TEMPFILE $HEADHASH{"$headstr"}, "\t", "$headstr", "\t",
                "$dnastr", "\n";
        }
        elsif ( $first_line =~ /^\@(.*)(?:\n|\r)/ ) {

            #$qualstr = trimall($qualstr);
            $qualstr = ""
                ; # need to escape this field properly, currently causing problems
            if ( $strip454 eq "1" && $stripped ) {
                print $TEMPFILE $HEADHASH{"$headstr"}, "\t", "$headstr",
                    "\t", "$dnastr", "\t", substr( $qualstr, 4 ), "\n";
            }
            else {
                print $TEMPFILE $HEADHASH{"$headstr"}, "\t", "$headstr",
                    "\t", "$dnastr", "\t", $qualstr, "\n";
            }
        }

    }

    print STDERR " (processed: $processed)\n";

}

# cleanup
close($TEMPFILE);
$count = $sth->execute();
$inserted += $count;
unlink("$TEMPDIR/fastareads_$DBNAME.txt");
$sth->finish;

AAA:

# reenable indices
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

$sth = $dbh->prepare('ALTER TABLE fasta_reads ENABLE KEYS;')
    or die "Couldn't prepare statement: " . $dbh->errstr;
$sth->execute or die "Cannot execute: " . $sth->errstr();
$sth->finish;

$sth = $dbh->prepare('ALTER TABLE replnk ENABLE KEYS;')
    or die "Couldn't prepare statement: " . $dbh->errstr;
$sth->execute or die "Cannot execute: " . $sth->errstr();
$sth->finish;

# disconnect
$dbh->disconnect();

SetStatistics( "NUMBER_READS", $totalReads );

# check
if ( $inserted == keys(%HEADHASH) ) {
    print STDERR "\n\n..."
        . $inserted
        . " reads inserted into database  ("
        . ( time() - $timestart )
        . ") secs..\n\n";
}
elsif ( $inserted != $processed ) {
    die "\n\nERROR: inserted into database "
        . $inserted
        . " entries, while gzipped fasta files have $processed matching entries. Aborting!\n\n";
}
else {
    die "\n\nERROR: hash contains " .
        keys(%HEADHASH)
        . " entries, while gzipped fasta files only have $inserted matching entries. Aborting!\n\n";
}

print STDERR "\n\nProcessing complete (insert_reads.pl).\n";

( $sec, $min, $hour, $mday, $mon, $year, $wday, $yday, $isdst )
    = localtime(time);
printf STDERR "\n\nend: %4d-%02d-%02d %02d:%02d:%02d\n", $year + 1900,
    $mon + 1, $mday, $hour, $min, $sec;

1;

