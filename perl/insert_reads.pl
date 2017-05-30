#!/usr/bin/perl

use strict;
use warnings;
use 5.010;
use Cwd;
use DBI;
use POSIX qw(strftime);
use FindBin;
use File::Basename;

use lib "$FindBin::Bin/vntr";
require "vutil.pm";
use lib "$FindBin::RealBin/lib";
use ProcInputReads qw(get_reader formats_regexs compressed_formats_regexs);

use vutil qw(get_credentials write_mysql stats_set);

my $RECORDS_PER_INFILE_INSERT = 100000;

my $timestart;

my $count;

################ main ####################
say STDERR strftime( "\n\nstart: %F %T\n\n", localtime );

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

my $dbh = DBI->connect( "DBI:mysql:$DBNAME;mysql_local_infile=1;host=$HOST",
    "$LOGIN", "$PASS" )
    || die "Could not connect to database: $DBI::errstr";
my %RHASH    = ();
my %SHASH    = ();
my %HEADHASH = ();
my $negcount = 0;

$timestart = time();

system("cp $clusterfile $rotatedfolder/allwithdups.clusters");

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
                $first   = $3;
                $last    = $4;
                $copy    = $5;
                $pat     = $6;
                $pattern = $7;
                $i++;

                $head = trim($head);

                unless ( exists $HEADHASH{"$head"} ) {
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

# get the read files
$totalReads = 0;
$inserted   = 0;
$processed  = 0;
$timestart  = time();
print STDERR
    "\nreading input read files and storing relevant entries in database..."
    . "\n\n";
opendir( my $dirhandle, $fastafolder );
my @dircontents = readdir($dirhandle);
closedir($dirhandle);

if ( $ENV{DEBUG} ) {
    warn join( ",", @dircontents ) . "\n";
}

# Get all supported files. See note above on priority of input formats
my @filenames;
my ( $input_format, $compression );

# Determine sequence format
my %supported_formats_regexs = formats_regexs();
while ( my ( $sf, $pat_re ) = each %supported_formats_regexs ) {
    if (@filenames = sort
        grep( /^${pat_re}.*$/, @dircontents )
        )
    {
        $input_format = $sf;

        # Determine compression format
        my %cmp_formats_regexs = compressed_formats_regexs();
        while ( my ( $cf, $cf_re ) = each %cmp_formats_regexs ) {
            if ( $filenames[0] =~ /.*\.(?:${cf_re})/ ) {
                $compression = $cf;
                last;
            }
        }
        last;
    }
}

my $files_to_process = @filenames;
unless (@filenames > 0) {
    die "Error: no supported files found in $fastafolder. Exiting...\n";
}

$sth
    = $dbh->prepare(
    "LOAD DATA LOCAL INFILE '$TEMPDIR/fastareads_$DBNAME.txt' INTO TABLE fasta_reads FIELDS TERMINATED BY '\t' LINES TERMINATED BY '\n';"
    ) or die "Couldn't prepare statement: " . $dbh->errstr;

my $headstr       = "";
my $dnastr        = "";
my $qualstr       = "";
my $HEADER_SUFFIX = "";

open( $TEMPFILE, ">$TEMPDIR/fastareads_$DBNAME.txt" ) or die $!;

my $files_processed = 0;
while ($files_processed < $files_to_process) {
    say STDERR "Reading file/file split ", ($files_processed + 1);

    my $reader = get_reader($fastafolder, $input_format, $compression, $files_processed, \$files_to_process, \@filenames);
    while (my ( $headstr, $dnastr ) = $reader->()) {
        $headstr = trim($headstr);
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

            say $TEMPFILE $HEADHASH{"$headstr"}, "\t", "$headstr",
                "\t", "$dnastr";

            if ( $processed % $RECORDS_PER_INFILE_INSERT == 0 ) {
                close($TEMPFILE);
                $count = $sth->execute();
                $inserted += $count;
                open( $TEMPFILE, ">$TEMPDIR/fastareads_$DBNAME.txt" )
                    or die $!;
            }

        }
    }
    $files_processed++;
    say STDERR " (processed: $processed)";
}

# cleanup
close($TEMPFILE);
$count = $sth->execute();
$inserted += $count;
unlink("$TEMPDIR/fastareads_$DBNAME.txt");
$sth->finish;

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

say STDERR "\n\nProcessing complete (insert_reads.pl).";

say STDERR strftime( "\n\nend: %F %T\n\n", localtime );

1;

############# subroutines ##############
sub read_file_line {
    my $fh = shift;

    if ( $fh and my $line = <$fh> ) {
        chomp $line;
        return $line;
    }
    return;
}

# Perl function to remove whitespace from the start and end of the string
sub trim {
    my $string = shift;
    # Remove leading ">", if any
    $string =~ s/^>//;
    $string =~ s/^\s+//;
    $string =~ s/\s+$//;
    return $string;
}

# Perl function to remove all whitespace
sub trimall {
    my $string = shift;
    $string =~ s/\s+//g;
    return $string;
}

# Perl flipc function to reverse direction of repeats
sub flipc {
    my $string = shift;

    (my $ret_str = $string) =~ tr/\'"/"\'/;
}

sub dummyquals {
    my $dna = shift;
    my @arr = split( //, $dna );
    my $len = scalar @arr;
    for ( my $i = 0; $i < $len; $i++ ) {
        $arr[$i] = 'a';
    }
    return join( "", @arr );
}

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