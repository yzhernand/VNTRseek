#!/usr/bin/perl

use strict;
use warnings;
use 5.010;
use Cwd;
use DBI;
use POSIX qw(strftime);
use FindBin;
use File::Basename;

use lib "$FindBin::RealBin/lib";
require "vutil.pm";
use ProcInputReads qw(get_reader init_bam formats_regexs compressed_formats_regexs);

use vutil qw(get_config get_dbh set_statistics get_trunc_query);

my $RECORDS_PER_INFILE_INSERT = 100000;

my $timestart;

my $count;

################ main ####################
warn strftime( "\n\nstart: %F %T\n\n", localtime );

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
my $DBSUFFIX        = $ARGV[6];
my $MSDIR           = $ARGV[7];
my $TEMPDIR         = $ARGV[8];
my $IS_PAIRED_READS = $ARGV[9];

# set these mysql credentials in vs.cnf (in installation directory)
my %run_conf = get_config( $MSDIR . "vs.cnf" );
my ( $LOGIN, $PASS, $HOST ) = @run_conf{qw(LOGIN PASS HOST)};

my $totalReads = 0;

# TODO for all files needing this function, maybe run get_config first
# to eliminate need for second arg
my $dbh = get_dbh( $DBSUFFIX, $MSDIR . "vs.cnf" )
    or die "Could not connect to database: $DBI::errstr";
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
my ( $trunc_query, $sth );
$trunc_query = get_trunc_query( $run_conf{BACKEND}, "replnk" );
$sth = $dbh->do($trunc_query)
    or die "Couldn't do statement: " . $dbh->errstr;
$trunc_query = get_trunc_query( $run_conf{BACKEND}, "fasta_reads" );
$sth = $dbh->do($trunc_query)
    or die "Couldn't do statement: " . $dbh->errstr;

if ( $run_conf{BACKEND} eq "mysql" ) {
    $sth = $dbh->do('ALTER TABLE replnk DISABLE KEYS;')
        or die "Couldn't do statement: " . $dbh->errstr;

    $sth = $dbh->do('ALTER TABLE fasta_reads DISABLE KEYS;')
        or die "Couldn't do statement: " . $dbh->errstr;

    $sth = $dbh->do('SET AUTOCOMMIT = 0;')
        or die "Couldn't do statement: " . $dbh->errstr;

    $sth = $dbh->do('SET FOREIGN_KEY_CHECKS = 0;')
        or die "Couldn't do statement: " . $dbh->errstr;

    $sth = $dbh->do('SET UNIQUE_CHECKS = 0;')
        or die "Couldn't do statement: " . $dbh->errstr;

    # now only insert the sequences that we saw in clusters
    #my $sth0 = $dbh->prepare('INSERT INTO fasta_reads (head) VALUES (?)')
    #                or die "Couldn't prepare statement: " . $dbh->errstr;

#$sth = $dbh->prepare('INSERT INTO replnk(rid,sid,first,last,copynum,patsize,pattern,profile,profilerc,profsize) VALUES (?,?,?,?,?,?,?,?,?,?)')
#                or die "Couldn't prepare statement: " . $dbh->errstr;
    $sth
        = $dbh->prepare(
        "LOAD DATA LOCAL INFILE '$TEMPDIR/replnk_$DBSUFFIX.txt' INTO TABLE replnk FIELDS TERMINATED BY ',' LINES TERMINATED BY '\n';"
        ) or die "Couldn't prepare statement: " . $dbh->errstr;
}
elsif ( $run_conf{BACKEND} eq "sqlite" ) {
    $dbh->do("PRAGMA foreign_keys = OFF");
    $dbh->{AutoCommit} = 0;

    # Insert all ref TRs from replnk file
    $sth
        = $dbh->prepare(
        qq{INSERT INTO replnk(rid,sid,first,last,patsize,copynum,pattern,profile,profilerc,profsize) VALUES (?,?,?,?,?,?,?,?,?,?)}
        );
}

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
my @replnk_rows;

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

                if ( $run_conf{BACKEND} eq "mysql" ) {
                    push @replnk_rows,
                        join( ",",
                        $id,      $HEADHASH{"$head"}, $first,
                        $last,    $pat,               $copy,
                        $pattern, $profile,           $profilerc,
                        $proflen );

                    if ( $processed % $RECORDS_PER_INFILE_INSERT == 0 )
                    {
                        open( my $TEMPFILE, ">", "$TEMPDIR/replnk_$DBSUFFIX.txt" )
                            or die $!;
                        for my $row (@replnk_rows) {
                            say $TEMPFILE $row;
                        }
                        close($TEMPFILE);
                        $count = $sth->execute();
                        $inserted += $count;
                        @replnk_rows = ();
                    }
                }
                elsif ( $run_conf{BACKEND} eq "sqlite" ) {
                    $sth->execute(
                        $id,      $HEADHASH{"$head"}, $first,
                        $last,    $pat,               $copy,
                        $pattern, $profile,           $profilerc,
                        $proflen
                    );
                    $inserted++;
                }

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

# Remaining rows
if ( ( $run_conf{BACKEND} eq "mysql" ) && @replnk_rows ) {
    open( my $TEMPFILE, ">", "$TEMPDIR/replnk_$DBSUFFIX.txt" ) or die $!;
    for my $row (@replnk_rows) {
        say $TEMPFILE $row;
    }
    close($TEMPFILE);
    $count = $sth->execute();
    $inserted += $count;
    @replnk_rows = ();
}

# $sth->finish;
$dbh->commit;
unlink("$TEMPDIR/replnk_$DBSUFFIX.txt");

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

unless ( @filenames > 0 ) {
    die "Error: no supported files found in $fastafolder. Exiting...\n";
}

# If BAM, init files list
if ($input_format eq "bam") {
    @filenames = init_bam($fastafolder, $IS_PAIRED_READS, \@filenames);
}

my $files_to_process = @filenames;

if ( $run_conf{BACKEND} eq "mysql" ) {
    $sth
        = $dbh->prepare(
        "LOAD DATA LOCAL INFILE '$TEMPDIR/fastareads_$DBSUFFIX.txt' INTO TABLE fasta_reads FIELDS TERMINATED BY '\t' LINES TERMINATED BY '\n';"
        ) or die "Couldn't prepare statement: " . $dbh->errstr;
}
elsif ( $run_conf{BACKEND} eq "sqlite" ) {
    $sth = $dbh->prepare( qq{INSERT INTO fasta_reads (sid, head, dna) VALUES(?, ?, ?)} )
        or die "Couldn't prepare statement: " . $dbh->errstr;
}

my $headstr       = "";
my $dnastr        = "";
my $qualstr       = "";
my $HEADER_SUFFIX = "";
my @fasta_reads_rows;

# TODO Can be made more efficient for BAM files?
my $files_processed = 0;
while ( $files_processed < $files_to_process ) {
    say STDERR "Reading file/file split ", ( $files_processed + 1 );

    my $reader = get_reader( $fastafolder, $input_format, $compression,
        $files_processed, \$files_to_process, \@filenames );
    while ( my ( $headstr, $dnastr ) = $reader->() ) {
        $headstr = trim($headstr);
        $dnastr  = trimall($dnastr);
        my $dnabak = $dnastr;
        $totalReads++;
        if ( exists $HEADHASH{"$headstr"} ) {

            $processed++;

            if ( $strip454 eq "1" && $dnastr !~ s/^TCAG//i ) {
                warn "Read does not start with keyseq TCAG : "
                    . $dnabak . " ("
                    . $headstr . ")\n";
            }


            if ( $run_conf{BACKEND} eq "mysql" ) {
                push @fasta_reads_rows,
                    join( "\t", $HEADHASH{"$headstr"}, "$headstr", "$dnastr" );

                if ( $processed % $RECORDS_PER_INFILE_INSERT == 0 )
                {
                    open( my $TEMPFILE, ">", "$TEMPDIR/fastareads_$DBSUFFIX.txt" )
                        or die $!;
                    for my $row (@fasta_reads_rows) {
                        say $TEMPFILE $row;
                    }
                    close($TEMPFILE);
                    $count = $sth->execute();
                    $inserted += $count;
                    @fasta_reads_rows = ();
                }
            }
            elsif ( $run_conf{BACKEND} eq "sqlite" ) {
                $sth->execute( $HEADHASH{"$headstr"}, "$headstr", "$dnastr" );
                $inserted++;
            }

        }
    }
    $files_processed++;
    say STDERR " (processed: $processed)";
}

# cleanup
if ( ( $run_conf{BACKEND} eq "mysql" ) && @fasta_reads_rows ) {
    open( my $TEMPFILE, ">", "$TEMPDIR/fastareads_$DBSUFFIX.txt" ) or die $!;
    for my $row (@fasta_reads_rows) {
        say $TEMPFILE $row;
    }
    close($TEMPFILE);
    $count = $sth->execute();
    $inserted += $count;
    @fasta_reads_rows = ();
}
# $sth->finish;
$dbh->commit;
unlink("$TEMPDIR/fastareads_$DBSUFFIX.txt");

# reenable indices
if ($run_conf{BACKEND} eq "mysql") {
    $sth = $dbh->do('SET AUTOCOMMIT = 1;')
        or die "Couldn't do statement: " . $dbh->errstr;

    $sth = $dbh->do('SET FOREIGN_KEY_CHECKS = 1;')
        or die "Couldn't do statement: " . $dbh->errstr;

    $sth = $dbh->do('SET UNIQUE_CHECKS = 1;')
        or die "Couldn't do statement: " . $dbh->errstr;

    $sth = $dbh->do('ALTER TABLE fasta_reads ENABLE KEYS;')
        or die "Couldn't do statement: " . $dbh->errstr;

    $sth = $dbh->do('ALTER TABLE replnk ENABLE KEYS;')
        or die "Couldn't do statement: " . $dbh->errstr;
}
elsif ($run_conf{BACKEND} eq "sqlite") {
    $dbh->{AutoCommit} = 1;
    $dbh->do("PRAGMA foreign_keys = ON");
}

# disconnect
$dbh->disconnect();

set_statistics( $DBSUFFIX, "NUMBER_READS", $totalReads );

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

warn strftime( "\n\nend: %F %T\n\n", localtime );

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

    ( my $ret_str = $string ) =~ tr/\'"/"\'/;
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
