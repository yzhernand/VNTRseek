#!/usr/bin/env perl

# command line usage example:
#  ./produce_indist.pl full.leb36 filtered.leb36
# where 'full.leb36' is the full path to the unfiltered or larger
# reference set and 'filtered.leb36' is the full path to the filtered
# or smaller reference set.
#
# For determining indistinguishables, this script automates the process
# described in (TODO: cite), where multiple runs of profile clustering
# are performed with varying flank lengths (specifically 10, 20, and 50).
# Indistinguishables at each flank length are determined simply by
# checking if the filtered reference TR mapped to any other TR in the
# full, unfiltered set. The union of indistinguishables at all flank
# lengths is taken and written out to a file.
#
# This script allows configurable flank lengths (given as a comma
# separated list as an argument), and also collects data such as the
# number of TRs to which a TR was mapped. The latter is never used
# to make any determination, but may be necessary should the method be
# modified.
#

use strict;
use warnings;
use feature 'say';
use Getopt::Long;
use IO::Handle;
use FindBin;
use DBI;
use Cwd;
use File::Copy;
use File::Basename;
use File::Temp qw/ tempfile tempdir /;
use lib "$FindBin::RealBin/lib";
use vutil qw(make_refseq_db load_profiles_if_not_exists run_redund);
use Data::Dumper;

my %opts = (
    perform_redund_full => 1,
    perform_redund_filt => 1,
    perform_psearch     => 1,
);
GetOptions( \%opts, "redo|r", );
die
    "Usage: $0 [options] full.leb36 filtered.leb36 flanklength1[,flanklength2,...] cutoff\n"
    . "\n"
    . "-r, --redo       Force recreating indist set and profile databases\n"
    unless @ARGV == 4;

my ( $full_set, $filtered_set, $flanklengths, $cutoff ) = @ARGV;
my $install_dir       = $FindBin::Bin;
my $full_sqlitedb     = $full_set =~ s/.leb36$/.db/r;
my $filtered_sqlitedb = $filtered_set =~ s/.leb36$/.db/r;

my $redund_executable = "$install_dir/redund.exe";

my $proclu_executable = "$install_dir/psearch.exe";

#=<<Prepare SQLite dbs for reference set>>
#=<<First the filtered set>>
# Prepares the db with the ref set sequences, if not already done.
# Only need this for the filtered set.
# MUST do before the DB connection.
warn "Initializing filtered set database (if needed).\n";
make_refseq_db($filtered_set, $opts{redo});

# check if the table exists in the database (do this for both reference sets)
warn "Loading filtered set profiles (if needed).\n";
if (load_profiles_if_not_exists( $filtered_set, $opts{redo} ) )
{

    #=<<Run redund.exe on the filtered file>>
    # Since this is the filtered set, we're not interested
    # in keeping the files it creates, at this stage
    warn "Running redundancy elimination on filtered set.\n";
    run_redund( $filtered_set, "filtered.leb36", 0 );
}

#=<<Use databases to fetch rotindex files>>
# First the filtered set
my $filtered_dbh = DBI->connect(
    "DBI:SQLite:dbname=$filtered_sqlitedb",
    undef, undef,
    {   AutoCommit                 => 1,
        RaiseError                 => 1,
        sqlite_see_if_its_a_number => 1,
    }
) or die "Could not connect to database: $DBI::errstr";

$filtered_dbh->sqlite_create_function( 'mkflank', 2, sub {
    my ($lflank, $rflank) = @_;
    return substr($lflank, -60) . "|" . substr($rflank, 0, 60);
    } );

# Get rotindex saved in db
my $get_rotindex = q{SELECT rotindex
    FROM files};
# For filtered set, fetch from the database, negate rid
my $get_filtered_set_profiles = q{SELECT -rid, length(pattern) AS patsize,
    round(copynum, 2), proflen, proflenrc, profile, profilerc, nA, nC, nG, nT,
    mkflank(flankleft, flankright) AS flanks
    FROM fasta_ref_reps JOIN ref_profiles USING (rid)
        JOIN minreporder USING (rid)
    ORDER BY minreporder.idx ASC};
# Create a temporary directory and write out filtered set
# and rotindex files there
my $tmpdir      = File::Temp->newdir();
my $tmpdir_name = $tmpdir->dirname;
my $tmp_filt_file = File::Temp->new( SUFFIX => ".leb36", DIR => $tmpdir_name );

# Filtered file, ordered by min representation as given by redund
my $get_filtered_set_profiles_sth = $filtered_dbh->prepare($get_filtered_set_profiles);
$get_filtered_set_profiles_sth->execute;
open my $tmp_filt_file_fh, ">", $tmp_filt_file;
while (my @fields = $get_filtered_set_profiles_sth->fetchrow_array) {
    say $tmp_filt_file_fh join(" ", @fields)
};
close $tmp_filt_file_fh;

# Filtered rotindex
my ($rotindex_str) = $filtered_dbh->selectrow_array($get_rotindex);
die "Error getting rotindex file from filtered set. Try rerunning with --redo option\n"
    unless ($rotindex_str);
open my $tmp_rotindex, ">", "${tmp_filt_file}.rotindex";
# Need to negate all indices
$rotindex_str =~ s/(\d+)/-$1/g;
print $tmp_rotindex $rotindex_str;
close $tmp_rotindex;

#=<<Then for the full file>>
my $tmp_full_dir;
warn "Loading full set profiles (if needed).\n";
if ( load_profiles_if_not_exists( $full_set, $opts{redo} ) ) {

    #=<<Run redund.exe on the full file>>
    warn "Running redundancy elimination on full set.\n";
    # Tell run_redund that we want to keep the files it creates
    $tmp_full_dir = run_redund( $full_set, "full.leb36", 1 );
}

# FUll file location
my $tmp_full_file = $tmp_full_dir->dirname . "/full.leb36";
# Full rotindex
my $tmp_full_rotindex = "${tmp_full_file}.rotindex";

print "Press enter to continue...";
my $dummy = <STDIN>;

#=<<Now, run psearch.exe for each flanklength given>>
my @flens = split /,/, $flanklengths;
my ( %indist, %files_seen_in );
for my $fl (@flens) {
    my $out_prefix = "reads${fl}_F_${cutoff}";
    my $proclu_cmd
        = qq($proclu_executable $tmp_filt_file $tmp_full_file $install_dir/eucledian.dst $cutoff 0 0 -r $fl -m);
    warn $proclu_cmd;
    system($proclu_cmd);
    if ( $? == -1 ) { die "command failed: $!\n"; }
    else {
        my $rc = ( $? >> 8 );
        if ( 0 != $rc ) { die "command exited with value $rc"; }
    }

    # Process each map file, get all the TRs in the filtered set which mapped
    # to a TR other than itself. Filtered set TRs are negative.
    open my $mapping, "<", "$tmp_full_file.map";
    my %num_links_per_tr;
    while ( my $line = <$mapping> ) {
        $line =~ /(-\d+)(['"])=>(\d+):.*/;

        # if ( ( exists( $trlist{$1} ) ) && ( -$1 != $3 ) ) {
        if ( -$1 != $3 ) {
            if ( $ENV{DEBUG} ) {
                warn
                    "$1 mapped to TR other than itself at flank length $fl: $3\n";
            }

            # Check if seen as indist for this param set. Only
            # increment global counter if we haven't seen a mapping
            # to another TR for this parameter set yet. Keeps track
            # of how many flank lengths agreed this was indist.
            # (Counts are not currently used.)
            $indist{$1}++ unless ( exists $num_links_per_tr{$1} );

            # Mark that we've seen as indist this param set and count
            # how many map links seen so far per TR.
            $num_links_per_tr{$1}++;
        }
    }
}

# Find all elements of %indist where the count matches the length of the
# @flens array. These are the indists.
$filtered_dbh->begin_work();
$filtered_dbh->do(
    q{CREATE TEMPORARY TABLE indist (
    rid integer PRIMARY KEY)}
) or die "Couldn't do statement: $DBI::errstr\n";
my $insert_indist_sth = $filtered_dbh->prepare(
    q{INSERT INTO indist
        (rid) VALUES (?)}
);

# NB Can get rid of this after we eliminate the need for indist files.
my $bname = basename( $filtered_set, ".leb36" );
open my $indist_out, ">", getcwd . "/$bname.indist";
for my $tr ( sort keys(%indist) ) {
    say $indist_out $tr;
    $insert_indist_sth->execute(-$tr);
}
close $indist_out;
# End NB

warn "Updating indistinguishable TRs...\n";
my $update_indist = q{UPDATE fasta_ref_reps
    SET is_indist = 1
    WHERE EXISTS (
    SELECT rid FROM indist t2
    WHERE fasta_ref_reps.rid = t2.rid
)};
my $indist_count = $filtered_dbh->do($update_indist)
    or die "Couldn't do statement: " . $filtered_dbh->errstr;
warn "$indist_count indistinguishable TRs...\n";
warn "Updating singleton TRs...\n";
my $update_singleton = q{UPDATE fasta_ref_reps
    SET is_singleton = 1,
        is_dist = 1
    WHERE is_indist = 0};
my $singleton_count = $filtered_dbh->do($update_singleton)
    or die "Couldn't do statement: " . $filtered_dbh->errstr;
warn "$singleton_count singleton TRs...\n";
$filtered_dbh->commit;
$filtered_dbh->disconnect;

