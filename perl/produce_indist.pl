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
use Getopt::Std;
use IO::Handle;
use FindBin;
use Cwd;
use File::Copy;
use File::Basename;
use File::Temp qw/ tempfile tempdir /;

die
    "Usage: $0 full.leb36 filtered.leb36 flanklength1[,flanklength2,...] cutoff"
    unless @ARGV == 4;

my ( $full_set, $filtered_set, $flanklengths, $cutoff ) = @ARGV;
my $install_dir = $FindBin::Bin;

my $redund_executable = "$install_dir/redund.exe";

# my $proclu_executable = "$install_dir/psearch.exe";
# Need a particular build of psearch (for now)
my $proclu_executable = "$install_dir/psearch_reftoref.exe";

#=<<Run renumber on the input filtered.leb36 file>>
# First, create a temporary directory and copy the filtered set there
my $tmpdir      = File::Temp->newdir();
my $tmpdir_name = $tmpdir->dirname;
copy( $filtered_set, "$tmpdir_name/filtered.leb36" ) or die "Copy failed: $!";

# Then renumber
system( $install_dir . "/renumber.pl -r -n " . $tmpdir_name );

# move($tmpdir->dirname . "/filtered.leb36", "./filtered.leb36");

#=<<Run redund.exe on the input leb36 files.>>
my $tmp_full = File::Temp->new( SUFFIX => ".leb36", DIR => $tmpdir_name );
my $tmp_filt = File::Temp->new( SUFFIX => ".leb36", DIR => $tmpdir_name );

# First, on the filtered file
system("$redund_executable $tmpdir_name/filtered.leb36 $tmp_filt -s -i");
if ( $? == -1 ) { die "command failed: $!\n"; }
else {
    my $rc = ( $? >> 8 );
    if ( 0 != $rc ) { die "command exited with value $rc"; }
}

system("$redund_executable $tmp_filt $tmpdir_name/filtered.leb36 -n -i");
if ( $? == -1 ) { die "command failed: $!\n"; }
else {
    my $rc = ( $? >> 8 );
    if ( 0 != $rc ) { die "command exited with value $rc"; }
}

#=<<Save list of TRs in filtered file>>
open my $filtered_set_fh, "<", "$tmpdir_name/filtered.leb36";
my %trlist;
while ( my $entry = <$filtered_set_fh> ) {

    # Read in each TR. Sign is switched by this point,
    # so trids are negative
    my ($trid) = split /\s/, $entry;
    $trlist{$trid} = 1;
}
close $filtered_set_fh;

#=<<Then run redund.exe on the full file>>
system("$redund_executable $full_set $tmp_full -s -i");
if ( $? == -1 ) { die "command failed: $!\n"; }
else {
    my $rc = ( $? >> 8 );
    if ( 0 != $rc ) { die "command exited with value $rc"; }
}

system("$redund_executable $tmp_full $tmpdir_name/full.leb36 -n -i");
if ( $? == -1 ) { die "command failed: $!\n"; }
else {
    my $rc = ( $? >> 8 );
    if ( 0 != $rc ) { die "command exited with value $rc"; }
}

#=<<Now, run psearch.exe for each flanklength given>>
my @flens = split /,/, $flanklengths;
my ( %indist, %files_seen_in );
for my $fl (@flens) {
    my $out_prefix = "reads${fl}_F_${cutoff}";
    my $proclu_cmd
        = qq($proclu_executable $tmpdir_name/filtered.leb36 $tmpdir_name/full.leb36 $install_dir/eucledian.dst $cutoff 0 0 -r $fl -m);
    warn $proclu_cmd;
    system($proclu_cmd);
    if ( $? == -1 ) { die "command failed: $!\n"; }
    else {
        my $rc = ( $? >> 8 );
        if ( 0 != $rc ) { die "command exited with value $rc"; }
    }

    # Once we verify that the map files are the same as before,
    # simply comment this line and process the .map files
    # move("$tmpdir_name/full.leb36.map", getcwd . "/$out_prefix.map");
    # TODO Only TRs indist in common to ALL flank lengths
    open my $mapping, "<", "$tmpdir_name/full.leb36.map";
    my %num_links_per_tr;
    while ( my $line = <$mapping> ) {
        $line =~ /(-\d+)(['"])=>(\d+):.*/;
        if ( ( exists( $trlist{$1} ) ) && ( -$1 != $3 ) ) {
            if ( $ENV{DEBUG} ) {
                warn
                    "$1 mapped to TR other than itself at flank length $fl: $3";
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
my $bname = basename( $filtered_set, ".leb36" );
open my $indist_out, ">", getcwd . "/$bname.indist";
for my $tr ( sort keys(%indist) ) {
    say $indist_out $tr;
}
close $indist_out;
