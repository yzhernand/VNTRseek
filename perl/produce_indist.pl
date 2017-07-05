#!/usr/bin/env perl

# command line usage example:
#  ./produce_indist.pl full.leb36 filtered.leb36
# where 'full.leb36' is the full path to the unfiltered or larger
# reference set and 'filtered.leb36' is the full path to the filtered
# or smaller reference set.
#

use strict;
use warnings;
use feature 'say';
use Getopt::Std;
use IO::Handle;
use FindBin;
use File::Copy;
use File::Temp qw/ tempfile tempdir /;

die "Usage: $0 full.leb36 filtered.leb36 flanklength1[,flanklength2,...] cutoff"
	unless @ARGV==4;

my ($full_set, $filtered_set, $flanklengths, $cutoff) = @ARGV;
my $install_dir = $FindBin::Bin;

my $temp_full_fh = File::Temp->new(SUFFIX => ".leb36");
my $temp_full_name = $temp_full_fh->filename;
my $temp_filt_fh = File::Temp->new(SUFFIX => ".leb36");
my $temp_filt_name = $temp_filt_fh->filename;

#=<<Run renumber on the input filtered.leb36 file>>
# First, create a temporary directory and copy the filtered set there
my $tmpdir = File::Temp->newdir();
copy($filtered_set, $tmpdir->dirname . "filtered.leb36") or die "Copy failed: $!";
# Then renumber
system($install_dir . "/renumber.pl -r -n " . $tmpdir->dirname);
move();

#=<<Run redund.exe on the input leb36 files. First, 
system("./redund.exe $full_set $temp_full_name -s -i");
if ( $? == -1 ) { die "command failed: $!\n"; }
else {
    my $rc = ($? >> 8);
    if ( 0 != $rc ) {  die "command exited with value $rc"; }
}

system("./redund.exe $temp_full_name full.leb36 -s -i");
if ( $? == -1 ) { die "command failed: $!\n"; }
else {
    my $rc = ($? >> 8);
    if ( 0 != $rc ) {  die "command exited with value $rc"; }
}

system("./redund.exe temp.leb36  filtered.leb36 -n -i");
if ( $? == -1 ) { die "command failed: $!\n"; }
else {
    my $rc = ($? >> 8);
    if ( 0 != $rc ) {  die "command exited with value $rc"; }
}
