#!/usr/bin/env perl

# command line usage example:
#  ./remove_sw_fpvntrs.pl <TODO>
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
use DBI;
# use File::Copy;
use File::Basename;
# use File::Temp qw/ tempfile tempdir /;
use lib "$FindBin::Bin/vntr";
require "vutil.pm";
use vutil qw(get_credentials get_config);

# Get input files
die "Usage: $0 <dbsuffix> [new reference set name]"
	unless (@ARGV == 1);
my ($dbsuffix, $refsetname) = @ARGV;
my $MSDIR = $ENV{HOME} . "/${dbsuffix}.";

# Connect to DB
my ( $login, $pass, $host ) = get_credentials($MSDIR);
my $dbh = DBI->connect( "DBI:mysql:VNTRPIPE_$dbsuffix;mysql_local_infile=1;host=$host",
    "$login", "$pass" )
    || die "Could not connect to database: $DBI::errstr";
# Retrieve from the database the list of all reads mapping to a reftr
# my $all_reads_mapped_query = q{SELECT DISTINCT map.refid AS reftrid, SUBSTRING_INDEX(fasta_reads.head, "_", -1) AS origintrid
# 	FROM map INNER JOIN replnk ON replnk.rid=map.readid
# 		INNER JOIN fasta_reads on fasta_reads.sid=replnk.sid
# 	WHERE map.bbb = 1 ORDER BY reftrid,origintrid};

# Retrieve from the database the list of all reads mapping to a reftr called as a VNTR with at least one spanning read
my $vntr_reads_mapped_query = q{SELECT DISTINCT map.refid AS reftrid, SUBSTRING_INDEX(fasta_reads.head, "_", -1) AS origintrid
	FROM map INNER JOIN replnk ON replnk.rid=map.readid
		INNER JOIN fasta_reads on fasta_reads.sid=replnk.sid
		INNER JOIN fasta_ref_reps ON map.refid = fasta_ref_reps.rid
	WHERE map.bbb = 1 AND support_vntr_span1 > 0 ORDER BY reftrid,origintrid;};
my $get_vntr_mapped_reads_sth = $dbh->prepare($vntr_reads_mapped_query);
$get_vntr_mapped_reads_sth->execute
	or die "Error executing query for all mapped reads: " . $get_vntr_mapped_reads_sth->errstr();
my ($vntr, $mapped_to, %discard);
$get_vntr_mapped_reads_sth->bind_columns(\($vntr, $mapped_to));

# Iterate through this list and remove all VNTRs and other mapped loci
while ($get_vntr_mapped_reads_sth->fetch) {
	$discard{$vntr} = 1;
	$discard{$mapped_to} = 1;
}

for my $d (sort keys %discard) {
	say $d;
}

my @ref_files = (get_config($MSDIR))[5..7];
for my $r (@ref_files) {
	my @lines;
	my ($name, $path, $suffix) = fileparse( $r, ".leb36", ".seq", ".indist" );
	# If not an absolute path, assume file is in install dir.
	# TODO maybe file IS in pwd!
	$path = "$FindBin::Bin" if ($path eq "./");
	open my $ref_in_fh, "<", "$path/$name.$suffix";
	while (my $line = <$ref_in_fh>) {
		my ($trid) = split /[,\s]/, $line;
		# If this is the .seq header, we keep that.
		if ($trid eq "Repeatid") {
			push @lines, $line;
			next;
		}

		push @lines, $line unless (exists $discard{abs($trid)});
	}
	close $ref_in_fh;

	open my $ref_out_fh, ">", ( ($refsetname) ? $refsetname : @lines ) . ".$suffix";
	print $ref_out_fh, @lines;
	close $ref_out_fh;
}

# my $span1_vntrs_query = q{SELECT rid FROM fasta_ref_reps WHERE support_vntr_span1 > 0 ORDER BY rid};

# TODO Replace this in perl
# awk '{if(ARGIND==1){vntrs[$1]=1}else if((ARGIND==2)&&(-$1 in vntrs)){print $0}}' vntrs_VNTRPIPE_228486_HG38_sw100_set.span1.txt 228486_sw100_set_all_reads_mapped.txt > 228486_sw100_set_vntr_reads_mapped.txt