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
# use File::Basename;
# use File::Temp qw/ tempfile tempdir /;
use lib "$FindBin::Bin/vntr";
require "vutil.pm";
use vutil qw(get_credentials stats_get stats_set);

# Step one, retrieve from the database the list of all reads mapping to a reftr
#TODO Test speed of DISTINCT (and syntax)
my $all_reads_mapped_query = q{SELECT DISTINCT map.refid AS reftrid, SUBSTRING_INDEX(fasta_reads.head, "_", -1) AS origintrid FROM map INNER JOIN rank ON rank.refid=map.refid AND rank.readid=map.readid INNER JOIN rankflank ON rankflank.refid=map.refid AND rankflank.readid=map.readid INNER JOIN replnk ON replnk.rid=map.readid INNER JOIN fasta_reads on fasta_reads.sid=replnk.sid INNER JOIN fasta_ref_reps ON map.refid=fasta_ref_reps.rid WHERE map.bbb = 1 ORDER BY reftrid,origintrid};
my $span1_vntrs_query = q{SELECT rid FROM fasta_ref_reps WHERE support_vntr_span1 > 0 ORDER BY rid};

# TODO Replace this in perl
# awk '{if(ARGIND==1){vntrs[$1]=1}else if((ARGIND==2)&&(-$1 in vntrs)){print $0}}' vntrs_VNTRPIPE_228486_HG38_sw100_set.span1.txt 228486_sw100_set_all_reads_mapped.txt > 228486_sw100_set_vntr_reads_mapped.txt