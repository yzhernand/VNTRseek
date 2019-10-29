#!/usr/bin/env perl

# command line usage example:
#  ./run_trf.pl 6 fasta
# where 6 is the number of files to process in one batch
# and fasta is the input directory containing zipped files
#

use v5.24;
use warnings;
use FindBin;
use lib ("$FindBin::RealBin/lib", "$FindBin::RealBin/local/lib/perl5");
use Getopt::Std;
use IO::Handle;
use Carp;
use Parallel::ForkManager;
use VNTRseek::SeqReader;

my $files_processed  = 0;    # files processed
my $files_to_process = 0;    # Either: the actual number of files to process
          # OR the number of splits of a BAM file (and others?)
my %p;    # associates forked pids with output pipe pids

my %opts;
getopts( 'p:t:u:swr', \%opts );

if ( !scalar(@ARGV) || !defined( $opts{'t'} ) || !defined( $opts{'u'} ) ) {
    die "Usage:
        perl $0 [-p <processes>] -t <trf> -u <trf2proclu> -s (strip_454_TCAG) -w (warn_454_TCAG) -r (IS_PAIRED_END) <input_dir> <output_dir>\n";
}

my $trf_param        = $opts{'t'};
my $trf2proclu_param = $opts{'u'};
my $strip_454_TCAG   = ( exists $opts{'s'} && $opts{'s'} ) ? 1 : 0;
my $warn_454_TCAG    = ( exists $opts{'w'} && $opts{'w'} ) ? 1 : 0;
my $is_paired_end    = ( exists $opts{'r'} && $opts{'r'} ) ? 1 : 0;
my $max_processes    = ( exists $opts{'p'} && $opts{'p'} ) ? $opts{'p'} : 2;

# Input directory
my $input_dir = shift || die "$0: Need to provide input directory\n";

# Output directory
my $output_dir = shift || die "$0: Need to provide output directory\n";

#die "Output directory already exists -- please delete it before running this program\n" if -e $output_dir;
mkdir $output_dir;

$trf_param =~ s/['"]([^'"]+)['"] //;
my $trf_bin = $1;
$trf2proclu_param =~ s/['"]([^'"]+)['"] //;
my $trf2proclu_bin = $1;

my $seq_reader = VNTRseek::SeqReader->new(
    input_dir        => $input_dir,
    output_dir       => $output_dir,
    is_paired_end    => $is_paired_end,
    reads_split      => 1e6,
    trf_param        => [ $trf_bin, split /\s+/, $trf_param ],
    trf2proclu_param => [ $trf2proclu_bin, split /\s+/, $trf2proclu_param ],
);

warn "Will use $max_processes processes\n";

my $totreads = 0;

# Asynchronously collect input from input reader
# As soon as we have 1e6 reads, launch a new TRF process, or wait if
# not enough workers available.
my $pm = Parallel::ForkManager->new($max_processes);
$pm->run_on_finish(
    sub {
        my ( $pid, $exit_code, $ident, $exit_signal, $core_dump, $res ) = @_;
        $totreads += $res->{reads} // 0;
        ( $ENV{DEBUG} && defined $res )
            && warn "Process $ident read $res->{reads} reads.\n";
    }
);
$pm->set_waitpid_blocking_sleep(0);

# $seq_reader->start_reading( max_processes => $max_processes );
my $split_index = 0;
READS:
while ( my $reads = $seq_reader->get_reads() ) {

    # say(@{$reads}+0 . " reads, $split_index index\n");
    $pm->start($split_index) and ( $split_index++, next READS );

    # Child code
    $pm->finish( 0, { reads => 0 } ) unless $reads->%*;
    my $start_id = ( $split_index * $seq_reader->{reads_split} ) + 1;
    say "Queuing worker: $split_index, starting id: $start_id";
    my $res = $seq_reader->run_trf(
        output_prefix => "$seq_reader->{output_dir}/$split_index",
        index         => $split_index,
        input         => $reads,
        start_id      => $start_id,
    );
    $pm->finish( 0, $res );
}

say "Finished reading. Waiting for TRF processes...";
$pm->wait_all_children;

warn
    "Processing complete -- processed $split_index file(s) and $totreads reads.\n";
