#!/usr/bin/perl

# command line usage example:
#  ./run_trf.pl 6 fasta
# where 6 is the number of files to process in one batch
# and fasta is the input directory containing zipped files
#

use strict;
use warnings;
use feature 'state';
use Getopt::Std;
use IO::Handle;
use Carp;

# TODO fewer globals. All functions get passed the variables they need (or refs to vars)
# TODO Write a generic function for generating sequences to be passed to a TRF instance/pipe
#     This function will use different functions as backends depending on the input format.
#     Should this function use IO::* in some way?
#     API should look like:
#     my @seq_splits = this_function($format, $compression, $max_processes, \@file_list);
#     for (my $i=0; $i < $max_processes; ++$i) {
#         $p{ fork_trf($seq_splits->())} } = 1;
#     }
#
#     @seq_splits is an array of function references. These are closures with predetermined
#     splits of the input data, either a subset of the files, as in FASTA/FASTQ files, or
#     a subset of ranges as in BAM files.

my $files_to_process = 1;    # number of files to process in one batch
my $files_processed  = 0;    # files processed
my %p;                       # associates forked pids with output pipe pids
my $max_processes = 0;

# List all supported file extensions and input formats here. Order of
# input formats is in priority order: first format if found is used.
# For new formats, simply add the name of the format here and require
# that input file names have it as a prefix (eg, fasta files should
# have prefix "fasta_")
my @supported_format_names = qw(fasta fastq bam);
my %supported_formats;
@supported_formats{@supported_format_names} = @supported_format_names;

# Dispatch table to call appropriate function for each format.
# Uses the list @supported_format_names above for hash keys. When
# adding a new format, add the name first to @supported_format_names,
# and then add the reference to the reader function in the correct
# place in the values list.
my %reader_table;
@reader_table{@supported_format_names}
    = ( \&read_fasta, \&read_fastq, \&read_bam );

# Similarly for compression formats, add a new one to the list, then
# add the regex and command needed to extract in the correct position
# in the values list
my @compressed_format_names = qw(targz gzip tarbzip bzip tarxz xz);
my %compressed_formats;
@compressed_formats{@compressed_format_names} = (
    qr/tgz|tar\.gz/,               qr/gz/,
    qr/tbz|tbz2|tar\.bz|tar\.bz2/, qr/bz|bz2/,
    qr/txz|tar\.xz/,               qr/xz/
);

# All decompression commands end with a space, for convenience
my %decompress_cmds;
@decompress_cmds{@compressed_format_names} = (
    "tar xzfmO ",
    "gunzip -c ",
    "tar xjfmO ",
    "bzip2 -c ",
    "tar xJfmO ",
    "xzcat "
);

my %opts;
getopts( 'p:t:u:swr', \%opts );

if ( !scalar(@ARGV) || !defined( $opts{'t'} ) || !defined( $opts{'u'} ) ) {
    die "Usage:
        perl $0 [-p <processes>] -t <trf> -u <trf2proclu> -s (strip_454_TCAG) -w (warn_454_TCAG) -r (IS_PAIRED_READS) <input_dir> <output_dir>\n";
}

my $TRF_PARAM        = $opts{'t'};
my $TRF2PROCLU_PARAM = $opts{'u'};
my $strip_454_TCAG   = ( defined $opts{'s'} && $opts{'s'} ) ? 1 : 0;
my $warn_454_TCAG    = ( defined $opts{'w'} && $opts{'w'} ) ? 1 : 0;
my $IS_PAIRED_READS  = ( defined $opts{'r'} && $opts{'r'} ) ? 1 : 0;

my $reverse_read = 1; # if 1, each read will be reversed and processed as well

my $HEADER_SUFFIX = ""
    ;   # set by program automatically for paired reads to distinguish headers

$max_processes = $opts{'p'} if defined $opts{'p'};

# Input directory
my $input_dir = $ARGV[0];
die "Need to provide input directory\n" if !defined $input_dir;
die "Input directory does not exist\n"  if !-d $input_dir;

# Output directory
my $output_dir = $ARGV[1];
die "Need to provide output directory\n" if !defined $output_dir;

#die "Output directory already exists -- please delete it before running this program\n" if -e $output_dir;
mkdir $output_dir;

# get a list of input files
opendir( my $dirhandle, $input_dir );
my @dircontents = readdir($dirhandle);
closedir($dirhandle);

if ( $ENV{DEBUG} ) {
    warn join( ",", @dircontents ) . "\n";
}

# Get all supported files. See note above on priority of input formats
my @filenames;
my ( $input_format, $compression );

# Determine sequence format
for my $sf (@supported_format_names) {
    my $pat_re = $supported_formats{$sf};
    if (@filenames = sort
        grep( /^${pat_re}.*$/, @dircontents )
        )
    {
        $input_format = $sf;

        # Determine compression format
        for my $cf (@compressed_format_names) {
            my $cf_re = $compressed_formats{$cf};
            if ( $filenames[0] =~ /.*\.(?:${cf_re})/ ) {
                $compression = $cf;
                last;
            }
        }
        last;
    }
}

my $file_count = @filenames;
my $compression_msg
    = ($compression)
    ? "compressed as $compression"
    : "unknown compression (assuming uncompressed)";
die "0 supported files found in $input_dir. Exiting\n" if $file_count == 0;
warn
    "$file_count supported files ($input_format format, $compression) found in $input_dir\n";

# /proc/cpuinfo is found on Linux/cygwin only: use other methods for other platforms
if ( $max_processes == 0 ) {
    open my $cpuinfo, "<", "/proc/cpuinfo";
    if ($cpuinfo) {
        $max_processes = scalar( map /^processor/, <$cpuinfo> );
        warn "$max_processes CPU core(s) detected\n";
        if ( $max_processes == 0 ) {
            warn "Unknown formatting in /proc/cpuinfo: Assuming single CPU\n";
            $max_processes = 1;
        }
    }
    else {
        warn "Could not open /proc/cpuinfo: Assuming single CPU\n";
        $max_processes = 1;
    }
}

warn "Will use $max_processes processes\n";

# fork as many new processes as there are CPUs
for ( my $i = 0; $i < $max_processes; $i++ ) {
    $p{ fork_proc( $input_format, $compression, $files_processed,
            \@filenames )
    } = 1;
    $files_processed++;
}

# wait for processes to finish and then fork new ones
while ( ( my $pid = wait ) != -1 ) {

    # check return value
    my ( $rc, $sig, $core ) = ( $? >> 8, $? & 127, $? & 128 );
    if ($core) {
        warn "run_trf process $pid dumped core\n";
        exit(1000);
    }
    elsif ( $sig == 9 ) {
        warn "run_trf process $pid was murdered!\n";
        exit(1001);
    }
    elsif ( $rc != 0 ) {
        warn "run_trf process $pid has returned $rc!\n";
        exit($rc);
    }

    if ( $p{$pid} ) {

        # one instance has finished processing -- start a new one
        delete $p{$pid};

     # For BAM files (and maybe other formats?) the parent does not know
     # how many total processes are needed to run the whole file. For these
     # formats, reader functions need to write out a temporary file signalling
     # the last processes that needs to run has already begun.
        if ( -e "$output_dir/trf_alldone" ) {
            unlink("$output_dir/trf_alldone");
            last;
        }
        else {
            $p{ fork_proc(
                    $input_format,    $compression,
                    $files_processed, \@filenames
                )
            } = 1;
            $files_processed++;
        }
    }
    else {
        die
            "ERROR: Process wait()ed on was not found in our process list. PID=$pid\n";
    }
}

warn "Processing complete -- processed $files_processed file(s).\n";

############################ Procedures ###############################################################

=item I<read_fasta()>

Given a list of files, return a sub which knows how to open each file and
where it is in the list of files.

=cut

=item I<make_file_streams()>

Takes file format and compression format names, and a list of full
paths to files in that format. Then returns a function which, when
called, returns a stream to sequences which can be used to run
TRF/TRF2PROCLU.

=cut

sub fork_proc {
    my ( $format, $compression, $files_processed, $filelist ) = @_;
    defined( my $pid = fork() )
        or die "Unable to fork: $!\n";
    if ( $pid == 0 ) {    #Child

        my $reader
            = $reader_table{$format}
            ->( $compression, $files_processed, $filelist );
        my $output_prefix = "$output_dir/$files_processed";
        exit unless $reader;
        warn "Running child, files_processed = $files_processed...\n";

        # TODO Error checking if TRF, in the start of the pipe, breaks down
        local $SIG{PIPE} = sub { die "Error in trf+trf2proclu pipe: $?\n" };
        open my $trf_pipe,
            "|$TRF_PARAM | $TRF2PROCLU_PARAM -o '$output_prefix.index' > '$output_prefix.leb36'"
            or die "Cannot start TRF+trf2proclu pipe: $!\n";

        # TODO Need way of logging TRF output?
        # open my $logfile, ">", "$output_prefix.log"
        #     or die "Error opening file $output_prefix.log: $!\n";
        # $logfile->autoflush;
        my $debug_trs_found = 0;

        while ( my ( $header, $body ) = $reader->() ) {

            # say $logfile $data[0] . "\n" . $data[1];
            pipe_to_trf( $trf_pipe, $header, $body );
        }

  # Normally, close() returns false for failure of a pipe. If the only problem
  # was that the exit status of the pipe was non-0, then $! == 0.
  # Important because trf2proclu returns non-0 on success.
        if ( !close $trf_pipe ) {

    # Here the process trf2proclu finished with non-0 AND there was some other
    # problem, since $! is not 0.
            if ($!) {
                warn "Error closing trf+trf2proclu pipe: $!\n";
                exit(1002);
            }

       # Here trf2proclu exited with non-0 status but that was the only issue,
       # so just report that value.
            elsif ( $? < -2 ) {
                warn "trf+trf2proclu pipe has returned $?\n";
                exit($?);
            }

            # TODO Doesn't account for TRF's negative return values
        }

        # Check exit error

        exit;
    }
    else {    # Parent
              # warn "Running parent, out_counter = $out_counter...\n";
        return $pid;
    }
}

sub reverse_complement {
    my $dna = shift;

    # reverse the DNA sequence
    my $revcomp = reverse($dna);

    # complement the reversed DNA sequence
    $revcomp =~ tr/ACGTacgt/TGCAtgca/;
    return $revcomp;
}

# Processes input FASTA records to reverse them or remove 454 tags
# Prints processed output directly to given file handle to TRF
# process.
sub pipe_to_trf {
    my ( $trf_fh, $header, $body ) = @_;

    # warn "Processing header $header";

    if ( $reverse_read && $header ne "" ) {
        say $trf_fh $header . "_" . length($body) . "_RCYES";
        say $trf_fh reverse_complement($body);
    }

    # FASTA header
    say $trf_fh $header;
    if ( $strip_454_TCAG && ( $body !~ s/^TCAG//i ) ) {
        if ($warn_454_TCAG) {
            warn
                "Read does not start with keyseq TCAG. Full sequence: $body\n";
        }
        else {
            die
                "Read does not start with keyseq TCAG. Full sequence: $body\n";
        }
    }
    say $trf_fh $body;
}

sub read_fasta {
    my ( $compression, $files_processed, $filelist ) = @_;

# warn "Need to process " . scalar(@$filelist) . " files, working on $files_processed\n";

    # Don't do more if we've exhausted the file list
    if ( $files_processed >= @$filelist ) {
        system("touch '$output_dir/trf_alldone'");
        return undef;
    }

 # If file uncompressed, simply read from file. Else open a pipe to a command.
    my $openmode = ($compression) ? "-|" : "<";

    warn "Processing file " . $filelist->[$files_processed] . "\n";
    my $filename
        = ( ($compression) ? $decompress_cmds{$compression} : "" ) . '"'
        . "$input_dir/"
        . $filelist->[$files_processed] . '"';

    # $files_processed contains how many files processed so far.
    # Use to index into filelist
    # warn $filename;
    local $/ = ">";

    # warn "Filename/command = '$filename'\n";
    open my $fasta_fh, $openmode, $filename
        or die "Error opening file " . $filename;

   # Consume first empty record because of the way $/ splits the FASTA format.
    <$fasta_fh>;
    return sub {
        local $/ = ">";
        my $fasta_rec = <$fasta_fh>;
        return () unless ($fasta_rec);
        my ( $header, @seqlines ) = split( /\n+/, $fasta_rec );
        chomp $header;
        chomp @seqlines;

        # warn "header: '$header'";
        my $seq = join( "", @seqlines );

        # warn "seq: '$seq'";
        return ( ">" . $header, $seq );
    };
}

sub read_fastq {

    # Code modified from https://www.biostars.org/p/11599/#11657
    my ( $compression, $files_processed, $filelist ) = @_;

    # Don't do more if we've exhausted the file list
    if ( $files_processed >= @$filelist ) {
        system("touch '$output_dir/trf_alldone'");
        return undef;
    }

 # If file uncompressed, simply read from file. Else open a pipe to a command.
    my $openmode = ($compression) ? "-|" : "<";

    warn "Processing file " . $filelist->[$files_processed] . "\n";
    my $filename
        = ( ($compression) ? $decompress_cmds{$compression} : "" ) . '"'
        . "$input_dir/"
        . $filelist->[$files_processed] . '"';

    # $files_processed contains how many files processed so far.
    # Use to index into filelist
    # warn "Filename/command = '$filename'\n";
    open my $fastq_fh, $openmode, $filename
        or die "Error opening file " . $filename;

    my $aux         = undef;
    my $seq_counter = 0;

    return sub {

        # warn "Seq number: " . $seq_counter++ . "\n";
        $aux = [ undef, 0 ] unless ( defined $aux );
        return () if ( $aux->[1] );
        if ( !defined( $aux->[0] ) ) {
            while (<$fastq_fh>) {
                chomp;
                if ( substr( $_, 0, 1 ) eq '>' || substr( $_, 0, 1 ) eq '@' )
                {
                    $aux->[0] = $_;
                    last;
                }
            }
            if ( !defined( $aux->[0] ) ) {
                $aux->[1] = 1;
                return ();
            }
        }

        # warn "After header: $.\n";
        my $name = $aux->[0] =~ /^.(\S+)/ ? $1 : '';
        my $seq = '';
        my $c;
        $aux->[0] = undef;
        while (<$fastq_fh>) {
            chomp;
            $c = substr( $_, 0, 1 );
            last if ( $c eq '>' || $c eq '@' || $c eq '+' );
            $seq .= $_;
        }

        # warn "After seq: $.\n";
        $aux->[0] = $_;
        $aux->[1] = 1 if ( !defined( $aux->[0] ) );
        croak "Next line is not qual header: c = $c" if ( $c ne '+' );

        # return (">" . $name, $seq) if ($c ne '+');
        my $qual = '';
        while (<$fastq_fh>) {
            chomp;
            $qual .= $_;
            if ( length($qual) >= length($seq) ) {
                $aux->[0] = undef;
                last;
            }
        }

        # warn "After qual: $.\n";
        $aux->[1] = 0;
        return ( ">" . $name, $seq );
    };
}

# Requires samtools
sub read_bam {
    my ( $compression, $files_processed, $filelist ) = @_;
    my $bamfile = "$input_dir/" . $filelist->[0];

    # warn "$bamfile\n";
    my $baifile = $bamfile . ".bai";

    # warn "$baifile\n";
    my $unmapped_template = "*";
    my @samcmds;

    # Check if .bai file exists and then run MakeBedFiles.jar
    die
        "Error reading bam file: corresponding bai file required for processing bam files."
        unless ( -e -r $baifile );

    # Requires samtools to be installed/available
    # Get all regions in the bam file
    my @regions = qx(
            samtools idxstats "$bamfile"
        );

    # Don't do more if we've exhausted the file list
    # warn "Regions: " . scalar(@regions) . "\n";

# We need to read the regions in the bam file and construct samtools commands
# These are all saved in an array which are interated through like the FASTA files before.
# Then we catch the output of these through a file handle, and process into FASTA.
    for my $r (@regions) {
        my ( $chr, $end, $num_aln, $num_unaln ) = split /\s+/, $r;
        my $unmapped = ( $chr eq $unmapped_template );

        # Don't save sequence with 0 reads
        next if ( ( $num_aln + $num_unaln ) == 0 );

        # $start is always 1
        my $region = "$chr:1-$end";
        my $scmd
            = "samtools view "
            . ( ($unmapped) ? "-f 4 " : "" )
            . $bamfile
            . ( ($unmapped) ? "" : " $region" );
        push @samcmds, $scmd;
    }

    if ( $files_processed >= @samcmds ) {
        system("touch '$output_dir/trf_alldone'");
        return undef;
    }

    # warn "$files_processed\n";
    warn "Processing bam chunk using: " . $samcmds[$files_processed] . "\n";

    local $SIG{PIPE} = sub { die "Error in samtools pipe: $?\n" };
    open my $samout, "-|", $samcmds[ $files_processed++ ]
        or die "Error opening samtools pipe: $!\n";
    return sub {
        my $bam_rec = <$samout>;
        return () unless ($bam_rec);

        # print ">" $1 "\n" $10
        my ($header, undef, undef, undef, undef,
            undef,   undef, undef, undef, $seq
        ) = split( /\s+/, $bam_rec );
        return ( ">" . $header, $seq );
        }
}
