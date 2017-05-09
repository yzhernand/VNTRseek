#!/usr/bin/perl

# command line usage example:
#  ./run_trf.pl 6 fasta
# where 6 is the number of files to process in one batch
# and fasta is the input directory containing zipped files
#

use strict;
use warnings;

use Getopt::Std;
use IO::Handle;

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

# Dispatch table for use by sequence reader
my %reader_table = (
    fasta => \&read_fasta,
    fastq => \&read_fastq,
    bam   => \&read_bam
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

# List all supported file extensions and input formats here. Order of
# input formats is in priority order: first format if found is used.
my %supported_formats = (
    fasta => "fasta",
    fastq => "fastq",
    bam   => "bam"
);
my @supported_format_names = qw(fasta fastq bam);
my %compressed_formats     = (
    targz   => qr/tgz|tar\.gz/,
    gzip    => qr/gz/,
    tarbzip => qr/tbz|tbz2|tar\.bz|tar\.bz2/,
    bzip    => qr/bz|bz2/,
    tarxz   => qr/txz|tar\.xz/,
    xz      => qr/xz/
);
my @compressed_format_types = qw(targz gzip tarbzip bzip tarxz xz);
# All decompression commands end with a space, for convenience
my %decompress_cmds         = (
    targz   => "tar xzfmO ",
    gzip    => "gunzip -c ",
    tarbzip => "tar xjfmO ",
    bzip    => "bzip2 -c ",
    tarxz   => "tar xJfmO ",
    xz      => "xzcat "
);

# Get all supported files. See note above on priority of input formats
my @tarballs;
my ( $input_format, $compression );
for my $sf (@supported_format_names) {
    my $pat_re = $supported_formats{$sf};
    if (@tarballs = sort
        grep( /^${pat_re}.*$/, @dircontents )
        )
    {
        $input_format = $sf;
        for my $cf (@compressed_format_types) {
            my $cf_re = $compressed_formats{$cf};
            if ( $tarballs[0] =~ /.*\.(?:${cf_re})/ ) {
                $compression = $cf;
                last;
            }
        }
        last;
    }
}

my $tarball_count = @tarballs;
my $compression_msg
    = ($compression)
    ? "compressed as $compression"
    : "unknown compression (assuming uncompressed)";
die "0 supported files found in $input_dir. Exiting\n" if $tarball_count == 0;
warn
    "$tarball_count supported files ($input_format format, $compression) found in $input_dir\n";

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

# TODO Now figure out how to best call Marzie's script for BAM files before anything else if input
# is BAM format.
if ( $input_format == "bam" ) {

# Will use samtools for reading instead. Pass on execution to the BAM file reading script
    my $bamfile = shift @tarballs;
    system(
        "./fromBam.sh $max_processes \"$bamfile\" \"$output_dir\" \"$TRF_PARAM\" \"$TRF2PROCLU_PARAM\""
    );
    return 0;
}

# because of the limit of number of open files at the same time, let's keep number of files to under 200
# redundancy step (3) opens them all at the same time, and if multiple pipelines are run it might be possible
# to reach the limit on some system (usually 1024)
while ( int( $tarball_count / $files_to_process ) > 200 ) {
    $files_to_process++;
}
warn "Will use $max_processes processes\n";
warn "Will process $files_to_process file(s) per batch\n";

# fork as many new processes as there are CPUs
for ( my $i = 0; $i < $max_processes; $i++ ) { $p{ fork_trf() } = 1 }

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
        $p{ fork_trf() } = 1;
    }
    else {
        die "ERROR: Do not remember process PID=$pid\n";
    }
}
warn "Processing complete -- processed $files_processed file(s).\n";
0;

############################ Procedures ###############################################################

=item I<read_fasta()>

Given a list of files, return a sub which knows how to open each file and
where it is in the list of files.

=cut

sub read_fasta {
    my ( $max_processes, $compression, $filelist ) = @_;

 # If file uncompressed, simply read from file. Else open a pipe to a command.
    my $openmode = ($compression) ? "-|" : "<";
    my $decom_cmd
        = ($compression) ? $decompress_cmds{$compression} : "";
    my $files_processed = 0;
    return sub {

        # $files_processed contains how many files processed so far.
        # Use to index into filelist
        open my $fasta_fh, $openmode,
              $decom_cmd
            . "'$input_dir/"
            . $filelist->[ $files_processed++ ] . "'";
        return $fasta_fh;
    };
}

=item I<make_file_streams()>

Takes file format and compression format names, and a list of full
paths to files in that format. Then returns a function which, when
called, returns a stream to sequences which can be used to run
TRF/TRF2PROCLU.

=cut

sub make_file_streams {

    # First, take in the format, and the list of items
    my ( $format, $compression, @files ) = @_;

    # S

    # Now, return the appropriate function
}

sub fork_trf {
    if ( $files_processed >= $tarball_count ) {
        return 0;
    }

    # unzip a predefined number of files
    my $until = $files_processed + $files_to_process - 1;
    $until = $tarball_count - 1 if $until > ( $tarball_count - 1 );
    warn 'Processing files '
        . ( $files_processed + 1 ) . ' to '
        . ( $until + 1 ) . "\n";
    my $output_prefix    = "$output_dir/$files_processed-$until";
    my @file_slice       = @tarballs[ ($files_processed) .. ($until) ];
    my $file_slice_count = @file_slice;
    $files_processed += $files_to_process;

    # Try to fork
    defined( my $pid = fork )
        or die "Unable to fork: $!\n";

    # Child process. This process use open() to fork again. Its child,
    # "grandchild" will run TRF. The parent process, "child", will have
    # a filehandle to the *output* of TRF and will use that as input
    if ( $pid == 0 ) {

        #warn "This is child\n";
        # Open a filehandle to the input of a child process. Forks a child.
        defined( my $grandchild_pid = open my $grandchild, '-|' )
            or die "Unable to open grandchild: $!\n";

        # This branch is that child process, "grandchild". Runs TRF.
        if ( $grandchild_pid == 0 ) {

            #warn "Starting TRF: $TRF_PARAM";
            defined( my $trf_pid = open( my $trf_fh, "| $TRF_PARAM" ) )
                or die "Cannot start TRF: $!\n";
            foreach (@file_slice) {
                next if not defined $_;

                # Read from sequence files. Detect if files are compressed,
                # and decompress using appropriate decom binary + options
                # If file has none of the known compression extensions, assume
                # decompressed and read from file.
                # TODO Test
                warn "Reading from sequence files $input_dir/$_\n";
                my $input_seq_fh;
                if ($compression) {
                    open( $input_seq_fh, "-|",
                        $decompress_cmds{$compression} . " '$input_dir/$_'" );
                }
                else {
                    warn "Assuming file $_ uncompressed...\n";
                    open( $input_seq_fh, "<", "$input_dir/$_" );
                }

                # if this is a paired file, add .1 and .2 to all headers
                if ( $IS_PAIRED_READS && $_ =~ /s_\d+_1/ ) {
                    $HEADER_SUFFIX = ".1";
                }
                elsif ( $IS_PAIRED_READS && $_ =~ /s_\d+_2/ ) {
                    $HEADER_SUFFIX = ".2";
                }
                else { $HEADER_SUFFIX = ""; }

                # Parse the input file. Currently supports FASTA and FASTQ
                parse_readfile( $trf_fh, $input_seq_fh );
                close $input_seq_fh;
            }

            #close TRF;

            # check return value
            #my ($rc, $sig, $core) = ($? >> 8, $? & 127, $? & 128);
            #if ($core){
            #    warn "trf process $trf_pid dumped core\n";
            #    exit (1000);
            #}elsif($sig == 9){
            #    print  STDERR "trf process $trf_pid was murdered!\n";
            #    exit (1001);
            #}elsif ($rc != 0){
            #    warn  "trf process $trf_pid has returned $rc!\n";
            #    exit ($rc);
            #}

            if ( !close $trf_fh ) {
                if ($!) {
                    warn "Error closing trf process $trf_pid pipe: $!\n";
                    exit(1002);
                }
                elsif ( $? != 0 ) {
                    warn "trf process $trf_pid has returned $?!\n";
                    exit($?);
                }
            }

            #warn "Exiting grandchild\n";
            exit 0;
        }

# This is the parent process, "child". Has the filehandle to the output of TRF.
# Uses TRF output to pipe into trf2proclu.
        else {
            #warn "This is parent of grandchild $grandchild_pid\n";
            # Open a pipe to trf2proclu process
            defined(
                my $trf2proclu_pid = open( TRF2PROCLU,
                    "| $TRF2PROCLU_PARAM -o '$output_prefix.index' > '$output_prefix.leb36'"
                )
            ) or die "Cannot start trf2proclu: $!\n";

            # DEBUG
            open( my $LOGFILE, ">${output_prefix}.log" )
                or die
                "Can't open for writing logfile ('${output_prefix}.log')!";
            $LOGFILE->autoflush;

            my $debug_trs_found = 0;

   # while input comes in from TRF, pipe into the input of trf2proclu process.
            while (<$grandchild>) {
                print TRF2PROCLU $_;

                #DEBUG
                #print $LOGFILE $_;
                $debug_trs_found++;
                if ( ( $debug_trs_found % 100 ) == 0 ) {
                    print $LOGFILE
                        "TRF output lines processed: $debug_trs_found\n";
                }

            }

            #DEBUG
            close $LOGFILE;

            if ( !close TRF2PROCLU ) {
                if ($!) {
                    warn
                        "Error closing trf2proclu process $trf2proclu_pid pipe: $!\n";
                    exit(1002);
                }
                elsif ( $? < -2 ) {
                    warn
                        "trf2proclu process $trf2proclu_pid has returned $?!\n";
                    exit($?);
                }
            }

            #close TRF2PROCLU;

      # check return value
      #my ($rc, $sig, $core) = ($? >> 8, $? & 127, $? & 128);
      #if ($core){
      #    warn "trf2proclu process $trf2proclu_pid dumped core\n";
      #    exit (1000);
      #}elsif($sig == 9){
      #    print  STDERR "trf2proclu process $trf2proclu_pid was murdered!\n";
      #    exit (1001);
      #}elsif ($rc < -1){
      #    warn  "trf2proclu process $trf2proclu_pid has returned $rc!\n";
      #    exit ($rc);
      #}

        }

        close $grandchild;

        # check return value
        my ( $rc, $sig, $core ) = ( $? >> 8, $? & 127, $? & 128 );
        if ($core) {
            warn "run_trf process $grandchild_pid dumped core\n";
            exit(1000);
        }
        elsif ( $sig == 9 ) {
            warn "run_trf process $grandchild_pid was murdered!\n";
            exit(1001);
        }
        elsif ( $rc != 0 ) {
            warn "run_trf process $grandchild_pid has returned $rc!\n";
            exit($rc);
        }

        #if (!close $grandchild) {
        # if ($!)  {
        #   warn "Error closing run_trf process $grandchild_pid pipe: $!\n";
        #   exit(1002);
        # } elsif ($? != 0) {
        #   warn "run_trf process $grandchild_pid has returned $?!\n";
        #   exit($?);
        # }
        #}

        #warn "Exiting child\n";
        exit 0;    # child must never return
    }
    else {
        # parent process -- do nothing
        #warn "This is parent of child $pid\n";
        return $pid;
    }

    die "************ Should never get here\n";
}

sub reverse_complement {
    my $dna = shift;

    # reverse the DNA sequence
    my $revcomp = reverse($dna);

    # complement the reversed DNA sequence
    $revcomp =~ tr/ACGTacgt/TGCAtgca/;
    return $revcomp;
}

sub parse_readfile {
    my ( $trf_fh, $input_seq_fh ) = @_;

    # Grab first line and save it so parser can use it
    my $line = <$input_seq_fh>;

    # Determine file format and run correct parser. Very simplistic.
    my $last_state;
    if ( $line =~ /^>/ ) {

        # File is a FASTA file
        $last_state = parse_fasta( $line, $trf_fh, $input_seq_fh );
    }
    elsif ( $line =~ /^@/ ) {

        # File is a FASTQ file
        $last_state = parse_fastq( $line, $trf_fh, $input_seq_fh );
    }
    else {
        warn
            "File $_ is not of a supported format (FASTA or FASTQ). Skipping this file\n";
        return;
    }

    # Use 0 to indicate a header in parser functions for this to work.
    # Check if the next state the parser was going to read is a header,
    # meaning it last read an entire record.
    die "Read must follow a header\n" if $last_state;
}

# Simple FSM to read a FASTA file
sub parse_fasta {
    my ( $first_line, $trf_fh, $input_seq_fh ) = @_;
    my $header = "";
    my $body   = "";
    chomp($first_line);
    print $trf_fh $first_line . $HEADER_SUFFIX . "\n";

    # Start in state 1 because calling function read header first
    my $read_state = 1;    # read state 1: read, 0: header
    while (<$input_seq_fh>) {
        if ( !$read_state && /^>.*(?:\n|\r)/ ) {

            # previous reversed read
            if ( $reverse_read && $header ne "" ) {

                #warn $header."_RC\n";
                print $trf_fh $header
                    . $HEADER_SUFFIX . "_"
                    . length($body)
                    . "_RCYES\n";
                print $trf_fh reverse_complement($body) . "\n";
            }

            # FASTA header
            chomp;
            print $trf_fh $_ . $HEADER_SUFFIX . "\n";
            $header     = $_;
            $read_state = 1;
        }
        elsif ($read_state) {

            # First line of the read
            if ( $strip_454_TCAG && !s/^TCAG//i ) {
                if ($warn_454_TCAG) {
                    warn
                        "Read does not start with keyseq TCAG. Full sequence: $_\n";
                }
                else {
                    die
                        "Read does not start with keyseq TCAG. Full sequence: $_\n";
                }
            }
            print $trf_fh $_;
            chomp;
            $body       = $_;
            $read_state = 0;
        }
        else {

            # Subsequent lines of the read
            print $trf_fh $_;
            chomp;
            $body .= $_;
        }
    }

    # last reversed read
    if ( $reverse_read && $header ne "" ) {
        print $trf_fh $header
            . $HEADER_SUFFIX . "_"
            . length($body)
            . "_RCYES\n";
        print $trf_fh reverse_complement($body) . "\n";
    }

    return $read_state;
}

# Simple FSM to read a FASTQ file
sub parse_fastq {
    my ( $first_line, $trf_fh, $input_seq_fh ) = @_;

    my $header = "";
    my $body   = "";
    $first_line =~ s/^@/>/;    # Change to FASTA header
    chomp($first_line);
    print $trf_fh $first_line . $HEADER_SUFFIX . "\n";

    # read state 0: header, 1: read, 2: quality header, 3: qualities
    # Start in state 1 since calling function read first line already.
    my $read_state      = 1;
    my $read_line_count = 0;
    while (<$input_seq_fh>) {

        #warn "Read state: $read_state\nLine: $_\n";
        if ( !$read_state && /^@/ ) {

            #warn "In header branch\n";

            # previous reversed read
            if ( $reverse_read && $header ne "" ) {
                print $trf_fh $header
                    . $HEADER_SUFFIX . "_"
                    . length($body)
                    . "_RCYES\n";
                print $trf_fh reverse_complement($body) . "\n";
            }

         # This should be a header. There are some ways to check this is not a
         # quality line, but they depend on the source.
            $_ =~ s/^@/>/;
            chomp;
            print $trf_fh $_ . $HEADER_SUFFIX . "\n";
            $header          = $_;
            $read_state      = 1;
            $read_line_count = 0;
        }
        elsif ( ( $read_state == 1 ) && !/^\+/ ) {

            #warn "In first read line branch\n";

            # First line of the read
            if ( $strip_454_TCAG && !s/^TCAG//i ) {

                if ($warn_454_TCAG) {
                    warn
                        "Read does not start with keyseq TCAG. Full sequence: $_\n";
                }
                else {
                    die
                        "Read does not start with keyseq TCAG. Full sequence: $_\n";
                }
            }
            print $trf_fh $_;
            chomp;
            $body            = $_;
            $read_state      = 2;
            $read_line_count = 1;
        }
        elsif ( ( $read_state == 2 ) && /^\+/ ) {

            #warn "In quality header branch\n";
            $read_state = 3;    #Ignore quality header
        }
        elsif ( $read_state == 3 ) {

#warn "In quality scores branch\nWe read $read_line_count lines of sequence\n";

            # Discard remaining score lines
            # Already consumed one line, so stop at 1
            while ( $read_line_count > 1 ) {
                my $ignore = <$input_seq_fh>;
                $read_line_count--;
            }
            $read_state = 0;    #Next line should be next header or EOF
        }
        elsif ( $read_state == 2 ) {

            #warn "In further read lines branch\n";

            # Should only get here on read sequences spanning multiple lines
            # Subsequent lines of the read
            print $trf_fh $_;
            chomp;
            $body .= $_;
            $read_line_count++;
        }
        else {
            die "Error: bad format (not a FASTQ?). Exiting...";
        }
    }

    # last reversed read
    if ( $reverse_read && $header ne "" ) {
        print $trf_fh $header
            . $HEADER_SUFFIX . "_"
            . length($body)
            . "_RCYES\n";
        print $trf_fh reverse_complement($body) . "\n";
    }

    return $read_state;
}
