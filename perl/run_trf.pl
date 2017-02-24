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

my $files_to_process = 1;    # number of files to process in one batch
my $files_processed  = 0;    # files processed
my %p;                       # associates forked pids with output pipe pids
my $max_processes = 0;

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
my $tgz_dir = $ARGV[0];
die "Need to provide input directory\n" if !defined $tgz_dir;
die "Input directory does not exist\n"  if !-d $tgz_dir;

# Output directory
my $output_dir = $ARGV[1];
die "Need to provide output directory\n" if !defined $output_dir;

#die "Output directory already exists -- please delete it before running this program\n" if -e $output_dir;
mkdir $output_dir;

# get a list of input files
opendir( DIR, $tgz_dir );

# the only extensions are .tgz, .tar.gz, and .gz
my @tarballs = sort grep( /^fasta.*\.(?:tgz|gz)$/, readdir(DIR) );

unless (@tarballs) {    # Look for fastq files if no fasta files found
    rewinddir(DIR);
    @tarballs = sort grep( /^fastq.*\.(?:tgz|gz)$/, readdir(DIR) );
}
closedir(DIR);

my $tarball_count = @tarballs;
print STDERR "$tarball_count supported files found in $tgz_dir\n";
die "Exiting\n" if $tarball_count == 0;

# /proc/cpuinfo is found on Linux only: use other methods for other platforms
if ( $max_processes == 0 ) {
    if ( open( CPUINFO, "/proc/cpuinfo" ) ) {
        while (<CPUINFO>) {
            $max_processes++ if /^processor\b/;
        }
        close(CPUINFO);
        print STDERR "$max_processes CPU core(s) detected\n";
    }
    else {
        warn "Could not open /proc/cpuinfo: Assuming single CPU\n";
        $max_processes = 1;
    }
    if ( $max_processes == 0 ) {
        warn "Unknown formatting in /proc/cpuinfo: Assuming single CPU\n";
        $max_processes = 1;
    }
}

# because of the limit of number of open files at the same time, let's keep number of files to under 200
# redundancy step (3) opens them all at the same time, and if multiple pipelines are run it might be possible
# to reach the limit on some system (usually 1024)
while ( int( $tarball_count / $files_to_process ) > 200 ) {
    $files_to_process++;
}
print STDERR "Will use $max_processes processes\n";
print STDERR "Will process $files_to_process file(s) per batch\n";

# fork as many new processes as there are CPUs
for ( my $i = 0; $i < $max_processes; $i++ ) { $p{ fork_trf() } = 1 }

# wait for processes to finish and then fork new ones
while ( ( my $pid = wait ) != -1 ) {

    # check return value
    my ( $rc, $sig, $core ) = ( $? >> 8, $? & 127, $? & 128 );
    if ($core) {
        print STDERR "run_trf process $pid dumped core\n";
        exit(1000);
    }
    elsif ( $sig == 9 ) {
        print STDERR "run_trf process $pid was murdered!\n";
        exit(1001);
    }
    elsif ( $rc != 0 ) {
        print STDERR "run_trf process $pid has returned $rc!\n";
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
print STDERR "Processing complete -- processed $files_processed file(s).\n";
0;

############################ Procedures ###############################################################

sub fork_trf {
    if ( $files_processed >= $tarball_count ) {
        return 0;
    }

    # unzip a predefined number of files
    my $until = $files_processed + $files_to_process - 1;
    $until = $tarball_count - 1 if $until > ( $tarball_count - 1 );
    print STDERR 'Processing files '
        . ( $files_processed + 1 ) . ' to '
        . ( $until + 1 ) . "\n";
    my $output_prefix    = "$output_dir/$files_processed-$until";
    my @file_slice       = @tarballs[ ($files_processed) .. ($until) ];
    my $file_slice_count = @file_slice;
    $files_processed += $files_to_process;

    defined( my $pid = fork )
        or die "Unable to fork: $!\n";
    if ( $pid == 0 ) {

        #print STDERR "This is child\n";
        defined( my $grandchild_pid = open GRANDCHILD, '-|' )
            or die "Unable to open grandchild: $!\n";
        if ( $grandchild_pid == 0 ) {

            #print STDERR "Starting TRF: $TRF_PARAM";
            defined( my $trf_pid = open( TRF, "| $TRF_PARAM" ) )
                or die "Cannot start TRF: $!\n";
        FILE:
            foreach (@file_slice) {
                next FILE if not defined $_;

                # TODO: not parsing FASTA files at the moment:
                # assuming all extracted files are in FASTA format
                print STDERR "Unzipping $tgz_dir/$_\n";
                if ( $_ =~ /\.(?:tgz|tar\.gz)$/ ) {
                    open( FASTA_IN, "tar xzfmoO '$tgz_dir/$_' |" );
                }
                elsif ( $_ =~ /\.gz$/ ) {
                    open( FASTA_IN, "gunzip -c '$tgz_dir/$_' |" );
                }
                else {
                    warn "File $_ has wrong extension. Skipping this file\n";
                    next FILE;
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
                parse_readfile();
                close FASTA_IN;
            }

            #close TRF;

            # check return value
            #my ($rc, $sig, $core) = ($? >> 8, $? & 127, $? & 128);
            #if ($core){
            #    print STDERR "trf process $trf_pid dumped core\n";
            #    exit (1000);
            #}elsif($sig == 9){
            #    print  STDERR "trf process $trf_pid was murdered!\n";
            #    exit (1001);
            #}elsif ($rc != 0){
            #    print STDERR  "trf process $trf_pid has returned $rc!\n";
            #    exit ($rc);
            #}

            if ( !close TRF ) {
                if ($!) {
                    print STDERR
                        "Error closing trf process $trf_pid pipe: $!\n";
                    exit(1002);
                }
                elsif ( $? != 0 ) {
                    print STDERR "trf process $trf_pid has returned $?!\n";
                    exit($?);
                }
            }

            #print STDERR "Exiting grandchild\n";
            exit 0;
        }
        else {
            #print STDERR "This is parent of grandchild $grandchild_pid\n";
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

            while (<GRANDCHILD>) {
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
                    print STDERR
                        "Error closing trf2proclu process $trf2proclu_pid pipe: $!\n";
                    exit(1002);
                }
                elsif ( $? < -2 ) {
                    print STDERR
                        "trf2proclu process $trf2proclu_pid has returned $?!\n";
                    exit($?);
                }
            }

            #close TRF2PROCLU;

  # check return value
  #my ($rc, $sig, $core) = ($? >> 8, $? & 127, $? & 128);
  #if ($core){
  #    print STDERR "trf2proclu process $trf2proclu_pid dumped core\n";
  #    exit (1000);
  #}elsif($sig == 9){
  #    print  STDERR "trf2proclu process $trf2proclu_pid was murdered!\n";
  #    exit (1001);
  #}elsif ($rc < -1){
  #    print STDERR  "trf2proclu process $trf2proclu_pid has returned $rc!\n";
  #    exit ($rc);
  #}

        }

        close GRANDCHILD;

        # check return value
        my ( $rc, $sig, $core ) = ( $? >> 8, $? & 127, $? & 128 );
        if ($core) {
            print STDERR "run_trf process $grandchild_pid dumped core\n";
            exit(1000);
        }
        elsif ( $sig == 9 ) {
            print STDERR "run_trf process $grandchild_pid was murdered!\n";
            exit(1001);
        }
        elsif ( $rc != 0 ) {
            print STDERR
                "run_trf process $grandchild_pid has returned $rc!\n";
            exit($rc);
        }

  #if (!close GRANDCHILD) {
  # if ($!)  {
  #   print STDERR "Error closing run_trf process $grandchild_pid pipe: $!\n";
  #   exit(1002);
  # } elsif ($? != 0) {
  #   print STDERR "run_trf process $grandchild_pid has returned $?!\n";
  #   exit($?);
  # }
  #}

        #print STDERR "Exiting child\n";
        exit 0;    # child must never return
    }
    else {
        # parent process -- do nothing
        #print STDERR "This is parent of child $pid\n";
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

    # Grab first line and save it so parser can use it
    my $line = <FASTA_IN>;

    # Determine file format and run correct parser. Very simplistic.
    my $last_state;
    if ( $line =~ /^>/ ) {

        # File is a FASTA file
        $last_state = parse_fasta($line);
    }
    elsif ( $line =~ /^@/ ) {

        # File is a FASTQ file
        $last_state = parse_fastq($line);
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

sub parse_fasta {

    # Simple FSM to read a FASTA file
    my $header     = "";
    my $body       = "";
    my $first_line = shift;
    chomp($first_line);
    print TRF $first_line . $HEADER_SUFFIX . "\n";

    # Start in state 1 because calling function read header first
    my $read_state = 1;    # read state 1: read, 0: header
    while (<FASTA_IN>) {
        if ( !$read_state && /^>.*(?:\n|\r)/ ) {

            # previous reversed read
            if ( $reverse_read && $header ne "" ) {

                #print STDERR $header."_RC\n";
                print TRF $header
                    . $HEADER_SUFFIX . "_"
                    . length($body)
                    . "_RCYES\n";
                print TRF reverse_complement($body) . "\n";
            }

            # FASTA header
            chomp;
            print TRF $_ . $HEADER_SUFFIX . "\n";
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
            print TRF $_;
            chomp;
            $body       = $_;
            $read_state = 0;
        }
        else {

            # Subsequent lines of the read
            print TRF $_;
            chomp;
            $body .= $_;
        }
    }

    # last reversed read
    if ( $reverse_read && $header ne "" ) {
        print TRF $header . $HEADER_SUFFIX . "_" . length($body) . "_RCYES\n";
        print TRF reverse_complement($body) . "\n";
    }

    return $read_state;
}

sub parse_fastq {

    # Simple FSM to read a FASTQ file
    my $header     = "";
    my $body       = "";
    my $first_line = shift;
    $first_line =~ s/^@/>/;    # Change to FASTA header
    chomp($first_line);
    print TRF $first_line . $HEADER_SUFFIX . "\n";

    # read state 0: header, 1: read, 2: quality header, 3: qualities
    # Start in state 1 since calling function read first line already.
    my $read_state      = 1;
    my $read_line_count = 0;
    while (<FASTA_IN>) {

        #warn "Read state: $read_state\nLine: $_\n";
        if ( !$read_state && /^@/ ) {

            #warn "In header branch\n";

            # previous reversed read
            if ( $reverse_read && $header ne "" ) {
                print TRF $header
                    . $HEADER_SUFFIX . "_"
                    . length($body)
                    . "_RCYES\n";
                print TRF reverse_complement($body) . "\n";
            }

         # This should be a header. There are some ways to check this is not a
         # quality line, but they depend on the source.
            $_ =~ s/^@/>/;
            chomp;
            print TRF $_ . $HEADER_SUFFIX . "\n";
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
            print TRF $_;
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
                my $ignore = <FASTA_IN>;
                $read_line_count--;
            }
            $read_state = 0;    #Next line should be next header or EOF
        }
        elsif ( $read_state == 2 ) {

            #warn "In further read lines branch\n";

            # Should only get here on read sequences spanning multiple lines
            # Subsequent lines of the read
            print TRF $_;
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
        print TRF $header . $HEADER_SUFFIX . "_" . length($body) . "_RCYES\n";
        print TRF reverse_complement($body) . "\n";
    }

    return $read_state;
}
