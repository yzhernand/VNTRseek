
=head1 NAME

ProcInputReads - We use this module to process input read files for
TRF. The functions in this module take in file names and return
functions which act like reader objects, eg as in BioPerl.

=head1 SYNOPSIS

    use ProcInputReads;
    my 

=head1 DESCRIPTION

This module does not really exist, it
was made for the sole purpose of
demonstrating how POD works.

=cut

package ProcInputReads;
use strict;
use warnings;
use Carp;
use feature 'say';
use Exporter qw(import);

our @EXPORT_OK = qw(fork_proc get_reader formats_regexs compressed_formats_regexs);

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
    = ( \&read_fastaq, \&read_fastaq, \&read_bam );

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

=head2 Methods

=over 12

=item C<formats_regexs>

Returns a hash of sequence format names and the regular expressions
used to determine it from a file name.

=cut

sub formats_regexs {
    return %supported_formats;
}

=item C<compressed_formats_regexs>

Returns a hash of sequence format names and the regular expressions
used to determine it from a file name.

=cut

sub compressed_formats_regexs {
    return %compressed_formats;
}

=item I<fork_proc()>

Takes input/output dirs, the TRF and TRF2proclu command calls with
parameters, flags for reversing reads and 454 tag handling, file
format and compression format names, the number of files processed
so far, and a list of names of files in that format. Then
returns a function which, when called, returns a stream to sequences
which can be used to run TRF/TRF2PROCLU.

Relies on a global, reader_table, which links sequence format names
to functions which can read those formats. Functions in that table
must return a sub which returns a list of exactly two values each
time it is called, and undef when there are no more sequences. The
two values are a FASTA header and a sequence string, in that order.

=cut

sub fork_proc {
    my ($input_dir,        $output_dir,       $trf_param,
        $trf2proclu_param, $reverse_read,     $strip_454_TCAG,
        $warn_454_TCAG,    $format,           $compression,
        $files_processed,  $files_to_process, $filelist
    ) = @_;
    defined( my $pid = fork() )
        or die "Unable to fork: $!\n";
    if ( $pid == 0 ) {    #Child

        my $reader = $reader_table{$format}->(
            $input_dir, $compression, $files_processed,
            $files_to_process, $filelist
        );
        my $output_prefix = "$output_dir/$files_processed";
        exit unless $reader;
        warn "Running child, files_processed = $files_processed...\n";

        # TODO Error checking if TRF, in the start of the pipe, breaks down
        local $SIG{PIPE} = sub { die "Error in trf+trf2proclu pipe: $?\n" };
        open my $trf_pipe,
            "|$trf_param | $trf2proclu_param -o '$output_prefix.index' > '$output_prefix.leb36'"
            or die "Cannot start TRF+trf2proclu pipe: $!\n";

        # TODO Need way of logging TRF output?
        # open my $logfile, ">", "$output_prefix.log"
        #     or die "Error opening file $output_prefix.log: $!\n";
        # $logfile->autoflush;
        # my $debug_reads_processed = 0;

        while ( my ( $header, $body ) = $reader->() ) {

            # say $logfile $data[0] . "\n" . $data[1];
            pipe_to_trf( $reverse_read, $strip_454_TCAG, $warn_454_TCAG,
                $trf_pipe, $header, $body );
                # $trf_pipe, $header, $body, $debug_reads_processed );
            # $debug_reads_processed++;
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

=item I<get_reader()>

Takes input dir, file format and compression format names,
the number of files processed so far, and a list of names of files
in that format. Then returns a function which, when called, returns
a stream to sequences from the read file.

=cut

sub get_reader {
    my ( $input_dir, $format, $compression, $files_processed,
        $files_to_process, $filelist )
        = @_;

    my $reader = $reader_table{$format}->(
        $input_dir, $compression, $files_processed,
        $files_to_process, $filelist
    );

    exit unless $reader;
    return $reader;
}

=item I<reverse_complement()>

Takes a DNA sequence string and returns the reverse complement.

=cut

sub reverse_complement {
    my $dna = shift;

    # reverse the DNA sequence
    my $revcomp = reverse($dna);

    # complement the reversed DNA sequence
    $revcomp =~ tr/ACGTacgt/TGCAtgca/;
    return $revcomp;
}

=item I<pipe_to_trf()>

Takes flags related to reversing reads and handling 454 tags, a file
handle to a TRF process pipe to TRF2proclu, a FASTA header, and a
FASTA sequence (plain sequence string) as input.

Processes input FASTA records to reverse them or remove 454 tags,
as dictated by user options. Prints processed output directly to
given file handle to TRF pipeline.

=cut

sub pipe_to_trf {
    my ( $reverse_read, $strip_454_TCAG, $warn_454_TCAG, $trf_fh, $header,
        $body )
        # $body, $reads_processed )
        = @_;

    # warn "Processing header $header";

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

    # NOTE: Must print reverse complement AFTER the forward read
    say $trf_fh "$header\n$body";

    if ( $reverse_read && $header ne "" ) {
    # if ( ($reads_processed > 1) && $reverse_read && $header ne "" ) {
        say $trf_fh $header . "_"
            . length($body)
            . "_RCYES\n"
            . reverse_complement($body);
    }
}

=item I<read_fasta()>

FASTA file reader.

Given the input dir, file and compression formats, the number
of files processed so far, a reference to a file count, and a list of
all files, return a sub which knows which file it is responsible for
and how to open it.

This returned sub itself returns a list comprising a FASTA header and
sequence each time it is called. When there are no more sequences to
read, it returns an empty list.

This particular sub does NOT modify the $files_to_process value.

=cut

sub read_fasta {
    my ( $input_dir, $compression, $files_processed,
        $files_to_process, $filelist )
        = @_;

# warn "Need to process " . scalar(@$filelist) . " files, working on $files_processed\n";

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

=item I<read_fastq()>

FASTQ file reader.

Given the input dir, file and compression formats, the number
of files processed so far, a reference to a file count, and a list of
all files, return a sub which knows which file it is responsible for
and how to open it.

This returned sub itself returns a list comprising a FASTA header and
sequence each time it is called. When there are no more sequences to
read, it returns an empty list.

This particular sub does NOT modify the $files_to_process value.

=cut

sub read_fastq {

    # Code modified from https://www.biostars.org/p/11599/#11657
    my ( $input_dir, $compression, $files_processed,
        $files_to_process, $filelist )
        = @_;

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

=item I<read_fastaq()>

FASTA/Q file reader. Requires seqtk.

Given the input dir, file and compression formats, the number
of files processed so far, a reference to a file count, and a list of
all files, return a sub which knows which file it is responsible for
and how to open it.

This returned sub itself returns a list comprising a FASTA header and
sequence each time it is called. When there are no more sequences to
read, it returns an empty list.

This particular sub does NOT modify the $files_to_process value.

=cut

sub read_fastaq {
    my ( $input_dir, $compression, $files_processed,
        $files_to_process, $filelist )
        = @_;

# warn "Need to process " . scalar(@$filelist) . " files, working on $files_processed\n";

    # Since we are using seqtk, use pipe open mode
    my $openmode = "-|";

    warn "Processing file " . $filelist->[$files_processed] . "\n";
    my $filename
        = ( ($compression) ? $decompress_cmds{$compression} : "" ) . '"'
        . "$input_dir/"
        . $filelist->[$files_processed] . '"'
        . "| seqtk seq -a -S";

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
        my ( $header, $seq ) = split( /\n+/, $fasta_rec );
        chomp $header;
        chomp $seq;
        # warn "header: '$header'";
        # my $seq = join( "", @seqlines );

        # warn "seq: '$seq'";
        return ( ">". $header, $seq );
    };
}

=item I<read_bam()>

BAM file reader.

Given the input dir, file and compression formats, the number
of files processed so far, a reference to a file count, and a list of
all files, return a sub which knows which file it is responsible for
and how to open it.

This returned sub itself returns a list comprising a FASTA header and
sequence each time it is called. When there are no more sequences to
read, it returns an empty list.

Requires samtools as an external dependency. This sub DOES modify the
value of $files_to_process.

=cut

sub read_bam {
    my ( $input_dir, $compression, $files_processed,
        $files_to_process, $filelist )
        = @_;
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

    # Save the number of samtools commands we'll need to use
    $$files_to_process = scalar @samcmds;

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

=item I<read_FORMAT() stub>

Example file reader.

Given the input dir, file and compression formats, the number
of files processed so far, a reference to a file count, and a list of
all files, return a sub which knows which file it is responsible for
and how to open it.

This returned sub itself returns a list comprising a FASTA header and
sequence each time it is called. When there are no more sequences to
read, it returns an empty list.

This is an example of how to write a sub which reads a given file
in a certain format and returns a subroutine which works as expected
by fork_proc().

sub read_FORMAT {
    my ( $input_dir, $compression, $files_processed, $files_to_process, $filelist )
        = @_;

    # Here make some decisions. For instance, if we expect to process multiple
    # files, such as FASTA or FASTQ files and unlike BAM files, then choose
    # the file to read using $files_processed. Otherwise, logically split
    # the input and decide which split this call will work on. You may have to
    # write the next section before this one in the case of multiple files.

    # In the case of files like BAM files, where there is only one input
    # file but the file is indexed in regions which can be read directly,
    # save the number of splits needed into the variable pointed to by
    # $files_to_process. NOT needed if there are multiple input files,
    # one for each child process.

    # Open a file handle to the file or pipe we will read from.
    # If file might be compressed, use the $compression value and
    # the decompress_cmds hash to construct the input pipe command.
    my $input_fh = open...

    # Handle any errors, and die if there is a fatal problem so this file
    # or split isn't considered further and another process can jump in on
    # the next item in the list.

    # Return an anonymous subroutine which will use that file handle.
    return sub {
        # Read from input file handle
        my $record = <$input_fh>;
        # Return empty list if nothing else to read
        return () unless ($record);

        # Process $record as needed, read more lines, etc (see other functions
        # for examples).

        # Return a list consisting of the header, prefixed by a ">", and the
        # sequence string.
        return(">" . $header, $seq);
    }
}

=back

=head1 LICENSE

This is released under the Artistic 
License. See L<perlartistic>.

=head1 AUTHOR

Yozen Hernandez
Yevgeniy Gelfand

=head1 SEE ALSO

L<perl>

=cut
