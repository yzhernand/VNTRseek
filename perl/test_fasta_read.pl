#!/usr/bin/env perl

use 5.010;
use strict;
use warnings;

# Test reading FASTA files and outputting a sequence and header.
die "Usage: $0 <format> <compression> <fasta file 1> ... <fasta file N>"
    unless @ARGV >= 3;
my $max_processes = 2;
my ( $format, $compression, @filenames ) = @ARGV;
my %decompress_cmds = (
    targz   => "tar xzfmO ",
    gzip    => "gunzip -c ",
    tarbzip => "tar xjfmO ",
    bzip    => "bzip2 -c ",
    tarxz   => "tar xJfmO ",
    xz      => "xzcat "
);
my %reader_table = (
    fasta => \&read_fasta,
    fastq => \&read_fastq,
    bam   => \&read_bam
);

my $make_bed_file_loc = "~/src/VNTRseek/java/MakeBedFiles.jar";

# Setting here because needed in full program
my $input_dir   = "";
my $output_dir  = ".";
my $out_counter = 0;
my %p;

# FIXME Need to return 0 from fork parent to stop spawning procs
for ( my $i = 0; $i < $max_processes; $i++ ) {
    my $reader = $reader_table{$format};
    my $pid    = fork();
    die "Unable to fork: $!\n" unless defined($pid);
    if ( $pid != 0 ) {    #Parent
        $p{$pid} = 1;
    }
    else {
        my $reader = $reader_table{$format}->("", \@filenames);
        open my $outfile, ">", "$output_dir/$out_counter.out"
            or die "Error opening file $output_dir/$out_counter.out: $!\n";
        while ( my @data = $reader->() ) {
            say $outfile $data[0] . "\n" . $data[1];
        }
        say "";
        $out_counter++;
    }
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
        my $reader = $reader_table{$format};
        my $newpid = fork();
        die "Unable to fork: $!\n" unless defined($newpid);
        if ( $newpid != 0 ) {    #Parent
            $p{$newpid} = 1;
        }
        else {
            my $reader = $reader_table{$format}->("", \@filenames);
            open my $outfile, ">", "$output_dir/$out_counter.out"
                or die
                "Error opening file $output_dir/$out_counter.out: $!\n";
            while ( my @data = $reader->() ) {
                say $outfile $data[0] . "\n" . $data[1];
            }
            say "";
            $out_counter++;
        }
    }
    else {
        die "ERROR: Do not remember process PID=$pid\n";
    }
}

# while ( my $fasta_reader = read_bam( "", \@filenames ) ) {
#     last if $out_counter == 4;

#     # my $fasta_reader = read_fasta( 8, "", \@filenames );
#     open my $outfile, ">", "$output_dir/$out_counter.out"
#         or die "Error opening file $output_dir/$out_counter.out: $!\n";
#     while ( my @data = $fasta_reader->() ) {
#         say $outfile $data[0] . "\n" . $data[1];
#     }
#     say "";
#     $out_counter++;
# }

sub read_fasta {
    my ( $compression, $filelist ) = @_;
    state $files_processed = 0;
    warn $files_processed;

    # Don't do more if we've exhausted the file list
    return undef if $files_processed == @$filelist;

 # If file uncompressed, simply read from file. Else open a pipe to a command.
    my $openmode = ($compression) ? "-|" : "<";

    say "Processing file: " . $filelist->[$files_processed];

    # warn $openmode;
    my $filename
        = ($compression)
        ? $decompress_cmds{$compression} . '"'
        . "$input_dir/"
        . $filelist->[$files_processed] . '"'
        : $filelist->[$files_processed];
    $files_processed++;

    # $files_processed contains how many files processed so far.
    # Use to index into filelist
    # warn $filename;
    local $/ = ">";
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
    #
}

# Requires samtools
sub read_bam {
    my ( $compression, $filelist ) = @_;
    state $files_processed = 0;
    my $bamfile = "$input_dir/" . $filelist->[0];
    # warn "$bamfile\n";
    my $baifile = $bamfile . ".bai";
    # warn "$baifile\n";
    my $unmapped_template = "*";
    state @samcmds;

    if ( $files_processed == 0 ) {

        # Check if .bai file exists and then run MakeBedFiles.jar
        die
            "Error reading bam file: corresponding bai file required for processing bam files."
            unless ( -e -r $baifile );

        # Requires samtools and java to be installed/available
        # Make only one bed file
        system(
            "samtools idxstats \"$bamfile\" | java -jar $make_bed_file_loc 1 $output_dir"
        );

# Produces $max_processes bed files. These are all named "bedN" where N is 0 to $max_processes-1
# Alternative: produce ONE bed file, read in all regions, split these into samtools procs over each TRF proc (like split FASTA files)
        my $filename = "$output_dir/bed0";
        open my $bed_fh, "<", $filename
            or die "Error opening file " . $filename;

# We need to read the bed file, line by line, and construct samtools commands
# These are all saved in an array which are interated through like the FASTA files before.
# Then we catch the output of these through a file handle, and process into FASTA.
        while ( my $bedline = <$bed_fh> ) {
            my ( $chr, $start, $end ) = split /\t/, $bedline;
            my $unmapped = ( $chr eq $unmapped_template );
            my $region   = "$chr:$start-$end";
            my $scmd
                = "samtools view "
                . ( ($unmapped) ? "-f 4 " : "" )
                . $bamfile
                . ( ($unmapped) ? " $region" : "" );
            push @samcmds, $scmd;
        }
    }

    # Don't do more if we've exhausted the file list
    return undef if $files_processed == @samcmds;
    warn "Processing bam chunk using: " . $samcmds[ $files_processed ] . "\n";

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
