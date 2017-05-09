#!/usr/bin/env perl

use 5.010;
use strict;
use warnings;

# Test reading FASTA files and outputting a sequence and header.
die "Usage: $0 <format> <compression> <fasta file 1> ... <fasta file N>"
    unless @ARGV >= 3;
my $max_processes = 8;
my ( $format, $compression, @filenames ) = @ARGV;
my %decompress_cmds = (
    targz   => "tar xzfmO ",
    gzip    => "gunzip -c ",
    tarbzip => "tar xjfmO ",
    bzip    => "bzip2 -c ",
    tarxz   => "tar xJfmO ",
    xz      => "xzcat "
);

my $make_bed_file_loc = "~/src/VNTRseek/perl";

# Setting here because needed in full program
my $input_dir = "";
my $output_dir = ".";

for my $file (@filenames) {
	say "Processing file: $file";
    my $fasta_reader = read_fasta( 8, "", \@filenames );
    open my $outfile, ">", $file . ".out"
    	or die "Error opening file: " . $file . ".out";
    while ( my @data = $fasta_reader->() ) {
        say $outfile $data[0] . "\n" . $data[1];
    }
    say "";
}

sub read_fasta {
    my ( $max_processes, $compression, $filelist ) = @_;
    state $files_processed = 0;

 # If file uncompressed, simply read from file. Else open a pipe to a command.
    my $openmode = ($compression) ? "-|" : "<";
    # warn $openmode;
    my $filename
        = ($compression)
        ? $decompress_cmds{$compression} . '"'
        . "'$input_dir/"
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


sub read_bam {
    my ( $max_processes, $compression, $filelist ) = @_;
    state $files_processed = 0;
    my $bamfile = $filelist->[0];
    my $baifile = $bamfile . ".bai";
	
	# Check if .bai file exists and then run MakeBedFiles.jar
	die "Error reading bam file: corresponding bai file required for processing bam files."
		unless (-e -r $baifile);

	# Requires samtools and java to be installed/available
	system("samtools idxstats $bamfile | java -jar $make_bed_file_loc $max_processes $output_dir");
	# Produces $max_processes bed files. These are all named "bedN" where N is 0 to $max_processes-1
	my $filename = "$output_dir/bed" . $files_processed++;
	open my $bed_fh, "<", $filename
        or die "Error opening file " . $filename;
    return sub {
    	
    }
}