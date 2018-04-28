#!/usr/bin/perl

use strict;
use warnings;
use Cwd;
use DBI;
use List::Util qw[min max];
use POSIX qw(strftime);
use Carp qw(croak carp);
use FindBin;
use File::Basename;
use Data::Dumper;
use lib "$FindBin::RealBin/lib";

use vutil
    qw(get_config get_dbh get_ref_dbh get_trunc_query set_statistics get_statistics);

#use GD::Graph::linespoints;

# Argument parsing & set up
my $argc = @ARGV;

if ( $argc < 7 ) {
    die
        "Usage: updaterefs.pl read_profiles_folder read_profiles_folder_clean dbname run_dir file_representative latexfile VERSION\n";
}

my $curdir        = getcwd;
my $readpf        = $ARGV[0];
my $rpfc          = $ARGV[1];
my $DBSUFFIX      = $ARGV[2];
my $run_dir       = $ARGV[3];
my $filerep       = $ARGV[4];
my $result_prefix = $ARGV[5];
my $VERSION       = $ARGV[6];

# set these mysql credentials in vs.cnf (in installation directory)
my %run_conf = get_config( $DBSUFFIX, $run_dir );
my ( $HTTPSERVER, $MIN_SUPPORT_REQUIRED, $TEMPDIR )
    = @run_conf{qw(SERVER MIN_SUPPORT_REQUIRED TMPDIR)};

############################ Procedures ###############################################################
sub RC {

    # complement reversed DNA sequence
    my $seq = shift;

    $seq = reverse $seq;

    $seq =~ tr/ACGT/TGCA/;

    return $seq;
}

sub write_vcf_rec {
    croak "write_vcf_rec requires three arguments"
        unless @_ == 3;
    my ( $spanN_vcffile, $allwithsupport_vcffile, $supported_tr ) = @_;

    # warn ">", $supported_tr->{rid}, "\n" if $ENV{DEBUG};
    warn "Supported TR hash (before joins):\n", Dumper( \$supported_tr ), "\n"
        if ( $ENV{DEBUG} );

    # extra code for single allele, v 1.07
    if ( $supported_tr->{alleles_supported} == 1 ) {
        $supported_tr->{gt_string}
            = ( $supported_tr->{subSameAsRef1} == 1 ) ? "0/0" : "1/1";
    }
    else {
        $supported_tr->{gt_string} = join "/",
            @{ $supported_tr->{gt_string} };
    }

    my $qual = ".";
    if ( "" eq $supported_tr->{seq} ) { $supported_tr->{seq} = "."; }
    if ( "" eq $supported_tr->{alt} ) { $supported_tr->{alt} = "."; }

    my $filter = ( $supported_tr->{singleton} == 1 ) ? "PASS" : "SC";

    $supported_tr->{read_support} = join ",",
        @{ $supported_tr->{read_support} };
    $supported_tr->{copy_diff} = join ",", @{ $supported_tr->{copy_diff} };
    $supported_tr->{alt}       = join ",", @{ $supported_tr->{alt} };

    warn "Supported TR hash (after joins):\n", Dumper( \$supported_tr ), "\n"
        if ( $ENV{DEBUG} );

    my $info = sprintf(
        "RC=%.2lf;RPL=%d;RAL=%d;RCP=%s;ALGNURL=http://%s/index.php?db=VNTRPIPE_%s&ref=-%d&isref=1&istab=1&ispng=1&rank=3",
        $supported_tr->{copiesfloat},
        length( $supported_tr->{consenuspat} ),
        $supported_tr->{arlen},
        $supported_tr->{consenuspat},
        $HTTPSERVER,
        $DBSUFFIX,
        $supported_tr->{rid}
    );
    my $format = "GT:SP:CGL";

    # Commented out since this should always be true
    # if ( $supported_tr->{alleles_supported} > 0 ) {
    my $vcf_rec = join(
        "\t",
        $supported_tr->{head},
        ( $supported_tr->{pos1} - 1 ),
        "td" . $supported_tr->{rid},
        $supported_tr->{seq},
        $supported_tr->{alt},
        $qual, $filter, $info, $format,
        join( ":",
            $supported_tr->{gt_string}, $supported_tr->{read_support},
            $supported_tr->{copy_diff} )
    ) . "\n";

    $supported_tr->{is_called_vntr} && print $spanN_vcffile $vcf_rec;
    print $allwithsupport_vcffile $vcf_rec;

    # }
}

####################################
sub update_ref_table {
    croak "update_ref_table requires two arguments"
        unless @_ == 2;
    my ( $ref, $update_sth ) = @_;
    return
        unless exists $ref->{refid};

    if ( $ref->{readsum} >= 1 ) {
        $ref->{span1} = 1;
    }
    if ( $ref->{readsum} >= $MIN_SUPPORT_REQUIRED ) {
        $ref->{spanN} = 1;
    }

    # hetez_multi is true if there is support
    # for more than one allele
    # Other genotype classes can only be true if hetez_multi
    # is false.
    # TODO Add global config for ploidy. Replace '2'
    # below with the value for the ploidy (eg, haploid == 1)
    # Modify last branch to deal with special haploid case
    if ( $ref->{nsupport} > 2 ) {
        $ref->{hetez_multi} = 1;
    }
    elsif ( $ref->{nsupport} == 2 ) {
        $ref->{hetez_same} = 1 * ( $ref->{nsameasref} );
        $ref->{hetez_diff} = 1 * !( $ref->{nsameasref} );
    }
    elsif ( $ref->{nsupport} == 1 ) {
        $ref->{homez_same} = 1 * ( $ref->{nsameasref} );
        $ref->{homez_diff} = 1 * !( $ref->{nsameasref} );
    }

    $ref->{support_vntr}
        = 1 * ($ref->{homez_diff}
            || $ref->{hetez_same}
            || $ref->{hetez_diff}
            || $ref->{hetez_multi} );

    warn "Updating fasta_ref_reps table with: ", Dumper($ref), "\n"
        if ( $ENV{DEBUG} );
    $update_sth->execute(
        $ref->{refid},        $ref->{nsupport},
        $ref->{nsameasref},   $ref->{entr},
        $ref->{has_support},  $ref->{span1},
        $ref->{spanN},        $ref->{homez_same},
        $ref->{homez_diff},   $ref->{hetez_same},
        $ref->{hetez_diff},   $ref->{hetez_multi},
        $ref->{support_vntr}, $ref->{support_vntr_span1}
    ) or die "Cannot execute: " . $update_sth->errstr;
}

#open (VCFFILE2, ">", "${latex}.span2.vcf") or die "Can't open for reading ${latex}.vcf.";

####################################
sub nowhitespace($) {
    my $string = shift;
    $string =~ s/\s+//g;
    return $string;
}

####################################
sub commify {
    my $input = shift;
    carp "Undef input" unless defined $input;
    $input = reverse $input;
    $input =~ s/(\d\d\d)(?=\d)(?!\d*\.)/$1,/g;
    return reverse $input;
}

####################################
# Takes a boolean as an argument. If the boolean is anything perl considers false,
# then this function will only produce a VCF file for supported VNTRs. If true, all
# supported TRs are included in the file.
sub print_vcf {

    # my $allwithsupport = shift;

    # Get needed stats
    my @stats = qw(PARAM_TRF
        FILE_REFERENCE_SEQ
        FILE_REFERENCE_LEB
        NUMBER_REF_TRS
        FOLDER_FASTA
        FOLDER_PROFILES
        NUMBER_READS
        NUMBER_TRS_IN_READS);
    my $stat_hash = get_statistics(@stats);
    my $dbh = get_dbh( { readonly => 1, userefdb => 1 } );

    # Get total number of TRs supported
    my $numsup_sth
        = $dbh->prepare(
        "select count(distinct vntr_support.refid) FROM vntr_support WHERE support >= $MIN_SUPPORT_REQUIRED;"
        ) or die "Couldn't prepare statement: " . $dbh->errstr;

    my $numsup = 0;
    $numsup_sth->execute() or die "Cannot execute: " . $numsup_sth->errstr();
    $numsup_sth->bind_columns( \$numsup );
    $numsup_sth->fetch;
    if ( !defined $numsup ) {
        die "Error getting number of supported TRs: " . $dbh->errstr;
    }
    $numsup_sth->finish;

    # Get number of VNTRs supported
    my $numvntrs_sth
        = $dbh->prepare(
        "select count(*) FROM fasta_ref_reps WHERE support_vntr>0;")
        or die "Couldn't prepare statement: " . $dbh->errstr;

    my $numvntrs = 0;
    $numvntrs_sth->execute()
        or die "Cannot execute: " . $numvntrs_sth->errstr();
    $numvntrs_sth->bind_columns( \$numvntrs );
    $numvntrs_sth->fetch;
    $numvntrs_sth->finish;
    if ( !defined $numsup ) {
        die "Error getting number of supported VNTRs: " . $dbh->errstr;
    }

    # $update_spanN_sth->finish;
    my $spanN_fn = "${result_prefix}.span${MIN_SUPPORT_REQUIRED}.vcf";
    open my $spanN_vcffile, ">", $spanN_fn
        or die "Can't open for writing $spanN_fn\n\n";
    my $allwithsupport_fn
        = "${result_prefix}.allwithsupport.span${MIN_SUPPORT_REQUIRED}.vcf";
    open my $allwithsupport_vcffile, ">", $allwithsupport_fn
        or die "Can't open for writing $allwithsupport_fn\n\n";

    my $vcf_header
        = "##fileformat=VCFv4.1\n"
        . strftime( "##fileDate=\"%Y%m%d\"\n", localtime )
        . qq[##source="Vntrseek ver. $VERSION"
##TRFParameters="$stat_hash->{PARAM_TRF}"
##referenceseq="$stat_hash->{FILE_REFERENCE_SEQ}"
##referenceprofile="$stat_hash->{FILE_REFERENCE_LEB}"
##numrefTRs="$stat_hash->{NUMBER_REF_TRS}"
##readseqfolder="$stat_hash->{FOLDER_FASTA}"
##readprofilefolder="$stat_hash->{FOLDER_PROFILES}"
##numreads="$stat_hash->{NUMBER_READS}"
##numreadTRs="$stat_hash->{NUMBER_TRS_IN_READS}"
##numVNTRs="$numvntrs"
##numTRsWithSupport="$numsup"
##database="VNTRPIPE_$DBSUFFIX"
##databaseurl="http://${HTTPSERVER}/result.php?db=VNTRPIPE_${DBSUFFIX}"
##INFO=<ID=RC,Number=1,Type=Float,Description="Reference Copies">
##INFO=<ID=RPL,Number=1,Type=Integer,Description="Reference Pattern Length">
##INFO=<ID=RAL,Number=1,Type=Integer,Description="Reference Tandem Array Length">
##INFO=<ID=RCP,Number=1,Type=String,Description="Reference Consensus Pattern">
##INFO=<ID=ALGNURL,Number=1,Type=String,Description="Alignment URL">
##FILTER=<ID=SC,Description="Reference is Singleton">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=SP,Number=A,Type=Integer,Description="Number of Spanning Reads">
##FORMAT=<ID=CGL,Number=A,Type=Integer,Description="Copies Gained or Lost with respect to reference">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t$DBSUFFIX
];

    print $spanN_vcffile $vcf_header;
    print $allwithsupport_vcffile $vcf_header;

# Get information on all VNTRs
# "SELECT rid,alleles_sup,allele_sup_same_as_ref,is_singleton,is_dist,is_indist,firstindex,lastindex,copynum,pattern,clusterid,reserved,reserved2,head,sequence,flankleft,direction FROM fasta_ref_reps INNER JOIN clusterlnk ON fasta_ref_reps.rid=-clusterlnk.repeatid INNER JOIN clusters ON clusters.cid=clusterlnk.clusterid WHERE support_vntr>0 ORDER BY head, firstindex;"
# $dbh->sqlite_create_function(
#     'mkflank',
#     2,
#     sub {
#         my ( $lflank, $rflank ) = @_;
#         return substr( $lflank, -60 ) . "|" . substr( $rflank, 0, 60 );
#     }
# );
    my $get_supported_reftrs_query
        = qq{SELECT reftab.rid, is_singleton, is_dist, firstindex,
            (lastindex - firstindex) + 1 AS arlen, reftab.copynum, reftab.pattern, reftab.head,
            sequence, c1.direction AS refdir, copies, sameasref,
            support, first, last, dna, c2.direction AS readdir, support_vntr
        FROM refdb.fasta_ref_reps reftab INNER JOIN vntr_support ON reftab.rid=-vntr_support.refid
        INNER JOIN main.fasta_ref_reps mainreftab ON mainreftab.rid=-vntr_support.refid
        INNER JOIN clusterlnk c1 ON vntr_support.refid=c1.repeatid
        INNER JOIN replnk ON vntr_support.representative=replnk.rid
        INNER JOIN clusterlnk c2 ON c2.repeatid=replnk.rid
        INNER JOIN fasta_reads ON replnk.sid=fasta_reads.sid
        ORDER BY reftab.head ASC, reftab.firstindex ASC, sameasref DESC};
    my $get_supported_reftrs_sth = $dbh->prepare($get_supported_reftrs_query);
    $get_supported_reftrs_sth->execute()
        or die "Cannot execute: " . $get_supported_reftrs_sth->errstr();
    my ($rid,     $singleton,   $disting,     $pos1,
        $arlen,   $copiesfloat, $consenuspat, $head,
        $seq,     $refdir,      $copies,      $sameasref,
        $support, $readTRStart, $readTRStop,  $dna,
        $readdir, $support_vntr
    );
    $get_supported_reftrs_sth->bind_columns(
        \(  $rid,     $singleton,   $disting,     $pos1,
            $arlen,   $copiesfloat, $consenuspat, $head,
            $seq,     $refdir,      $copies,      $sameasref,
            $support, $readTRStart, $readTRStop,  $dna,
            $readdir, $support_vntr
        )
    );

    # Loop over all supported ref TRs
    my $supported_tr = { rid => -1 };

# my $oldrefid     = -1;
# my ( $alleleWithSupportFound, $first, $subSameAsRef1, $j, $al, $patlen )
#     = (0) x 6;
# my ( $gt_string, $read_support, $copy_diff, $read_support1, $copy_diff1, $alt )
#     = ("") x 6;
# TODO Doesn't work properly. Look at original code. Might need to rewrite completely
    while ( $get_supported_reftrs_sth->fetch() ) {
        if (   $supported_tr->{rid} != -1
            && $supported_tr->{rid} != $rid )
        {
            write_vcf_rec( $spanN_vcffile, $allwithsupport_vcffile,
                $supported_tr );
        }

        if ( $supported_tr->{rid} != $rid ) {
            $supported_tr = {
                rid   => $rid,
                first => $copies,

                # Start at 0 if the ref allele was supported, else start at 1
                al => ($sameasref) ? 0 : 1,
                arlen             => $arlen,
                pos1              => $pos1,
                copiesfloat       => $copiesfloat,
                singleton         => $singleton,
                head              => $head,
                seq               => ($seq) ? uc( nowhitespace($seq) ) : "",
                consenuspat       => $consenuspat,
                subSameAsRef1     => $sameasref,
                is_called_vntr    => $support_vntr,
                alleles_supported => 0,
                gt_string         => [],
                read_support      => [],
                copy_diff         => [],
                read_support1     => "",
                copy_diff1        => "",
                alt               => []
            };
        }

        # Loop over all allele supporting read TRs
        # else {
        $dna = ($dna) ? uc( nowhitespace($dna) ) : "";
        my $cdiff = $copies;
        ( $ENV{DEBUG} ) && warn "\tAllele: $cdiff, support $support, alt: ",
            Dumper( $supported_tr->{alt} ), "\n";

        # If an alternate allele
        if ( 0 == $sameasref ) {

            # If there is no DNA string for this read, exit with error
            if ( !$dna || $dna eq "" ) {
                die
                    "Error: read source sequence not found in database for ref ($rid) alternate allele",
                    $supported_tr->{al}, "! Stopped at";
            }

# If the stop position is beyond the length of the read sequence, exit with error
            if ( $readTRStop > ( length($dna) ) ) {
                die
                    "Error: last position outside of the range of the source sequence buffer for ref ($rid) alternate allele",
                    $supported_tr->{al}, "! Stopped at";
            }

            # flip read if opposite dirs
            if ( $refdir ne $readdir ) {
                ( $ENV{DEBUG} )
                    && warn "\t\tRead TR is on complement strand\n";
                my $tdlen = length($dna);
                $dna = RC($dna);

                #print STDERR "$dna\n";
                my $temp = $readTRStart;
                $readTRStart = ( $tdlen - $readTRStop ) + 1;
                $readTRStop  = ( $tdlen - $temp ) + 1;
            }

            ( $ENV{DEBUG} )
                && warn
                "\t\tSequence: $dna, TR start: $readTRStart, TR stop: $readTRStop\n";

            my $atemp = substr(
                $dna,
                $readTRStart - 1,
                ( ( $readTRStop - $readTRStart ) + 1 )
            );

            push @{ $supported_tr->{alt} }, $atemp;

            # warn "\t\tAlt: ", $supported_tr->{alt}, "\n" if $ENV{DEBUG};

        }

        push @{ $supported_tr->{gt_string} },    $supported_tr->{al}++;
        push @{ $supported_tr->{read_support} }, $support;
        push @{ $supported_tr->{copy_diff} },    $cdiff;

        $supported_tr->{alleles_supported}++;
    }

    # Print last record:
    write_vcf_rec( $spanN_vcffile, $allwithsupport_vcffile, $supported_tr );

    $get_supported_reftrs_sth->finish;

    close $spanN_vcffile;
    close $allwithsupport_vcffile;
    $dbh->disconnect;
}

####################################
sub print_distr {
    my $LARGEST_PSIZE = 122;
    my $LARGEST_ASIZE = 311;

    #my $LARGEST_PSIZE = 5000;
    #my $LARGEST_ASIZE = 5000;

    my $BUCKET_PAT = 1;
    my $BUCKET_AS  = 20;

    my %PHASH1      = ();
    my %PHASH2      = ();
    my %PHASH3      = ();
    my %PHASH4      = ();
    my %PHASH5      = ();
    my %PHASH6      = ();
    my %PHASHB      = ();
    my %PHASH_SN    = ();
    my %PHASH_SN_VN = ();

    open my $distrfh, ">", "${result_prefix}.span${MIN_SUPPORT_REQUIRED}.txt"
        or die
        "\nCan't open for reading ${result_prefix}.span${MIN_SUPPORT_REQUIRED}.txt\n";
    print $distrfh
        "\n\nPatternSize, All, Span1, PercentageS1, Span${MIN_SUPPORT_REQUIRED}, PercentageS${MIN_SUPPORT_REQUIRED}, VntrSpan${MIN_SUPPORT_REQUIRED}, PercentageVS${MIN_SUPPORT_REQUIRED}\n";

    my $dbh = get_dbh( { readonly => 1, userefdb => 1 } );

    # 1 patsize (all)
    my $sth = $dbh->prepare(
        q{SELECT length(pattern) as patsize, count(*)
        FROM refdb.fasta_ref_reps GROUP BY patsize ORDER BY patsize ASC}
    ) or die "Couldn't prepare statement: " . $dbh->errstr;
    $sth->execute() or die "Cannot execute: " . $sth->errstr();
    my $i = 0;
    while ( my @data = $sth->fetchrow_array() ) {
        if ( $data[0] < $LARGEST_PSIZE ) {
            $PHASH1{ $data[0] } = $data[1];
        }
        else {
            $PHASH1{$LARGEST_PSIZE} += $data[1];
        }
        $i++;
    }
    $sth->finish;

    # 2 patsize (span1)
    $sth = $dbh->prepare(
        q{SELECT length(pattern) as patsize, count(*)
        FROM main.fasta_ref_reps mainreftab JOIN refdb.fasta_ref_reps reftab USING (rid)
        WHERE mainreftab.span1>0
            GROUP BY patsize ORDER BY patsize ASC}
    ) or die "Couldn't prepare statement: " . $dbh->errstr;
    $sth->execute() or die "Cannot execute: " . $sth->errstr();
    while ( my @data = $sth->fetchrow_array() ) {
        if ( $data[0] < $LARGEST_PSIZE ) {
            $PHASH2{ $data[0] } = $data[1];
        }
        else {
            $PHASH2{$LARGEST_PSIZE} += $data[1];
        }

        $i++;
    }
    $sth->finish;

    # 3 patsize (spanN)
    $sth = $dbh->prepare(
        q{SELECT length(pattern) as patsize, count(*)
        FROM main.fasta_ref_reps mainreftab JOIN refdb.fasta_ref_reps reftab USING (rid)
        WHERE mainreftab.spanN>0 GROUP BY patsize ORDER BY patsize ASC}
    ) or die "Couldn't prepare statement: " . $dbh->errstr;
    $sth->execute() or die "Cannot execute: " . $sth->errstr();
    while ( my @data = $sth->fetchrow_array() ) {
        if ( $data[0] < $LARGEST_PSIZE ) {
            $PHASH_SN{ $data[0] } = $data[1];
        }
        else {
            $PHASH_SN{$LARGEST_PSIZE} += $data[1];
        }

        $i++;
    }
    $sth->finish;

    # 4 patsize (vntr spanN)
    $sth = $dbh->prepare(
        q{SELECT length(pattern) as patsize, count(*)
        FROM main.fasta_ref_reps mainreftab JOIN refdb.fasta_ref_reps reftab USING (rid)
        WHERE support_vntr>0 GROUP BY patsize ORDER BY patsize ASC}
    ) or die "Couldn't prepare statement: " . $dbh->errstr;
    $sth->execute() or die "Cannot execute: " . $sth->errstr();
    while ( my @data = $sth->fetchrow_array() ) {
        if ( $data[0] < $LARGEST_PSIZE ) {
            $PHASH_SN_VN{ $data[0] } = $data[1];
        }
        else {
            $PHASH_SN_VN{$LARGEST_PSIZE} += $data[1];
        }

        $i++;
    }
    $sth->finish;

    # print
    my $maxval     = 10;
    my @arX0       = ();
    my @arX        = ();
    my @arY1       = ();
    my @arY2       = ();
    my %KEYS12     = ();
    my $tall       = 0;
    my $tspan1     = 0;
    my $tspanN     = 0;
    my $tspanNvntr = 0;

    foreach my $key ( keys %PHASH1 )      { $KEYS12{$key} = 1; }
    foreach my $key ( keys %PHASH2 )      { $KEYS12{$key} = 1; }
    foreach my $key ( keys %PHASH_SN )    { $KEYS12{$key} = 1; }
    foreach my $key ( keys %PHASH_SN_VN ) { $KEYS12{$key} = 1; }
    foreach my $key ( sort { $a <=> $b } ( keys %KEYS12 ) ) {
        my $val1 = 0;
        if ( exists $PHASH1{$key} ) { $val1 = $PHASH1{$key}; }
        my $val2 = 0;
        if ( exists $PHASH2{$key} ) { $val2 = $PHASH2{$key}; }
        my $val3 = 0;
        if ( exists $PHASH_SN{$key} ) { $val3 = $PHASH_SN{$key}; }
        my $val4 = 0;
        if ( exists $PHASH_SN_VN{$key} ) { $val4 = $PHASH_SN_VN{$key}; }

        if ( $key < $LARGEST_PSIZE ) {
            print $distrfh $key . ", "
                . $val1 . ", "
                . $val2 . ", "
                . int( 100 * $val2 / $val1 ) . "\%, "
                . $val3 . ", "
                . int( 100 * $val3 / $val1 ) . "\%, "
                . $val4 . ", "
                . int( 100 * $val4 / $val1 ) . "\%\n";
            push( @arX0, "$key" );
        }
        else {
            print $distrfh $key . "+, "
                . $val1 . ", "
                . $val2 . ", "
                . int( 100 * $val2 / $val1 ) . "\%, "
                . $val3 . ", "
                . int( 100 * $val3 / $val1 ) . "\%, "
                . $val4 . ", "
                . int( 100 * $val4 / $val1 ) . "\%\n";
            push( @arX0, "$key+" );
        }
        push( @arX,  $key );
        push( @arY1, $val1 );
        push( @arY2, $val2 );

        $tall       += $val1;
        $tspan1     += $val2;
        $tspanN     += $val3;
        $tspanNvntr += $val4;

        $maxval = max( $maxval, $val1 );
    }
    warn "TOTAL, $tall, $tspan1, \n";
    print $distrfh "TOTAL, $tall, $tspan1, "
        . int( 100 * $tspan1 / $tall )
        . "\%, $tspanN, "
        . int( 100 * $tspanN / $tall )
        . "\%, $tspanNvntr, "
        . int( 100 * $tspanNvntr / $tall ) . "\%\n";

#my $temp = commify(GetStatistics("NUMBER_REF_TRS"));
#(my $day, my $month, my $year) = (localtime)[3,4,5];
#my $temp2 = sprintf("%02d\/%02d\/%04d", $month+1, $day, $year+1900);
#my $title = "Distribution of References Spanned by Pattern Size (total refs: $temp, $temp2)";

    #my $imwidth = 2000;
    #my $imheight = 2000;
    #my $graph = GD::Graph::linespoints->new($imwidth, $imheight);

    #  $graph->set(
    #      transparent  => 0,
    #      bgclr    => 'white',
    #      x_label           => 'pattern size',
    #      y_label           => 'number of references',
    #      title             => $title,
    #      x_min_value       => 0,
    #      y_max_value       => int($maxval + $maxval*.1),
    #      x_label_skip      => 10
    #  ) or die $graph->error;

#my @legend_keys = ('input reference','references with at least one spanning read');
#$graph->set_legend(@legend_keys);

    # @data = (
    #    [ @arX0 ],
    #    [ @arY1 ],
    #    [ @arY2 ]
    #  );

    #my $gd = $graph->plot(\@data) or die $graph->error;

    #open(IMG, ">${latex}.span${MIN_SUPPORT_REQUIRED}.psize.png") or die $!;
    #binmode IMG;
    #print IMG $gd->png;
    #close IMG;

    $tall = 0;
    my $tclustered = 0;
    my $tmapped    = 0;
    my $tbest      = 0;
    $tspan1 = 0;

    # 3
    print $distrfh
        "\n\nArraySize, All, Clustered, PercentageC, MappedFlanks, PercentageM, BBB, PercentageB, Span1, PercentageS\n";
    $sth = $dbh->prepare(
        q{SELECT (lastindex - firstindex + 1)  as arraysize, count(distinct rid)
        FROM refdb.fasta_ref_reps
        GROUP BY arraysize ORDER BY arraysize ASC}
    ) or die "Couldn't prepare statement: " . $dbh->errstr;
    $sth->execute() or die "Cannot execute: " . $sth->errstr();
    while ( my @data = $sth->fetchrow_array() ) {
        if ( $data[0] < $LARGEST_ASIZE ) {
            $PHASH3{ $data[0] } = $data[1];
        }
        else {
            $PHASH3{$LARGEST_ASIZE} += $data[1];
        }

        $i++;
    }
    $sth->finish;

    # 4
    $sth = $dbh->prepare(
        q{SELECT (lastindex - firstindex + 1)  as arraysize, count(distinct mainreftab.rid)
        FROM main.fasta_ref_reps mainreftab JOIN refdb.fasta_ref_reps reftab USING (rid)
        WHERE mainreftab.span1>0 GROUP BY arraysize ORDER BY arraysize ASC}
    ) or die "Couldn't prepare statement: " . $dbh->errstr;
    $sth->execute() or die "Cannot execute: " . $sth->errstr();
    while ( my @data = $sth->fetchrow_array() ) {
        if ( $data[0] < $LARGEST_ASIZE ) {
            $PHASH4{ $data[0] } = $data[1];
        }
        else {
            $PHASH4{$LARGEST_ASIZE} += $data[1];
        }

        $i++;
    }
    $sth->finish;

    # bbb
    $sth = $dbh->prepare(
        q{SELECT (lastindex - firstindex + 1)  as arraysize, count(distinct reftab.rid)
        FROM refdb.fasta_ref_reps reftab
        WHERE rid IN (SELECT DISTINCT map.refid FROM map WHERE bbb=1)
        GROUP BY arraysize ORDER BY arraysize ASC}
    ) or die "Couldn't prepare statement: " . $dbh->errstr;
    $sth->execute() or die "Cannot execute: " . $sth->errstr();
    while ( my @data = $sth->fetchrow_array() ) {
        if ( $data[0] < $LARGEST_ASIZE ) {
            $PHASHB{ $data[0] } = $data[1];
        }
        else {
            $PHASHB{$LARGEST_ASIZE} += $data[1];
        }

        $i++;
    }
    $sth->finish;

    # percent clustered
    $sth = $dbh->prepare(
        q{SELECT (lastindex - firstindex + 1)  as arraysize, count(distinct reftab.rid)
        FROM refdb.fasta_ref_reps reftab INNER JOIN clusterlnk ON reftab.rid=-clusterlnk.repeatid
        GROUP BY arraysize ORDER BY arraysize ASC}
    ) or die "Couldn't prepare statement: " . $dbh->errstr;
    $sth->execute() or die "Cannot execute: " . $sth->errstr();
    while ( my @data = $sth->fetchrow_array() ) {
        if ( $data[0] < $LARGEST_ASIZE ) {
            $PHASH5{ $data[0] } = $data[1];
        }
        else {
            $PHASH5{$LARGEST_ASIZE} += $data[1];
        }

        $i++;
    }
    $sth->finish;

    # percent mapped
    $sth = $dbh->prepare(
        q{SELECT (lastindex - firstindex + 1)  as arraysize, count(distinct reftab.rid)
        FROM refdb.fasta_ref_reps reftab
        WHERE rid IN (SELECT DISTINCT map.refid FROM map)
        GROUP BY arraysize ORDER BY arraysize ASC}
    ) or die "Couldn't prepare statement: " . $dbh->errstr;
    $sth->execute() or die "Cannot execute: " . $sth->errstr();
    while ( my @data = $sth->fetchrow_array() ) {
        if ( $data[0] < $LARGEST_ASIZE ) {
            $PHASH6{ $data[0] } = $data[1];
        }
        else {
            $PHASH6{$LARGEST_ASIZE} += $data[1];
        }

        $i++;
    }
    $sth->finish;

    # print
    $maxval = 10;
    @arX0   = ();
    @arX    = ();
    @arY1   = ();
    @arY2   = ();
    foreach my $key ( sort { $a <=> $b } ( keys %PHASH3 ) ) {
        my $val1 = 0;
        if ( exists $PHASH3{$key} ) { $val1 = $PHASH3{$key}; }
        my $val2 = 0;
        if ( exists $PHASH4{$key} ) { $val2 = $PHASH4{$key}; }
        my $val4 = 0;
        if ( exists $PHASH5{$key} ) { $val4 = $PHASH5{$key}; }
        my $val3 = 0;
        if ( exists $PHASH6{$key} ) { $val3 = $PHASH6{$key}; }
        my $valb = 0;
        if ( exists $PHASHB{$key} ) { $valb = $PHASHB{$key}; }

        if ( $key < $LARGEST_ASIZE ) {
            print $distrfh $key . ", "
                . $val1 . ", "
                . $val4 . ", "
                . int( 100 * $val4 / $val1 ) . "%, "
                . $val3 . ", "
                . int( 100 * $val3 / $val1 ) . "\%, "
                . $valb . ", "
                . int( 100 * $valb / $val1 ) . "%, "
                . $val2 . ", "
                . int( 100 * $val2 / $val1 ) . "\%\n";
            push( @arX0, "$key" );
        }
        else {
            print $distrfh $key . "+, "
                . $val1 . ", "
                . $val4 . ", "
                . int( 100 * $val4 / $val1 ) . "%, "
                . $val3 . ", "
                . int( 100 * $val3 / $val1 ) . "\%, "
                . $valb . ", "
                . int( 100 * $valb / $val1 ) . "%, "
                . $val2 . ", "
                . int( 100 * $val2 / $val1 ) . "\%\n";
            push( @arX0, "$key+" );
        }
        push( @arX,  $key );
        push( @arY1, $val1 );
        push( @arY2, $val2 );

        $tall       += $val1;
        $tclustered += $val4;
        $tmapped    += $val3;
        $tbest      += $valb;
        $tspan1     += $val2;

        $maxval = max( $maxval, $val1 );
    }

    print $distrfh "TOTAL, $tall, $tclustered, "
        . int( 100 * $tclustered / $tall )
        . "\%, $tmapped, "
        . int( 100 * $tmapped / $tall )
        . "\%, $tbest, "
        . int( 100 * $tbest / $tall )
        . "\%, $tspan1, "
        . int( 100 * $tspan1 / $tall ) . "\%\n";

    # copies gained/lost
    my $total = 0;
    print $distrfh
        "\n\n(vntr support>=$MIN_SUPPORT_REQUIRED) Copies Gained, Frequency\n";
    $sth = $dbh->prepare(
        q{SELECT copies, COUNT(*)
        FROM vntr_support
        WHERE copies!=0 AND support>=$MIN_SUPPORT_REQUIRED GROUP BY copies ORDER BY copies ASC}
    ) or die "Couldn't prepare statement: " . $dbh->errstr;
    $sth->execute() or die "Cannot execute: " . $sth->errstr();

    while ( my @data = $sth->fetchrow_array() ) {
        print $distrfh $data[0] . "," . $data[1] . "\n";
        $i++;
        $total += $data[1];
    }
    $sth->finish;
    print $distrfh "TOTAL: $total\n";

    # copies by patsize
    $total = 0;
    print $distrfh
        "\n\n(vntr support>=$MIN_SUPPORT_REQUIRED) PatternSize, Copies Gained, Frequency\n";
    $sth = $dbh->prepare(
        q{SELECT length(pattern) AS patsize, copies, count(*)
        FROM vntr_support INNER JOIN refdb.fasta_ref_reps reftab ON reftab.rid = -vntr_support.refid
        WHERE copies!=0 AND support>=$MIN_SUPPORT_REQUIRED
        GROUP BY patsize, copies
        ORDER BY patsize ASC, copies ASC}
    ) or die "Couldn't prepare statement: " . $dbh->errstr;
    $sth->execute() or die "Cannot execute: " . $sth->errstr();
    while ( my @data = $sth->fetchrow_array() ) {
        print $distrfh $data[0] . "," . $data[1] . "," . $data[2] . "\n";
        $i++;
        $total += $data[2];
    }
    $sth->finish;
    print $distrfh "TOTAL: $total\n";

    # copies by array size
    $total = 0;
    print $distrfh
        "\n\n(vntr support>=$MIN_SUPPORT_REQUIRED) ArraySize, Copies Gained, Frequency\n";
    $sth = $dbh->prepare(
        q{SELECT (lastindex-firstindex+1) AS arraysize, copies, count(*)
        FROM vntr_support INNER JOIN refdb.fasta_ref_reps reftab ON reftab.rid = -vntr_support.refid
        WHERE copies!=0 AND support>=$MIN_SUPPORT_REQUIRED
        GROUP BY arraysize, copies ORDER BY arraysize ASC, copies ASC}
    ) or die "Couldn't prepare statement: " . $dbh->errstr;
    $sth->execute() or die "Cannot execute: " . $sth->errstr();

    while ( my @data = $sth->fetchrow_array() ) {
        print $distrfh $data[0] . "," . $data[1] . "," . $data[2] . "\n";
        $i++;
        $total += $data[2];
    }
    $sth->finish;
    print $distrfh "TOTAL: $total\n";

    # spanning reads

    my %SpanningSing   = ();
    my %SpanningDist   = ();
    my %SpanningIndist = ();
    my %SpanningBBB    = ();

    $sth = $dbh->prepare(
        q{SELECT refid,sum(support)
        FROM vntr_support INNER JOIN main.fasta_ref_reps mainreftab ON reftab.rid = -vntr_support.refid
            INNER JOIN refdb.fasta_ref_reps reftab USING (rid)
        WHERE is_singleton=1 GROUP BY refid}
    ) or die "Couldn't prepare statement: " . $dbh->errstr;
    $sth->execute() or die "Cannot execute: " . $sth->errstr();
    while ( my @data = $sth->fetchrow_array() ) {
        $SpanningSing{ $data[1] }++;
        $i++;
    }
    $sth->finish;

    $sth = $dbh->prepare(
        q{SELECT refid,sum(support)
        FROM vntr_support INNER JOIN main.fasta_ref_reps mainreftab ON reftab.rid = -vntr_support.refid
            INNER JOIN refdb.fasta_ref_reps reftab USING (rid)
        WHERE is_dist=1 GROUP BY refid}
    ) or die "Couldn't prepare statement: " . $dbh->errstr;
    $sth->execute() or die "Cannot execute: " . $sth->errstr();
    my %SpanninDist;

    while ( my @data = $sth->fetchrow_array() ) {
        $SpanninDist{ $data[1] }++;
        $i++;
    }
    $sth->finish;

    $sth = $dbh->prepare(
        q{SELECT refid,sum(support)
        FROM vntr_support INNER JOIN main.fasta_ref_reps mainreftab ON reftab.rid = -vntr_support.refid
            INNER JOIN refdb.fasta_ref_reps reftab USING (rid)
        WHERE is_indist=1 GROUP BY refid}
    ) or die "Couldn't prepare statement: " . $dbh->errstr;
    $sth->execute() or die "Cannot execute: " . $sth->errstr();
    while ( my @data = $sth->fetchrow_array() ) {
        $SpanningIndist{ $data[1] }++;
        $i++;
    }
    $sth->finish;

    $sth = $dbh->prepare(
        q{SELECT refid,sum(support)
        FROM vntr_support INNER JOIN refdb.fasta_ref_reps reftab ON reftab.rid = -vntr_support.refid
        GROUP BY refid}
    ) or die "Couldn't prepare statement: " . $dbh->errstr;
    $sth->execute() or die "Cannot execute: " . $sth->errstr();
    while ( my @data = $sth->fetchrow_array() ) {
        $SpanningBBB{ $data[1] }++;
        $i++;
    }
    $sth->finish;

    # SPANNING
    print $distrfh "\n\nSpanning reads per locus\n\nRefClass     ";
    my $maxkey = 0;
    foreach my $key ( sort { $a <=> $b } ( keys %SpanningBBB ) ) {
        $maxkey = $key;
    }

    for ( my $key = 0; $key <= $maxkey; $key++ ) {
        print $distrfh "\t$key";
    }
    print $distrfh "\tTOTAL";

    $total = 0;
    print $distrfh "\nSINGLETON";
    for ( my $key = 0; $key <= $maxkey; $key++ ) {
        my $val1 = 0;
        if ( exists $SpanningSing{$key} ) { $val1 = $SpanningSing{$key}; }
        print $distrfh "\t$val1";
        $total += $val1;
    }
    print $distrfh "\t$total";

    $total = 0;
    print $distrfh "\nDISTING   ";
    for ( my $key = 0; $key <= $maxkey; $key++ ) {
        my $val1 = 0;
        if ( exists $SpanninDist{$key} ) { $val1 = $SpanninDist{$key}; }
        print $distrfh "\t$val1";
        $total += $val1;
    }
    print $distrfh "\t$total";

    $total = 0;
    print $distrfh "\nINDIST      ";
    for ( my $key = 0; $key <= $maxkey; $key++ ) {
        my $val1 = 0;
        if ( exists $SpanningIndist{$key} ) { $val1 = $SpanningIndist{$key}; }
        print $distrfh "\t$val1";
        $total += $val1;
    }
    print $distrfh "\t$total";

    $total = 0;
    print $distrfh "\nBBB         ";
    for ( my $key = 0; $key <= $maxkey; $key++ ) {
        my $val1 = 0;
        if ( exists $SpanningBBB{$key} ) { $val1 = $SpanningBBB{$key}; }
        print $distrfh "\t$val1";
        $total += $val1;
    }
    print $distrfh "\t$total";

    # alleles spanning reads

    my %AllelesSing   = ();
    my %AllelesDist   = ();
    my %AllelesIndist = ();
    my %AllelesBBB    = ();

    $sth = $dbh->prepare(
        q{SELECT refid,copies,support
        FROM vntr_support INNER JOIN main.fasta_ref_reps mainreftab ON reftab.rid = -vntr_support.refid
            INNER JOIN refdb.fasta_ref_reps reftab USING (rid)
        WHERE is_singleton=1}
    ) or die "Couldn't prepare statement: " . $dbh->errstr;
    $sth->execute() or die "Cannot execute: " . $sth->errstr();
    while ( my @data = $sth->fetchrow_array() ) {
        $AllelesSing{ $data[2] }++;
        $i++;
    }
    $sth->finish;

    $sth = $dbh->prepare(
        q{SELECT refid,copies,support
        FROM vntr_support INNER JOIN main.fasta_ref_reps mainreftab ON reftab.rid = -vntr_support.refid
            INNER JOIN refdb.fasta_ref_reps reftab USING (rid)
        WHERE is_dist=1}
    ) or die "Couldn't prepare statement: " . $dbh->errstr;
    $sth->execute() or die "Cannot execute: " . $sth->errstr();
    while ( my @data = $sth->fetchrow_array() ) {
        $AllelesDist{ $data[2] }++;
        $i++;
    }
    $sth->finish;

    $sth = $dbh->prepare(
        q{SELECT refid,copies,support
        FROM vntr_support INNER JOIN main.fasta_ref_reps mainreftab ON reftab.rid = -vntr_support.refid
            INNER JOIN refdb.fasta_ref_reps reftab USING (rid)
        WHERE is_indist=1}
    ) or die "Couldn't prepare statement: " . $dbh->errstr;
    $sth->execute() or die "Cannot execute: " . $sth->errstr();
    while ( my @data = $sth->fetchrow_array() ) {
        $AllelesIndist{ $data[2] }++;
        $i++;
    }
    $sth->finish;

    $sth = $dbh->prepare(
        q{SELECT refid,copies,support
        FROM vntr_support INNER JOIN main.fasta_ref_reps mainreftab ON mainreftab.rid = -vntr_support.refid}
    ) or die "Couldn't prepare statement: " . $dbh->errstr;
    $sth->execute() or die "Cannot execute: " . $sth->errstr();
    while ( my @data = $sth->fetchrow_array() ) {
        $AllelesBBB{ $data[2] }++;
        $i++;
    }
    $sth->finish;

    # ALLELES
    print $distrfh "\n\nSpanning reads per allele\n\nRefClass     ";
    $maxkey = 0;
    foreach my $key ( sort { $a <=> $b } ( keys %AllelesBBB ) ) {
        $maxkey = $key;
    }

    for ( my $key = 0; $key <= $maxkey; $key++ ) {
        print $distrfh "\t$key";
    }
    print $distrfh "\tTOTAL";

    $total = 0;
    print $distrfh "\nSINGLETON";
    for ( my $key = 0; $key <= $maxkey; $key++ ) {
        my $val1 = 0;
        if ( exists $AllelesSing{$key} ) { $val1 = $AllelesSing{$key}; }
        print $distrfh "\t$val1";
        $total += $val1;
    }
    print $distrfh "\t$total";

    $total = 0;
    print $distrfh "\nDISTING   ";
    for ( my $key = 0; $key <= $maxkey; $key++ ) {
        my $val1 = 0;
        if ( exists $AllelesDist{$key} ) { $val1 = $AllelesDist{$key}; }
        print $distrfh "\t$val1";
        $total += $val1;
    }
    print $distrfh "\t$total";

    $total = 0;
    print $distrfh "\nINDIST      ";
    for ( my $key = 0; $key <= $maxkey; $key++ ) {
        my $val1 = 0;
        if ( exists $AllelesIndist{$key} ) { $val1 = $AllelesIndist{$key}; }
        print $distrfh "\t$val1";
        $total += $val1;
    }
    print $distrfh "\t$total";

    $total = 0;
    print $distrfh "\nBBB         ";
    for ( my $key = 0; $key <= $maxkey; $key++ ) {
        my $val1 = 0;
        if ( exists $AllelesBBB{$key} ) { $val1 = $AllelesBBB{$key}; }
        print $distrfh "\t$val1";
        $total += $val1;
    }
    print $distrfh "\t$total";

    # alleles spanning reads for VNTRs

    %AllelesSing   = ();
    %AllelesDist   = ();
    %AllelesIndist = ();
    %AllelesBBB    = ();

    $sth = $dbh->prepare(
        q{SELECT refid,copies,support
        FROM vntr_support INNER JOIN main.fasta_ref_reps mainreftab ON mainreftab.rid = -vntr_support.refid
            INNER JOIN refdb.fasta_ref_reps reftab USING(rid)
        WHERE  support_vntr>0 AND is_singleton=1}
    ) or die "Couldn't prepare statement: " . $dbh->errstr;
    $sth->execute() or die "Cannot execute: " . $sth->errstr();
    while ( my @data = $sth->fetchrow_array() ) {
        $AllelesSing{ $data[2] }++;
        $i++;
    }
    $sth->finish;

    $sth = $dbh->prepare(
        q{SELECT refid,copies,support
        FROM vntr_support INNER JOIN main.fasta_ref_reps mainreftab ON mainreftab.rid = -vntr_support.refid
            INNER JOIN refdb.fasta_ref_reps reftab USING(rid)
        WHERE support_vntr>0 AND is_dist=1}
    ) or die "Couldn't prepare statement: " . $dbh->errstr;
    $sth->execute() or die "Cannot execute: " . $sth->errstr();
    while ( my @data = $sth->fetchrow_array() ) {
        $AllelesDist{ $data[2] }++;
        $i++;
    }
    $sth->finish;

    $sth = $dbh->prepare(
        q{SELECT refid,copies,support
        FROM vntr_support INNER JOIN main.fasta_ref_reps mainreftab ON mainreftab.rid = -vntr_support.refid
            INNER JOIN refdb.fasta_ref_reps reftab USING(rid)
        WHERE support_vntr>0 AND is_indist=1}
    ) or die "Couldn't prepare statement: " . $dbh->errstr;
    $sth->execute() or die "Cannot execute: " . $sth->errstr();
    while ( my @data = $sth->fetchrow_array() ) {
        $AllelesIndist{ $data[2] }++;
        $i++;
    }
    $sth->finish;

    $sth = $dbh->prepare(
        q{SELECT refid,copies,support
        FROM vntr_support INNER JOIN main.fasta_ref_reps mainreftab ON mainreftab.rid = -vntr_support.refid
        WHERE support_vntr>0}
    ) or die "Couldn't prepare statement: " . $dbh->errstr;
    $sth->execute() or die "Cannot execute: " . $sth->errstr();
    while ( my @data = $sth->fetchrow_array() ) {
        $AllelesBBB{ $data[2] }++;
        $i++;
    }
    $sth->finish;

    # ALLELES for VNTRs
    print $distrfh "\n\nSpanning reads per allele for VNTRs\n\nRefClass     ";
    $maxkey = 0;
    foreach my $key ( sort { $a <=> $b } ( keys %AllelesBBB ) ) {
        $maxkey = $key;
    }

    for ( my $key = 0; $key <= $maxkey; $key++ ) {
        print $distrfh "\t$key";
    }
    print $distrfh "\tTOTAL";

    $total = 0;
    print $distrfh "\nSINGLETON";
    for ( my $key = 0; $key <= $maxkey; $key++ ) {
        my $val1 = 0;
        if ( exists $AllelesSing{$key} ) { $val1 = $AllelesSing{$key}; }
        print $distrfh "\t$val1";
        $total += $val1;
    }
    print $distrfh "\t$total";

    $total = 0;
    print $distrfh "\nDISTING   ";
    for ( my $key = 0; $key <= $maxkey; $key++ ) {
        my $val1 = 0;
        if ( exists $AllelesDist{$key} ) { $val1 = $AllelesDist{$key}; }
        print $distrfh "\t$val1";
        $total += $val1;
    }
    print $distrfh "\t$total";

    $total = 0;
    print $distrfh "\nINDIST      ";
    for ( my $key = 0; $key <= $maxkey; $key++ ) {
        my $val1 = 0;
        if ( exists $AllelesIndist{$key} ) { $val1 = $AllelesIndist{$key}; }
        print $distrfh "\t$val1";
        $total += $val1;
    }
    print $distrfh "\t$total";

    $total = 0;
    print $distrfh "\nBBB         ";
    for ( my $key = 0; $key <= $maxkey; $key++ ) {
        my $val1 = 0;
        if ( exists $AllelesBBB{$key} ) { $val1 = $AllelesBBB{$key}; }
        print $distrfh "\t$val1";
        $total += $val1;
    }
    print $distrfh "\t$total\n";

    # alleles spanning reads for VNTRs by VNTR class

    my %AllelesBBBInf   = ();
    my %AllelesBBBObs   = ();
    my %AllelesBBBTotal = ();

    $sth = $dbh->prepare(
        q{SELECT refid,copies,support
        FROM vntr_support INNER JOIN main.fasta_ref_reps mainreftab ON mainreftab.rid = -vntr_support.refid
        WHERE support_vntr>0 AND mainreftab.homez_diff=1}
    ) or die "Couldn't prepare statement: " . $dbh->errstr;
    $sth->execute() or die "Cannot execute: " . $sth->errstr();
    while ( my @data = $sth->fetchrow_array() ) {
        $AllelesBBBInf{ $data[2] }++;
        $i++;
    }
    $sth->finish;

    $sth = $dbh->prepare(
        q{SELECT refid,copies,support
        FROM vntr_support INNER JOIN main.fasta_ref_reps mainreftab ON mainreftab.rid = -vntr_support.refid
        WHERE support_vntr>0 AND (mainreftab.hetez_same=1 OR mainreftab.hetez_diff=1 OR mainreftab.hetez_multi=1)}
    ) or die "Couldn't prepare statement: " . $dbh->errstr;
    $sth->execute() or die "Cannot execute: " . $sth->errstr();
    while ( my @data = $sth->fetchrow_array() ) {
        $AllelesBBBObs{ $data[2] }++;
        $i++;
    }
    $sth->finish;

    print $distrfh
        "\n\nSpanning reads per allele for VNTRs by VNTR class\n\nVNTRClass     ";
    $maxkey = 0;
    foreach my $key ( sort { $a <=> $b } ( keys %AllelesBBBInf ) ) {
        $maxkey = $key;
        $AllelesBBBTotal{$key} = $AllelesBBBInf{$key};
    }
    foreach my $key ( sort { $a <=> $b } ( keys %AllelesBBBObs ) ) {
        $maxkey = $key;
        $AllelesBBBTotal{$key} += $AllelesBBBObs{$key};
    }

    for ( my $key = 0; $key <= $maxkey; $key++ ) {
        print $distrfh "\t$key";
    }
    print $distrfh "\tTOTAL";

    $total = 0;
    print $distrfh "\nInferred";
    for ( my $key = 0; $key <= $maxkey; $key++ ) {
        my $val1 = 0;
        if ( exists $AllelesBBBInf{$key} ) { $val1 = $AllelesBBBInf{$key}; }
        print $distrfh "\t$val1";
        $total += $val1;
    }
    print $distrfh "\t$total";

    $total = 0;
    print $distrfh "\nObserved";
    for ( my $key = 0; $key <= $maxkey; $key++ ) {
        my $val1 = 0;
        if ( exists $AllelesBBBObs{$key} ) { $val1 = $AllelesBBBObs{$key}; }
        print $distrfh "\t$val1";
        $total += $val1;
    }
    print $distrfh "\t$total";

    $total = 0;
    print $distrfh "\nTotal";
    for ( my $key = 0; $key <= $maxkey; $key++ ) {
        my $val1 = 0;
        if ( exists $AllelesBBBTotal{$key} ) {
            $val1 = $AllelesBBBTotal{$key};
        }
        print $distrfh "\t$val1";
        $total += $val1;
    }
    print $distrfh "\t$total\n";

    # counts for support 1 alleles
    print $distrfh "\n1A: ";
    $sth = $dbh->prepare(
        q{SELECT count(distinct refid)
        FROM vntr_support INNER JOIN main.fasta_ref_reps mainreftab ON mainreftab.rid = -vntr_support.refid
        WHERE copies!=0 AND mainreftab.has_support=0 AND support_vntr=0 AND support=1}
    );
    $sth->execute() or die "Cannot execute: " . $sth->errstr();
    if ( my @data = $sth->fetchrow_array() ) {
        print $distrfh $data[0];
    }
    $sth->finish;
    print $distrfh "\n1B: ";
    $sth = $dbh->prepare(
        q{SELECT count(distinct refid)
        FROM vntr_support INNER JOIN main.fasta_ref_reps mainreftab ON mainreftab.rid = -vntr_support.refid
        WHERE copies!=0 AND has_support=1 AND support_vntr=0 AND support=1}
    );
    $sth->execute() or die "Cannot execute: " . $sth->errstr();
    if ( my @data = $sth->fetchrow_array() ) {
        print $distrfh $data[0];
    }
    $sth->finish;
    print $distrfh "\n2A: ";
    $sth = $dbh->prepare(
        q{SELECT count(distinct refid)
        FROM vntr_support INNER JOIN main.fasta_ref_reps mainreftab ON mainreftab.rid = -vntr_support.refid
        WHERE homez_diff=1 and support=1}
    );
    $sth->execute() or die "Cannot execute: " . $sth->errstr();
    if ( my @data = $sth->fetchrow_array() ) {
        print $distrfh $data[0];
    }
    $sth->finish;
    print $distrfh "\n2B: ";
    $sth = $dbh->prepare(
        q{SELECT count(distinct refid)
        FROM vntr_support INNER JOIN main.fasta_ref_reps mainreftab ON mainreftab.rid = -vntr_support.refid
        WHERE (hetez_same=1 OR hetez_diff=1 OR hetez_multi=1) AND support=1}
    );
    $sth->execute() or die "Cannot execute: " . $sth->errstr();
    if ( my @data = $sth->fetchrow_array() ) {
        print $distrfh $data[0];
    }
    $sth->finish;
    $dbh->disconnect;
    close($distrfh);

    # images

    #$imwidth = 2000;
    #$imheight = 2000;
    #$graph = GD::Graph::linespoints->new($imwidth, $imheight);

#$title = "Distribution of References Spanned by Array Size (total refs: $temp, $temp2)";

    #  $graph->set(
    #      x_label           => 'array size',
    #      y_label           => 'number of references',
    #      title             => $title,
    #      x_min_value       => 0,
    #      y_max_value       => int($maxval + $maxval*.1),
    #      x_label_skip      => 10
    #  ) or die $graph->error;

#@legend_keys = ('input reference','references with at least one spanning read');
#$graph->set_legend(@legend_keys);

    # @data = (
    #    [ @arX0 ],
    #    [ @arY1 ],
    #    [ @arY2 ]
    #  );

    #$gd = $graph->plot(\@data) or die $graph->error;

    #open(IMG, ">${latex}.span${MIN_SUPPORT_REQUIRED}.asize.png") or die $!;
    #binmode IMG;
    #print IMG $gd->png;
    #close IMG;

}

###################
sub print_latex {
    my $ReadTRsSupport = shift;

    # warn "Read TRs supported: $ReadTRsSupport\n"
    #     if ($ENV{DEBUG});

    my $sum_has_support  = 0;
    my $sum_span1        = 0;
    my $sum_spanN        = 0;
    my $sum_homez_same   = 0;
    my $sum_homez_diff   = 0;
    my $sum_hetez_same   = 0;
    my $sum_hetez_diff   = 0;
    my $sum_hetez_multi  = 0;
    my $sum_support_vntr = 0;

    # Get needed stats
    my @stats = qw(NUMBER_READS
        NUMBER_TRS_IN_READS
        NUMBER_TRS_IN_READS_GE7
        NUMBER_READS_WITHTRS_GE7
        NUMBER_READS_WITHTRS
        NUMBER_READS_WITHTRS_GE7_AFTER_REDUND
        CLUST_NUMBER_OF_READ_REPS_IN_CLUSTERS
        NUMBER_REF_TRS
        NUMBER_REFS_TRS_AFTER_REDUND
        CLUST_NUMBER_OF_REF_REPS_IN_CLUSTERS
        CLUST_NUMBER_OF_REFS_WITH_PREDICTED_VNTR);
    my $stat_hash = get_statistics(@stats);

    my $dbh = get_dbh( { readonly => 1, userefdb => 1 } );
    my $sth = $dbh->prepare(
        q{SELECT sum(has_support),sum(span1),sum(spanN),sum(homez_same),sum(homez_diff),
            sum(hetez_same),sum(hetez_diff),sum(hetez_multi),sum(support_vntr)
            FROM main.fasta_ref_reps}
    ) or die "Couldn't prepare statement: " . $dbh->errstr;
    $sth->execute() or die "Cannot execute: " . $sth->errstr();
    if ( my @data = $sth->fetchrow_array() ) {

        $sum_has_support  = $data[0];
        $sum_span1        = $data[1];
        $sum_spanN        = $data[2];
        $sum_homez_same   = $data[3];
        $sum_homez_diff   = $data[4];
        $sum_hetez_same   = $data[5];
        $sum_hetez_diff   = $data[6];
        $sum_hetez_multi  = $data[7];
        $sum_support_vntr = $data[8];
    }
    $sth->finish();

# Definitions
# READS -- Reads from the subject data set
# READ-TRs -- read TRs from the subject data set
# REF-TRs -- reference TRs from the hg19 data set
# INITIAL -- total started with in the original data set
# GE7 -- only patterns >= 7
# RDE (Redundancy Elimination) -- left after redundancy elimination
# CRDE (Cyclic Redundancy Elimination) -- left after cyclic redundancy elimination
# PC (Profile Clustered) -- PROCLU clustered into a cluster with at least one reference and at least one read
# MAP (Mapped) -- flank TTT clustered into a cluster with one reference and at least one read
# TIE-OK -- applies to a read which is MAP, but may appear in more than one flank TTT cluster (due to ties or no flanking sequence)
# ADDBACK -- includes those eliminated by CRDE, but added back for TTT flank alignment
# SINGLETON -- PC into a cluster with EXACTLY ONE REFERENCE and MAP
# DISTINGUISHABLE -- PC into a cluster with MORE THAN ONE REFERENCE but flank distinguishable and MAP
# INDISTINGUISHABLE -- PC into a cluster with MORE THAN ONE REFERENCE, NOT flank indistinguishable and MAP
# SPAN1 -- MAP with at least 1 read that is NOT assembly required
# SPANN -- MAP with at least N reads that are NOT assembly required
# SUPPORT -- MAP with support equal to 2 or more for at least one copy number
# HOMOZYGOUS -- MAP only one copy number has SUPPORT
# HETEROZYGOUS -- MAP two or more copy numbers have SUPPORT
# MULTI -- MAP three or more copy numbers have support
# SAME -- one copy number with SUPPORT matches reference copy number
# DIFF -- NO copy number with SUPPORT matches reference copy number

    # ***************************************************
    # Yevgeniy, assign real numbers to these variables.

    my $totalReads = $stat_hash->{NUMBER_READS};    # INITIAL READS

#$totalReadsWithTRs = $data[0];       # INITIAL READS which contain TRs
#$totalReadsWithTRsPatternGE7  = $data[0];        # INITIAL READS which contain TRs GE7
#$readTRsWithPatternGE7  = $data[0];       # INITIAL READ-TRs GE7
#$readTRsWPGE7AfterCyclicRedundancyElimination = $data[0];       # INITIAL READ-TRs GE7 CRDE

    my $totalReadTRs = $stat_hash->{NUMBER_TRS_IN_READS};

    my $readTRsWithPatternGE7       = $stat_hash->{NUMBER_TRS_IN_READS_GE7};
    my $totalReadsWithTRsPatternGE7 = $stat_hash->{NUMBER_READS_WITHTRS_GE7};
    my $totalReadsWithTRs           = $stat_hash->{NUMBER_READS_WITHTRS};
    my $readTRsWPGE7AfterCyclicRedundancyElimination
        = $stat_hash->{NUMBER_READS_WITHTRS_GE7_AFTER_REDUND};

    my $readTRsProfileClustered
        = $stat_hash->{CLUST_NUMBER_OF_READ_REPS_IN_CLUSTERS}
        ;    # INITIAL READ-TRs GE7 PC ADDBACK

    $sth
        = $dbh->prepare(
        "SELECT count(distinct map.readid) FROM map WHERE bbb=1;")
        or die "Couldn't prepare statement: " . $dbh->errstr;
    my $readTRsMapped = 0;
    $sth->execute() or die "Cannot execute: " . $sth->errstr();
    {
        my @data = $sth->fetchrow_array();
        $readTRsMapped
            = $data[0];    # INITIAL READ-TRs RDE GE7 PC ADDBACK MAP TIE-OK
    }
    $sth->finish();

    my $refTRsAREWithPatternGE7
        = $stat_hash->{NUMBER_REF_TRS};    # INITIAL REF-TRs RDE GE7
    my $refTRsAREWPGE7AfterCyclicRedundancyElimination
        = $stat_hash->{NUMBER_REFS_TRS_AFTER_REDUND}
        ;                                  # INITIAL REF-TRs RDE GE7 CRDE
    my $refTRsProfileClustered
        = $stat_hash->{CLUST_NUMBER_OF_REF_REPS_IN_CLUSTERS}
        ;    # INITIAL REF-TRs RDE GE7 PC ADDBACK

    $sth
        = $dbh->prepare(
        "SELECT count(distinct map.refid) FROM map WHERE bbb=1;")
        or die "Couldn't prepare statement: " . $dbh->errstr;
    my $refTRsMapped = 0;
    $sth->execute() or die "Cannot execute: " . $sth->errstr();
    {
        my @data = $sth->fetchrow_array();
        $refTRsMapped = $data[0];    # INITIAL REF-TRs RDE GE7 PC ADDBACK MAP
    }
    $sth->finish();

    my $refTRsMappedSpan1
        = $sum_span1;    # INITIAL REF-TRs RDE GE7 PC ADDBACK MAP SPAN1
    my $refTRsMappedSpanN
        = $sum_spanN;    # INITIAL REF-TRs RDE GE7 PC ADDBACK MAP SPANN
     #my $refTRsMappedSingleton = GetStatistics("NUMBER_REFS_SINGLE_REF_CLUSTER_WITH_READS_MAPPED");  # INITIAL REF-TRs RDE GE7 PC ADDBACK MAP SINGLETON
    $sth
        = $dbh->prepare(
        "SELECT count(distinct rid) FROM refdb.fasta_ref_reps reftab INNER JOIN map ON reftab.rid=map.refid WHERE bbb=1 AND is_singleton=1;"
        ) or die "Couldn't prepare statement: " . $dbh->errstr;
    $sth->execute() or die "Cannot execute: " . $sth->errstr();

    my $refTRsMappedSingleton;
    if ( my @data = $sth->fetchrow_array() ) {
        $refTRsMappedSingleton = $data[0];
    }
    $sth->finish();

    my $refTRsMappedDistinguishable
        = -999;    # INITIAL REF-TRs RDE GE7 PC ADDBACK MAP DISTINGUISHABLE
    my $refTRsMappedIndistinguishable
        = -999;    # INITIAL REF-TRs RDE GE7 PC ADDBACK MAP INDISTINGUISHABLE

 # this would be 0 now that we use psearch to get singleton/indist information
    $sth
        = $dbh->prepare(
        "SELECT count(distinct rid) FROM refdb.fasta_ref_reps reftab INNER JOIN map ON reftab.rid=map.refid WHERE bbb=1 AND is_dist=1 and is_singleton=0;"
        ) or die "Couldn't prepare statement: " . $dbh->errstr;
    $sth->execute() or die "Cannot execute: " . $sth->errstr();
    if ( my @data = $sth->fetchrow_array() ) {
        $refTRsMappedDistinguishable = $data[0];
    }
    $sth->finish();

    $sth
        = $dbh->prepare(
        "SELECT count(distinct rid) FROM refdb.fasta_ref_reps reftab INNER JOIN map ON reftab.rid=map.refid WHERE bbb=1 AND is_indist=1;"
        )

        or die "Couldn't prepare statement: " . $dbh->errstr;
    $sth->execute() or die "Cannot execute: " . $sth->errstr();
    if ( my @data = $sth->fetchrow_array() ) {
        $refTRsMappedIndistinguishable = $data[0];
    }
    $sth->finish();

    my $readTRsMappedToSingleton
        = -777;    # INITIAL READ-TRs RDE GE7 PC ADDBACK MAP SINGLETON
    my $readTRsMappedToDistinguishable = -777
        ;    # INITIAL READ-TRs RDE GE7 PC ADDBACK MAP TIE-OK DISTINGUISHABLE
    my $readTRsMappedToIndistinguishable = -777
        ;   # INITIAL READ-TRs RDE GE7 PC ADDBACK MAP TIE-OK INDISTINGUISHABLE

    $sth
        = $dbh->prepare(
        "SELECT count(distinct map.readid) FROM refdb.fasta_ref_reps reftab INNER JOIN map ON reftab.rid=map.refid WHERE bbb=1 AND is_singleton=1;"
        ) or die "Couldn't prepare statement: " . $dbh->errstr;
    $sth->execute() or die "Cannot execute: " . $sth->errstr();
    if ( my @data = $sth->fetchrow_array() ) {
        $readTRsMappedToSingleton = $data[0];
    }
    $sth->finish();

    $sth
        = $dbh->prepare(
        "SELECT count(distinct map.readid) FROM refdb.fasta_ref_reps reftab INNER JOIN map ON reftab.rid=map.refid WHERE bbb=1 AND is_dist=1 AND is_singleton=0;"
        ) or die "Couldn't prepare statement: " . $dbh->errstr;
    $sth->execute() or die "Cannot execute: " . $sth->errstr();
    if ( my @data = $sth->fetchrow_array() ) {
        $readTRsMappedToDistinguishable = $data[0];
    }
    $sth->finish();

    $sth
        = $dbh->prepare(
        "SELECT count(distinct map.readid) FROM refdb.fasta_ref_reps reftab INNER JOIN map ON reftab.rid=map.refid WHERE bbb=1 AND is_indist=1;"
        )

        or die "Couldn't prepare statement: " . $dbh->errstr;
    $sth->execute() or die "Cannot execute: " . $sth->errstr();
    if ( my @data = $sth->fetchrow_array() ) {
        $readTRsMappedToIndistinguishable = $data[0];
    }
    $sth->finish();

    my $VNTRasSingleton
        = -666;    # INITIAL REF-TRs RDE GE7 PC ADDBACK MAP SINGLETON VNTR
    my $VNTRasDistinguishable
        = -666;  # INITIAL REF-TRs RDE GE7 PC ADDBACK MAP DISTINGUISHABLE VNTR
    my $VNTRasIndistinguishable = -666
        ;    # INITIAL REF-TRs RDE GE7 PC ADDBACK MAP INDISTINGUISHABLE VNTR

    $sth
        = $dbh->prepare(
        "SELECT count(distinct rid) FROM refdb.fasta_ref_reps reftab INNER JOIN main.fasta_ref_reps USING (rid) INNER JOIN map ON reftab.rid=map.refid WHERE bbb=1 AND is_singleton=1 AND support_vntr=1;"
        ) or die "Couldn't prepare statement: " . $dbh->errstr;
    $sth->execute() or die "Cannot execute: " . $sth->errstr();
    if ( my @data = $sth->fetchrow_array() ) {
        $VNTRasSingleton = $data[0];
    }
    $sth->finish();

    $sth
        = $dbh->prepare(
        "SELECT count(distinct rid) FROM refdb.fasta_ref_reps reftab INNER JOIN main.fasta_ref_reps USING (rid) INNER JOIN map ON reftab.rid=map.refid WHERE bbb=1 AND is_dist=1 AND is_singleton=0 AND support_vntr=1;"
        ) or die "Couldn't prepare statement: " . $dbh->errstr;
    $sth->execute() or die "Cannot execute: " . $sth->errstr();
    if ( my @data = $sth->fetchrow_array() ) {
        $VNTRasDistinguishable = $data[0];
    }
    $sth->finish();

    $sth
        = $dbh->prepare(
        "SELECT count(distinct rid) FROM refdb.fasta_ref_reps reftab INNER JOIN main.fasta_ref_reps USING (rid) INNER JOIN map ON reftab.rid=map.refid WHERE bbb=1 AND is_indist=1 AND support_vntr=1;"
        )

        or die "Couldn't prepare statement: " . $dbh->errstr;
    $sth->execute() or die "Cannot execute: " . $sth->errstr();
    if ( my @data = $sth->fetchrow_array() ) {
        $VNTRasIndistinguishable = $data[0];
    }
    $sth->finish();

    my $refWithSupport
        = $sum_has_support;   # INITIAL REF-TRs RDE GE7 PC ADDBACK MAP SUPPORT
    my $refAsHomozygousSame = $sum_homez_same
        ;    # INITIAL REF-TRs RDE GE7 PC ADDBACK MAP SUPPORT HOMZYGOUS SAME
    my $refAsHomozygousDiff = $sum_homez_diff
        ;    # INITIAL REF-TRs RDE GE7 PC ADDBACK MAP SUPPORT HOMZYGOUS DIFF
    my $refAsHeterozygousSame = $sum_hetez_same
        ;   # INITIAL REF-TRs RDE GE7 PC ADDBACK MAP SUPPORT HETEROZYGOUS SAME
    my $refAsHeterozygousDiff = $sum_hetez_diff
        ;   # INITIAL REF-TRs RDE GE7 PC ADDBACK MAP SUPPORT HETEROZYGOUS DIFF
    my $refAsMulti = $sum_hetez_multi
        ;    # INITIAL REF-TRs RDE GE7 PC ADDBACK MAP SUPPORT MULTI

    # ***************************************************

    # formatted to show commas at thousands
    my $form_totalReads        = commify($totalReads);
    my $form_totalReadsWithTRs = commify($totalReadsWithTRs);
    my $form_totalReadsWithTRsPatternGE7
        = commify($totalReadsWithTRsPatternGE7);
    my $form_totalReadTRs          = commify($totalReadTRs);
    my $form_readTRsWithPatternGE7 = commify($readTRsWithPatternGE7);
    my $form_readTRsWPGE7AfterCyclicRedundancyElimination
        = commify($readTRsWPGE7AfterCyclicRedundancyElimination);
    my $form_readTRsProfileClustered  = commify($readTRsProfileClustered);
    my $form_readTRsMapped            = commify($readTRsMapped);
    my $form_readTRsMappedToSingleton = commify($readTRsMappedToSingleton);
    my $form_readTRsMappedToDistinguishable
        = commify($readTRsMappedToDistinguishable);
    my $form_readTRsMappedToIndistinguishable
        = commify($readTRsMappedToIndistinguishable);

    my $form_refTRsAREWithPatternGE7 = commify($refTRsAREWithPatternGE7);
    my $form_refTRsAREWPGE7AfterCyclicRedundancyElimination
        = commify($refTRsAREWPGE7AfterCyclicRedundancyElimination);
    my $form_refTRsProfileClustered = commify($refTRsProfileClustered);
    my $form_refTRsMapped           = commify($refTRsMapped);
    my $form_refTRsMappedSpan1      = commify($refTRsMappedSpan1);
    my $form_refTRsMappedSpanN      = commify($refTRsMappedSpanN);
    my $form_refTRsMappedSingleton  = commify($refTRsMappedSingleton);
    my $form_refTRsMappedDistinguishable
        = commify($refTRsMappedDistinguishable);
    my $form_refTRsMappedIndistinguishable
        = commify($refTRsMappedIndistinguishable);

    my $form_VNTRasSingleton         = commify($VNTRasSingleton);
    my $form_VNTRasDistinguishable   = commify($VNTRasDistinguishable);
    my $form_VNTRasIndistinguishable = commify($VNTRasIndistinguishable);

    my $form_refWithSupport        = commify($refWithSupport);
    my $form_refAsHomozygousSame   = commify($refAsHomozygousSame);
    my $form_refAsHomozygousDiff   = commify($refAsHomozygousDiff);
    my $form_refAsHeterozygousSame = commify($refAsHeterozygousSame);
    my $form_refAsHeterozygousDiff = commify($refAsHeterozygousDiff);
    my $form_refAsMulti            = commify($refAsMulti);
    my $form_ReadTRsSupport        = commify($ReadTRsSupport);

    # DO NOT COMPUTE
    my $percentReadTRsProfileClustered;
    my $percentReadTRsMapped;
    my $percentRefTRsProfileClustered;
    my $percentRefTRsMapped;
    my $percentRefTRsMappedSingleton;
    my $percentRefTRsMappedDistinguishable;
    my $percentRefTRsMappedIndistinguishable;
    my $percentReadTRsMappedToSingleton;
    my $percentReadTRsMappedToDistinguishable;
    my $percentReadTRsMappedToIndistinguishable;
    my $percentVNTRasSingleton;
    my $percentVNTRasDistinguishable;
    my $percentVNTRasIndistinguishable;
    my $percentRefTRsMappedSpan1;
    my $percentRefTRsMappedSpanN;
    my $percentRefWithSupport;
    my $percentRefAsHomozygousSame;
    my $percentRefAsHomozygousDiff;
    my $percentRefAsHeterozygousSame;
    my $percentRefAsHeterozygousDiff;
    my $percentRefAsMulti;
    my $percentVNTRTotalByGeneticClass;
    my $VNTRTotalByRefClass
        = commify( $stat_hash->{CLUST_NUMBER_OF_REFS_WITH_PREDICTED_VNTR} );
    my $VNTRTotalByGeneticClass;

    # latex header
    open( my $texfh, ">", "${result_prefix}.span${MIN_SUPPORT_REQUIRED}.tex" )
        or die
        "\nCan't open for writing ${result_prefix}.span${MIN_SUPPORT_REQUIRED}.tex\n\n";
    printf( $texfh "\n\\documentclass[twoside,10pt]{article}"
            . "\n\\usepackage{underscore}"
            . "\n\\setlength{\\topmargin}{-1cm}"
            . "\n\\setlength{\\oddsidemargin}{-0.5cm}"
            . "\n\\setlength{\\evensidemargin}{-0.5cm}"
            . "\n\\setlength{\\textwidth}{6.5in}"
            . "\n\\setlength{\\textheight}{9in}"
            . "\n\\setlength{\\parskip}{0.1 in}"
            . "\n\\setlength{\\parindent}{0.0 in}"
            . "\n\\begin{document}"

            . "\n\\begin{center}"
            . "\n{\\Large{\\bf VNTRseek analysis of %s}}"
            . "\n\\end{center}", $DBSUFFIX
    );

    # TR Clustering and Mapping Outcomes
    $percentRefTRsProfileClustered
        = $refTRsAREWithPatternGE7
        ? int( 100 * $refTRsProfileClustered / $refTRsAREWithPatternGE7 )
        : 0;
    $percentRefTRsMapped
        = $refTRsAREWithPatternGE7
        ? int( 100 * $refTRsMapped / $refTRsAREWithPatternGE7 )
        : 0;
    $percentReadTRsProfileClustered
        = $readTRsWithPatternGE7
        ? int( 100 * $readTRsProfileClustered / $readTRsWithPatternGE7 )
        : 0;
    $percentReadTRsMapped
        = $readTRsWithPatternGE7
        ? int( 100 * $readTRsMapped / $readTRsWithPatternGE7 )
        : 0;
    $percentRefTRsMappedSingleton
        = $refTRsMapped
        ? int( 100 * $refTRsMappedSingleton / $refTRsMapped )
        : 0;
    $percentRefTRsMappedDistinguishable
        = $refTRsMapped
        ? int( 100 * $refTRsMappedDistinguishable / $refTRsMapped )
        : 0;
    $percentRefTRsMappedIndistinguishable
        = $refTRsMapped
        ? int( 100 * $refTRsMappedIndistinguishable / $refTRsMapped )
        : 0;
    $percentReadTRsMappedToSingleton
        = $readTRsMapped
        ? int( 100 * $readTRsMappedToSingleton / $readTRsMapped )
        : 0;
    $percentReadTRsMappedToDistinguishable
        = $readTRsMapped
        ? int( 100 * $readTRsMappedToDistinguishable / $readTRsMapped )
        : 0;
    $percentReadTRsMappedToIndistinguishable
        = $readTRsMapped
        ? int( 100 * $readTRsMappedToIndistinguishable / $readTRsMapped )
        : 0;
    $VNTRTotalByRefClass
        = $VNTRasSingleton
        + $VNTRasDistinguishable
        + $VNTRasIndistinguishable;
    my $form_VNTRTotalByRefClass = commify($VNTRTotalByRefClass);
    $percentVNTRasSingleton
        = $VNTRTotalByRefClass
        ? int( 100 * $VNTRasSingleton / $VNTRTotalByRefClass )
        : 0;
    $percentVNTRasDistinguishable
        = $VNTRTotalByRefClass
        ? int( 100 * $VNTRasDistinguishable / $VNTRTotalByRefClass )
        : 0;
    $percentVNTRasIndistinguishable
        = $VNTRTotalByRefClass
        ? int( 100 * $VNTRasIndistinguishable / $VNTRTotalByRefClass )
        : 0;

    # VNTR calling
    $percentRefTRsMappedSpan1
        = $refTRsAREWithPatternGE7
        ? int( 100 * $refTRsMappedSpan1 / $refTRsAREWithPatternGE7 )
        : 0;
    $percentRefTRsMappedSpanN
        = $refTRsAREWithPatternGE7
        ? int( 100 * $refTRsMappedSpanN / $refTRsAREWithPatternGE7 )
        : 0;
    $percentRefWithSupport
        = $refTRsAREWithPatternGE7
        ? int( 100 * $refWithSupport / $refTRsAREWithPatternGE7 )
        : 0;
    $VNTRTotalByGeneticClass
        = $refAsHomozygousDiff
        + $refAsHeterozygousSame
        + $refAsHeterozygousDiff
        + $refAsMulti;

    #$form_VNTRTotalByGeneticClass = commify($VNTRTotalByGeneticClass);
    $percentRefAsHomozygousSame
        = $refWithSupport
        ? int( 100 * $refAsHomozygousSame / $refWithSupport )
        : 0;
    $percentRefAsHomozygousDiff
        = $refWithSupport
        ? int( 100 * $refAsHomozygousDiff / $refWithSupport )
        : 0;
    $percentRefAsHeterozygousSame
        = $refWithSupport
        ? int( 100 * $refAsHeterozygousSame / $refWithSupport )
        : 0;
    $percentRefAsHeterozygousDiff
        = $refWithSupport
        ? int( 100 * $refAsHeterozygousDiff / $refWithSupport )
        : 0;
    $percentRefAsMulti
        = $refWithSupport ? int( 100 * $refAsMulti / $refWithSupport ) : 0;
    $percentVNTRTotalByGeneticClass
        = $refWithSupport
        ? int( 100 * $VNTRTotalByGeneticClass / $refWithSupport )
        : 0;

#$percentReadTRsSupport = $readTRsWithPatternGE7 ? int(100*$ReadTRsSupport/$readTRsWithPatternGE7) : 0;

    my $form_readsMapped   = "";
    my $percentReadsMapped = $form_readsMapped;

    $sth
        = $dbh->prepare(
        "SELECT count(distinct head) FROM map INNER JOIN replnk ON map.readid=replnk.rid INNER JOIN fasta_reads ON fasta_reads.sid=replnk.sid  WHERE bbb=1;"
        ) or die "Couldn't prepare statement: " . $dbh->errstr;
    $sth->execute() or die "Cannot execute: " . $sth->errstr();
    {
        my @data = $sth->fetchrow_array();
        $form_readsMapped = commify( int( $data[0] ) );
        $percentReadsMapped
            = $totalReadsWithTRsPatternGE7
            ? int( 100 * int( $data[0] ) / $totalReadsWithTRsPatternGE7 )
            : 0;
    }
    $sth->finish();

    #New Table Formats
    my $percentRefTRsMappedSpan1byMapped;
    my $percentRefTRsMappedSpanNbyMapped;
    my $percentRefWithSupportbyMapped;
    my $form_refOneAlleleDiff;
    my $form_refTwoAllelesSame;
    my $form_refTwoAllelesDiff;
    my $form_refMultiAlleles;
    my $percentRefOneAlleleDiff;
    my $percentRefTwoAllelesSame;
    my $percentRefTwoAllelesDiff;
    my $percentRefMultiAlleles;

    $percentRefTRsMappedSpan1byMapped
        = $refTRsMapped ? int( 100 * $refTRsMappedSpan1 / $refTRsMapped ) : 0;
    $percentRefTRsMappedSpanNbyMapped
        = $refTRsMapped ? int( 100 * $refTRsMappedSpanN / $refTRsMapped ) : 0;
    $percentRefWithSupportbyMapped
        = $refTRsMapped ? int( 100 * $refWithSupport / $refTRsMapped ) : 0;
    $form_refOneAlleleDiff  = commify($refAsHomozygousDiff);
    $form_refTwoAllelesSame = commify($refAsHeterozygousSame);
    $form_refTwoAllelesDiff = commify($refAsHeterozygousDiff);
    $form_refMultiAlleles   = commify($refAsMulti);
    $percentRefOneAlleleDiff
        = $VNTRTotalByGeneticClass
        ? int( 100 * $refAsHomozygousDiff / $VNTRTotalByGeneticClass )
        : 0;
    $percentRefTwoAllelesSame
        = $VNTRTotalByGeneticClass
        ? int( 100 * $refAsHeterozygousSame / $VNTRTotalByGeneticClass )
        : 0;
    $percentRefTwoAllelesDiff
        = $VNTRTotalByGeneticClass
        ? int( 100 * $refAsHeterozygousDiff / $VNTRTotalByGeneticClass )
        : 0;
    $percentRefMultiAlleles
        = $VNTRTotalByGeneticClass
        ? int( 100 * $refAsMulti / $VNTRTotalByGeneticClass )
        : 0;

    printf( $texfh "\n\\begin{table}[htdp]"
            . "\n\\begin{center}"
            . "\nA. Mapping\\\\"
            . "\n\\vspace{.1in}"
            . "\n\\begin{tabular}{|c|c||c|c|c||}"
            . "\n\\hline"
            . "\n\\multicolumn{5}{|c||}{Mapping}\\\\"
            . "\n\\hline"
            . "\n&&Input&&\\\\"
            . "\n&Total&(After Filters)&Mapped&\\%%\\\\"
            . "\n\\hline"
            . "\nReference TRs&---&%s&%s&%s\\\\"
            . "\n\\hline"
            . "\nRead TRs&%s&%s&%s&%s\\\\"
            . "\n\\hline"
            . "\nReads&%s&%s&%s&%s\\\\"
            . "\n\\hline"
            . "\n\\end{tabular}\\\\"
            . "\n\\vspace{.2in}" . "\n"
            . "\nB. Reference Results\\\\"
            . "\n\\vspace{.1in}"
            . "\n\\begin{tabular}{|c|c||c||c|c||}"
            . "\n\\hline"
            . "\n\\multicolumn{5}{|c||}{Mapped Reference TRs}\\\\"
            . "\n\\hline"
            . "\n\\multicolumn{2}{|c||}{}&At Least&\\multicolumn{2}{c||}{}\\\\"
            . "\n\\multicolumn{2}{|c||}{Mapped by at least}&One Allele&\\multicolumn{2}{c||}{By Reference Category}\\\\"
            . "\n\\cline{1-2}"
            . "\n\\cline{4-5}"
            . "\nOne Read&Two Reads& Supported&Singleton&Indistinguishable\\\\\\hline"
            . "\n%s&%s&%s&%s&%s\\\\"
            . "\n\\hline"
            . "\n%s\\\%%&%s\\%%&%s\\%%&%s\\%%&%s\\%%\\\\"
            . "\n\\hline"
            . "\n\\end{tabular}"
            . "\n\\vspace{.2in}" . "\n" . "\n"
            . "\nC. VNTR Results\\\\"
            . "\n\\vspace{.1in}"
            . "\n\\begin{tabular}{|c||c|c|c|c||c|c||c||}"
            . "\n\\hline"
            . "\n\\multicolumn{7}{|c||}{VNTRs}\\\\"
            . "\n\\hline"
            . "\n&\\multicolumn{4}{c||}{Alleles Supported}&\\multicolumn{2}{c||}{}\\\\"
            . "\n\\cline{2-5}"
            . "\n&One&\\multicolumn{3}{c||}{Two or More}&\\multicolumn{2}{c||}{By Reference Category}\\\\"
            . "\n\\cline{2-5}"
            . "\n&\$\\star\$&\$\\bullet\$&\$\\bullet\$&\$\\bullet\$&\\multicolumn{2}{c||}{}\\\\"
            . "\n\\cline{6-7}"
            . "\nTotal&Diff&Same/Diff&Diff/Diff &Multi&Singleton&Indistinguishable\\\\"
            . "\n\\hline"
            . "\n%s&%s&%s&%s&%s&%s&%s\\\\"
            . "\n\\hline"
            . "\n100\\%%&%s\\%%&%s\\%%&%s\\%%&%s\\%%&%s\\%%&%s\\%%\\\\"
            . "\n\\hline"
            . "\n\\end{tabular}\\\\"
            . "\n\\ \\\\"
            . "\n\$\\star\$ Inferred VNTR \\quad \$\\bullet\$ Observed VNTR"
            . "\n\\vspace{.2in}"
            . "\n\\end{center}"
            . "\n\\caption{{\\bf VNTRseek Results.} }"
            . "\n\\label{table:mapping and mapped by reference category}"
            . "\n\\end{table}\%%" . "\n",
        $form_refTRsAREWithPatternGE7,
        $form_refTRsMapped,
        $percentRefTRsMapped,
        $form_totalReadTRs,
        $form_readTRsWithPatternGE7,
        $form_readTRsMapped,
        $percentReadTRsMapped,
        $form_totalReads,
        $form_totalReadsWithTRsPatternGE7,
        $form_readsMapped,
        $percentReadsMapped,
        $form_refTRsMappedSpan1,
        $form_refTRsMappedSpanN,
        $form_refWithSupport,
        $form_refTRsMappedSingleton,
        $form_refTRsMappedIndistinguishable,
        ,
        $percentRefTRsMappedSpan1byMapped,
        $percentRefTRsMappedSpanNbyMapped,
        $percentRefWithSupportbyMapped,
        $percentRefTRsMappedSingleton,
        $percentRefTRsMappedIndistinguishable,
        $form_VNTRTotalByRefClass,
        $form_refOneAlleleDiff,
        $form_refTwoAllelesSame,
        $form_refTwoAllelesDiff,
        $form_refMultiAlleles,
        $form_VNTRasSingleton,
        $form_VNTRasIndistinguishable,
        ,
        $percentRefOneAlleleDiff,
        $percentRefTwoAllelesSame,
        $percentRefTwoAllelesDiff,
        $percentRefMultiAlleles,
        $percentVNTRasSingleton,
        $percentVNTRasIndistinguishable
    );
    printf( $texfh "\n\\end{document}" );
    $dbh->disconnect;
    close($texfh);
}

####################################
sub calc_entropy {
    my @ACGTcount = (0) x 4;
    my @diversity = ();
    my $in        = uc(shift);
    my $len       = length($in);
    my $i         = 0;
    my $count     = 0;

    while ( $i < $len ) {
        my $chr = substr( $in, $i, 1 );
        if ( $chr eq 'A' ) { $ACGTcount[0]++; $count++; }
        if ( $chr eq 'C' ) { $ACGTcount[1]++; $count++; }
        if ( $chr eq 'G' ) { $ACGTcount[2]++; $count++; }
        if ( $chr eq 'T' ) { $ACGTcount[3]++; $count++; }
        $i++;
    }

    if ( $count == 0 ) { return 0.0; }

    $diversity[0] = $ACGTcount[0] / $count;
    $diversity[1] = $ACGTcount[1] / $count;
    $diversity[2] = $ACGTcount[2] / $count;
    $diversity[3] = $ACGTcount[3] / $count;

    for ( my $e = 2; $e >= 0; $e-- ) {
        for ( my $f = 0; $f <= $e; $f++ ) {
            if ( $diversity[$f] < $diversity[ $f + 1 ] ) {
                my $temp = $diversity[ $f + 1 ];
                $diversity[ $f + 1 ] = $diversity[$f];
                $diversity[$f] = $temp;
            }
        }
    }

    my $entropy = (
        (     ( $diversity[0] == 0 ) ? 0
            : ( $diversity[0] * ( log( $diversity[0] ) / log(2) ) )
        ) + (
            ( $diversity[1] == 0 ) ? 0
            : ( $diversity[1] * ( log( $diversity[1] ) / log(2) ) )
        ) + (
            ( $diversity[2] == 0 ) ? 0
            : ( $diversity[2] * ( log( $diversity[2] ) / log(2) ) )
        ) + (
            ( $diversity[3] == 0 ) ? 0
            : ( $diversity[3] * ( log( $diversity[3] ) / log(2) ) )
        )
    );

    if ( $entropy < 0 ) { $entropy = -$entropy; }

    return $entropy;
}

############ Main #############

my $num;
my $i;

# read representative file
my %DHASH   = ();
my %REPHASH = ();

#open (MYFILE, "$filerep") or die "Can't open for reading $filerep.";
#while (<MYFILE>) {
#
#  if ($_ =~ /(\d+)_(\d+)\W(.+)\W(.+)/i) {
#    $id=$1;
#    $fd=$3;
#    $pt=$4;
#print $id. " " . $fd. " ". $pt . "\n";
#    $DHASH{$id} = $fd;
#    $REPHASH{$id} = $pt;
#  }
#}

set_statistics( { N_MIN_SUPPORT => $MIN_SUPPORT_REQUIRED } );

my $dbh = get_dbh()
    or die "Could not connect to database: $DBI::errstr";

#goto AAA;

# Select all refids from vntr_support table.
# Must negate refid to match.
our $supported_alleles = $dbh->selectall_arrayref(
    q{SELECT -refid,copies,sameasref,support
    FROM vntr_support}
) or die "Couldn't select from supported VNTRs: " . $dbh->errstr;
my ($supported_vntr_count) = $dbh->selectrow_array(q{SELECT COUNT(DISTINCT refid) FROM vntr_support});

( $ENV{DEBUG} )
    && warn "Supported ref alleles:\n" . Dumper($supported_alleles) . "\n";

# Connect to refdb
my $ref_dbh
    = get_ref_dbh(
    @run_conf{qw( REFERENCE_SEQ REFERENCE_FILE REFERENCE_INDIST REDO_REFDB )}
    );

# register the module and declare the virtual table
$ref_dbh->sqlite_create_module(
    perl => "DBD::SQLite::VirtualTable::PerlData" );
$ref_dbh->do(
    q{CREATE VIRTUAL TABLE temp.vntr_support
    USING perl(refid INT, copies INT, sameasref INT, support INT, arrayrefs="main::supported_alleles")}
);

# warn "\nTurning off AutoCommit\n";
$dbh->do("PRAGMA foreign_keys = OFF");
$dbh->do("PRAGMA synchronous = OFF");
$dbh->begin_work;
$dbh->do( get_trunc_query( $run_conf{BACKEND}, "main.fasta_ref_reps" ) )
    or die "Couldn't do statement: " . $dbh->errstr;

# now we can SELECT from reference TR table, using the virtual table
# to JOIN.
my $get_supported_reftrs_sth = $ref_dbh->prepare(
    q{SELECT rid,pattern,copies,sameasref,support
    FROM fasta_ref_reps INNER JOIN temp.vntr_support ON (rid=refid)
    ORDER BY rid ASC}
) or die "Couldn't prepare statement: " . $dbh->errstr;

# my $sth6
#     = $dbh->prepare(
#     "SELECT copies,sameasref,support FROM vntr_support WHERE refid=?")
#     or die "Couldn't prepare statement: " . $dbh->errstr;

warn "\n\nUpdating fasta_ref_reps table...\n";

my $query = qq{INSERT INTO main.fasta_ref_reps
        (rid, alleles_sup, allele_sup_same_as_ref, entropy,
            has_support, span1, spanN, homez_same, homez_diff,
            hetez_same, hetez_diff, hetez_multi, support_vntr,
            support_vntr_span1)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
    };

my $update_ref_table_sth = $dbh->prepare($query)
    or die "Couldn't prepare statement: " . $dbh->errstr;

my $ReadTRsSupport = 0;
my @supported_refTRs;

# my $refid = -1;
$get_supported_reftrs_sth->execute()
    or die "Cannot execute: " . $get_supported_reftrs_sth->errstr();
$i = 0;
my $ref = { refid => -1 };
while ( my @data = $get_supported_reftrs_sth->fetchrow_array() ) {
    warn "Row $i: ", join( ", ", @data ), "\n"
        if ( $ENV{DEBUG} );
    if ( $i > 0
        && ( $ref->{refid} != $data[0] && ( @supported_refTRs % 1e4 == 0 ) ) )
    {
        while ( my $r = shift @supported_refTRs ) {
            update_ref_table( $r, $update_ref_table_sth );
        }
    }

    if ( $ref->{refid} != $data[0] ) {
        ( $i > 0 ) && push( @supported_refTRs, $ref );

        # Note: refid is positive
        $ref = {
            refid              => $data[0],
            entr               => calc_entropy( $data[1] ),
            homez_same         => 0,
            homez_diff         => 0,
            hetez_same         => 0,
            hetez_diff         => 0,
            hetez_multi        => 0,
            support_vntr       => 0,
            span1              => 0,
            spanN              => 0,
            nsupport           => 0,
            has_support        => 0,
            nsameasref         => 0,
            readsum            => 0,
            support_vntr_span1 => 0
        };

        $i++;
    }

    my $copies    = $data[2];
    my $sameasref = $data[3];
    my $support   = $data[4];

    if ( $support >= $MIN_SUPPORT_REQUIRED ) {

        $ReadTRsSupport += $support;

        $ref->{has_support} = 1;
        $ref->{nsupport}++;

        if ( $sameasref > 0 ) {
            $ref->{nsameasref} = 1;
        }
    }

    if ( $support >= 1 && $sameasref == 0 ) {
        $ref->{support_vntr_span1} = 1;
    }

    $ref->{readsum} += $support;
}

my $updfromtable = $i;
$get_supported_reftrs_sth->finish;
$ref_dbh->disconnect;

# $sth6->finish;

# Insert last rows:
push @supported_refTRs, $ref;
while ( my $r = shift @supported_refTRs ) {
    update_ref_table( $r, $update_ref_table_sth );
}
$dbh->commit;

if ( $updfromtable != $supported_vntr_count ) {
    die
        "Updated number of entries($updfromtable) not equal to the number of references, aborting!";
}

# cleanup temp file
unlink("$TEMPDIR/updaterefs_$DBSUFFIX.txt");

# updating stats table
$dbh->begin_work;
print STDERR "Updating stats table...\n";
my $update_stats_sth = $dbh->prepare(
    q{UPDATE stats SET NUMBER_REFS_SINGLE_REF_CLUSTER_WITH_READS_MAPPED=(
    select count(*)
        FROM (
            select count(*) as thecount
            from clusterlnk
            where repeatid<0 group by clusterid having thecount=1
        ) f
    )}
) or die "Couldn't prepare statement: " . $dbh->errstr;
$update_stats_sth->execute()
    or die "Cannot execute: " . $update_stats_sth->errstr();
$update_stats_sth->finish;

$dbh->do('CREATE TABLE t1 (c1 INT PRIMARY KEY NOT NULL);')
    or die "Couldn't do statement: " . $dbh->errstr;

my $sth = $dbh->prepare(
    q{insert into t1 select urefid
    from (
        select count(*) as thecount,max(repeatid) as urefid
        from clusterlnk
        where repeatid<0
        group by clusterid
        having thecount=1
    ) f}
) or die "Couldn't prepare statement: " . $dbh->errstr;
$sth->execute() or die "Cannot execute: " . $sth->errstr();

$sth = $dbh->prepare(
    q{UPDATE stats SET NUMBER_REFS_SINGLE_REF_CLUSTER = (
        select count(distinct repeatid)
        FROM clusterlnk 
        WHERE repeatid IN (select c1 from t1)
    )}
) or die "Couldn't prepare statement: " . $dbh->errstr;
$sth->execute() or die "Cannot execute: " . $sth->errstr();

$sth
    = $dbh->prepare(
    'UPDATE stats SET NUMBER_REFS_SINGLE_REF_CLUSTER_WITH_READS_MAPPED = (select count(distinct map.refid) FROM clusterlnk INNER JOIN map ON map.refid=-clusterlnk.repeatid WHERE repeatid IN (select c1 from t1));'
    ) or die "Couldn't prepare statement: " . $dbh->errstr;
$sth->execute() or die "Cannot execute: " . $sth->errstr();

$sth = $dbh->prepare('drop table t1;')
    or die "Couldn't prepare statement: " . $dbh->errstr;
$sth->execute() or die "Cannot execute: " . $sth->errstr();

$sth
    = $dbh->prepare(
    "UPDATE stats SET NUMBER_REFS_SINGLE_REF_CLUSTER_WITH_NO_READS_MAPPED = NUMBER_REFS_SINGLE_REF_CLUSTER - NUMBER_REFS_SINGLE_REF_CLUSTER_WITH_READS_MAPPED;"
    ) or die "Couldn't prepare statement: " . $dbh->errstr;
$sth->execute() or die "Cannot execute: " . $sth->errstr();

# Update last few stats:
my ( $mapped, $rank, $rankflank, $num_spanN );

($mapped) = $dbh->selectrow_array(q{SELECT count(*) FROM map})
    or die "Couldn't select map count: " . $dbh->errstr;

($rank) = $dbh->selectrow_array(q{SELECT count(*) FROM rank})
    or die "Couldn't select rank count: " . $dbh->errstr;

($rankflank) = $dbh->selectrow_array(q{SELECT count(*) FROM rankflank})
    or die "Couldn't select rankflank count: " . $dbh->errstr;

# update spanN number on stats
($num_spanN) = $dbh->selectrow_array(
    q{SELECT count(*) FROM fasta_ref_reps 
    WHERE support_vntr > 0}
) or die "Couldn't select span N count: " . $dbh->errstr;

$dbh->do(
    qq{UPDATE stats SET
        NUMBER_MAPPED=$mapped,
        NUMBER_RANK=$rank,
        NUMBER_RANKFLANK = $rankflank,
        NUMBER_REFS_VNTR_SPAN_N = $num_spanN}
) or die "Couldn't do statement: " . $dbh->errstr;

$dbh->commit;

# set old db settings
$dbh->do("PRAGMA foreign_keys = ON");
$dbh->do("PRAGMA synchronous = ON");
$dbh->disconnect;

warn "Producing output LaTeX file...\n";
print_latex($ReadTRsSupport);

warn "Producing distribution file...\n";
print_distr();

# Print VCF
warn "Producing VCF...\n";
print_vcf();

# All with support
# warn "Producing VCF (all with support)...\n";
# print_vcf(1);

warn "Done!\n";
