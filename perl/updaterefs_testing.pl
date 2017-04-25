#!/usr/bin/env perl

use 5.010;
use strict;
use warnings;
use DBI;
use POSIX qw(strftime);
use FindBin;
use lib "$FindBin::Bin/vntr";
use vutil qw(get_credentials stats_get stats_set);

my $MIN_SUPPORT_REQUIRED = 2;
my $VERSION              = 1.08;
my $DBNAME               = shift;
my $HTTPSERVER           = shift;
my $BASENAME             = $DBNAME;
$BASENAME =~ s/VNTRPIPE_//;
my $MSDIR = $ENV{HOME} . "/${BASENAME}.";

my ( $LOGIN, $PASS, $HOST ) = get_credentials($MSDIR);

############################ Procedures ###############################################################
sub RC {

    # complement reversed DNA sequence
    my $seq = shift;

    $seq = reverse $seq;

    $seq =~ tr/ACGT/TGCA/;

    return $seq;
}

####################################
sub nowhitespace($) {
    my $string = shift;
    $string =~ s/\s+//g;
    return $string;
}

####################################
sub GetStatistics {

    my $argc = @_;
    if ( $argc < 1 ) {
        die "stats_set: expects 1 parameter, passed $argc !\n";
    }

    my $NAME = $_[0];

    return stats_get( $DBNAME, $LOGIN, $PASS, $HOST, $NAME );
}

####################################
sub commify {
    my $input = shift;
    $input = reverse $input;
    $input =~ s/(\d\d\d)(?=\d)(?!\d*\.)/$1,/g;
    return reverse $input;
}

sub print_vcf {
    my $dbh = shift;

    # Get total number of TRs supported
    my $numsup_sth
        = $dbh->prepare(
        "select count(distinct vntr_support.refid) FROM vntr_support WHERE support >= $MIN_SUPPORT_REQUIRED;"
        ) or die "Couldn't prepare statement: " . $dbh->errstr;

    my $numsup = 0;
    $numsup_sth->execute() or die "Cannot execute: " . $numsup_sth->errstr();
    $numsup_sth->bind_columns( \$numsup );
    $numsup_sth->fetch;

    # $numsup_sth->finish;

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

    # update spanN number on stats
    my $update_spanN_sth
        = $dbh->prepare('UPDATE stats SET NUMBER_REFS_VNTR_SPAN_N = ?;')
        or die "Couldn't prepare statement: " . $dbh->errstr;
    $update_spanN_sth->execute($numvntrs)
        or die "Cannot execute: " . $update_spanN_sth->errstr();

    # $update_spanN_sth->finish;
    open my $vcffile, ">", "${DBNAME}.span${MIN_SUPPORT_REQUIRED}.vcf"
        or die
        "\nCan't open for writing ${DBNAME}.span${MIN_SUPPORT_REQUIRED}.vcf\n\n";
    print $vcffile "##fileformat=VCFv4.1\n";
    print $vcffile strftime( "##fileDate=\"%Y%m%d\"\n", localtime );
    print $vcffile "##source=\"Vntrseek ver. $VERSION\"\n";
    print $vcffile "##TRFParameters=\"", GetStatistics("PARAM_TRF"), "\"\n";
    print $vcffile "##referenceseq=\"", GetStatistics("FILE_REFERENCE_SEQ"),
        "\"\n";
    print $vcffile "##referenceprofile=\"",
        GetStatistics("FILE_REFERENCE_LEB"), "\"\n";
    print $vcffile "##numrefTRs=\"", GetStatistics("NUMBER_REF_TRS"), "\"\n";
    print $vcffile "##readseqfolder=\"", GetStatistics("FOLDER_FASTA"), "\"\n";
    print $vcffile "##readprofilefolder=\"", GetStatistics("FOLDER_PROFILES"),
        "\"\n";
    print $vcffile "##numreads=\"", GetStatistics("NUMBER_READS"), "\"\n";
    print $vcffile "##numreadTRs=\"", GetStatistics("NUMBER_TRS_IN_READS"),
        "\"\n";
    print $vcffile "##numVNTRs=\"",          $numvntrs, "\"\n";
    print $vcffile "##numTRsWithSupport=\"", $numsup,   "\"\n";
    print $vcffile "##database=\"",          $DBNAME,   "\"\n";
    print $vcffile
        "##databaseurl=\"http://${HTTPSERVER}/result.php?db=${DBNAME}\" \n";
    print $vcffile
        "##INFO=<ID=RC,Number=1,Type=Float,Description=\"Reference Copies\">\n";
    print $vcffile
        "##INFO=<ID=RPL,Number=1,Type=Integer,Description=\"Reference Pattern Length\">\n";
    print $vcffile
        "##INFO=<ID=RAL,Number=1,Type=Integer,Description=\"Reference Tandem Array Length\">\n";
    print $vcffile
        "##INFO=<ID=RCP,Number=1,Type=String,Description=\"Reference Consensus Pattern\">\n";
    print $vcffile
        "##INFO=<ID=ALGNURL,Number=1,Type=String,Description=\"Alignment URL\">\n";
    print $vcffile "##FILTER=<ID=SC,Description=\"Reference is Singleton\">\n";

    print $vcffile
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
    print $vcffile
        "##FORMAT=<ID=SP,Number=A,Type=Integer,Description=\"Number of Spanning Reads\">\n";
    print $vcffile
        "##FORMAT=<ID=CGL,Number=A,Type=Integer,Description=\"Copies Gained or Lost with respect to reference\">\n";

    print $vcffile
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t$BASENAME\n";

# Get information on all VNTRs
# "SELECT rid,alleles_sup,allele_sup_same_as_ref,is_singleton,is_dist,is_indist,firstindex,lastindex,copynum,pattern,clusterid,reserved,reserved2,head,sequence,flankleft,direction FROM fasta_ref_reps INNER JOIN clusterlnk ON fasta_ref_reps.rid=-clusterlnk.repeatid INNER JOIN clusters ON clusters.cid=clusterlnk.clusterid WHERE support_vntr>0 ORDER BY head, firstindex;"
    my $get_vntrs_sth
        = $dbh->prepare(
        "SELECT rid,is_singleton,is_dist,firstindex,(lastindex - firstindex) + 1 AS arlen,copynum,pattern,head,sequence,flankleft,direction FROM fasta_ref_reps INNER JOIN clusterlnk ON fasta_ref_reps.rid=-clusterlnk.repeatid INNER JOIN clusters ON clusters.cid=clusterlnk.clusterid WHERE support_vntr>0 ORDER BY head, firstindex;"
        ) or die "Couldn't prepare statement: " . $dbh->errstr;
    $get_vntrs_sth->execute()
        or die "Cannot execute: " . $get_vntrs_sth->errstr();
    my ($rid,   $singleton,   $disting,     $pos1,
        $arlen, $copiesfloat, $consenuspat, $head,
        $seq,   $leftflank,   $refdir
    );
    $get_vntrs_sth->bind_columns(
        \(  $rid,   $singleton,   $disting,     $pos1,
            $arlen, $copiesfloat, $consenuspat, $head,
            $seq,   $leftflank,   $refdir
        )
    );

    # Given a TRID, get all read TRs supporting VNTR call
    my $vntr_support_sth
        = $dbh->prepare(
        "SELECT copies,sameasref,support,first,last,dna,direction FROM vntr_support LEFT OUTER JOIN replnk ON vntr_support.representative=replnk.rid LEFT OUTER JOIN clusterlnk ON replnk.rid=clusterlnk.repeatid LEFT OUTER JOIN fasta_reads ON replnk.sid=fasta_reads.sid  WHERE refid=-? ORDER BY sameasref DESC;"
        ) or die "Couldn't prepare statement: " . $dbh->errstr;

    # Loop over all VNTRs
    while ( $get_vntrs_sth->fetch() ) {
        $seq = ($seq) ? uc( nowhitespace($seq) ) : "";
        $leftflank
            = ($leftflank) ? uc( nowhitespace($leftflank) ) : "";

        warn ">$rid\n" if $ENV{DEBUG};

        my ( $first, $subSameAsRef1, $j, $al ) = (0) x 4;
        my ( $subjectA, $subjectB, $subjectC, $subjectB1, $subjectC1, $alt )
            = ("") x 6;

        my $patlen = length($consenuspat);

        # start 1 position before, else put an N there
        $seq = ( ($leftflank) ? substr( $leftflank, -1 ) : "N" ) . $seq;

# get vntr support
# "SELECT copies,sameasref,support,first,last,dna,direction FROM vntr_support LEFT OUTER JOIN replnk ON vntr_support.representative=replnk.rid LEFT OUTER JOIN clusterlnk ON replnk.rid=clusterlnk.repeatid LEFT OUTER JOIN fasta_reads ON replnk.sid=fasta_reads.sid  WHERE refid=-? ORDER BY sameasref DESC;"
        $vntr_support_sth->execute($rid);
        my ( $copies, $sameasref, $support, $readTRStart, $readTRStop, $dna,
            $readdir );
        $vntr_support_sth->bind_columns(
            \(  $copies,     $sameasref, $support, $readTRStart,
                $readTRStop, $dna,       $readdir
            )
        );
        my $jnum = $vntr_support_sth->rows;
        if ( 0 == $jnum ) {
            print STDERR
                "Error: could not get any VNTR_SUPPORT records for rid=$rid!";
            exit(1);
        }

        #print STDERR "\n";

        my $alleleWithSupportFound = 0;
        warn "\tAlt (before read TR loop): $alt\n" if $ENV{DEBUG};

        # Loop over all allele supporting read TRs
        while ( $vntr_support_sth->fetch ) {
            $dna = ($dna) ? uc( nowhitespace($dna) ) : "";

            #if (0 != $j) { print STDERR "ref: $refdir = read: $readdir\n"; }

            if ( $j == 0 ) {
                $first = $copies;
            }

            my $cdiff = $copies - $first;

            warn "\tAllele: $cdiff, support $support, alt: $alt\n"
                if $ENV{DEBUG};

            if ( $support >= $MIN_SUPPORT_REQUIRED ) {
                $alleleWithSupportFound++;
                $subjectB1     = "$support";
                $subjectC1     = "$cdiff";
                $subSameAsRef1 = $sameasref;
            }

            if ( $support >= $MIN_SUPPORT_REQUIRED ) {

                if ( $subjectA ne "" ) {
                    $subjectA .= "/";
                    $subjectB .= ",";
                    $subjectC .= ",";
                }

                if ( $alt ne "" ) {
                    $alt .= ",";
                }
                warn "\t\tAlt (after empty check): $alt\n" if $ENV{DEBUG};

                $subjectA .= "$al";
                $subjectB .= "$support";
                $subjectC .= "$cdiff";

                # If an alternate allele
                if ( 0 == $sameasref ) {

                    # If there is no DNA string for this read, exit with error
                    if ( !$dna || $dna eq "" ) {
                        print STDERR
                            "Error: read source sequence not found in database for ref ($rid) alternate allele $al!";
                        exit(1);
                    }

# If the stop position is beyond the length of the read sequence, exit with error
                    if ( $readTRStop > ( length($dna) ) ) {
                        print STDERR
                            "Error: last position outside of the range of the source sequence buffer for ref ($rid) alternate allele $al!";
                        exit(1);
                    }

                    # flip read if opposite dirs
                    if ( $refdir ne $readdir ) {
                        warn "\t\tRead TR is on complement strand\n"
                            if $ENV{DEBUG};
                        my $tdlen = length($dna);
                        $dna = RC($dna);

                        #print STDERR "$dna\n";
                        my $temp = $readTRStart;
                        $readTRStart = ( $tdlen - $readTRStop ) + 1;
                        $readTRStop  = ( $tdlen - $temp ) + 1;
                    }

                    warn
                        "\t\tSequence: $dna, TR start: $readTRStart, TR stop: $readTRStop\n"
                        if $ENV{DEBUG};

                    # find the first similar character
                    my $fl = substr( $seq, 0, 1 );

                    warn "\t\tFirst similar character in refseq: $fl\n"
                        if $ENV{DEBUG};

                    my $atemp = "N";
                    $atemp .= substr(
                        $dna,
                        $readTRStart - 1,
                        ( $readTRStop - $readTRStart + 1 )
                    );    # in case character is not encountered

                    my $extra = 1;

                 # Proposed fix: Use start index factoring in extra nucleotide
                 # by subtracting 2 from the start (one for 1 extra base, one
                 # more since substr indexes start at 0). Decrement start by 1
                 # each iter, and increment extra by one.
                    for (
                        my $start = $readTRStart - 2;
                        $start >= 1;
                        --$start, ++$extra
                        )
                    {
                        my $ll = substr( $dna, $start, 1 );

         # Check if backtracked character matches last from left flank in ref.
                        if ( $fl eq $ll ) {
                            my $substr_len
                                = ( $readTRStop - $readTRStart + 1 + $extra );
                            $atemp = substr( $dna, $start, $substr_len );
                            warn
                                "\t\t\tfl == ll; extra: $extra, start: $start, ll: $ll, atemp: $atemp\n"
                                if $ENV{DEBUG};
                            warn
                                "\t\t\t\tWill take substring at (0-indexed) position "
                                . $start
                                . " with length "
                                . $substr_len . "\n"
                                if $ENV{DEBUG};
                            last;
                        }
                        warn
                            "\t\t\textra: $extra, ll: $ll, start: $start, atemp: $atemp\n"
                            if $ENV{DEBUG};
                    }

                    $alt .= $atemp;
                    warn "\t\tAlt: $alt\n" if $ENV{DEBUG};

                }

                $al++;
            }
            elsif ($sameasref)
            {    # if nothing for reference, we want to increment counter
                $al++;
            }

            $j++;
        }

        # extra code for single allele, v 1.07
        if ( 1 == $alleleWithSupportFound ) {

            if ($subSameAsRef1) {
                $subjectA = "0/0";
            }
            else {
                $subjectA = "1/1";
            }
            $subjectB = $subjectB1;
            $subjectC = $subjectC1;
        }

        my $qual = ".";
        if ( "" eq $seq ) { $seq = "."; }
        if ( "" eq $alt ) { $alt = "."; }

        my $filter = ( $singleton == 1 ) ? "PASS" : "SC";

        my $info
            = sprintf(
            "RC=%.2lf;RPL=$patlen;RAL=$arlen;RCP=%s;ALGNURL=http://${HTTPSERVER}/index.php?db=${DBNAME}&ref=-$rid&isref=1&istab=1&ispng=1&rank=3",
            $copiesfloat, $consenuspat );
        my $format = "GT:SP:CGL";

        if ($alleleWithSupportFound) {
            print $vcffile "$head\t"
                . ( $pos1 - 1 )
                . "\ttd$rid\t$seq\t$alt\t$qual\t$filter\t$info\t$format\t$subjectA:$subjectB:$subjectC\n";
        }
    }

    # $get_vntrs_sth->finish;

    # $vntr_support_sth->finish;
    close $vcffile;
}

my $dbh = DBI->connect( "DBI:mysql:$DBNAME;mysql_local_infile=1;host=$HOST",
    "$LOGIN", "$PASS" )
    || die "Could not connect to database: $DBI::errstr";

print_vcf($dbh);
