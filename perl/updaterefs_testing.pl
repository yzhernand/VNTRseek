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

    # TODO Add code for "all with support" functionality
    my $dbh            = shift;
    my $allwithsupport = shift;

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
    my $filename
        = "${DBNAME}."
        . ( ($allwithsupport) ? "allwithsupport." : "" )
        . "span${MIN_SUPPORT_REQUIRED}" . ".vcf";
    open my $vcffile, ">", $filename
        or die "\nCan't open for writing $filename\n\n";

    print $vcffile "##fileformat=VCFv4.1\n"
        . strftime( "##fileDate=\"%Y%m%d\"\n", localtime )
        . qq[##source="Vntrseek ver. $VERSION"
##TRFParameters="] . GetStatistics("PARAM_TRF") . qq["
##referenceseq="] . GetStatistics("FILE_REFERENCE_SEQ") . qq["
##referenceprofile="] . GetStatistics("FILE_REFERENCE_LEB") . qq["
##numrefTRs="] . GetStatistics("NUMBER_REF_TRS") . qq["
##readseqfolder="] . GetStatistics("FOLDER_FASTA") . qq["
##readprofilefolder="] . GetStatistics("FOLDER_PROFILES") . qq["
##numreads="] . GetStatistics("NUMBER_READS") . qq["
##numreadTRs="] . GetStatistics("NUMBER_TRS_IN_READS") . qq["
##numVNTRs="$numvntrs"
##numTRsWithSupport="$numsup"
##database="$DBNAME"
##databaseurl="http://${HTTPSERVER}/result.php?db=${DBNAME}"
##INFO=<ID=RC,Number=1,Type=Float,Description="Reference Copies">
##INFO=<ID=RPL,Number=1,Type=Integer,Description="Reference Pattern Length">
##INFO=<ID=RAL,Number=1,Type=Integer,Description="Reference Tandem Array Length">
##INFO=<ID=RCP,Number=1,Type=String,Description="Reference Consensus Pattern">
##INFO=<ID=ALGNURL,Number=1,Type=String,Description="Alignment URL">
##FILTER=<ID=SC,Description="Reference is Singleton">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=SP,Number=A,Type=Integer,Description="Number of Spanning Reads">
##FORMAT=<ID=CGL,Number=A,Type=Integer,Description="Copies Gained or Lost with respect to reference">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t$BASENAME
];

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

                    my $atemp = substr(
                        $dna,
                        $readTRStart - 1,
                        ( ( $readTRStop - $readTRStart ) + 1 )
                    );

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

# Span N
print_vcf($dbh);
# all with support
print_vcf($dbh, 1);
