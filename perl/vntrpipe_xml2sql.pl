#!/usr/bin/env perl

use strict;
use warnings;

use DBI;
#use XML::LibXML;
use Getopt::Std;

use FindBin;
use File::Basename;

use lib "$FindBin::RealBin/lib"; 
require "vutil.pm";

use vutil qw(get_config get_dbh get_trunc_query);

sub trim($)
{
        my $string = shift;
        $string =~ s/^\s+//;
        $string =~ s/\s+$//;
        return $string;
}

sub enter_new_paramset {
	my ($dbh, $backend, $flength, $ferrors) = @_;
	my $sth = $dbh->prepare('INSERT INTO flank_params (ferrors, flength) VALUES (?, ?)');
	$sth->execute($flength, $ferrors) or die $sth->errstr;
	$dbh->commit;
	if ($backend eq "mysql") {
		$sth = $dbh->prepare('SELECT LAST_INSERT_ID()');
	}
	elsif ($backend eq "sqlite") {
		$sth = $dbh->prepare('SELECT last_insert_rowid()');
	}
	$sth->execute() or die $sth->errstr;
	my @ret = $sth->fetchrow_array();
	my $ret = $ret[0];
	$sth->finish;
	return $ret;
}

sub set_indist {
    my ($dbh, $refid) = @_;
    my $sth = $dbh->prepare('UPDATE fasta_ref_reps SET is_singleton=0,is_indist=1 WHERE rid=?');
    $sth->execute($refid) or die $sth->errstr;
}

######################################################

my %opts;
getopts('rk:t:d:u:', \%opts);
if (!defined($opts{'k'}) || 
	!defined($opts{'t'}) || 
	!defined($opts{'d'}) || 
	!defined($opts{'u'})) 
{
        die "Usage:
        perl $0 [-r] -k <num> -t <num> -d <str> -u <str>

REQUIRED PARAMETERS
	-k <num>	 max_errors parameter used to generate flank comparison
	-t <num>	 trim_to parameter used to generate flank comparison

	-d <string>      Database suffix
	-u <string>      Folder where master file is located

OPTIONS
	-r	recreate MySQL tables from scratch (WARNING: may lead to loss of data)
";
}

# set these mysql credentials in vs.cnf (in installation directory)
my %run_conf = get_config( $opts{u} . "vs.cnf" );
my ( $LOGIN, $PASS, $HOST, $indistfile ) = @run_conf{qw(LOGIN PASS HOST REFERENCE_INDIST)};

my $dbh = get_dbh($opts{d}, $opts{u} . "vs.cnf");
if ($opts{r}) {
	warn "Clearing old tables\n";
	$dbh->do( get_trunc_query( $run_conf{BACKEND}, "flank_params" ) )
		or die "Couldn't do statement: " . $dbh->errstr;
	$dbh->do( get_trunc_query( $run_conf{BACKEND}, "flank_connection" ) )
		or die "Couldn't do statement: " . $dbh->errstr;
	# Assumes fasta_ref_reps is populated, sets values for ALL rows.
	warn "Resetting singleton/indist flags\n";
	$dbh->do("UPDATE fasta_ref_reps SET is_dist=0,is_singleton=1,is_indist=0;") 
        or die "Couldn't do statement: " . $dbh->errstr;
}
my $paramsetid = enter_new_paramset($dbh, $run_conf{BACKEND},  $opts{'k'}, $opts{'t'});

open FILE, "<$indistfile" or die $!;
my $count = 0;
while (<FILE>) {

 $count++;

 my @values = split(' ', $_);
 my $repcount = @values - 1;

 # if ($ENV{DEBUG}) {
 # 	warn $count . ". " . @values . "\n";
 # }

 my $i=0;
 foreach my $val (@values) {

   $i++;

   $val = trim($val);

   if ($val < 0) {
	set_indist($dbh, -$val); 
   }
 }
}
close(FILE);

$dbh->disconnect();
