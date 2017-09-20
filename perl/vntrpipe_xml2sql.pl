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

use vutil ('get_credentials');

sub trim($)
{
        my $string = shift;
        $string =~ s/^\s+//;
        $string =~ s/\s+$//;
        return $string;
}

sub my_connect {
	my ($DBNAME, $LOGIN, $PASS, $HOST) = @_;
	my $dbh = DBI->connect("DBI:mysql:$DBNAME;host=$HOST", "$LOGIN", "$PASS") || die $DBI::errstr;
	return $dbh;
}
sub my_disconnect {
	my $dbh = shift;
	$dbh->disconnect();
}
sub create_tables {
	my $dbh = shift;
	my $sth;

	$sth = $dbh->prepare("DROP TABLE IF EXISTS flank_connection");
	$sth->execute() or die $sth->errstr;
	$sth->finish;

        $sth = $dbh->prepare("DROP TABLE IF EXISTS flank_params");
	$sth->execute() or die $sth->errstr;
	$sth->finish;

	$sth = $dbh->prepare("
CREATE TABLE flank_params (
`paramsetid` INT(11) NOT NULL PRIMARY KEY AUTO_INCREMENT,
`flength` INT(11) NOT NULL,
`ferrors` INT(11) NOT NULL
) ENGINE=INNODB;
") or die $dbh->errstr;
        $sth->execute() or die $sth->errstr;
        $sth->finish;

	$sth = $dbh->prepare("
CREATE TABLE flank_connection (
`refid` INT(11) NOT NULL,
`clusterid` INT(11) NOT NULL,
`fcomponentsize` INT(11) NOT NULL,
`fcomponentid` INT(11) NOT NULL,
`paramsetid` INT(11) NOT NULL,
PRIMARY KEY (refid, paramsetid),
FOREIGN KEY (paramsetid) REFERENCES flank_params (paramsetid)
) ENGINE=INNODB;
") or die $dbh->errstr;
	$sth->execute() or die $sth->errstr;
	$sth->finish;


  # set default 
  $sth = $dbh->prepare("UPDATE fasta_ref_reps SET is_dist=0,is_singleton=1,is_indist=0;") 
                or die "Couldn't prepare statement: " . $dbh->errstr;
  $sth->execute() or die "Cannot execute: " . $sth->errstr();;
  $sth->finish;

}

sub enter_new_paramset {
	my ($dbh, $flength, $ferrors) = @_;
	my $sth = $dbh->prepare('INSERT INTO flank_params (ferrors, flength) VALUES (?, ?)');
	$sth->execute($flength, $ferrors) or die $sth->errstr;
	$sth->finish;
	$sth = $dbh->prepare('SELECT LAST_INSERT_ID()');
	$sth->execute() or die $sth->errstr;
	my @ret = $sth->fetchrow_array();
	my $ret = $ret[0];
	$sth->finish;
	return $ret;
}
#sub enter_new_ref {
#        my ($dbh, $refid, $clusterid, $fcomponentsize, $fcomponentid, $paramset) = @_;
#        my $sth = $dbh->prepare('INSERT IGNORE INTO flank_connection (refid, clusterid, fcomponentsize, fcomponentid, paramsetid) VALUES (?, ?, ?, ?, ?)');
#        $sth->execute($refid, $clusterid, $fcomponentsize, $fcomponentid, $paramset) or die $sth->errstr;
#        $sth->finish;
#}
#sub set_sing {
#        my ($dbh, $refid) = @_;
#        my $sth = $dbh->prepare('UPDATE fasta_ref_reps SET is_singleton=1 WHERE rid=?;');
#        $sth->execute($refid) or die $sth->errstr;
#        $sth->finish;
#}
sub set_indist {
        my ($dbh, $refid) = @_;
        my $sth = $dbh->prepare('UPDATE fasta_ref_reps SET is_singleton=0,is_indist=1 WHERE rid=?;');
        $sth->execute($refid) or die $sth->errstr;
        $sth->finish;
}

######################################################

my %opts;
getopts('rk:t:d:u:', \%opts);
if (!scalar(@ARGV) || 
	!defined($opts{'k'}) || 
	!defined($opts{'t'}) || 
	!defined($opts{'d'}) || 
	!defined($opts{'u'})) 
{
        die "Usage:
        perl $0 [-r] -k <num> -t <num> -d <str> -u <str>  <dir>

REQUIRED PARAMETERS
	-k <num>	 max_errors parameter used to generate flank comparison
	-t <num>	 trim_to parameter used to generate flank comparison

	-d <string>      Database name
	-u <string>      Folder where master file is located

OPTIONS
	-r	recreate MySQL tables from scratch (WARNING: may lead to loss of data)
";
}

# set these mysql credentials in vs.cnf (in installation directory)
my ($LOGIN,$PASS,$HOST) = get_credentials($opts{'u'});

my $dbh = my_connect($opts{'d'}, $LOGIN, $PASS, $HOST);
if ($opts{'r'}) {
	print STDERR "Dropping old tables and creating new ones\n";
	create_tables($dbh);
}
my $paramsetid = enter_new_paramset($dbh, $opts{'k'}, $opts{'t'});

my $refclusfile = $ARGV[0];


open FILE, "<$refclusfile" or die $!;
my $count = 0;
while (<FILE>) {

 $count++;

 my @values = split(' ', $_);
 my $repcount = @values - 1;

 print STDERR $count . ". " . @values;

 my $i=0;
 foreach my $val (@values) {

   $i++;

   $val = trim($val);

   if ($val < 0) {
	set_indist($dbh, -$val); 
   }
 }

 print STDERR "\n";

}
close(FILE);




my_disconnect($dbh);
