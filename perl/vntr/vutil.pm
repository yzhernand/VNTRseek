use strict;
use Cwd;
use DBI;

use FindBin;

package vutil;
use base 'Exporter';
our @EXPORT_OK = ('read_global_config_file','get_config','get_config_vars','get_credentials','set_config','set_config_vars','set_credentials','write_mysql','stats_set','stats_get','set_datetime','print_config');

# vutil.pm
# author: Yevgeniy Gelfand
# create date: Oct 30, 2010
# function: create mysql database for vntr pipleline, 
# provide functions for database management

my %VSCNF_FILE = ();
my $VSREAD = 0;

################################################################
sub trim($)
{
        my $string = shift;
        $string =~ s/^\s+//;
        $string =~ s/\s+$//;
        return $string;
}

################################################################
sub read_config_file {

 my $argc = @_;
 if ($argc <1) { die "read_config_file: expects 1 parameters, passed $argc !\n"; }
 my $startdir = $_[0];
 

 print STDERR "Using config file '${startdir}vs.cnf'...\n";


 # read local config
 if (open(MFREAD,"${startdir}vs.cnf")) {

  while (<MFREAD>) {
   chomp;
   
   # skip start comments
   if (/^\#/) { next; }
   
   if ( /(.+)=(.*)/ ) {
     my $key = trim($1);
     my $val = trim($2);
     $val =~ s/\s*\#.*//; # strip end comments
     $VSCNF_FILE{uc($key)}=$val; 
     #print $key."=".$val."\n";
   }
  }

  close(MFREAD);

 } else {

  #die "\nread_config_file: Cannot find config file (${startdir}vs.cnf). Aborting!\n";

  # create config file
  print_config($startdir);

 }

 $VSREAD = 1;

}

################################################################
sub read_global_config_file {

 # install dir
 my $argc = @_;
 if ($argc <1) { die "read_global_config_file: expects 1 parameters, passed $argc !\n"; }
 my $installdir = $_[0];
 
 # read global config
 if (open(MFREAD,"$installdir/vs.cnf")) {

  while (<MFREAD>) {
   chomp;
   
   # skip start comments
   if (/^\#/) { next; }
   
   if ( /(.+)=(.*)/ ) {
     my $key = trim($1);
     my $val = trim($2);
     $val =~ s/\s*\#.*//; # strip end comments
     $VSCNF_FILE{uc($key)}=$val; 
     #print $key."=".$val."\n";
   }
  }

  close(MFREAD);

 } else {

  warn "read_global_config_file: can't read global config file '$installdir/vs.cnf'!\n";
 }

}

################################################################
sub get_config {

 my $SERVER=""; 
 my $TMPDIR=""; 
 my $html_dir=""; 
 my $fasta_folder=""; 
 my $output_root=""; 
 my $reference_file=""; 
 my $reference_seq=""; 
 my $reference_indist="";

 my $argc = @_;
 if ($argc <1) { die "get_config: expects 1 parameters, passed $argc !\n"; }

 my $startdir = $_[0];

 if (!$VSREAD) { read_config_file($startdir); }

 if (defined $VSCNF_FILE{"SERVER"}) { $SERVER = $VSCNF_FILE{"SERVER"}; } 
 if (defined $VSCNF_FILE{"TMPDIR"}) { $TMPDIR = $VSCNF_FILE{"TMPDIR"}; } 
 if (defined $VSCNF_FILE{"HTML_DIR"}) { $html_dir = $VSCNF_FILE{"HTML_DIR"}; } 
 if (defined $VSCNF_FILE{"FASTA_DIR"}) { $fasta_folder = $VSCNF_FILE{"FASTA_DIR"}; } 
 if (defined $VSCNF_FILE{"OUTPUT_ROOT"}) { $output_root = $VSCNF_FILE{"OUTPUT_ROOT"}; } 
 if (defined $VSCNF_FILE{"REFERENCE_FILE"}) { $reference_file = $VSCNF_FILE{"REFERENCE_FILE"}; } 
 if (defined $VSCNF_FILE{"REFERENCE_SEQ"}) { $reference_seq = $VSCNF_FILE{"REFERENCE_SEQ"}; } 
 if (defined $VSCNF_FILE{"REFERENCE_INDIST"}) { $reference_indist = $VSCNF_FILE{"REFERENCE_INDIST"}; } 
 
 return ($SERVER,$TMPDIR,$html_dir,$fasta_folder,$output_root,$reference_file,$reference_seq,$reference_indist);
}

################################################################
sub set_config {

 my $argc = @_;
 if ($argc <8) { die "set_config: expects 9 parameters, passed $argc !\n"; }

 $VSCNF_FILE{"SERVER"}=$_[0];
 $VSCNF_FILE{"TMPDIR"}=$_[1]; 
 $VSCNF_FILE{"HTML_DIR"}=$_[2]; 
 $VSCNF_FILE{"FASTA_DIR"}=$_[3]; 
 $VSCNF_FILE{"OUTPUT_ROOT"}=$_[4]; 
 $VSCNF_FILE{"REFERENCE_FILE"}=$_[5]; 
 $VSCNF_FILE{"REFERENCE_SEQ"}=$_[6]; 
 $VSCNF_FILE{"REFERENCE_INDIST"}=$_[7]; 
}

################################################################
sub get_config_vars {

 my $NPROCESSES=-1; 
 my $strip_454_keytags=-1; 
 my $is_paired_reads=-1; 
 my $reference_indist_produce=-1; 
 my $MIN_SUPPORT_REQUIRED=-1; 
 my $MIN_FLANK_REQUIRED=-1;
 my $MAX_FLANK_CONSIDERED=-1;
 my $REFS_TOTAL=-1;

 my $argc = @_;
 if ($argc <1) { die "get_config_vars: expects 1 parameters, passed $argc !\n"; }

 my $startdir = $_[0];

 if (!$VSREAD) { read_config_file($startdir); }

 if (defined $VSCNF_FILE{"NPROCESSES"}) { $NPROCESSES = $VSCNF_FILE{"NPROCESSES"}; } 
 if (defined $VSCNF_FILE{"STRIP_454_KEYTAGS"}) { $strip_454_keytags = $VSCNF_FILE{"STRIP_454_KEYTAGS"}; } 
 if (defined $VSCNF_FILE{"IS_PAIRED_READS"}) { $is_paired_reads = $VSCNF_FILE{"IS_PAIRED_READS"}; } 
 if (defined $VSCNF_FILE{"REFERENCE_INDIST_PRODUCE"}) { $reference_indist_produce = $VSCNF_FILE{"REFERENCE_INDIST_PRODUCE"}; } 
 if (defined $VSCNF_FILE{"MIN_SUPPORT_REQUIRED"}) { $MIN_SUPPORT_REQUIRED = $VSCNF_FILE{"MIN_SUPPORT_REQUIRED"}; } 
 if (defined $VSCNF_FILE{"MIN_FLANK_REQUIRED"}) { $MIN_FLANK_REQUIRED = $VSCNF_FILE{"MIN_FLANK_REQUIRED"}; } 
 if (defined $VSCNF_FILE{"MAX_FLANK_CONSIDERED"}) { $MAX_FLANK_CONSIDERED = $VSCNF_FILE{"MAX_FLANK_CONSIDERED"}; } 
 if (defined $VSCNF_FILE{"REFS_TOTAL"}) { $REFS_TOTAL = $VSCNF_FILE{"REFS_TOTAL"}; } 
 
 return ($NPROCESSES,$strip_454_keytags,$is_paired_reads,$reference_indist_produce,$MIN_SUPPORT_REQUIRED,$MIN_FLANK_REQUIRED,$MAX_FLANK_CONSIDERED,$REFS_TOTAL);
}

################################################################
sub set_config_vars {

 my $argc = @_;
 if ($argc <8) { die "set_config_vars: expects 8 parameters, passed $argc !\n"; }

 $VSCNF_FILE{"NPROCESSES"}=$_[0]; 
 $VSCNF_FILE{"STRIP_454_KEYTAGS"}=$_[1]; 
 $VSCNF_FILE{"IS_PAIRED_READS"}=$_[2]; 
 $VSCNF_FILE{"REFERENCE_INDIST_PRODUCE"}=$_[3]; 
 $VSCNF_FILE{"MIN_SUPPORT_REQUIRED"}=$_[4]; 
 $VSCNF_FILE{"MIN_FLANK_REQUIRED"}=$_[5]; 
 $VSCNF_FILE{"MAX_FLANK_CONSIDERED"}=$_[6]; 
 $VSCNF_FILE{"REFS_TOTAL"}=$_[7]; 
}

################################################################
sub get_credentials {

 my $LOGIN = "";
 my $PASS = "";
 my $HOST = "";

 my $argc = @_;
 if ($argc <1) { die "get_credentials: expects 1 parameters, passed $argc !\n"; }

 my $startdir = $_[0];

 if (!$VSREAD) { read_config_file($startdir); }
 
 if (defined $VSCNF_FILE{"LOGIN"}) { $LOGIN = $VSCNF_FILE{"LOGIN"}; } 
 if (defined $VSCNF_FILE{"PASS"})  { $PASS = $VSCNF_FILE{"PASS"}; } 
 if (defined $VSCNF_FILE{"HOST"})  { $HOST = $VSCNF_FILE{"HOST"}; } 
 
 return ($LOGIN, $PASS, $HOST);
}

################################################################
sub set_credentials {

 my $argc = @_;
 if ($argc <3) { die "set_credentials: expects 3 parameters, passed $argc !\n"; }

 $VSCNF_FILE{"LOGIN"}=$_[0]; 
 $VSCNF_FILE{"PASS"}=$_[1]; 
 $VSCNF_FILE{"HOST"}=$_[2]; 
}

################################################################

sub set_datetime {

 my $argc = @_;
 if ($argc <5) { die "set_datetime: expects 5 parameters, passed $argc !\n"; }

 my $DBNAME = $_[0];
 my $LOGIN = $_[1];
 my $PASS = $_[2];
 my $HOST = $_[3];
 my $NAME = $_[4];

 my $dbh = DBI->connect("DBI:mysql:$DBNAME;host=$HOST", "$LOGIN", "$PASS"
                   ) || die "Could not connect to database: $DBI::errstr"; 

 my $sth = $dbh->prepare("UPDATE stats SET $NAME=now()")
                or die "Couldn't prepare statement: " . $dbh->errstr;

 $sth->execute()             # Execute the query
            or die "Couldn't execute statement: " . $sth->errstr;

 $sth->finish;
 $dbh->disconnect();

 return 0;

}

################################################################

sub stats_set {

 my $argc = @_;
 if ($argc <6) { die "stats_set: expects 6 parameters, passed $argc !\n"; }

 my $DBNAME = $_[0];
 my $LOGIN = $_[1];
 my $PASS = $_[2];
 my $HOST = $_[3];
 my $NAME = $_[4];
 my $VALUE = $_[5];

 my $dbh = DBI->connect("DBI:mysql:$DBNAME;host=$HOST", "$LOGIN", "$PASS"
                   ) || die "Could not connect to database: $DBI::errstr"; 

 my $sth = $dbh->prepare("UPDATE stats SET $NAME=?")
                or die "Couldn't prepare statement: " . $dbh->errstr;

 $sth->execute($VALUE)             # Execute the query
            or die "Couldn't execute statement: " . $sth->errstr;

 $sth->finish;
 $dbh->disconnect();

 return 0;
}

################################################################

sub stats_get {

 my $argc = @_;
 if ($argc <5) { die "stats_set: expects 5 parameters, passed $argc !\n"; }

 my $DBNAME = $_[0];
 my $LOGIN = $_[1];
 my $PASS = $_[2];
 my $HOST = $_[3];
 my $NAME = $_[4];
 my $VALUE = undef;

 my $dbh = DBI->connect("DBI:mysql:INFORMATION_SCHEMA;host=$HOST", "$LOGIN", "$PASS"
                   ) || die "Could not connect to database: $DBI::errstr"; 

 # check if database exists first, return undef if not
 my $sth = $dbh->prepare("SELECT SCHEMA_NAME FROM INFORMATION_SCHEMA.SCHEMATA WHERE SCHEMA_NAME = '$DBNAME';")
                or die "Couldn't prepare statement: " . $dbh->errstr;
 $sth->execute()             # Execute the query
            or die "Couldn't execute statement: " . $sth->errstr;
 if ($sth->rows == 0) {
   $sth->finish;
   $dbh->disconnect();
   return undef;
 }
 $sth->finish;

 # get the namve/value pair
 my $sth = $dbh->prepare("SELECT $NAME FROM ${DBNAME}.stats;")
                or die "Couldn't prepare statement: " . $dbh->errstr;

 $sth->execute()             # Execute the query
            or die "Couldn't execute statement: " . $sth->errstr;

 my @data = $sth->fetchrow_array();

 if ($sth->rows == 0) {
            print STDERR "No field in database  stats.`$NAME'. Aborting!\n\n";
            exit(1);
  }

 $VALUE = $data[0];
 if (!defined $VALUE) { $VALUE=""; }

 $sth->finish;
 $dbh->disconnect();

 return $VALUE;
}

################################################################

sub write_mysql {

 my $argc = @_;
 if ($argc <2) { die "stats_set: expects 2 parameter, passed $argc !\n"; }

 my $DBNAME = $_[0];
 my $TMP = $_[1];

open FILE, ">$TMP/${DBNAME}.sql" or die $!;


print FILE "CREATE database IF NOT EXISTS ${DBNAME};\n\n";
print FILE "USE ${DBNAME};\n\n";


print FILE <<TEST;

drop table IF EXISTS vntr_support;
CREATE TABLE
vntr_support (
`refid` INT(11) NOT NULL,
`copies` INT(11) NOT NULL,
`sameasref` INT(11) NOT NULL,
`support` INT(11) NOT NULL DEFAULT 0,
`copiesfloat` FLOAT(11) NOT NULL,
`representative` INT(11) NULL,
 PRIMARY KEY (refid,copies)
) ENGINE=INNODB;

drop table IF EXISTS fasta_ref_reps;
CREATE TABLE
fasta_ref_reps (
`rid` INT(11) NOT NULL PRIMARY KEY, 

`firstindex` INT(11) NOT NULL, 
`lastindex` INT(11) NOT NULL, 
`copynum` FLOAT(11) NOT NULL, 
`head` VARCHAR(100) NULL,

`flankleft` VARCHAR(8000) NULL,
`pattern` VARCHAR(5000) NOT NULL,
`sequence` TEXT NOT NULL,
`flankright` VARCHAR(8000) NULL,
`conserved` FLOAT(11) NULL,
`comment` VARCHAR(500) NULL,
`flank_disting` INT(11) NULL,
`entropy` FLOAT(11) NOT NULL, 
`has_support` INT(11) NULL,
`span1` INT(11) NULL,
`spanN` INT(11) NULL,
`is_singleton` INT(11) NOT NULL DEFAULT 0,
`is_dist` INT(11) NOT NULL DEFAULT 0,
`is_indist` INT(11) NOT NULL DEFAULT 0,
`homez_same` INT(11) NULL,
`homez_diff` INT(11) NULL,
`hetez_same` INT(11) NULL,
`hetez_diff` INT(11) NULL,
`hetez_multi` INT(11) NULL,
`support_vntr` INT(11) NULL,
`support_vntr_span1` INT(11) NULL,
`alleles_sup` INT(11) NULL,
`allele_sup_same_as_ref` INT(11) NULL,
INDEX(`head`)
) ENGINE=INNODB;
CREATE UNIQUE INDEX comment_index on fasta_ref_reps(rid,comment);

drop table IF EXISTS fasta_reads;
CREATE TABLE 
fasta_reads (
`sid` BIGINT(11) UNSIGNED NOT NULL AUTO_INCREMENT PRIMARY KEY, 
`head` CHAR(100) NOT NULL,
`dna` VARCHAR(8000) NULL,
`qual` VARCHAR(8000) NULL,
 UNIQUE (`head`)
) ENGINE=INNODB;


drop table IF EXISTS replnk;
CREATE TABLE 
replnk (
`rid` INT(11) NOT NULL PRIMARY KEY, 
`sid` BIGINT(11) UNSIGNED NOT NULL,
`first` INT(11) NOT NULL,
`last` INT(11) NOT NULL,
`patsize` INT(11) NOT NULL,
`copynum`  FLOAT(11) NOT NULL,
`pattern` VARCHAR(5000) NOT NULL,
`profile` VARCHAR(8000) NULL,
`profilerc` VARCHAR(8000) NULL,
`profsize` INT(11) NULL,
 INDEX(`sid`)
) ENGINE=INNODB;

drop table IF EXISTS clusters;
CREATE TABLE
clusters (
`cid` INT(11) NOT NULL PRIMARY KEY,
`minpat` INT(11) NOT NULL,
`maxpat` INT(11) NOT NULL,

`refs_flank_undist` INT(11) NULL,
`refs_flank_dist` INT(11) NULL,
`refs_flank_undist_l` INT(11) NULL,
`refs_flank_undist_r` INT(11) NULL,
`refs_flank_undist_lr` INT(11) NULL,

`repeatcount` INT(11) NOT NULL,
`refcount` INT(11) NOT NULL,
`variability` INT(11) NOT NULL DEFAULT 0,
`assemblyreq` INT(11) NULL,
`profdensity` FLOAT(11) NULL,
`flankdensity` FLOAT(11) NULL,
`mcpattern` VARCHAR(5000) NULL,
`aveentropy` FLOAT(11) NULL
) ENGINE=INNODB;


drop table IF EXISTS clusterlnk;
CREATE TABLE
clusterlnk (
`clusterid` INT(11) NOT NULL,
`repeatid` INT(11) NOT NULL,
`direction` CHAR(1) NOT NULL,
`reserved` INT(11) NOT NULL,
`reserved2` INT(11) NOT NULL,
 PRIMARY KEY (clusterid,repeatid)
) ENGINE=INNODB;
CREATE UNIQUE INDEX repeat_index on clusterlnk(repeatid);


drop table IF EXISTS map;
CREATE TABLE
map (
`refid` INT(11) NOT NULL,
`readid` INT(11) NOT NULL,
`reserved` INT(11) NOT NULL,
`reserved2` INT(11) NOT NULL,
`bbb` tinyint(3) NOT NULL DEFAULT 0,
 PRIMARY KEY (refid,readid)
) ENGINE=INNODB;
CREATE INDEX read_index ON map (readid);

drop table IF EXISTS rank;
CREATE TABLE
rank (
`refid` INT(11) NOT NULL,
`readid` INT(11) NOT NULL,
`score` FLOAT(11) NULL,
`ties` smallint(11) NOT NULL DEFAULT 0,
`refdir` CHAR(1) NOT NULL,
 PRIMARY KEY (refid,readid)
) ENGINE=INNODB;

drop table IF EXISTS rankflank;
CREATE TABLE
rankflank (
`refid` INT(11) NOT NULL,
`readid` INT(11) NOT NULL,
`score` FLOAT(11) NULL,
`ties` smallint(11) NOT NULL DEFAULT 0,
 PRIMARY KEY (refid,readid)
) ENGINE=INNODB;


drop table IF EXISTS flank_connection;
CREATE TABLE `flank_connection` (
  `refid` int(11) NOT NULL,
  `clusterid` int(11) NOT NULL,
  `fcomponentsize` int(11) NOT NULL,
  `fcomponentid` int(11) NOT NULL,
  `paramsetid` int(11) NOT NULL,
  PRIMARY KEY  (`refid`,`paramsetid`),
  INDEX `paramsetid` (`paramsetid`)
) ENGINE=InnoDB;


drop table IF EXISTS flank_params;
CREATE TABLE `flank_params` (
  `paramsetid` int(11) NOT NULL auto_increment,
  `flength` int(11) NOT NULL,
  `ferrors` int(11) NOT NULL,
  PRIMARY KEY  (`paramsetid`)
) ENGINE=InnoDB;

drop table IF EXISTS stats;
CREATE TABLE
stats (
`id` INT(11) NOT NULL PRIMARY KEY,

`MAP_ROOT` VARCHAR(500) NULL,

`PARAM_TRF` VARCHAR(500) NULL,
`PARAM_PROCLU` VARCHAR(500) NULL,

`FOLDER_FASTA` VARCHAR(500) NULL,
`FOLDER_PROFILES` VARCHAR(500) NULL,
`FOLDER_PROFILES_CLEAN` VARCHAR(500) NULL,
`FOLDER_REFERENCE` VARCHAR(500) NULL,

`FILE_REFERENCE_LEB` VARCHAR(500) NULL,
`FILE_REFERENCE_SEQ` VARCHAR(500) NULL,

`NUMBER_READS` BIGINT NULL,
`NUMBER_TRS_IN_READS` BIGINT NULL,

`NUMBER_TRS_IN_READS_GE7` BIGINT NULL,
`NUMBER_READS_WITHTRS` BIGINT NULL,
`NUMBER_READS_WITHTRS_GE7` BIGINT NULL,
`NUMBER_READS_WITHTRS_GE7_AFTER_REDUND` BIGINT NULL,

`NUMBER_TRS_IN_READS_AFTER_REDUND` BIGINT NULL,
`NUMBER_REF_TRS` BIGINT NULL,
`NUMBER_REFS_TRS_AFTER_REDUND` BIGINT NULL,

`CLUST_NUMBER_OF_PROCLU_CLUSTERS` BIGINT NULL,
`CLUST_NUMBER_OF_PROCLU_CLUSTERS_BEFORE_REJOIN` BIGINT NULL,
`CLUST_NUMBER_OF_EXACTPAT_CLUSTERS` BIGINT NULL,
`CLUST_NUMBER_OF_REF_REPS_IN_CLUSTERS` BIGINT NULL,
`CLUST_NUMBER_OF_READ_REPS_IN_CLUSTERS` BIGINT NULL,


`CLUST_LARGEST_NUMBER_OF_TRS_IN_PROCLU_CLUSTER` BIGINT NULL,
`CLUST_LARGEST_NUMBER_OF_REFS_IN_PROCLU_CLUSTER` BIGINT NULL,
`CLUST_LARGEST_PATRANGE_IN_PROCLU_CLUSTER` BIGINT NULL,
`CLUST_LARGEST_NUMBER_OF_TRS_IN_EXACTPAT_CLUSTER` BIGINT NULL,
`CLUST_LARGEST_NUMBER_OF_REFS_IN_EXACTPAT_CLUSTER` BIGINT NULL,

`CLUST_NUMBER_OF_REFS_WITH_PREDICTED_VNTR` BIGINT NULL,
`CLUST_NUMBER_OF_CLUSTERS_WITH_PREDICTED_VNTR` BIGINT NULL,
`NUMBER_REFS_VNTR_SPAN_N` BIGINT NULL,

`NUMBER_REFS_SINGLE_REF_CLUSTER` BIGINT NULL,
`NUMBER_REFS_SINGLE_REF_CLUSTER_WITH_READS_MAPPED` BIGINT NULL,
`NUMBER_REFS_SINGLE_REF_CLUSTER_WITH_NO_READS_MAPPED` BIGINT NULL,

`NUMBER_MAPPED` BIGINT NULL,
`NUMBER_RANK` BIGINT NULL,
`NUMBER_RANKFLANK` BIGINT NULL,
`INTERSECT_RANK_AND_RANKFLANK` BIGINT NULL,
`INTERSECT_RANK_AND_RANKFLANK_BEFORE_PCR` BIGINT NULL,

`BBB_WITH_MAP_DUPS` BIGINT NULL,
`BBB` BIGINT NULL,

`RANK_EDGES_OVERCUTOFF` BIGINT NULL,
`RANK_REMOVED_SAMEREF` BIGINT NULL,
`RANK_REMOVED_SAMESEQ` BIGINT NULL,
`RANK_REMOVED_PCRDUP` BIGINT NULL,

`RANKFLANK_EDGES_INSERTED` BIGINT NULL,
`RANKFLANK_REMOVED_SAMEREF` BIGINT NULL,
`RANKFLANK_REMOVED_SAMESEQ` BIGINT NULL,
`RANKFLANK_REMOVED_PCRDUP` BIGINT NULL,


`TIME_MYSQLCREATE` INT(11) NULL,
`TIME_TRF` INT(11) NULL,
`TIME_RENUMB` INT(11) NULL,
`TIME_REDUND` INT(11) NULL,
`TIME_PROCLU` INT(11) NULL,
`TIME_JOINCLUST` INT(11) NULL,
`TIME_DB_INSTERT_REFS` INT(11) NULL,
`TIME_DB_INSTERT_READS` INT(11) NULL,
`TIME_WRITE_FLANKS` INT(11) NULL,
`TIME_MAP_FLANKS` INT(11) NULL,
`TIME_MAP_REFFLANKS` INT(11) NULL,
`TIME_MAP_INSERT` INT(11) NULL,
`TIME_EDGES` INT(11) NULL,
`TIME_INDEX_PCR` INT(11) NULL,
`TIME_PCR_DUP` INT(11) NULL,
`TIME_MAP_DUP` INT(11) NULL,
`TIME_VNTR_PREDICT` INT(11) NULL,
`TIME_ASSEMBLYREQ` INT(11) NULL,
`TIME_REPORTS` INT(11) NULL,

`DATE_MYSQLCREATE` DATETIME NULL,
`DATE_TRF` DATETIME NULL,
`DATE_RENUMB` DATETIME NULL,
`DATE_REDUND` DATETIME NULL,
`DATE_PROCLU` DATETIME NULL,
`DATE_JOINCLUST` DATETIME NULL,
`DATE_DB_INSTERT_REFS` DATETIME NULL,
`DATE_DB_INSTERT_READS` DATETIME NULL,
`DATE_WRITE_FLANKS` DATETIME NULL,
`DATE_MAP_FLANKS` DATETIME NULL,
`DATE_MAP_REFFLANKS` DATETIME NULL,
`DATE_MAP_INSERT` DATETIME NULL,
`DATE_EDGES` DATETIME NULL,
`DATE_INDEX_PCR` DATETIME NULL,
`DATE_PCR_DUP` DATETIME NULL,
`DATE_MAP_DUP` DATETIME NULL,
`DATE_VNTR_PREDICT` DATETIME NULL,
`DATE_ASSEMBLYREQ` DATETIME NULL,
`DATE_REPORTS` DATETIME NULL,

`ERROR_STEP` INT(11) NOT NULL DEFAULT 0,
`ERROR_DESC` VARCHAR(500) NOT NULL DEFAULT '',
`ERROR_CODE` INT(11) NOT NULL DEFAULT 0,

`N_MIN_SUPPORT` INT(11) NOT NULL DEFAULT 0,
`MIN_FLANK_REQUIRED` INT(11) NOT NULL DEFAULT 0,
`MAX_FLANK_CONSIDERED` INT(11) NOT NULL DEFAULT 0
);

INSERT INTO stats (id) VALUES(1);


TEST


 close(FILE);

 return 0;
}
################################################################

sub print_config {

 my $argc = @_;
 if ($argc <1) { die "print_config: expects 1 parameters, passed $argc!\n"; }

 my $startdir = $_[0];

 if (!defined $VSCNF_FILE{"LOGIN"}) { $VSCNF_FILE{"LOGIN"}=""; } 
 if (!defined $VSCNF_FILE{"PASS"})  {  $VSCNF_FILE{"PASS"}=""; } 
 if (!defined $VSCNF_FILE{"HOST"})  {  $VSCNF_FILE{"HOST"}="localhost"; } 
 if (!defined $VSCNF_FILE{"NPROCESSES"})  {  $VSCNF_FILE{"NPROCESSES"}=-1; } 
 if (!defined $VSCNF_FILE{"MIN_FLANK_REQUIRED"})  {  $VSCNF_FILE{"MIN_FLANK_REQUIRED"}=-1; } 
 if (!defined $VSCNF_FILE{"MAX_FLANK_CONSIDERED"})  {  $VSCNF_FILE{"MAX_FLANK_CONSIDERED"}=-1; } 
 if (!defined $VSCNF_FILE{"MIN_SUPPORT_REQUIRED"})  {  $VSCNF_FILE{"MIN_SUPPORT_REQUIRED"}=-1; } 

 if (!defined $VSCNF_FILE{"SERVER"})  {  $VSCNF_FILE{"SERVER"}=""; } 
 if (!defined $VSCNF_FILE{"STRIP_454_KEYTAGS"})  {  $VSCNF_FILE{"STRIP_454_KEYTAGS"}=-1; } 
 if (!defined $VSCNF_FILE{"IS_PAIRED_READS"})  {  $VSCNF_FILE{"IS_PAIRED_READS"}=-1; } 
 if (!defined $VSCNF_FILE{"HTML_DIR"})  {  $VSCNF_FILE{"HTML_DIR"}=""; } 
 if (!defined $VSCNF_FILE{"FASTA_DIR"})  {  $VSCNF_FILE{"FASTA_DIR"}=""; } 
 if (!defined $VSCNF_FILE{"OUTPUT_ROOT"})  {  $VSCNF_FILE{"OUTPUT_ROOT"}=""; } 
 if (!defined $VSCNF_FILE{"TMPDIR"})  {  $VSCNF_FILE{"TMPDIR"}=""; } 
 if (!defined $VSCNF_FILE{"REFERENCE_FILE"})  {  $VSCNF_FILE{"REFERENCE_FILE"}=""; } 
 if (!defined $VSCNF_FILE{"REFERENCE_SEQ"})  {  $VSCNF_FILE{"REFERENCE_SEQ"}=""; } 
 if (!defined $VSCNF_FILE{"REFERENCE_INDIST"})  {  $VSCNF_FILE{"REFERENCE_INDIST"}=""; } 
 if (!defined $VSCNF_FILE{"REFERENCE_INDIST_PRODUCE"})  {  $VSCNF_FILE{"REFERENCE_INDIST_PRODUCE"}=-1; } 
 if (!defined $VSCNF_FILE{"REFS_TOTAL"})  {  $VSCNF_FILE{"REFS_TOTAL"}=-1; } 


 # look in the directory the script was started in
 if (open(MFREAD,">${startdir}vs.cnf")) {

print MFREAD <<CNF;
# COPY THIS FILE TO A DIFFERENT LOCATION AND SET ALL VARIABLES. 
# DO NOT FORGET TO CHMOD THIS FILE TO PREVENT OTHER PEOPLE ON 
# THE SYSTEM FROM LEARNING YOUR MYSQL CREDENTIALS.

# mysql credentials
LOGIN=$VSCNF_FILE{"LOGIN"}
PASS=$VSCNF_FILE{"PASS"}
HOST=$VSCNF_FILE{"HOST"}

# set this to the number of processors on your system 
# (or less if sharing the system with others or RAM is limited)
# eg, 8
NPROCESSES=$VSCNF_FILE{"NPROCESSES"}

# minimum required flank on both sides for a read TR to be considered
# eg, 10
MIN_FLANK_REQUIRED=$VSCNF_FILE{"MIN_FLANK_REQUIRED"}

# maximum flank length used in flank alignments
# set to big number to use all
# if read flanks are long with a lot of errors, 
# it might be useful to set this to something like 50
# max number of errors per flank is currently set to 8 (can be changed in main script only)
# eg, 1000
MAX_FLANK_CONSIDERED=$VSCNF_FILE{"MAX_FLANK_CONSIDERED"}

# minimum number of mapped reads which agree on copy number to call an allele
# eg, 2
MIN_SUPPORT_REQUIRED=$VSCNF_FILE{"MIN_SUPPORT_REQUIRED"}

# server name, used for html generating links
# eg, orca.bu.edu
SERVER=$VSCNF_FILE{"SERVER"}

# for 454 platform, strip leading 'TCAG' 
# eg, 1 - yes
# eg, 0 - no (use no for all other platforms)
STRIP_454_KEYTAGS=$VSCNF_FILE{"STRIP_454_KEYTAGS"}

# data is paired reads
# eg, 0 = no 
# eg, 1 - yes
IS_PAIRED_READS=$VSCNF_FILE{"IS_PAIRED_READS"}

# html directory (must be writable and executable!)
# eg, /var/www/html/vntrview
HTML_DIR=$VSCNF_FILE{"HTML_DIR"}

# input data directory 
# (plain or gzipped fasta/fastq files)
# eg, /input
FASTA_DIR=$VSCNF_FILE{"FASTA_DIR"}

# output directory (must be writable and executable!)
# eg, /output
OUTPUT_ROOT=$VSCNF_FILE{"OUTPUT_ROOT"}

# temp (scratch) directory (must be executable!)
# eg, /tmp
TMPDIR=$VSCNF_FILE{"TMPDIR"}

# names for the reference files 

# (leb36 file, sequence plus flank data file, indistinguishable references file) 
# files must be in install directory

# eg, reference.leb36
REFERENCE_FILE=$VSCNF_FILE{"REFERENCE_FILE"} 

# eg, reference.seq
REFERENCE_SEQ=$VSCNF_FILE{"REFERENCE_SEQ"} 

# this file can be generated bu setting reference_indist_produce to 1
# eg, reference.indist 
REFERENCE_INDIST=$VSCNF_FILE{"REFERENCE_INDIST"}

# generate a file of indistinguishable references, 
# necessary only if a file is not already available for the reference set
# eg, 1- generate
# eg, 0 - don't generate
REFERENCE_INDIST_PRODUCE=$VSCNF_FILE{"REFERENCE_INDIST_PRODUCE"}

# total number of reference TRs prior to filtering
# this is a fixed number to be printed in the latex file
# set to 0 if it is not applicable
# eg, 1188939 - human
REFS_TOTAL=$VSCNF_FILE{"REFS_TOTAL"} 

CNF

    close(MFREAD);

chmod 0600, "${startdir}vs.cnf";

 } else {

   die "print_config: can't open '${startdir}vs.cnf' for writing!\n";
 }

}




####################################################################################



1;


