package vutil;
use strict;
use Cwd;
use DBI;
use Carp;
use Config::Simple;

use FindBin;

use base 'Exporter';
our @EXPORT_OK = qw(read_config_file get_config get_credentials set_config set_credentials get_dbh write_mysql write_sqlite stats_set set_statistics stats_get set_datetime print_config trim create_blank_file get_trunc_query);

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

sub create_blank_file {
  my $filename = shift or croak("Error: filename is a required argument\n");
  open my $blank_file, ">", $filename or croak("Error creating blank file $filename: $!");
  close $blank_file or croak("Error closing blank file $filename: $!");
}

################################################################
sub read_config_file_configsimple {
  my $output_folder = shift
    or croak "Error: function expects 1 parameter (got none)\n";
  my $config_ref;
  carp "$output_folder/vs.cnf";
  my $cfg = new Config::Simple( syntax => "ini" );
  $cfg->read("$output_folder/vs.cnf");
  $config_ref = $cfg->get_block('default');
  $cfg->write("$output_folder/vs.cnf.bak");
  return $config_ref;
}

################################################################
sub read_config_file {

 # Get file location
 unless (@_) { die "read_config_file: expects 1 parameters.\n"; }
 my $file_loc = shift;
 
 # read global config
 if (open(my $cnf, "<", "$file_loc")) {

  while (<$cnf>) {
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

  close($cnf);
  return 1;

 } else {
  return 0;
 }

}

################################################################
sub get_config {
  croak "Error: function expects 1 parameter (got none)\n"
    unless (@_ == 1);
  my $installdir = "$FindBin::RealBin";
  unless ($VSREAD) {
    # Must read global file first. Sets up the defaults.
    warn "Could not read global config\n"
      unless read_config_file("$installdir/vs.cnf");
    my $config_loc = shift;
    warn "Could not read run config (harmless if this is a new run)\n"
      unless read_config_file($config_loc);
    # Set VSREAD;
    $VSREAD = 1;
  }

  return %VSCNF_FILE;
}

################################################################
sub set_config {
  # TODO Option validation (trusting caller sends hash with valid options)
  my %in_hash = @_;

  %VSCNF_FILE = %in_hash;
}

################################################################
sub get_credentials {

 my $LOGIN = "";
 my $PASS = "";
 my $HOST = "";

 my $argc = @_;
 if ($argc <1) { die "get_credentials: expects 1 parameters, passed $argc !\n"; }

 my $startdir = $_[0];

 if (!$VSREAD) { read_config_file($startdir . "vs.cnf"); }
 
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

 my $dbh = get_dbh($DBNAME, $ENV{HOME} . "/". $VSCNF_FILE{DBSUFFIX} . ".vs.cnf");

 my $sth = $dbh->prepare("UPDATE stats SET $NAME=?")
                or die "Couldn't prepare statement: " . $dbh->errstr;

 $sth->execute($VALUE)             # Execute the query
            or die "Couldn't execute statement: " . $sth->errstr;

 $sth->finish;
 $dbh->disconnect();

 return 0;
}

####################################
sub set_statistics {

  my $argc = @_;
  if ( $argc < 3 ) {
      die "stats_set: expects 3 parameters, passed $argc !\n";
  }

  my ($DBSUFFIX, $NAME, $VALUE)  = @_;
  my $dbh = get_dbh($DBSUFFIX, $ENV{HOME} . "/". $DBSUFFIX . ".vs.cnf");

  if ($ENV{DEBUG}) {
    warn "Setting stat: $NAME to $VALUE\n";
  }

  my $sth = $dbh->prepare("UPDATE stats SET $NAME=?")
    or croak "Couldn't prepare statement: " . $dbh->errstr;

  $sth->execute($VALUE)             # Execute the query
    or croak "Couldn't execute statement: " . $sth->errstr;

  $dbh->disconnect();
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

sub get_dbh {
  unless (@_ == 2) {
    carp "get_dbh requires 2 parameters: database suffix and the config file location.";
    return undef;
  }

  my ($dbsuffix, $config_loc) = @_;
  read_config_file($config_loc)
    unless ($VSREAD);

  my $dbh;
  if ($VSCNF_FILE{BACKEND} eq "sqlite") {
    my $dbfile = "$VSCNF_FILE{OUTPUT_ROOT}/vntr_$dbsuffix/$dbsuffix.db";
    if ($ENV{DEBUG}) {
      warn "Using SQLite db at: $dbfile\n";
    }
    $dbh = DBI->connect("DBI:SQLite:dbname=$dbfile", undef, undef, {
      AutoCommit => 1,
      RaiseError => 1,
      sqlite_see_if_its_a_number => 1,
    })
    or croak "Could not connect to database: $DBI::errstr";
  }
  else {
    $dbh = DBI->connect( "DBI:mysql:VNTRPIPE_$dbsuffix;mysql_local_infile=1;host=$VSCNF_FILE{HOST}",
    "$VSCNF_FILE{LOGIN}", "$VSCNF_FILE{PASS}" )
      or croak "Could not connect to database: $DBI::errstr";
  }

  return $dbh;
}

####################################
sub get_trunc_query {
    die "Error: need table name for truncate query\n"
        unless @_ == 2;
    my ($backend, $table) = @_;
    my $trunc_query;
    if ( $backend eq "sqlite" ) {
        $trunc_query = qq{DELETE FROM $table};
    }
    elsif ( $backend eq "mysql" ) {
        $trunc_query = qq{TRUNCATE TABLE $table};
    }

    return $trunc_query;
}

################################################################

sub write_sqlite {
  unless (@_ == 3) {
    carp "write_sqlite requires 3 parameters: database suffix, the run output directory, and the path to the run configuration file.";
    return undef;
  }

  my ($dbsuffix, $output_dir, $config_loc) = @_;
  my $installdir = "$FindBin::RealBin";
  my $exestring = "sqlite3 $output_dir/vntr_$dbsuffix/$dbsuffix.db < $installdir/sqlite_schema.sql";
  warn "Executing: $exestring\n";
  system($exestring);
  # read_config_file($config_loc)
  #   unless ($VSREAD);
  # open my $schema_fh, "<", "$installdir/sqlite_schema.sql"
  #   or croak "Error opening SQLite schema file '$installdir/sqlite_schema.sql': $?";
  # my $sqlite_schema;
  # while (<$schema_fh>) {
  #   # chomp;
  #   $sqlite_schema .= $_;
  # }
  # close $schema_fh;

  # my $dbh = get_dbh($dbsuffix, $config_loc)
  #   or croak "Could not connect to database: $DBI::errstr";
  # $dbh->do($sqlite_schema)
  #   or croak "Error creating tables in database: $DBI::errstr";
  # $dbh->commit;
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
CREATE INDEX read_index ON vntr_support (representative);

drop table IF EXISTS fasta_ref_reps;
CREATE TABLE
fasta_ref_reps (
`rid` INT(11) NOT NULL PRIMARY KEY, 

`firstindex` INT(11) NOT NULL, 
`lastindex` INT(11) NOT NULL, 
`copynum` FLOAT(11) NOT NULL, 
`head` VARCHAR(100) NULL,

`flankleft` TEXT(8000) NULL,
`pattern` TEXT(5000) NOT NULL,
`sequence` TEXT NOT NULL,
`flankright` TEXT(8000) NULL,
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
CREATE UNIQUE INDEX comment_index on fasta_ref_reps(rid,comment(100));


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
`pattern` TEXT(5000) NOT NULL,
`profile` TEXT(8000) NULL,
`profilerc` TEXT(8000) NULL,
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

# Database backend
BACKEND=$VSCNF_FILE{"BACKEND"}

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


