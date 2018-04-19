package vutil;
use strict;
use Cwd;
use DBI;
use Carp;
use FindBin;
use File::Temp;
use POSIX qw(strftime);

if ( $ENV{DEBUG} ) {
    use Data::Dumper;
}

use base 'Exporter';
our @EXPORT_OK
    = qw(read_config_file get_config set_config set_credentials get_dbh get_ref_dbh write_mysql make_refseq_db load_refprofiles_db run_redund write_sqlite set_statistics get_statistics stats_get set_datetime print_config trim create_blank_file get_trunc_query);

# vutil.pm
# author: Yevgeniy Gelfand
# create date: Oct 30, 2010
# function: create mysql database for vntr pipleline,
# provide functions for database management

my %VSCNF_FILE = ();
my $VSREAD     = 0;

################################################################
sub trim {
    my $string = shift;
    $string =~ s/^\s+//;
    $string =~ s/\s+$//;
    return $string;
}

sub create_blank_file {
    my $filename = shift or croak("Error: filename is a required argument\n");
    open my $blank_file, ">", $filename
        or croak("Error creating blank file $filename: $!");
    close $blank_file or croak("Error closing blank file $filename: $!");
}

################################################################
sub read_config_file {

    # Get file location
    unless (@_) { die "read_config_file: expects 1 parameters.\n"; }
    my $file_loc = shift;

    ( $ENV{DEBUG} ) && warn "Reading configuration file: $file_loc\n";

    # read global config
    if ( open( my $cnf, "<", "$file_loc" ) ) {

        while (<$cnf>) {
            chomp;

            # skip start comments
            if (/^\#/) { next; }

            if (/(.+)=(.*)/) {
                my $key = trim($1);
                my $val = trim($2);
                $val =~ s/\s*\#.*//;    # strip end comments
                $VSCNF_FILE{ uc($key) } = $val;

                #print $key."=".$val."\n";
            }
        }

        close($cnf);
        return 1;

    }
    else {
        return 0;
    }

}

################################################################
sub get_config {
    croak "Error: function expects 2 parameters\n"
        unless ( @_ == 2 );

    my ( $dbsuffix, $config_loc ) = @_;
    my $installdir = "$FindBin::RealBin";
    unless ($VSREAD) {

        # Must read global file first. Sets up the defaults.
        warn "Could not read global config\n"
            unless read_config_file("$installdir/vs.cnf");
        warn "Could not read run config (harmless if this is a new run)\n"
            unless read_config_file($config_loc);

        # Set VSREAD;
        $VSREAD = 1;
        $VSCNF_FILE{DBSUFFIX} = $dbsuffix;
        ( $VSCNF_FILE{REFERENCE_DB} = $VSCNF_FILE{REFERENCE_SEQ} )
            =~ s/.seq$/.db/;
    }

    return %VSCNF_FILE;
}

################################################################
sub set_config {
    my %in_hash = @_;

    # Validation
    unless ( $in_hash{SERVER} ) {
        croak(
            "Please set machine name (SERVER) variable on the command line or in the configuration file.\n"
        );
    }

    if ( $in_hash{NPROCESSES} <= 0 ) {
        croak(
            "Please set number of processes to be used by the pipeline (NPROCESSES) variable on the command line or in the configuration file.\n "
        );
    }

    if ( $in_hash{MIN_FLANK_REQUIRED} <= 0 ) {
        croak(
            "Please set min flank required to be used by the pipeline (MIN_FLANK_REQUIRED) variable on the command line or in the configuration file.\n"
        );
    }

    if ( $in_hash{MAX_FLANK_CONSIDERED} <= 0 ) {
        croak(
            "Please set max flank required to be used by the pipeline (MAX_FLANK_CONSIDERED) variable on the command line or in the configuration file.\n"
        );
    }

    if ( $in_hash{MIN_SUPPORT_REQUIRED} <= 0 ) {
        croak(
            "Please set min support required to be used by the pipeline (MIN_SUPPORT_REQUIRED) variable on the command line or in the configuration file.\n"
        );
    }

    if ( $in_hash{STRIP_454_KEYTAGS} < 0 ) {
        croak(
            "Please set strip_454_keytags (strip_454_keytags) variable on the command line or in the configuration file.\n"
        );
    }

    if ( $in_hash{REFERENCE_INDIST_PRODUCE} < 0 ) {
        croak(
            "Please set reference_indist_produce to be used by the pipeline (reference_indist_produce) variable on the command line or in the configuration file.\n"
        );
    }

    if ( $in_hash{IS_PAIRED_READS} < 0 ) {
        croak(
            "Please set is_paired_reads (is_paired_reads) variable on the command line or in the configuration file.\n"
        );
    }

    if ( "" eq $in_hash{HTML_DIR} ) {
        croak(
            "Please set html directory (html_dir) variable on the command line or in the configuration file. ($in_hash{HTML_DIR})\n"
        );
    }
    if ( !( -e $in_hash{HTML_DIR} ) && !mkdir("$in_hash{HTML_DIR}") ) {
        croak("Could not create html_dir directory ($in_hash{HTML_DIR}).\n");
    }

    unless ( exists $in_hash{FASTA_DIR} && $in_hash{FASTA_DIR} ) {
        croak(
            "Please set fasta directory (fasta_dir) variable on the command line or in the configuration file.\n"
        );
    }

    if ( "" eq $in_hash{OUTPUT_ROOT} ) {
        croak(
            "Please set output directory (output_root) variable on the command line or in the configuration file.\n"
        );
    }

    if ( "" eq $in_hash{REFERENCE_FILE} ) {
        croak(
            "Please set reference .leb36 file (reference_file) variable on the command line or in the configuration file.\n"
        );
    }

    if ( "" eq $in_hash{REFERENCE_SEQ} ) {
        croak(
            "Please set reference .seq file (reference_seq) variable on the command line or in the configuration file.\n"
        );
    }

    if ( "" eq $in_hash{REFERENCE_INDIST} ) {
        croak(
            "Please set reference .indist file (reference_indist) variable on the command line or in the configuration file.\n"
        );
    }

    if ( "" eq $in_hash{TMPDIR} ) {
        croak(
            "Please set temporary directory (TMPDIR) variable on the command line or in the configuration file.\n"
        );
    }

    unless ( -e $in_hash{REFERENCE_FILE} ) {
        die("File '$in_hash{REFERENCE_FILE}' not found!");
    }
    unless ( -e $in_hash{REFERENCE_SEQ} ) {
        die("File '$in_hash{REFERENCE_SEQ}' not found!");
    }
    if (   !( -e $in_hash{REFERENCE_INDIST} )
        && !$in_hash{REFERENCE_INDIST_PRODUCE} )
    {
        die("File '$in_hash{REFERENCE_INDIST}' not found!");
    }
    unless ( -e $in_hash{FASTA_DIR} ) {
        die("Directory '$in_hash{FASTA_DIR}' not found!");
    }
    if ( !-e $in_hash{OUTPUT_ROOT} && !mkdir( $in_hash{OUTPUT_ROOT} ) ) {
        die("Error creating output root ($in_hash{OUTPUT_ROOT}). Please check the path or manually create this directory and try again."
        );
    }
    unless ( -e $in_hash{OUTPUT_ROOT} ) {
        die("Directory '$in_hash{OUTPUT_ROOT}' not found!");
    }
    unless ( -e $in_hash{TMPDIR} ) {
        die("Temporary directory '$in_hash{TMPDIR}' not found!");
    }
    unless ( -x $in_hash{HTML_DIR} ) {
        die("Directory '$in_hash{HTML_DIR}' not executable!");
    }
    unless ( -x $in_hash{OUTPUT_ROOT} ) {
        die("Directory '$in_hash{OUTPUT_ROOT}' not executable!");
    }
    unless ( -x $in_hash{TMPDIR} ) {
        die("Directory '$in_hash{TMPDIR}' not executable!");
    }

    unless ( -w $in_hash{HTML_DIR} ) {
        die("Directory '$in_hash{HTML_DIR}' not writable!");
    }
    unless ( -w $in_hash{OUTPUT_ROOT} ) {
        die("Directory '$in_hash{OUTPUT_ROOT}' not writable!");
    }
    unless ( -w $in_hash{TMPDIR} ) {
        die("Directory '$in_hash{TMPDIR}' not writable!");
    }

    %VSCNF_FILE = %in_hash;
    $VSREAD     = 1;
}

################################################################
sub set_credentials {

    my $argc = @_;
    if ( $argc < 3 ) {
        die "set_credentials: expects 3 parameters, passed $argc !\n";
    }

    $VSCNF_FILE{"LOGIN"} = $_[0];
    $VSCNF_FILE{"PASS"}  = $_[1];
    $VSCNF_FILE{"HOST"}  = $_[2];
}

####################################
sub set_statistics {

# my $argc = @_;
# if ( $argc < 2 ) {
#     croak "set_statistics: expects at least 2 parameters, passed $argc !\n";
# }

    my $stats = shift;
    my $dbh   = get_dbh();
    my ( $sql_clause, @sql_qual, @sql_bind );
    ( $ENV{DEBUG} ) && warn Dumper($stats) . "\n";

    while ( my ( $key, $val ) = each %$stats ) {
        ( $ENV{DEBUG} ) && warn "Setting stat: $key to $val\n";
        push @sql_qual, "$key=?";
        push @sql_bind, $val;
    }

    $sql_clause = join ",", @sql_qual;
    my $sth = $dbh->prepare("UPDATE stats SET $sql_clause")
        or croak "Couldn't prepare statement: " . $dbh->errstr;

    $sth->execute(@sql_bind)    # Execute the query
        or croak "Couldn't execute statement: " . $sth->errstr;

    $dbh->disconnect();
}

################################################################

sub set_datetime {

    my $argc = @_;
    if ( $argc < 1 ) {
        croak "set_datetime: expects 1 parameter, passed $argc !\n";
    }

    my $NAME = shift;
    my $VALUE = strftime( "%F %T", localtime );

    return set_statistics( { $NAME, $VALUE } );
}

####################################
sub get_statistics {

  # my $argc = @_;
  # if ( $argc < 2 ) {
  #     die "get_statistics: expects at least 2 parameters, passed $argc !\n";
  # }
    unless ($VSREAD) {
        carp
            "Error: must call `get_config(dbsuffix, config_loc)` before this function.\n";
    }

    # my $DBSUFFIX = shift;
    my @stats = @_;
    my $dbh   = get_dbh();
    my $sql_clause;
    ( $ENV{DEBUG} ) && warn Dumper( \@stats ) . "\n";

    $sql_clause = join ", ", @stats;
    ( $ENV{DEBUG} ) && warn "Getting stats: " . $sql_clause . "\n";

    my $sql_res = $dbh->selectrow_hashref("SELECT $sql_clause FROM stats")
        or croak "Couldn't execute statement: " . $dbh->errstr;

    $dbh->disconnect();
    return $sql_res;
}

################################################################

sub stats_get {
    croak 'Cannot use "stats_get" in list context' if wantarray;

    my $argc = @_;
    if ( $argc < 5 ) {
        die "stats_set: expects 5 parameters, passed $argc !\n";
    }

    my $DBSUFFIX = $_[0];
    my $LOGIN    = $_[1];
    my $PASS     = $_[2];
    my $HOST     = $_[3];
    my $NAME     = $_[4];
    my $VALUE    = undef;

    my $dbh = get_dbh();

    # check if database exists first, return undef if not
    unless ($dbh) {
        return;
    }

    # get the namve/value pair
    my $sth = $dbh->prepare("SELECT $NAME FROM stats;")
        or die "Couldn't prepare statement: " . $dbh->errstr;

    $sth->execute()    # Execute the query
        or die "Couldn't execute statement: " . $sth->errstr;

    my @data = $sth->fetchrow_array();

    unless (@data) {
        print STDERR "No field in database  stats.`$NAME'. Aborting!\n\n";
        exit(1);
    }

    $VALUE = $data[0];
    $sth->finish;

    if ( !defined $VALUE ) { $VALUE = ""; }
    ( $ENV{DEBUG} ) && warn "$NAME stat is $VALUE\n";

    $dbh->disconnect();

    return $VALUE;
}

################################################################
sub _init_ref_dbh {
    my ( $refseq, $refleb36, $refindist, $redo ) = @_;
    ( my $refdbfile = $refseq ) =~ s/.seq$/.db/;

    # TODO Maybe first connect to a temp location, backup database
    # to that location then return handle to that location. This
    # is primarily for running on clusters/over NFS.
    my $dbh = DBI->connect(
        "DBI:SQLite:dbname=" . $refdbfile,
        undef, undef,
        {   AutoCommit                 => 1,
            RaiseError                 => 1,
            sqlite_see_if_its_a_number => 1,
        }
        )
        or die "Could not connect to database "
        . $refdbfile
        . ": $DBI::errstr";

    # Set some pragmas to make this part faster
    $dbh->do(q{PRAGMA synchronous = OFF});
    $dbh->do(q{PRAGMA journal_mode = WAL});

    make_refseq_db( $dbh, $refseq, $refleb36, $redo );
    load_refprofiles_db( $dbh, $redo );
    set_indist( $dbh, $refindist, $redo );
    if ( $redo && $VSREAD ) {

# # TODO This was grafted in here. Need to make sure caller tries to catch this to set error.
#         print STDERR "\n\n(updating database with dist/indist info)...";
#         my $installdir = "$FindBin::RealBin";
#         my $exstring   = "$installdir/update_indist.pl";
#         my @exargs     = (
#             qw(-r -k5 -t50 -d),
#             $VSCNF_FILE{DBSUFFIX}, "-u",
#             $ENV{HOME} . "/" . $VSCNF_FILE{DBSUFFIX} . "."
#         );
#         system( $exstring, @exargs );
#         if ( $? == -1 ) {
#             die "Error adding indistinguishable information: $!\n";
#         }
#         else {
#             my $rc = ( $? >> 8 );
#             if ( 0 != $rc ) {
#                 die
#                     "Error: updating database with dist/undist info failed with exit status $rc";
#             }
#         }
# Always set redo flag back to 0 when done redoing
        $VSCNF_FILE{REDO_REFDB} = 0;
        print_config( $ENV{HOME} . "/" . $VSCNF_FILE{DBSUFFIX} . "." );
    }

    # END TODO

    # Return the db handle
    return $dbh;
}
################################################################
sub get_ref_dbh {
    unless ( @_ == 4 ) {
        croak "Error: expecting 4 arguments, got " . @_ * 1 . "\n";
    }
    my ( $refseq, $refleb36, $refindist, $redo ) = @_;
    ( $ENV{DEBUG} ) && warn "get_ref_dbh arguments: " . Dumper( \@_ ) . "\n";
    my $dbh = _init_ref_dbh( $refseq, $refleb36, $refindist, $redo );

    return $dbh;
}

################################################################
sub get_dbh {
    unless ($VSREAD) {
        croak
            "Error: must call `get_config(dbsuffix, config_loc)` before this function.\n";
    }

    # Might use in the future, connection options
    my $opts = shift // {};

    my $dbh;
    if ( $VSCNF_FILE{BACKEND} eq "sqlite" ) {
        my $dbfile
            = "$VSCNF_FILE{OUTPUT_ROOT}/vntr_$VSCNF_FILE{DBSUFFIX}/$VSCNF_FILE{DBSUFFIX}.db";

# TODO First connect to a temp location, backup database to that location
# then return handle to that location. This is primarily for running on clusters.
        $dbh = DBI->connect(
            "DBI:SQLite:dbname=$dbfile",
            undef, undef,
            {   AutoCommit                 => 1,
                RaiseError                 => 1,
                sqlite_see_if_its_a_number => 1,
                ReadOnly => 1 * ( exists $opts->{readonly} ) || 0
            }
        ) or die "Could not connect to database $dbfile: $DBI::errstr";

        # Set default journal to write-ahead log
        $dbh->do(q{PRAGMA journal_mode = WAL})
            or die "Could not do statement on database $dbfile: $DBI::errstr";

        # 800MB cache
        $dbh->do("PRAGMA cache_size = 800000")
            or die "Could not do statement on database $dbfile: $DBI::errstr";

        # Attach reference set database
        if ( exists $opts->{userefdb} && $opts->{userefdb} ) {
            my $refdbfile = $VSCNF_FILE{REFERENCE_DB};
            ( $ENV{DEBUG} )
                && warn "Connecting to refseq db at $refdbfile\n";
            $dbh->do(qq{ATTACH DATABASE "$refdbfile" AS refdb})
                or die
                "Could not attach refseq db '$refdbfile': $DBI::errstr";
        }
    }

    return $dbh;
}

####################################
sub make_refseq_db {

    # TODO DBI error checking
    unless ( @_ == 4 ) {
        croak "Error: expecting 4 arguments, got " . @_ * 1 . "\n";
    }
    my ( $dbh, $refseq, $refdbfile, $redo ) = @_;

    ( $ENV{DEBUG} ) && warn "CWD = " . getcwd() . "\n";

    # Create the table in this new db.
    my $create_fasta_ref_reps_q
        = q{CREATE TABLE IF NOT EXISTS `fasta_ref_reps` (
    `rid` integer NOT NULL,
    `firstindex` integer NOT NULL,
    `lastindex` integer NOT NULL,
    `copynum` float NOT NULL,
    `head` varchar(100) DEFAULT NULL,
    `flankleft` text COLLATE BINARY,
    `pattern` text NOT NULL,
    `sequence` text NOT NULL,
    `flankright` text COLLATE BINARY,
    `conserved` float DEFAULT NULL,
    `comment` varchar(500) DEFAULT NULL,
    `is_singleton` integer NOT NULL DEFAULT '0',
    `is_dist` integer NOT NULL DEFAULT '0',
    `is_indist` integer NOT NULL DEFAULT '0',
    PRIMARY KEY (`rid`),
    UNIQUE (`rid`,`comment`))};
    my $fasta_ref_reps_index_q = q{CREATE INDEX IF NOT EXISTS
        "idx_fasta_ref_reps_head" ON "fasta_ref_reps" (`head`)};
    $dbh->do($create_fasta_ref_reps_q);
    $dbh->do($fasta_ref_reps_index_q);

    # If the table does not exist, or we are forcing a redo, load
    # the profiles into the db.
    my ($num_rows)
        = $dbh->selectrow_array(q{SELECT COUNT(*) FROM fasta_ref_reps});

    ( $ENV{DEBUG} )
        && warn "Rows in fasta_ref_reps == $num_rows, redo set to $redo\n";

    if ( $num_rows == 0 || $redo ) {
        warn "Creating reference sequence database...\n";
        $dbh->do(q{DROP TABLE IF EXISTS fasta_ref_reps});
        $dbh->do($create_fasta_ref_reps_q);
        $dbh->do($fasta_ref_reps_index_q);

        # Read in ref file and create a virtual table out of it
        $dbh->sqlite_create_module(
            perl => "DBD::SQLite::VirtualTable::PerlData" );

        our $seq_rows = [];

        # my $installdir = "$FindBin::RealBin";
        open my $refset, "<", $refseq
            or confess "Error opening reference sequences file "
            . $refseq
            . ": $!.\nStopped at";

        $dbh->begin_work;

        # Skip header
        <$refset>;

        # Insert all ref TRs from refset file
        while ( my $r = <$refset> ) {
            chomp $r;
            my @fields = split /,/, $r;
            push @$seq_rows, \@fields;
        }

        ( $ENV{DEBUG} )
            && warn "Read "
            . scalar(@$seq_rows)
            . " lines from refseq file\n";

        close $refset;

        # Create a virtual tables for the inputs
        $dbh->do(
            q{CREATE VIRTUAL TABLE temp.seqtab
            USING perl(rid integer,
                firstindex integer,
                lastindex integer,
                copynum float,
                head varchar(100),
                flankleft text,
                pattern text,
                sequence text,
                flankright text,
                conserved float,
            arrayrefs="vutil::seq_rows")}
        );

        my $cols = join(
            ",", qw(rid
                firstindex
                lastindex
                copynum
                head
                flankleft
                pattern
                sequence
                flankright
                conserved)
        );

        my $num_rows = $dbh->do(
            qq{INSERT INTO fasta_ref_reps ($cols)
            SELECT * FROM temp.seqtab}
        );

        if ( $num_rows != @$seq_rows ) {
            $dbh->rollback;
            $dbh->disconnect;

            # unlink($refdbfile);
            die "Error inserting reference sequences into "
                . $refdbfile
                . ": inserted $num_rows but read "
                . scalar(@$seq_rows)
                . " lines from file. Aborting...";
        }

        $dbh->commit;
    }
}

####################################
sub get_trunc_query {
    die "Error: need table name for truncate query\n"
        unless @_ == 2;
    my ( $backend, $table ) = @_;
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
# Use this to load reference set db with profiles (leb36 file)
sub load_refprofiles_db {
    my ( $dbh, $refleb36, $redo ) = @_;

    # By default, set redund flag to 1 so we can set
    # non-redundant trs to 0 later;
    my $create_ref_profiles_q = q{CREATE TABLE IF NOT EXISTS `ref_profiles` (
    `rid` integer NOT NULL,
    `proflen` integer NOT NULL,
    `proflenrc` integer NOT NULL,
    `profile` text NOT NULL,
    `profilerc` text NOT NULL,
    `nA` integer NOT NULL,
    `nC` integer NOT NULL,
    `nG` integer NOT NULL,
    `nT` integer NOT NULL,
    `redund` integer NOT NULL default 1,
    `minrepidx` integer default 0,
    PRIMARY KEY (`rid`))};
    $dbh->do($create_ref_profiles_q);

    # If the table does not exist, or we are forcing a redo, load
    # the profiles into the db.
    my ($num_rows)
        = $dbh->selectrow_array(q{SELECT COUNT(*) FROM ref_profiles});
    if ( $num_rows == 0 || $redo ) {
        $dbh->begin_work();
        $dbh->do(q{DROP TABLE IF EXISTS `ref_profiles`});
        $dbh->do($create_ref_profiles_q);

        # Read in ref leb36 file and create a virtual table out of it
        $dbh->sqlite_create_module(
            perl => "DBD::SQLite::VirtualTable::PerlData" );

        our $leb36_rows = [];
        open my $refleb36, "<", "$refleb36"
            or die "Error opening reference profiles file: $!\n";

        while ( my $r = <$refleb36> ) {
            chomp $r;
            my @fields = split /\s+/, $r;
            push @$leb36_rows, [ @fields[ 0, 3 .. 11 ] ];
        }

        ( $ENV{DEBUG} )
            && warn "Read "
            . scalar(@$leb36_rows)
            . " lines from refleb file\n";

        close $refleb36;

        $dbh->do(
            q{CREATE VIRTUAL TABLE temp.leb36tab
            USING perl(rid integer,
                proflen integer,
                proflenrc integer,
                profile text,
                profilerc text,
                nA integer,
                nC integer,
                nG integer,
                nT integer,
            arrayrefs="vutil::leb36_rows")}
        );

        my $cols = join(
            ",", qw(rid
                proflen
                proflenrc
                profile
                profilerc
                nA
                nC
                nG
                nT)
        );

        # Insert all ref TRs from refleb36 file
        my $num_rows = $dbh->do(
            qq{INSERT INTO ref_profiles ($cols)
            SELECT * FROM temp.leb36tab}
        );

        if ( $num_rows != @$leb36_rows ) {
            $dbh->rollback;
            $dbh->disconnect;
            die
                "Error inserting reference profiles into database: inserted $num_rows but read "
                . scalar(@$leb36_rows)
                . " lines from file. Aborting...";
        }

        $dbh->do(q{DROP TABLE temp.leb36tab});

        $dbh->commit;

        # Run redund
        run_redund( $refleb36, "reference.leb36", 0, 0 );
        return 1;
    }

    return 0;
}

################################################################
sub set_indist {
    my ( $dbh, $indist_file, $redo ) = @_;
    return
        unless $redo;

    ( $ENV{DEBUG} )
        && warn "Updating indistinguishable TRs...\n";
    open my $indist_fh, "<", $indist_file
        or die "Error opening indist file $indist_file: $!\n";
    my $count = 0;
    while ( my $line = <$indist_fh> ) {
        $count++;
        $line =~ /-(\d+)/;

        $dbh->do(
            qq{UPDATE refdb.fasta_ref_reps SET is_singleton=0,is_indist=1 WHERE rid=$1}
        ) or die $dbh->errstr;
    }
}

################################################################
# Use database to produce the profiles needed to run redund.exe
# and update database on redundant TRs.
sub run_redund {

    # First, create a temporary directory and copy the filtered set there
    my $tmpdir      = File::Temp->newdir();
    my $tmpdir_name = $tmpdir->dirname;
    my ( $dbh, $input_refset, $output_basename, $keep_files, $redo ) = @_;

    # # Add - sign if negating for ref set and set new input to same
    # # file path as the final output
    # if ($use_as_ref) {
    #     my $cmd = q(awk '{print "-"$0}')
    #         . qq{ $input_refset > $tmpdir_name/$output_basename};

    #     # warn "awk command: $cmd\n";
    #     system($cmd);
    #     $input_refset = "$tmpdir_name/$output_basename";

    #     # warn "New input: $input_refset\n";
    # }

    # Check if we need to run redund, or are forcing a rerun

    my ($num_rows)
        = $dbh->selectrow_array(
        q{SELECT COUNT(*) FROM ref_profiles WHERE redund = 0});

# Return if redund has already been run for this refset and were are not forcing a redo
    return
        unless ( $num_rows == 0 || $redo );

    #=<<Run redund.exe on the input leb36 files.>>
    # First run sorts by minimum representation
    my $tmp_file = File::Temp->new( SUFFIX => ".leb36", DIR => $tmpdir_name );
    my $installdir = "$FindBin::RealBin";
    system("$installdir/redund.exe $input_refset $tmp_file -s -i");
    if ( $? == -1 ) { die "command failed: $!\n"; }
    else {
        my $rc = ( $? >> 8 );
        if ( 0 != $rc ) { die "command exited with value $rc"; }
    }

    # Second run eliminates redundancies
    system(
        "$installdir/redund.exe $tmp_file $tmpdir_name/$output_basename -n -i"
    );
    if ( $? == -1 ) { die "command failed: $!\n"; }
    else {
        my $rc = ( $? >> 8 );
        if ( 0 != $rc ) { die "command exited with value $rc"; }
    }

    #=<<Mark all TRs in file as non-redundant>>
    my $rotindex_fn = "$tmpdir_name/$output_basename.rotindex";
    open my $set_fh, "<", $rotindex_fn;

    my @unique_trs;
    while ( my $entry = <$set_fh> ) {
        $entry =~ /(\d+)['"]/;
        my $trid = $1;
        push @unique_trs, $trid;
    }
    close $set_fh;

    # warn Dumper($unique_trs) ."\n";

    warn "Marking non-redundant TRs in database.\n";
    $dbh->begin_work();
    $dbh->do(
        q{CREATE TEMPORARY TABLE temp.unique_trs
        (`rid` integer PRIMARY KEY)}
    ) or die "Couldn't do statement: " . $dbh->errstr;

    my $insert_sth
        = $dbh->prepare(q{INSERT INTO temp.unique_trs (rid) VALUES (?)});
    for my $rid (@unique_trs) {
        $insert_sth->execute($rid);
    }

    # Update the ref_profiles table
    my $update_redund = q{UPDATE ref_profiles SET redund=0
        WHERE EXISTS (
        SELECT rid FROM temp.unique_trs t2
        WHERE ref_profiles.rid = t2.rid
    )};
    my $unique_tr_count = $dbh->do($update_redund)
        or die "Couldn't do statement: " . $dbh->errstr;
    $dbh->commit;

    $dbh->do(
        q{CREATE TABLE IF NOT EXISTS `files` (
        `rotindex` BLOB
        )}
    );

    #=<<Set minimum representation>>
    # Attach database created by redund.exe
    warn "Getting sort order of TRs by minimum representation.\n";
    my $minreporder_db = "$tmp_file.db";
    $dbh->do(qq{ATTACH DATABASE "$minreporder_db" AS minrep})
        or carp "Couldn't do statement: $DBI::errstr\n";
    ( $ENV{DEBUG} ) && warn "Connecting to SQLite db at $minreporder_db\n";

    # print "Press enter to continue...";
    # my $dummy = <STDIN>;

    $dbh->begin_work();

    # Copy the CREATE TABLE query for the minreporder table
    # and create the same table in the ref set db, copying
    # over the data
    my ($minrep_sql) = $dbh->selectrow_array(
        q{SELECT sql FROM minrep.sqlite_master
        WHERE name = 'minreporder'}
    ) or carp "Couldn't do statement: $DBI::errstr\n";
    $minrep_sql =~ s/minreporder/main.minreporder/;
    $dbh->do(q{DROP TABLE IF EXISTS main.minreporder})
        or carp "Couldn't do statement: $DBI::errstr\n";
    $dbh->do($minrep_sql)
        or carp "Couldn't do statement: $DBI::errstr\n";
    $dbh->do(
        q{INSERT INTO main.minreporder
            SELECT * FROM minrep.minreporder}
    ) or carp "Couldn't do statement: $DBI::errstr\n";
    $dbh->commit;
    $dbh->do(qq{DETACH minrep})
        or carp "Couldn't do statement: $DBI::errstr\n";

    # Save the rotindex file
    warn "Saving rotindex (redundancy index) into database\n";
    my $rotindex;
    {
        local $/;
        open my $fh, '<', $rotindex_fn or die "can't open $rotindex_fn: $!";
        $rotindex = <$fh>;
        close $fh;
    }
    my $load_rotindex_sth
        = $dbh->prepare(qq{INSERT INTO files (rotindex) VALUES (?)})
        or die "Couldn't do statement: $DBI::errstr";
    $load_rotindex_sth->bind_param( 1, $rotindex, DBI::SQL_BLOB );
    $load_rotindex_sth->execute;

    warn "$unique_tr_count unique TRs in filtered set\n";

    unlink($minreporder_db);

    if ($keep_files) {
        return $tmpdir;
    }
}

################################################################
# TODO Function to generate leb36 file from db after redund was
# run. Will execute run_redund unless 'minreporder' table exists
# SQL query: SELECT COUNT(*) FROM sqlite_master WHERE type = 'table' AND name='minreporder';
# Code to import from produce_indist

################################################################
sub write_sqlite {
    unless ($VSREAD) {
        carp
            "Error: must call `get_config(dbsuffix, config_loc)` before this function.\n";
    }

    my $installdir = "$FindBin::RealBin";
    my $exestring
        = "sqlite3 $VSCNF_FILE{OUTPUT_ROOT}/vntr_$VSCNF_FILE{DBSUFFIX}/$VSCNF_FILE{DBSUFFIX}.db < $installdir/sqlite_schema.sql";
    warn "Executing: $exestring\n";
    system($exestring);
}

################################################################

sub write_mysql {
    unless ($VSREAD) {
        carp
            "Error: must call `get_config(dbsuffix, config_loc)` before this function.\n";
    }

    open my $sql_file, ">",
        "$VSCNF_FILE{TMPDIR}/VNTRPIPE_$VSCNF_FILE{DBNAME}.sql"
        or die $!;

    print $sql_file
        "CREATE database IF NOT EXISTS VNTRPIPE_$VSCNF_FILE{DBNAME};\n\n";
    print $sql_file "USE VNTRPIPE_$VSCNF_FILE{DBNAME};\n\n";

    print $sql_file <<TEST;

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

    close($sql_file);

    return 0;
}

################################################################

sub print_config {

    my $argc = @_;
    if ( $argc < 1 ) {
        croak "print_config: expects 1 parameters, passed $argc!\n";
    }

    my $startdir = $_[0];

    if ( !defined $VSCNF_FILE{"NPROCESSES"} ) {
        $VSCNF_FILE{"NPROCESSES"} = -1;
    }
    if ( !defined $VSCNF_FILE{"MIN_FLANK_REQUIRED"} ) {
        $VSCNF_FILE{"MIN_FLANK_REQUIRED"} = -1;
    }
    if ( !defined $VSCNF_FILE{"MAX_FLANK_CONSIDERED"} ) {
        $VSCNF_FILE{"MAX_FLANK_CONSIDERED"} = -1;
    }
    if ( !defined $VSCNF_FILE{"MIN_SUPPORT_REQUIRED"} ) {
        $VSCNF_FILE{"MIN_SUPPORT_REQUIRED"} = -1;
    }

    if ( !defined $VSCNF_FILE{"SERVER"} ) { $VSCNF_FILE{"SERVER"} = ""; }
    if ( !defined $VSCNF_FILE{"STRIP_454_KEYTAGS"} ) {
        $VSCNF_FILE{"STRIP_454_KEYTAGS"} = -1;
    }
    if ( !defined $VSCNF_FILE{"IS_PAIRED_READS"} ) {
        $VSCNF_FILE{"IS_PAIRED_READS"} = -1;
    }
    if ( !defined $VSCNF_FILE{"HTML_DIR"} ) { $VSCNF_FILE{"HTML_DIR"} = ""; }
    if ( !defined $VSCNF_FILE{"FASTA_DIR"} ) {
        $VSCNF_FILE{"FASTA_DIR"} = "";
    }
    if ( !defined $VSCNF_FILE{"OUTPUT_ROOT"} ) {
        $VSCNF_FILE{"OUTPUT_ROOT"} = "";
    }
    if ( !defined $VSCNF_FILE{"TMPDIR"} ) { $VSCNF_FILE{"TMPDIR"} = ""; }
    if ( !defined $VSCNF_FILE{"REFERENCE_FILE"} ) {
        $VSCNF_FILE{"REFERENCE_FILE"} = "";
    }
    if ( !defined $VSCNF_FILE{"REFERENCE_SEQ"} ) {
        $VSCNF_FILE{"REFERENCE_SEQ"} = "";
    }
    if ( !defined $VSCNF_FILE{"REFERENCE_INDIST"} ) {
        $VSCNF_FILE{"REFERENCE_INDIST"} = "";
    }
    if ( !defined $VSCNF_FILE{"REFERENCE_INDIST_PRODUCE"} ) {
        $VSCNF_FILE{"REFERENCE_INDIST_PRODUCE"} = -1;
    }

    # look in the directory the script was started in
    if ( open( MFREAD, ">${startdir}vs.cnf" ) ) {

        print MFREAD <<CNF;
# COPY THIS FILE TO A DIFFERENT LOCATION AND SET ALL VARIABLES. 
# DO NOT FORGET TO CHMOD THIS FILE TO PREVENT OTHER PEOPLE ON 
# THE SYSTEM FROM LEARNING YOUR MYSQL CREDENTIALS.

# Database backend
BACKEND=$VSCNF_FILE{"BACKEND"}

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

# Rebuild reference database
# eg, 0 = no 
# eg, 1 - yes
REDO_REFDB=$VSCNF_FILE{"REDO_REFDB"}

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

CNF

        close(MFREAD);

        chmod 0600, "${startdir}vs.cnf";

    }
    else {

        die "print_config: can't open '${startdir}vs.cnf' for writing!\n";
    }

}

####################################################################################

1;

