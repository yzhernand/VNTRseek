#!/usr/bin/perl

# command line usage example:
#  ./checkl3b36.pl inputfolder referencefolder 

use strict;
use warnings;
use Cwd;
use POSIX qw(strftime);

if (@ARGV<2) {
 print STDERR "\n\ncheckl3b36.pl: NOT ENOUGH INPUT PARAMS!\n";
 exit 1;
}

say STDERR strftime( "\n\nstart: %F %T\n\n", localtime );

my $curdir =  getcwd;
my $tgz_dir = $ARGV[0];
my $reffolder = $ARGV[1];

# get a list of input files
opendir(DIR, $tgz_dir);
# the only extensions are .leb36
my @tarballs = grep(/\.(?:leb36)$/, readdir(DIR));
closedir(DIR);
my $tarball_count = @tarballs;
print STDERR "$tarball_count  supported files found in $tgz_dir\n";
die "Exiting\n" if $tarball_count == 0;

# enter dir
chdir($tgz_dir);


# make sure the ids are unique (renumbering was not ran correctly)
print STDERR "Checking reference leb36 file...\n";
my %uhash = ();
open FILE, "<$reffolder/reference.leb36" or die $!;
while (<FILE>) {
  if (/^(\d+)/) {
      if (exists $uhash{$1}) { die "Non-unique id ($1) detected in reads. Were steps 2 and 3 executed?\n"; }
      $uhash{$1} = 1;
  }
}
close FILE;

print STDERR "Checking read leb36 files...\n";
%uhash = ();
foreach my $ifile (@tarballs) {
  open FILE, "<$ifile" or die $!;
  while (<FILE>) {
    if (/^(\d+)/) {
        if (exists $uhash{$1}) { die "Non-unique id ($1) detected in reads. Were steps 2 and 3 executed?\n"; }
        $uhash{$1} = 1;
    }
  }
  close FILE;
}

%uhash = ();

say STDERR strftime( "\n\nend: %F %T\n\n", localtime );
