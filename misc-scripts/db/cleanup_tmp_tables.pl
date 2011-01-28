#!/usr/local/ensembl/bin/perl

=head1 NAME

cleanup_tmp_tables.pl - delete temporary and backup tables from a database

=head1 SYNOPSIS

cleanup_tmp_tables.pl [arguments]

Required arguments:

  --dbname, db_name=NAME              database name NAME
  --host, --dbhost, --db_host=HOST    database host HOST
  --port, --dbport, --db_port=PORT    database port PORT
  --user, --dbuser, --db_user=USER    database username USER
  --pass, --dbpass, --db_pass=PASS    database passwort PASS

Optional arguments:

  --conffile, --conf=FILE             read parameters from FILE
                                      (default: conf/Conversion.ini)

  --logfile, --log=FILE               log to FILE (default: *STDOUT)
  --logpath=PATH                      write logfile to PATH (default: .)
  --logappend, --log_append           append to logfile (default: truncate)

  -v, --verbose=0|1                   verbose logging (default: false)
  -i, --interactive=0|1               run script interactively (default: true)
  -n, --dry_run, --dry=0|1            don't write results to database
  -h, --help, -?                      print help (this message)

=head1 DESCRIPTION


=head1 LICENCE

This code is distributed under an Apache style licence:
Please see http://www.ensembl.org/code_licence.html for details

=head1 AUTHOR

Patrick Meidl <meidl@ebi.ac.uk>, Ensembl core API team

=head1 CONTACT

Please post comments/questions to the Ensembl development list
<dev@ensembl.org>

=cut

use strict;
use warnings;
no warnings 'uninitialized';

use FindBin qw($Bin);
use vars qw($SERVERROOT);

BEGIN {
    $SERVERROOT = "$Bin/../..";
    unshift(@INC, "$SERVERROOT/ensembl/modules");
    unshift(@INC, "$SERVERROOT/bioperl-live");
}

use Getopt::Long;
use Pod::Usage;
use Bio::EnsEMBL::Utils::ConversionSupport;

$| = 1;

my $support = new Bio::EnsEMBL::Utils::ConversionSupport($SERVERROOT);

# parse options
$support->parse_common_options(@_);
$support->allowed_params(
  $support->get_common_params,
);

if ($support->param('help') or $support->error) {
    warn $support->error if $support->error;
    pod2usage(1);
}

# ask user to confirm parameters to proceed
$support->confirm_params;

# get log filehandle and print heading and parameters to logfile
$support->init_log;

$support->check_required_params;

# connect to database and get adaptors
my $dba = $support->get_database('ensembl');
my $dbh = $dba->dbc->db_handle;

# find all temporary and backup tables
my @tables;

foreach my $pattern (qw(tmp temp bak backup)) {
  my $sql = qq(SHOW TABLES LIKE '\%$pattern\%');
  my $sth = $dbh->prepare($sql);
  $sth->execute;
  while (my ($table) = $sth->fetchrow_array) {
    push @tables, $table;
  }
}

if ($support->param('dry_run')) {
  # for a dry run, only show which databases would be deleted
  $support->log("Temporary and backup tables found:\n");
  foreach my $table (sort @tables) {
    $support->log("$table\n", 1);
  }

} else {
  # delete tables
  foreach my $table (sort @tables) {
    if ($support->user_proceed("Drop table $table?")) {
      $support->log("Dropping table $table...\n");
      $dbh->do("DROP TABLE $table");
      $support->log("Done.\n");
    }
  }
}

# finish logfile
$support->finish_log;

