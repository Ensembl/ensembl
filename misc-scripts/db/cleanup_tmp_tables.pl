#!/usr/bin/env perl
# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2018] EMBL-European Bioinformatics Institute
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#      http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


=head1 NAME

cleanup_tmp_tables.pl - delete temporary and backup tables from a database

=head1 SYNOPSIS

  ./cleanup_tmp_tables.pl [arguments]
  
  ./cleanup_tmp_tables.pl --nolog --host localhost --port 3306 --user user --dbname DB --dry_run
  
  ./cleanup_tmp_tables.pl --nolog --host localhost --port 3306 --user user --dbname '%mydbs%' --dry_run
  
  ./cleanup_tmp_tables.pl --nolog --host localhost --port 3306 --user user --dbname DB --interactive 0

Required arguments:

  --dbname, db_name=NAME              database name NAME (can be a pattern)
  --host, --dbhost, --db_host=HOST    database host HOST
  --port, --dbport, --db_port=PORT    database port PORT
  --user, --dbuser, --db_user=USER    database username USER
  --pass, --dbpass, --db_pass=PASS    database passwort PASS

Optional arguments:

  --mart                              Indicates we wish to search for mart
                                      temporary tables which are normally
                                      prefixed with MTMP_. 
                                      ONLY RUN IF YOU ARE A MEMBER OF 
                                      THE PRODUCTION TEAM

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

A script which looks for any table which we believe could be a temporary
table. This means any table which contains

=over 8

=item bak

=item backup

=item MTMP_ (only used when --mart is specified)

=back

You can run this over multiple DBs but caution is advised


=head1 AUTHOR

Ensembl core API team

=head1 CONTACT

Please post comments/questions to the Ensembl development list
<http://lists.ensembl.org/mailman/listinfo/dev>

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
$support->parse_extra_options(qw/mart!/);
$support->allowed_params(
  $support->get_common_params, 'mart'
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

my @databases;

# connect to database
my $dbh;
my $original_dbname = $support->param('dbname'); 
if($original_dbname =~ /%/) {
  $support->param('dbname', q{});
  $dbh = $support->get_dbconnection('');
  my $ref = $dbh->selectall_arrayref('show databases like ?', {}, $original_dbname);
  push(@databases, map {$_->[0]} @{$ref})
}
else {
  $dbh = $support->get_dbconnection('');
  push(@databases, $original_dbname);
}

my @patterns = map { '%\\_'.$_.'%' } qw/bak backup/;
if($support->param('mart')) {
  if($support->user_proceed('--mart was specified on the command line. Do not run this during a mart build. Do you wish to continue?')) {
    push(@patterns, 'MTMP\\_%');
  }
}

foreach my $db (@databases) {
  my %tables;
  $support->log('Switching to '.$db."\n");
  $dbh->do('use '.$db);
  foreach my $pattern (@patterns) {
    my $ref = $dbh->selectall_arrayref('show tables like ?', {}, $pattern);
    $tables{$_->[0]} = 1 for @{$ref};
  }
  
  my @tables = sort keys %tables;
  
  if ($support->param('dry_run')) {
    # for a dry run, only show which databases would be deleted
    if(scalar(@tables) > 0) {
      $support->log("Temporary and backup tables found:\n");
      foreach my $table (@tables) {
        $support->log("$table\n", 1);
      }
    }
  
  } else {
    # delete tables
    foreach my $table (@tables) {
      if ($support->user_proceed("Drop table $table?")) {
        $support->log("Dropping table $table...\n");
        $dbh->do("DROP TABLE $table");
        $support->log("Done.\n");
      }
    }
  }
}
# finish logfile
$support->finish_log;

