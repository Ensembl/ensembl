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

# Don't change the above line.
# Change the PATH in the myRun.ksh script if you want to use another perl.

=head1 NAME


=head1 SYNOPSIS

.pl [arguments]

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
  --loglevel=LEVEL                    define log level (default: INFO)

  -i, --interactive=0|1               run script interactively (default: true)
  -n, --dry_run, --dry=0|1            don't write results to database
  -h, --help, -?                      print help (this message)

=head1 DESCRIPTION



=head1 AUTHOR

Patrick Meidl <meidl@ebi.ac.uk>, Ensembl core API team

=head1 CONTACT

Please post comments/questions to the Ensembl development list
<http://lists.ensembl.org/mailman/listinfo/dev>

=cut

use strict;
use warnings;
no warnings 'uninitialized';

use FindBin qw($Bin);
use Bio::EnsEMBL::Utils::ConfParser;
use Bio::EnsEMBL::Utils::Logger;
use Bio::EnsEMBL::Utils::ScriptUtils qw(inject path_append);
use Bio::EnsEMBL::Utils::Exception qw(throw warning);

# parse configuration and commandline arguments
my $conf = new Bio::EnsEMBL::Utils::ConfParser(
  -SERVERROOT => "$Bin/../../..",
  -DEFAULT_CONF => "$Bin/default.conf"
);

$conf->parse_options(
  'sourcehost|source_host=s' => 1,
  'sourceport|source_port=n' => 1,
  'sourceuser|source_user=s' => 1,
  'sourcepass|source_pass=s' => 0,
  'sourcedbname|source_dbname=s' => 1,
  'targethost|target_host=s' => 1,
  'targetport|target_port=n' => 1,
  'targetuser|target_user=s' => 1,
  'targetpass|target_pass=s' => 0,
  'targetdbname|target_dbname=s' => 1,
  'basedir|basedir=s' => 1,
  'biotypes=s@' => 0,
  'biotypes_include=s@' => 0,
  'biotypes_exclude=s@' => 0,
  'dbtype=s' => 1,
  'cache_impl=s' => 1,
  'index|i=n' => 1,
);

# set default logpath
unless ($conf->param('logpath')) {
  $conf->param('logpath', path_append($conf->param('basedir'), 'log'));
}
# log to a subdirectory to prevent clutter
$conf->param('logpath', path_append($conf->param('logpath'),
  'dump_by_seq_region'));  

# append job index to logfile name
my $dbtype = $conf->param('dbtype');
my $index = $conf->param('index');
my $logautobase = ($conf->param('logautobase') || 'dump_by_seq_region').
  ".$dbtype";

# get log filehandle and print heading and parameters to logfile
my $logger = new Bio::EnsEMBL::Utils::Logger(
  -LOGFILE      => $conf->param('logfile'),
  -LOGAUTO      => $conf->param('logauto'),
  -LOGAUTOBASE  => $logautobase,
  -LOGAUTOID    => $conf->param('logautoid'),
  -LOGPATH      => $conf->param('logpath'),
  -LOGAPPEND    => 1,
  -LOGLEVEL     => $conf->param('loglevel'),
  -IS_COMPONENT => 1,
);

# build cache
my $cache_impl = $conf->param('cache_impl');
inject($cache_impl);

my $cache = $cache_impl->new(
  -LOGGER       => $logger,
  -CONF         => $conf,
);

# determine which slice to process. to do so, read the file containing the
# slices to be processed, and take the one at position $index
my $logpath = $conf->param('logpath');
my $filename = "$dbtype.dump_cache.slices.txt";
open(my $fh, '<', "$logpath/$filename") or
  throw("Unable to open $logpath/$filename for reading: $!");
my @slice_names = <$fh>;
my $slice_name = $slice_names[$index-1];
chomp($slice_name);
close($fh);

# no build the cache for this slice
$cache->build_cache_by_slice($dbtype, $slice_name);

# set flag to indicate everything went fine
my $success_file = $conf->param('logpath')."/$logautobase.$index.success";
open(TMPFILE, '>', $success_file) and close TMPFILE
  or throw "Can't open $success_file for writing: $!";

