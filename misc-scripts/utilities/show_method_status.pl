#!/usr/bin/env perl
# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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

show_method_status.pl - script to print "Status" of API methods

=head1 SYNOPSIS

show_method_status.pl [arguments]

Optional arguments:

  --path, --root=PATH                 directory root path to check (use absolute
                                      path or relative to cwd)
  --exclude=LIST                      don't show LISTed statuses
  --include=LIST                      only show LISTed statuses
  --show_empty, --show-empty, --empty print method name even if there's no
                                      status to report for it (only required 
                                      when using --exclude or --include)

  --conffile, --conf=FILE             read parameters from FILE

  --logfile, --log=FILE               log to FILE (default: *STDOUT)
  --logpath=PATH                      write logfile to PATH (default: .)
  --logappend, --log_append           append to logfile (default: truncate)

  -i, --interactive                   run script interactively (default: false)
  -h, --help, -?                      print help (this message)

=head1 DESCRIPTION

This script will print the "Status" documentation for each method in each perl
module found in the directory specified by --path (recursively). Output can be
limited to certain Statuses with the --exclude or --include options.

=head1 EXAMPLES

Show all methods which are 'At Risk':

  $ ./show_method_status.pl --path ../../modules/Bio/EnsEMBL --include 'at risk'

Show all methods except those that are 'Stable':

  $ ./show_method_status.pl --path ../../modules/Bio/EnsEMBL --exclude 'stable'


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
use File::Find qw(find);
use Cwd qw(getcwd abs_path);

# parse configuration and commandline arguments
my $conf = new Bio::EnsEMBL::Utils::ConfParser(
  -SERVERROOT => "$Bin/../../..",
  -DEFAULT_CONF => ""
);
warn $Bin;

$conf->parse_options(
  'path|root=s' => 0,
  'exclude=s@' => 0,
  'include=s@' => 0,
  'show_empty|show-empty|empty' => 0,
);

$conf->param('path', '.') unless $conf->param('path');

# get log filehandle and print heading and parameters to logfile
my $logger = new Bio::EnsEMBL::Utils::Logger(
  -LOGFILE    => $conf->param('logfile'),
  -LOGPATH    => $conf->param('logpath'),
  -LOGAPPEND  => $conf->param('logappend'),
  -LOGLEVEL   => $conf->param('loglevel'),
);

# initialise log
$logger->init_log($conf->list_param_values);

# recursively process all files
my $path = $conf->param('path');
$path = abs_path(getcwd."/$path") if ($path =~ /^\./);
find(\&parse_files, $path);

# finish logfile
$logger->finish_log;


### END main ###

sub parse_files {
  my $file = $_;

  # only read perl modules
  return unless ($file =~ /\.pm$/);

  # read file
  open(IN, $file) or die "Unable to open $file: $!\n";;
  
  my $pod_flag;
  my $method;
  my $status;
  my $result;
  
  LINE:
  while (my $line = <IN>) {
    chomp $line;
    
    # start of method pod
    if ($line =~ /=head2 (.*)$/) {
      
      $method = sprintf("%-40s", $1);
      $pod_flag = 1;
    
    # status
    } elsif ($line =~ /Status\s*:\s*(.+)/) {
    
      $status = $1;
    
    # end of method pod
    } elsif ($line =~ /=cut/ and $pod_flag) {
    
      # set status to unknown if not found
      $status ||= 'unknown';
      
      # exclude specified statuses from output
      foreach my $pattern ($conf->param('exclude')) {
        next LINE if ($status =~ /$pattern/i);
      }

      # only include specified statuses in output
      foreach my $pattern ($conf->param('include')) {
        next LINE if ($pattern and !($status =~ /$pattern/i));
      }

      $result .= "  $method $status\n";
      
      $status = undef;
      $pod_flag = 0;
    
    }
  }

  # log result for this module
  if ($result or $conf->param('show_empty')) {
    my $filepath = $File::Find::name;
    $filepath =~ s/$path\///;
    $logger->info("\n$filepath\n$result");
  }

}

