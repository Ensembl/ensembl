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


=head1 SYNOPSIS

.pl [arguments]

Required arguments:

  --basedir=PATH                     base directory of ID mapping results

Optional arguments:

  --conffile, --conf=FILE             read parameters from FILE
                                      (default: conf/Conversion.ini)

  --logfile, --log=FILE               log to FILE (default: *STDOUT)
  --logpath=PATH                      write logfile to PATH (default: .)
  --logappend, --log_append           append to logfile (default: truncate)
  --loglevel=LEVEL                    define log level (default: INFO)

  -i, --interactive=0|1               run script interactively (default: true)
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
use Bio::EnsEMBL::Utils::ScriptUtils qw(path_append);
use Bio::EnsEMBL::Utils::Exception qw(throw warning);

# parse configuration and commandline arguments
my $conf = new Bio::EnsEMBL::Utils::ConfParser(
  -SERVERROOT => "$Bin/../../../..",
  -DEFAULT_CONF => "$Bin/default.conf"
);

$conf->parse_options(
  'path1|p1=s' => 1,
  'path2|p2=s' => 1,
  'type|t=s@' => 0,
);

# set default logpath
unless ($conf->param('logpath')) {
  $conf->param('logpath', path_append($conf->param('basedir'), 'log'));
}

# get log filehandle and print heading and parameters to logfile
my $logger = new Bio::EnsEMBL::Utils::Logger(
  -LOGFILE      => $conf->param('logfile'),
  -LOGAUTO      => $conf->param('logauto'),
  -LOGAUTOBASE  => 'compare_scores',
  -LOGAUTOID    => $conf->param('logautoid'),
  -LOGPATH      => $conf->param('logpath'),
  -LOGAPPEND    => $conf->param('logappend'),
  -LOGLEVEL     => $conf->param('loglevel'),
);

# initialise log
$logger->init_log($conf->list_param_values);

my @types = $conf->param('type') || qw(exon transcript gene);

foreach my $type (@types) {
  &compare_scores($type);
}


# finish logfile
$logger->finish_log;


### END main ###


sub compare_scores {
  my $type = shift;

  # read scores from file
  $logger->info("Reading $type scores...\n", 0, 'stamped');
  
  my $scores1 = &parse_file($conf->param('path1')."/${type}_scores.txt");
  my $scores2 = &parse_file($conf->param('path2')."/${type}_scores.txt");

  $logger->info("Done.\n\n", 0, 'stamped');

  # look for pairs with scores in both result sets
  my %stats;
  my @both = ();
  my @only1 = ();
  my @only2 = ();
  
  foreach my $key (keys %$scores1) {
    if ($scores2->{$key}) {
      push @both, $key;
    } else {
      push @only1, $key;
    }
  }
  
  foreach my $key (keys %$scores2) {
    unless ($scores1->{$key}) {
      push @only2, $key;
    }
  }

  $stats{'TOT1'} = keys %$scores1;
  $stats{'TOT2'} = keys %$scores2;
  $stats{'BOTH'} = @both;
  $stats{'ONLY1'} = @only1;
  $stats{'ONLY2'} = @only2;

  $logger->info("Only in set 1 (first 10 shown):\n");
  my $i;
  foreach my $key (sort @only1) {
    $logger->info(sprintf("%-20s%-10s\n", $key, $scores1->{$key}), 1);
    last if ($i++ == 10);
  }

  $logger->info("\nOnly in set 2 (first 10 shown):\n");
  my $j;
  foreach my $key (sort @only2) {
    $logger->info(sprintf("%-20s%-10s\n", $key, $scores2->{$key}), 1);
    last if ($j++ == 10);
  }
  
  # compare scores which are present in both result sets
  $logger->info("\nScores different (first 10 shown):\n");
  
  foreach my $key (@both) {

    my $s1 = $scores1->{$key};
    my $s2 = $scores2->{$key};
    my $diff = $s1 - $s2;
    $diff = -$diff if ($diff < 0);

    unless ($diff < 0.000002) {
      $stats{'BOTH_DIFF'}++;
      if ($stats{'BOTH_DIFF'} <= 10) {
        $logger->info(sprintf("%-20s%-10s%-10s\n", $key, $s1, $s2), 1);
      }
    }
  }

  $logger->info("\nStats:\n");
  foreach my $t (qw(TOT1 TOT2 BOTH ONLY1 ONLY2 BOTH_DIFF)) {
    $logger->info(sprintf("%-10s%8d\n", $t, $stats{$t}), 1);
  }
}


sub parse_file {
  my $file = shift;

  open(my $fh, '<', $file) or
    throw("Unable to open $file for reading: $!");

  my %scores = ();

  while (my $line = <$fh>) {
    chomp $line;
    my ($old_id, $new_id, $score) = split(/\s+/, $line);
    $scores{"$old_id:$new_id"} = sprintf("%.6f", $score);
  }

  close($fh);

  return \%scores;
}


