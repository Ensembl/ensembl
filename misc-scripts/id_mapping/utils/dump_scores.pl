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

dump_scores.pl - dump scores from serialised ScoredMappingMatrix'es for debugging

=head1 SYNOPSIS

dump_scores.pl [arguments]

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

This script reads gene, transcript and exon scores from serialised
ScoredMappingMatrix files and dumps the data (old_internal_id, new_internal_id,
score) to a text file for debugging.

Note that if you ran the ID mapping with loglevel=DEBUG, these dummps are
generated automatically so you won't need to run this script.


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
use Bio::EnsEMBL::IdMapping::ScoredMappingMatrix;

# parse configuration and commandline arguments
my $conf = new Bio::EnsEMBL::Utils::ConfParser(
  -SERVERROOT => "$Bin/../../../..",
  -DEFAULT_CONF => "$Bin/default.conf"
);

$conf->parse_options(
  'basedir|basedir=s' => 1,
);

# set default logpath
unless ($conf->param('logpath')) {
  $conf->param('logpath', path_append($conf->param('basedir'), 'log'));
}

# get log filehandle and print heading and parameters to logfile
my $logger = new Bio::EnsEMBL::Utils::Logger(
  -LOGFILE      => $conf->param('logfile'),
  -LOGAUTO      => $conf->param('logauto'),
  -LOGAUTOBASE  => 'dump_scores',
  -LOGAUTOID    => $conf->param('logautoid'),
  -LOGPATH      => $conf->param('logpath'),
  -LOGAPPEND    => $conf->param('logappend'),
  -LOGLEVEL     => $conf->param('loglevel'),
);

# initialise log
$logger->init_log($conf->list_param_values);

my $dump_path = path_append($conf->param('basedir'), 'matrix');

# genes
my $gene_matrix = &read_matrix('gene');
&dump_scores('gene', $gene_matrix);

# transcripts
my $transcript_matrix = &read_matrix('transcript');
&dump_scores('transcript', $transcript_matrix);

# exons
my $exon_matrix = &read_matrix('exon_overlap');
my $exonerate_matrix = &read_matrix('exon_exonerate');
$exon_matrix->merge($exonerate_matrix);
&dump_scores('exon', $exon_matrix);


# finish logfile
$logger->finish_log;


### END main ###


sub read_matrix {
  my $type = shift;

  my $matrix = Bio::EnsEMBL::IdMapping::ScoredMappingMatrix->new(
    -DUMP_PATH   => $dump_path,
    -CACHE_FILE  => "${type}_matrix.ser",
  );

  my $cache = $matrix->cache_file;

  if (-s $cache) {
    
    # read from file
    $logger->info("Reading $type scoring matrix from file...\n", 0, 'stamped');
    $logger->debug("Cache file $cache.\n", 1);
    $matrix->read_from_file;
    $logger->info("Done.\n\n", 0, 'stamped');
    
  } else {
    $logger->warning("No cache file found at $cache.\n");
  }

  return $matrix;
}


sub dump_scores {
  my $type = shift;
  my $matrix = shift;

  $logger->info("Dumping  $type scores to file...\n", 0, 'stamped');

  my $debug_path = path_append($conf->param('basedir'), 'debug');
  my $logfile = "$debug_path/${type}_scores.txt";
  
  open(my $fh, '>', $logfile) or
    throw("Unable to open $logfile for writing: $!");

  #my $i = 0;
  foreach my $entry (@{ $matrix->get_all_Entries }) {
    #$logger->info($entry->to_string."\n");
    #last if (++$i == 10);
    print $fh ($entry->to_string."\n");
  }

  close($fh);
  
  $logger->info("Done.\n\n", 0, 'stamped');
}


