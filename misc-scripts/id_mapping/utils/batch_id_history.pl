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

batch_id_history.pl - find stable IDs in the archive

=head1 SYNOPSIS

find_stable_ids_in_archive.pl [arguments]

Required arguments:

  --host=hOST                 database host HOST
  --port=PORT                 database port PORT
  --user=USER                 database username USER
  --dbname=NAME               database name NAME
  --stable_id_file=FILE       read stable ID list from FILE

Optional arguments:

  --pass=PASS                 database passwort PASS
  --outfile=FILE              write output to FILE
  --pep_seq                   print peptide sequence


=head1 DESCRIPTION

This script reads a list of stable IDs from a file and sees if it can find them
in the stable ID archive. It will print the ID history for each of them and
optinally the peptide sequence found there as well. Note that this will not 
print the full history network, but rather branch out from your focus stable ID
only. If you are interested in the full network, have a look at
Bio::EnsEMBL::StableIdHistoryTree and related modules.


=head1 AUTHOR

Patrick Meidl <meidl@ebi.ac.uk>, Ensembl core API team

=head1 CONTACT

Please post comments/questions to the Ensembl development list
<http://lists.ensembl.org/mailman/listinfo/dev>

=cut

use strict;
use warnings;
no warnings 'uninitialized';

use Getopt::Long;
use Bio::EnsEMBL::SimpleFeature;
use Bio::EnsEMBL::Analysis;
use Bio::EnsEMBL::DBSQL::DBAdaptor;

$| = 1;

my ($host, $port, $user, $pass, $dbname, $stable_id_file, $outfile, $pep_seq);

GetOptions(
    "host=s",             \$host,
    "port=i",             \$port,
    "user=s",             \$user,
    "pass=s",             \$pass,
    "dbname=s",           \$dbname,
    "stable_id_file=s",   \$stable_id_file,
    "outfile=s",          \$outfile,
    "pep_seq",            \$pep_seq,
);

# check required params
unless ($host && $port && $user && $dbname && $stable_id_file) {
  die "ERROR: Unable to run script.\nNeed host, port, user, dbname and stable_id_file parameters.\n";
}

# connect to database and get adaptors
my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
    -HOST     => $host,
    -PORT     => $port,
    -USER     => $user,
    -PASS     => $pass,
    -DBNAME   => $dbname,
);

my $aa = $db->get_ArchiveStableIdAdaptor;

# read list of stable IDs from file
my $infh;
open($infh, "<", $stable_id_file) or
  die("Can't open $stable_id_file for reading: $!");

# get output filehandle
my $outfh;
if ($outfile) {
  open($outfh, ">", $outfile) or die("Can't open $outfile for writing: $!");
} else {
  $outfh = \*STDOUT;
}

while (my $sid = <$infh>) {
  
  # skip comments and empty lines
  next if (/^#/ or /^\s?\n$/);

  chomp($sid);
  print $outfh "\n$sid\n\n";

  my $archive_id = $aa->fetch_by_stable_id($sid);
  
  unless ($archive_id) {
    print $outfh "  Not found in database.\n";
    next;
  }

  my $history = $archive_id->get_history_tree;
  next unless $history;

  if ($history->is_incomplete) {
    print $outfh "  NOTE: History tree is incomplete.\n\n";
  }

  my $matrix = [];

  # get unique stable IDs (regardless of version)
  my @unique_ids = @{ $history->get_unique_stable_ids };

  my $i = 0;
  foreach my $id (@unique_ids) {
    $matrix->[$i++]->[0] = $id;
  }

  # get all releases for which we have nodes in this graph
  my @releases  = @{ $history->get_release_display_names };

  my $j = 1;
  foreach my $release (@releases) {
    $matrix->[scalar(@unique_ids)]->[$j++] = $release;
  }

  # print a "graphical" representation of the tree
  my $fmt = "  %-20s" . ("%-6s" x scalar(@releases)) . "\n";

  foreach my $a_id (@{ $history->get_all_ArchiveStableIds }) {
    my ($x, $y) = @{ $history->coords_by_ArchiveStableId($a_id) };
    $matrix->[$y]->[$x+1] = $a_id->version;
  }

  for (my $i = 0; $i < @$matrix; $i++) {
    print $outfh sprintf($fmt, @{ $matrix->[$i] });
  }

  # current versions in history
  print $outfh "\n  Current stable IDs in this tree:\n";
  my @current = @{ $history->get_all_current_ArchiveStableIds };
  if (@current) {
    map { print $outfh "    ".$_->stable_id.".".$_->version."\n" } @current;
  } else {
    print $outfh "    none\n";
  }

}

close($infh);
close($outfh);


