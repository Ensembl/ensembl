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

# New xref system maps xrefs to ensembl translations, transcripts etc
# Old mappings are always to translations, so need to convert new mappings
# to be all translations in order to facilitate comparison.

use strict;
use warnings;

use Getopt::Long;
use Bio::EnsEMBL::DBSQL::DBAdaptor;

my $host   = 'ecs2';
my $port   = 3365;
my $user   = 'ensro';
my $dbname = 'homo_sapiens_core_26_35';
my $file;

GetOptions ('file=s'  => \$file);

die "Must specify -file argument" if (!$file);

my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(-host   => $host,
                                            -user   => $user,
					    -port   => $port,
                                            -dbname => $dbname);

my $transcript_adaptor = $db->get_TranscriptAdaptor();
my $translation_adaptor = $db->get_TranslationAdaptor();

# loop over each mapping in the file and convert if necessary

open(FILE, $file) || die "Could not open $file\n";

#print <FILE>; # echo header

while (<FILE>) {

  my ($xref_id, $type, $ensembl_id) = split;

  # Just echo Translation mappings unchanged
  print if ($type =~ /Translation/i);

  # Convert transcripts to translations
  if ($type =~ /Transcript/i) {
    my $translation_id = get_translation_for_transcript($ensembl_id);
    print "$xref_id\tTranslation\t$translation_id\n" if ($translation_id);
  }

  # Other types here ...
}

close(FILE);

# ----------------------------------------------------------------------

sub get_translation_for_transcript {

  my $transcript_id = shift;

  my $transcript = $transcript_adaptor->fetch_by_dbID($transcript_id);
  my $translation = $translation_adaptor->fetch_by_Transcript($transcript);

  return undef if (!$translation);

  return $translation->dbID();

}

