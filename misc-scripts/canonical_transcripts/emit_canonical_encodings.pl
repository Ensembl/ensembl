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


use strict;
use warnings;

use Getopt::Long;

use Bio::EnsEMBL::Utils::Exception;
use Bio::EnsEMBL::Utils::TranscriptSelector;
use Bio::EnsEMBL::DBSQL::DBAdaptor;

my ($host, $port, $dbname, $user,$pass);
my ($dnahost, $dnaport, $dnadbname, $dnauser, $dnapass);
my ($ccds_host, $ccds_dbname, $ccds_user, $ccds_port, $ccds_pass);

my $transcript;
my $gene;

GetOptions( 'dbhost:s'            => \$host,
            'dbport:n'            => \$port,
            'dbname:s'            => \$dbname,
            'dbuser:s'            => \$user,
            'dbpass:s'            => \$pass,
            'dnahost:s'           => \$dnahost,
            'dnadbname:s'         => \$dnadbname,
            'dnauser:s'           => \$dnauser,
            'dnapass:s'           => \$dnapass,
            'dnaport:s'           => \$dnaport,
            'ccdshost:s'          => \$ccds_host,
            'ccdsdbname:s'        => \$ccds_dbname,
            'ccdsuser:s'          => \$ccds_user,
            'ccdsport:s'          => \$ccds_port,
            'ccdspass:s'          => \$ccds_pass,
            'transcript:s'        => \$transcript,
            'gene:s'              => \$gene,
            );


my $dba =
  new Bio::EnsEMBL::DBSQL::DBAdaptor( -host   => $host,
                                      -user   => $user,
                                      -port   => $port,
                                      -dbname => $dbname,
                                      -pass   => $pass,
                                      -species => 'default',
                                      );
                                      
if($dnadbname) {
  if(!$dnauser || !$dnahost) {
    throw ("You must provide user, host and dbname details to connect to DNA DB!");
  }
  my $dna_db =   new Bio::EnsEMBL::DBSQL::DBAdaptor( -host   => $dnahost,
                                      -user       => $dnauser,
                                      -port       => $dnaport,
                                      -dbname     => $dnadbname,
                                      -pass       => $dnapass,
                                      -species    => 'dna_'.$dba->species()
                                      );
  $dba->dnadb($dna_db);
}

my $ccds_dba;

if ($ccds_dbname) {
  if (!$ccds_user || !$ccds_host) {
    throw ("You must provide user, host and dbname details to connect to CCDS DB!");
  }
  $ccds_dba = 
  new Bio::EnsEMBL::DBSQL::DBAdaptor( -host   => $ccds_host,
                                      -user   => $ccds_user,
                                      -port   => $ccds_port,
                                      -pass   => $ccds_pass,
                                      -dbname => $ccds_dbname,
                                      -species => 'ccds_'.$dba->species() );
}

my $transcript_selector = Bio::EnsEMBL::Utils::TranscriptSelector->new($ccds_dba, 1);

my $gene_object;

if($gene) {
  printf 'Using "%s" as our Gene stable ID'."\n", $gene;
  $gene_object = $dba->get_GeneAdaptor()->fetch_by_stable_id($gene);
}
elsif($transcript) {
  printf 'Using "%s" as a Transcript to find a Gene'."\n", $transcript;
  my $t = $dba->get_TranscriptAdaptor()->fetch_by_stable_id($transcript);
  $gene_object = $t->get_Gene();
  printf 'Using "%s" as our Gene'."\n", $gene_object->stable_id();
}

$transcript_selector->select_canonical_transcript_for_Gene($gene_object);
print "Original: ".$gene_object->canonical_transcript->stable_id."\n";

sub usage {
print "
Example usage: perl emit_canonical_encodings.pl -dbhost host -dbuser user 
     -dbpass *** -dbname dbname -dbport 3306 -transcript ENST00019112
     
Example usage: perl emit_canonical_encodings.pl -dbhost host -dbuser user 
     -dbpass *** -dbname dbname -dbport 3306 -gene ENSG00019111

Script options:

    -dbname       Database name

    -dbhost       Database host

    -dbport       Database port

    -dbuser       Database user

    -dbpass       Database password

Optional DB connection arguments:

    -dnadbname    DNA Database name

    -dnadbhost    DNA Database host

    -dnadbuser    DNA Database user
    
    -dnadbport    DNA Database port
    
    -dnadbpass    DNA Database pass

    -ccdsdbname  CCDS database name

    -ccdshost    CCDS database host

    -ccdsuser    CCDS database user
    
    -ccdspass    CCDS database pass
    
    -ccdsport    CCDS database port

Search params:

    -transcript    The transcript stable ID to use to find the gene in question

    -gene          The gene stable ID to emit encoded transcripts for

A warning about not using CCDS is perfectly acceptible when not running on
Human, Mouse and Zebrafish.
";
    
}
