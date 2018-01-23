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

## no critic (RequireFilenameMatchesPackage)
package Test::SO::Term;
use strict;
use warnings;

sub new {
  my ($class) = @_;
  return bless({}, ref($class) || $class);
}

# default lookup and always returns region
sub name {
  my ($self) = @_;
  return 'region';
}

package Test::SO;

use base qw/Bio::EnsEMBL::DBSQL::OntologyTermAdaptor/;

sub new {
  my ($class) = @_;
  return bless({}, ref($class) || $class);
}

sub fetch_by_accession {
  my ($self) = @_;
  return Test::SO::Term->new();
}

package main;

use strict;
use warnings;
use Test::More;
use Test::Warnings;
use Bio::EnsEMBL::Test::MultiTestDB;
use Bio::EnsEMBL::Utils::IO::GFFSerializer;
use Bio::EnsEMBL::Feature;
use Bio::EnsEMBL::Slice;
use IO::String;
use Test::Differences;

my $db = Bio::EnsEMBL::Test::MultiTestDB->new();
my $dba = $db->get_DBAdaptor('core');
my $omulti = Bio::EnsEMBL::Test::MultiTestDB->new('ontology');
my $odb = $omulti->get_DBAdaptor('ontology');

my $id = 'ENSG00000131044';

my $ga = $dba->get_GeneAdaptor();

{
  my $gene = $ga->fetch_by_stable_id($id);
  delete $gene->{source};
  $gene->{description} = undef; #empty value means don't emit the key
  my $expected = <<'OUT';
##gff-version 3
##sequence-region   20 30274334 30300924
OUT
  #Have to do this outside of the HERETO thanks to tabs
  $expected .= join("\t", 
    qw/20  ensembl region 30274334  30300924  . + ./,
    'ID=gene:ENSG00000131044;Name=C20orf125;biotype=protein_coding;gene_id=ENSG00000131044;logic_name=ensembl;projection_parent_gene=ENSG_PARENT_GENE;version=1'
  );
  $expected .= "\n";

  assert_gff3($gene, $expected, 'Gene with no source serialises to GFF3 as expected. Source is ensembl');
}

{
  my $cs = $dba->get_CoordSystemAdaptor()->fetch_by_name('chromosome');
  my $feature = Bio::EnsEMBL::Feature->new(
    -SLICE => Bio::EnsEMBL::Slice->new(
      -COORD_SYSTEM => $cs,
      -SEQ => ('A'x10),
      -SEQ_REGION_NAME => 'wibble',
      -START => 1,
      -END => 10
    ),
    -START => 1,
    -END => 10,
    -STRAND => 1,
  );
  my $expected = <<'OUT';
##gff-version 3
##sequence-region   wibble 1 10
OUT
  #Have to do this outside of the HERETO thanks to tabs
  $expected .= join("\t", 
    qw/wibble  . region 1  10  . + ./,
    '' 
  );
  $expected .= "\n";

  assert_gff3($feature, $expected, 'Default feature should seralise without attributes but leave a trailing \t');
}


{
  my $gene = $ga->fetch_by_stable_id($id);
  $gene->source('wibble');
  my $expected = <<'OUT';
##gff-version 3
##sequence-region   20 30274334 30300924
OUT
  #Have to do this outside of the HERETO thanks to tabs
  $expected .= join("\t", 
    qw/20  wibble region 30274334  30300924  . + ./,
    'ID=gene:ENSG00000131044;Name=C20orf125;biotype=protein_coding;description=DJ310O13.1.2 (NOVEL PROTEIN SIMILAR DROSOPHILA PROTEIN CG7474%2C ISOFORM 2 ) (FRAGMENT). [Source:SPTREMBL%3BAcc:Q9BR18];gene_id=ENSG00000131044;logic_name=ensembl;projection_parent_gene=ENSG_PARENT_GENE;version=1'
  );
  $expected .= "\n";
  assert_gff3($gene, $expected, 'Gene with custom source serialises to GFF3 as expected. Source is wibble');
  $expected = <<'OUT';
##gff-version 3
##sequence-region   20 30274334 30298904
OUT
  $expected .= join("\t", 
  qw/20      ensembl region  30274334        30298904        .       +       ./,
  'ID=transcript:ENST00000310998;Name=C20orf125;Parent=gene:ENSG00000131044;biotype=protein_coding;logic_name=ensembl;projection_parent_transcript=ENSG_PARENT_TRANSCRIPT;transcript_id=ENST00000310998;version=1'
  );
  $expected .= "\n";
  assert_gff3($gene->canonical_transcript(), $expected, 'Transcript with custom source serialises to GFF3 as expected. Source is wibble');

  $expected = <<'OUT';
##gff-version 3
##sequence-region   20 30274334 30274425
OUT
  $expected .= join("\t",
  qw/20      ensembl region  30274334        30274425        .       +       0/,
  'ID=region:ENSP00000308980;Parent=transcript:ENST00000310998;protein_id=ENSP00000308980'
  );
  $expected .= "\n";
  my $cds = $gene->canonical_transcript->get_all_CDS();
  assert_gff3($cds->[0], $expected, 'CDS with custom source serialises to GFF3 as expected. Source is wibble');

  my $cds_expected = $expected;

  $expected = <<'OUT';
##gff-version 3
##sequence-region   20 30274334 30274425
OUT
  $expected .= join("\t",
  qw/20 ensembl       region  30274334        30274425        .       +       ./,
  'Name=ENSE00001155821;Parent=transcript:ENST00000310998;constitutive=0;ensembl_end_phase=2;ensembl_phase=0;exon_id=ENSE00001155821;rank=1;version=1'
  );
  $expected .= "\n";
  my $exon = $gene->canonical_transcript->get_all_ExonTranscripts();
  assert_gff3($exon->[0], $expected, 'Exon with custom source serialises to GFF3 as expected. Source is wibble');

  my $new_id = 'ENSG00000126003';
  my $utr_gene = $ga->fetch_by_stable_id($new_id);
  my $utrs = $utr_gene->canonical_transcript->get_all_five_prime_UTRs();
  $expected = <<'OUT';
##gff-version 3
##sequence-region   20 30583501 30583588
OUT
  $expected .= join("\t",
  qw/20 ensembl       region  30583501  30583588        .       -       ./,
  'Parent=transcript:ENST00000246229'
  );
  $expected .= "\n";
  assert_gff3($utrs->[0], $expected, 'UTR feature serialises to GFF3 as expected');

  $cds_expected .= join("\t",
  qw/20 ensembl       region  30274334        30274425        .       +       ./,
  'Name=ENSE00001155821;Parent=transcript:ENST00000310998;constitutive=0;ensembl_end_phase=2;ensembl_phase=0;exon_id=ENSE00001155821;rank=1;version=1'
  );
  $cds_expected .= "\n";
  my @list;
  push @list, $cds->[0];
  push @list, $exon->[0];
  assert_gff3_list(\@list, $cds_expected, "List of features serialises to GFF3 as expected");


}

{
  my $gene = $ga->fetch_by_stable_id($id);
  
  my $summary = $gene->summary_as_hash;
  $$summary{'Dbxref'} = ['bibble', 'fibble'];
  $$summary{'Ontology_term'} = 'GO:0001612';
  local undef &{Bio::EnsEMBL::Gene::summary_as_hash};
  local *{Bio::EnsEMBL::Gene::summary_as_hash} = sub {return $summary};
  
  my $expected = <<'OUT';
##gff-version 3
##sequence-region   20 30274334 30300924
OUT
  #Have to do this outside of the HERETO thanks to tabs
  $expected .= join("\t", 
    qw/20  ensembl region 30274334  30300924  . + ./,
    'ID=gene:ENSG00000131044;Name=C20orf125;Dbxref=bibble,fibble;Ontology_term=GO:0001612;biotype=protein_coding;description=DJ310O13.1.2 (NOVEL PROTEIN SIMILAR DROSOPHILA PROTEIN CG7474%2C ISOFORM 2 ) (FRAGMENT). [Source:SPTREMBL%3BAcc:Q9BR18];gene_id=ENSG00000131044;logic_name=ensembl;projection_parent_gene=ENSG_PARENT_GENE;version=1'
  );
  $expected .= "\n";

  assert_gff3($gene, $expected, 'Gene with array- and string-valued attributes as expected.');
}

sub assert_gff3 {
  my ($feature, $expected, $msg) = @_;
  my $ota = Test::SO->new();
  my $fh = IO::String->new();
  my $ser = Bio::EnsEMBL::Utils::IO::GFFSerializer->new($ota, $fh);
  $ser->print_main_header([$feature->feature_Slice()]);
  $ser->print_feature($feature);
  eq_or_diff(${$fh->string_ref()}, $expected, $msg);
}

sub assert_gff3_list {
  my ($features, $expected, $msg) = @_;
  my $ota = Test::SO->new();
  my $fh = IO::String->new();
  my $ser = Bio::EnsEMBL::Utils::IO::GFFSerializer->new($ota, $fh);
  $ser->print_main_header([$features->[0]->feature_Slice()]);
  $ser->print_feature_list($features);
  eq_or_diff(${$fh->string_ref()}, $expected, $msg);
}

done_testing();
