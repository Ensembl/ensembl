=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=cut
use strict;
use warnings;

use File::Spec;
use File::Basename qw/dirname/;
use File::Temp qw/tempdir/;
use Test::More;
use Test::Warnings;
use Bio::EnsEMBL::Test::MultiTestDB;
use Bio::EnsEMBL::Test::TestUtils;
use XrefParser::Database;

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new;
my $dba = $multi->get_DBAdaptor("xref");
plan skip_all => "xref database schema is mysql specific - won't work with a different driver"
   unless $dba->dbc->driver eq 'mysql';

my $database = XrefParser::Database->new( $dba->dbc);

$database->populate(
  File::Spec->catdir($multi->curr_dir, "../../misc-scripts/xref_mapping"),
  "with force please",
);

my %xref_tables_expected_empty_by_default = (
  checksum_xref=>0,
  coordinate_xref=>0,
  dependent_xref=>0,
  gene_direct_xref=>0,
  go_xref=>0,
  identity_xref=>0,
  object_xref=>0,
  primary_xref=>0,
  transcript_direct_xref=>0,
  translation_direct_xref=>0,
  xref=>0,
);
my $tmp_dir = tempdir(CLEANUP=>1);
sub store_in_temporary_file {
  my $path = "$tmp_dir/tmp";
  open(my $fh, ">", $path) or die $path;  
  print $fh @_;
  close($fh);
  return $path;
}
sub test_parser {
  my ($parser, $content, $source_id, $expected, $test_name) = @_;
  require_ok($parser);
  $parser->new($database)->run({
   files => [store_in_temporary_file($content)],
   source_id => $source_id,
   species_id => 1 #Happens to be right, but doesn't matter anyway - we are not testing the mapping
  });
  my $expected_table_counts = {%xref_tables_expected_empty_by_default, %$expected};
  subtest "$parser $test_name" => sub {
    plan tests => scalar(keys %$expected_table_counts);
    for my $table (keys %$expected_table_counts){
      my $actual_count = count_rows($dba, $table);
      $dba->dbc->prepare("delete from $table;")->execute() if $actual_count;
      my $expected_count = $expected_table_counts->{$table};
      is($actual_count, $expected_count, "$table has $expected_count rows") or diag "$table has $actual_count rows";
    }
  }
}
test_parser("XrefParser::WormbaseDirectParser", "", "source_id (unused)", {}, "null case");
my $wormbase_celegans_xrefs_head= <<EOF;
//
// WormBase Caenorhabditis elegans XREFs for WS265
//
// Columns (tab separated) are:
//    1. WormBase Gene sequence name
//    2. WormBase Gene accession
//    3. WormBase Gene CGC name
//    4. WormBase Transcript sequence name
//    5. WormPep protein accession
//    6. INSDC parent sequence accession
//    7. INSDC locus_tag id
//    8. INSDC protein_id
//    9. UniProt accession
//
// Missing or not applicable data (e.g. protein identifiers for non-coding RNAs) is denoted by a "."
//
2L52.1	WBGene00007063	.	2L52.1b	CE50569	BX284602	CELE_2L52.1	CTQ86426	A0A0K3AWR5
2L52.1	WBGene00007063	.	2L52.1a	CE32090	BX284602	CELE_2L52.1	CCD61130	A4F336
2L52.2	WBGene00200402	.	2L52.2	.	BX284602	CELE_2L52.2	.	.
EOF
test_parser("XrefParser::WormbaseDirectParser", $wormbase_celegans_xrefs_head, "source_id (unused)", {
xref=>9,
gene_direct_xref => 6,
transcript_direct_xref => 3,
translation_direct_xref => 2,
}, "Direct xrefs: genes: count currently off due to some questionable duplicates, transcripts: as in column 4, translations: as in column 5. At least one direct xref per xref (but should be one to one)");
done_testing();

