use strict;
use warnings;

use Test::More;
use Bio::EnsEMBL::Test::MultiTestDB;
use Bio::EnsEMBL::Utils::IO::GTFSerializer;
use IO::String;

my $db = Bio::EnsEMBL::Test::MultiTestDB->new();
my $dba = $db->get_DBAdaptor('core');

my $id = 'ENSG00000131044';

my $ga = $dba->get_GeneAdaptor();
my $gene = $ga->fetch_by_stable_id($id);

SKIP: {
  my $fh = IO::String->new();
  my $ser = Bio::EnsEMBL::Utils::IO::GTFSerializer->new($fh);
  # $ser->print_main_header([$gene->feature_Slice()]);
  $ser->print_feature($gene);
  
  my $expected = <<'OUT';
##gff-version 3
##sequence-region   20 30274334 30300924
OUT
  #Have to do this outside of the HERETO thanks to tabs
  $expected .= join("\t", 
    qw/20  EnsEMBL feature 30274334  30300924  . + ./,
    'ID=ENSG00000131044;logic_name=ensembl;external_name=C20orf125;description=DJ310O13.1.2 (NOVEL PROTEIN SIMILAR DROSOPHILA PROTEIN CG7474%2C ISOFORM 2 ) (FRAGMENT). [Source:SPTREMBL%3BAcc:Q9BR18];biotype=protein_coding' 
  );
  $expected .= "\n";

  is(${$fh->string_ref()}, $expected, 'Gene serialises to GFF3 as expected');
}

done_testing();