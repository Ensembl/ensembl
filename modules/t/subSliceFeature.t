use strict;
use warnings;
use Test::More;
use Bio::EnsEMBL::Test::MultiTestDB;
use Bio::EnsEMBL::SubSlicedFeature;

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new();
my $dba = $multi->get_DBAdaptor('core');

my $ga = $dba->get_GeneAdaptor();
my $gene = $ga->fetch_by_stable_id('ENSG00000125964');

diag($gene->stable_id);
diag($gene->start);
diag($gene->end);

my $transcript_list = $gene->get_all_Transcripts;

diag(scalar(@$transcript_list));
foreach (@$transcript_list) {
    diag($_->stable_id);
    diag($_->start);
    diag($_->end);
}
my $fake_gene = Bio::EnsEMBL::SubSlicedFeature->new(-feature => $gene, -start => 30840810, -end => 30859270);

$transcript_list = $fake_gene->get_all_Transcripts;
diag (scalar(@$transcript_list));
foreach (@$transcript_list) {
    diag($_->stable_id);
    diag($_->start);
    diag($_->end);
}

is ($transcript_list->[0]->stable_id,"ENST00000216932", "Only one transcript found in subsliced Gene");

my $exon_list = $fake_gene->get_all_Exons;
foreach (@$exon_list) {
    diag($_->stable_id);
    diag($_->start);
    diag($_->end);
}
ok(scalar(@$exon_list) == 4, "Correct Exons for subsliced Gene");
is($exon_list->[0]->stable_id,'ENSE00001048819', "Correct Exon returned for subsliced Gene");

$fake_gene = Bio::EnsEMBL::SubSlicedFeature->new(-feature => $gene, -start => 1, -end => 2);
$transcript_list = $fake_gene->get_all_Transcripts;
ok (scalar(@$transcript_list) == 0, "Out of bounds search for features");

# Check normal feature functions are unaffected.

ok($fake_gene->start == 30822726, "Normal feature function behaves through proxy");

done_testing;