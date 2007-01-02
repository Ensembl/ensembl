use strict;
use Bio::EnsEMBL::Test::TestUtils;

use Bio::EnsEMBL::Test::MultiTestDB;


BEGIN { $| = 1;
	use Test;
	plan tests => 6;
}

my $multi_db = Bio::EnsEMBL::Test::MultiTestDB->new;
my $db = $multi_db->get_DBAdaptor('core');

my $verbose = 0;

# Test Creation

my $gene = $db->get_GeneAdaptor->fetch_by_dbID(18256);
my $transcript = $db->get_TranscriptAdaptor()->fetch_by_dbID(21718);

my $utaa = $db->get_UnconventionalTranscriptAssociationAdaptor();

ok(ref($utaa) && $utaa->isa('Bio::EnsEMBL::DBSQL::UnconventionalTranscriptAssociationAdaptor'));

# test fetch all by interaction type
ok(@{$utaa->fetch_all_by_interaction_type('antisense')} == 3);

# test fetch all by gene, with and without type
ok(@{$utaa->fetch_all_by_gene($gene)} == 3);
ok(@{$utaa->fetch_all_by_gene($gene, 'antisense')} == 2);

# test fetch all by transcript, with and without type
ok(@{$utaa->fetch_all_by_transcript($transcript)} == 2);
ok(@{$utaa->fetch_all_by_transcript($transcript, 'antisense')} == 1);

# TODO - test functionality in gene.t
# Test store

$multi_db->save('core', 'unconventional_transcript_association');


my $uta = Bio::EnsEMBL::UnconventionalTranscriptAssociation->new($transcript, $gene, 'antisense');

$utaa->store($uta);

# TODO - check contents

$multi_db->restore('core');
