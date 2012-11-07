use strict;

use Bio::EnsEMBL::Test::TestUtils;

use Test::More;

use Bio::EnsEMBL::Test::MultiTestDB;
use Bio::EnsEMBL::UnconventionalTranscriptAssociation;


our $verbose = 0;

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new;

# get a core DBAdaptor
#
my $db = $multi->get_DBAdaptor("core");
my $utaa = $db->get_UnconventionalTranscriptAssociationAdaptor;

my $gene = $db->get_GeneAdaptor->fetch_by_stable_id('ENSG00000355555');
my $transcript = my $ta = $db->get_TranscriptAdaptor()->fetch_by_stable_id("ENST00000217347");

#
# 1 create a new UnconventionalTranscriptAssociation
#
my $uta = new Bio::EnsEMBL::UnconventionalTranscriptAssociation($transcript, $gene, 'antisense');
ok($uta);

#
# 2-4 test the basic getter and setters
#

# 2 gene
ok(test_getter_setter($uta, 'gene', $gene));

# 3 transcript
ok(test_getter_setter($uta, 'transcript', $transcript));

# 4 type
ok(test_getter_setter($uta, 'interaction_type', 'antisense'));


done_testing();