use strict;
use Bio::EnsEMBL::Test::TestUtils;
use Bio::EnsEMBL::Test::MultiTestDB;

use Bio::EnsEMBL::RegulatoryFactor;

BEGIN { $| = 1;
	use Test;
	plan tests => 11;
}

my $verbose = 0;

# get a core DBAdaptor and a RegualtoryFactorAdaptor
#
my $multi = Bio::EnsEMBL::Test::MultiTestDB->new;
my $dba = $multi->get_DBAdaptor("core");
my $factor_adaptor = $dba->get_RegulatoryFactorAdaptor();


#
# Test constructor
#
my $rm = Bio::EnsEMBL::RegulatoryFactor->new(-NAME => 'Joe',
					    -TYPE => 'promoter');

ok($rm);
ok(ref($rm));
ok($rm->isa('Bio::EnsEMBL::RegulatoryFactor'));

ok($rm->name eq 'Joe');
ok($rm->type eq 'promoter');

ok(test_getter_setter($rm,'name','Fred'));
ok(test_getter_setter($rm,'type','miRNA_target'));


# test coding_transcript & gene
my $factor = $factor_adaptor->fetch_by_dbID(1);
my $transcript = $factor->coding_transcript();
ok($transcript->dbID() == 21716);
$factor = $factor_adaptor->fetch_by_dbID(5);
my $gene = $factor->coding_gene();
ok($gene->dbID() == 18271);

# test getting factors coded for by gene and transcript
my $gene_adaptor = $dba->get_GeneAdaptor();
my $gene = $gene_adaptor->fetch_by_dbID(18271);
my @factors = @{$factor_adaptor->fetch_factors_coded_for_by_gene($gene)};
ok($factors[0]->dbID() == 5);

my $transcript_adaptor = $dba->get_TranscriptAdaptor();
my $transcript = $transcript_adaptor->fetch_by_dbID(21717);
my @factors = @{$factor_adaptor->fetch_factors_coded_for_by_transcript($transcript)};
ok($factors[0]->dbID() == 2);
