use strict;
use warnings;
BEGIN { $| = 1;
	use Test;
	plan tests => 17;
}

use Bio::EnsEMBL::Test::MultiTestDB;
use Bio::EnsEMBL::Test::TestUtils;
use Bio::EnsEMBL::DBSQL::OperonAdaptor;
use Bio::EnsEMBL::DBSQL::OperonTranscriptAdaptor;
debug( "Startup test" );
ok(1);

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new();

my $dba = $multi->get_DBAdaptor( "core" );

debug( "Test database instatiated" );
ok( $dba );
									  
my $operon_adaptor =  Bio::EnsEMBL::DBSQL::OperonAdaptor->new($dba);

# get a named operon
my $operon = $operon_adaptor->fetch_by_name("16152-16153-4840");
ok(defined $operon);
	debug("O ".$operon->dbID());
		ok(defined $operon->analysis());
# iterate over its transcripts
my $transcripts = $operon->get_all_OperonTranscripts();
ok(defined $transcripts);
ok(scalar(@$transcripts)>0);
for my $ot (@$transcripts) {
	ok(defined $ot->analysis());
	debug("OT ".$ot->dbID());
	for my $gene (@{$ot->get_all_Genes()}) {
	debug("G ".$gene->dbID());
	}
}

									  
# get a slice
my $slice = $dba->get_SliceAdaptor()->fetch_by_seq_region_id(469283);
my $operons = $operon_adaptor->fetch_all_by_Slice($slice);
ok(defined $operons);
ok(scalar(@$operons)>0);
for my $o (@$operons) {
		ok(defined $o->analysis());
	debug("O ".$o->dbID());
	for my $ot (@{$o->get_all_OperonTranscripts()}) {
	debug("OT ".$ot->dbID());
	for my $gene (@{$ot->get_all_Genes()}) {
	debug("G ".$gene->dbID());
	}
}
}

my $operon_transcript_adaptor =  Bio::EnsEMBL::DBSQL::OperonTranscriptAdaptor->new($dba);
# get a named operon
my $ot = $operon_transcript_adaptor->fetch_by_name("T16152-16153-4840");
ok(defined $ot);
	debug("OT ".$ot->dbID());
	debug("OTP ".$ot->operon()->dbID());
	for my $gene (@{$ot->get_all_Genes()}) {
	debug("G ".$gene->dbID());
	}
	
my $ots = $operon_transcript_adaptor->fetch_all_by_Slice($slice);
ok(defined $ots);
ok(scalar(@$ots)>0);	

# get an operon transcript by gene and then find the operons
my $gene_adaptor = Bio::EnsEMBL::DBSQL::GeneAdaptor->new($dba);
my ($gene) = @{$gene_adaptor->fetch_all_by_external_name('16152')};
ok(defined $gene);
	debug("GQ ".$gene->dbID());
$ots = $operon_transcript_adaptor->fetch_all_by_gene($gene);
ok(defined $ots && scalar(@$ots)>0);	
for my $ot (@$ots) {
		debug("OT ".$ot->dbID());
		for my $gene (@{$ot->get_all_Genes()}) {
	debug("G ".$gene->dbID());
	}
}

