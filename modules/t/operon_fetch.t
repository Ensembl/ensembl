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

use Test::More;

use Bio::EnsEMBL::Test::MultiTestDB;
use Bio::EnsEMBL::Test::TestUtils;
use Bio::EnsEMBL::DBSQL::OperonAdaptor;
use Bio::EnsEMBL::DBSQL::OperonTranscriptAdaptor;
note( "Startup test" );
ok(1);

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new();

my $dba = $multi->get_DBAdaptor( "core" );

note( "Test database instatiated" );
ok( $dba );
									  
my $operon_adaptor =  Bio::EnsEMBL::DBSQL::OperonAdaptor->new($dba);

# get a named operon
my $operon = $operon_adaptor->fetch_by_name("16152-16153-4840");
ok(defined $operon);
	note("O ".$operon->dbID());
		ok(defined $operon->analysis());
# iterate over its transcripts
my $transcripts = $operon->get_all_OperonTranscripts();
ok(defined $transcripts);
ok(scalar(@$transcripts)>0);
for my $ot (@$transcripts) {
	ok(defined $ot->analysis());
	note("OT ".$ot->dbID());
	for my $gene (@{$ot->get_all_Genes()}) {
	note("G ".$gene->dbID());
	}
}

									  
# get a slice
my $slice = $dba->get_SliceAdaptor()->fetch_by_seq_region_id(469283);
my $operons = $operon_adaptor->fetch_all_by_Slice($slice);
ok(defined $operons);
ok(scalar(@$operons)>0);
for my $o (@$operons) {
		ok(defined $o->analysis());
	note("O ".$o->dbID());
	for my $ot (@{$o->get_all_OperonTranscripts()}) {
	note("OT ".$ot->dbID());
	for my $gene (@{$ot->get_all_Genes()}) {
	note("G ".$gene->dbID());
	}
}
}

my $operon_transcript_adaptor =  Bio::EnsEMBL::DBSQL::OperonTranscriptAdaptor->new($dba);
# get a named operon
my $ot = $operon_transcript_adaptor->fetch_by_name("T16152-16153-4840");
ok(defined $ot);
	note("OT ".$ot->dbID());
	note("OTP ".$ot->operon()->dbID());
	for my $gene (@{$ot->get_all_Genes()}) {
	note("G ".$gene->dbID());
	}
	
my $ots = $operon_transcript_adaptor->fetch_all_by_Slice($slice);
ok(defined $ots);
ok(scalar(@$ots)>0);	

# get an operon transcript by gene and then find the operons
my $gene_adaptor = Bio::EnsEMBL::DBSQL::GeneAdaptor->new($dba);
my ($gene) = @{$gene_adaptor->fetch_all_by_external_name('16152')};
ok(defined $gene);
	note("GQ ".$gene->dbID());
$ots = $operon_transcript_adaptor->fetch_all_by_gene($gene);
ok(defined $ots && scalar(@$ots)>0);	
for my $ot (@$ots) {
		note("OT ".$ot->dbID());
		for my $gene (@{$ot->get_all_Genes()}) {
	note("G ".$gene->dbID());
	}
}

done_testing();
