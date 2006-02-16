
use strict;

BEGIN { $| = 1;
        use Test;
        plan tests => 10;
}

use Bio::EnsEMBL::Test::MultiTestDB;
use Bio::EnsEMBL::Test::TestUtils;
use Bio::EnsEMBL::UnmappedObject;
use Bio::EnsEMBL::DBSQL::UnmappedObjectAdaptor;

our $verbose = 0;

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new();
ok( $multi );

my $db = $multi->get_DBAdaptor( "core" );

my $uma = $db->get_UnmappedObjectAdaptor();

ok(ref($uma));

#
# test fetch operations
#

my @objects = @{$uma->fetch_all()};

ok(scalar(@objects) == 4);

@objects = @{$uma->fetch_all_by_type('xref')};

ok(scalar(@objects) == 2);

my $analysis_adaptor = $db->get_AnalysisAdaptor();

my $analysis = $analysis_adaptor->fetch_by_logic_name("Unigene");

if(!defined($analysis)){
   die "ARSE\n";
}		
@objects = @{$uma->fetch_all_by_analysis($analysis)};

ok(scalar(@objects) == 2);

@objects = @{$uma->fetch_all_by_analysis($analysis,"UniGene")};

ok(scalar(@objects) == 2);

@objects = @{$uma->fetch_all_by_analysis($analysis,"RFAM")};

ok(scalar(@objects) == 0);

@objects = @{$uma->fetch_by_identifier("X5678")};

ok(scalar(@objects) == 1);

@objects = @{$uma->fetch_by_identifier("X5678","UniGene")};

ok(scalar(@objects) == 1);

@objects = @{$uma->fetch_by_identifier("X5678","RFAM")};

ok(scalar(@objects) == 0);

