use strict;

BEGIN { $| = 1;  
	use Test ;
	plan tests => 7
}

use Bio::EnsEMBL::Map::MarkerSynonym;

use Bio::EnsEMBL::Test::TestUtils;

our $verbose = 0; #set to 1 to turn on debug printouts


######
# 1  #
######

#test constructor
my $source = 'genbank';
my $name = 'DS1234';
my $dbID = 10;


my $ms = Bio::EnsEMBL::Map::MarkerSynonym->new($dbID, $source, $name);

ok($ms && ref $ms && $ms->isa('Bio::EnsEMBL::Map::MarkerSynonym'));

#######
# 2-3 #
#######

#test source

ok($source eq $ms->source);
ok(&test_getter_setter($ms, 'source', 'uniSTS'));

#######
# 4-5 #
#######

#test name
ok($name eq $ms->name);
ok(&test_getter_setter($ms, 'name', 'new_name'));


#######
# 6-7 #
#######

#test dbID
ok($dbID == $ms->dbID);
ok(&test_getter_setter($ms, 'dbID', 123)); 

