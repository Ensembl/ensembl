use strict;
use warnings;

use lib 't';

BEGIN { $| = 1;  
	use Test;
	plan tests => 4;
}


use MultiTestDB;
use TestUtils qw(debug test_getter_setter);
use Bio::EnsEMBL::Chromosome;


#
#1 TEST - Chromosome Compiles
#
ok(1); 


my $CHR           = '20';
my $DBID          = 123;
my $LENGTH        = 250_000_000;

my $multi_db = MultiTestDB->new;
my $db = $multi_db->get_DBAdaptor('core');

my $ca = $db->get_ChromosomeAdaptor;

#
# Chromosome constructor
#
my $chromosome = Bio::EnsEMBL::Chromosome->new
  (-chr_name => $CHR,
   -dbID     => $DBID,
   -adaptor  => $ca,
   -length   => $LENGTH);

ok($chromosome->isa('Bio::EnsEMBL::Chromosome'));
ok($chromosome->adaptor == $ca);
ok($chromosome->length());
ok($chromosome->name());
