use strict;
use warnings;


BEGIN { $| = 1;  
	use Test;
	plan tests => 5;
}


use Bio::EnsEMBL::Test::MultiTestDB;
use Bio::EnsEMBL::Test::TestUtils;

use Bio::EnsEMBL::Utils::Exception qw(verbose);
use Bio::EnsEMBL::Chromosome;


######################################################################
# 
# Chromosome is a deprecated class but needed for backwards 
# compatibility.  These tests ensure that it actually works,
# but verbosity is turned off to avoid all of the deprecated warnings
#
#######################################################################

verbose(-1);

#
#1 TEST - Chromosome Compiles
#
ok(1); 


my $CHR           = '20';
my $DBID          = 123;
my $LENGTH        = 250_000_000;

my $multi_db = Bio::EnsEMBL::Test::MultiTestDB->new;
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


verbose(0);