use strict;
use warnings;

use lib 't';

BEGIN { $| = 1;  
	use Test;
	plan tests => 14;
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
my $KNOWN_GENES   = 1300;
my $UNKNOWN_GENES = 300;
my $SNPS          = 300_000;

my $multi_db = MultiTestDB->new;
my $db = $multi_db->get_DBAdaptor('core');

my $ca = $db->get_ChromosomeAdaptor;

#
#2-8 Chromosome constructor
#
my $chromosome = Bio::EnsEMBL::Chromosome->new
  (-chr_name => $CHR,
   -dbID     => $DBID,
   -adaptor  => $ca,
   -length   => $LENGTH,
   -known_genes => $KNOWN_GENES,
   -unknown_genes => $UNKNOWN_GENES,
   -snps     => $SNPS);

ok($chromosome->isa('Bio::EnsEMBL::Chromosome'));
ok($chromosome->length == $LENGTH);
ok($chromosome->dbID  == $DBID);
ok($chromosome->adaptor == $ca);
ok($chromosome->known_genes == $KNOWN_GENES);
ok($chromosome->unknown_genes == $UNKNOWN_GENES);
ok($chromosome->snps == $SNPS);

#
# 9-14 getter /setters
#
ok(test_getter_setter($chromosome, 'chr_name', 'Y'));
ok(test_getter_setter($chromosome, 'dbID', 321));
ok(test_getter_setter($chromosome, 'length', 123_000_000));
ok(test_getter_setter($chromosome, 'known_genes', 400));
ok(test_getter_setter($chromosome, 'unknown_genes', 1200));
ok(test_getter_setter($chromosome, 'snps', 12_345));






