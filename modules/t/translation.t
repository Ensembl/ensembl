use strict;
use warnings;

use lib 't';
use TestUtils qw(debug test_getter_setter);
use Bio::EnsEMBL::Translation;
use Bio::EnsEMBL::Exon;

BEGIN { $| = 1;  
	use Test;
	plan tests => 21;
}

my $loaded = 0;
END {print "not ok 1\n" unless $loaded;}

#turn on/off debug prints:
our $verbose = 0;

use MultiTestDB;

my $multi = MultiTestDB->new();

$loaded = 1;

ok(1);

my $db = $multi->get_DBAdaptor( 'core' );

my $t = Bio::EnsEMBL::Translation->new();

ok($t);

ok(test_getter_setter($t,'stable_id',1));

ok(test_getter_setter($t,'dbID',3));

ok(test_getter_setter($t,'start',42));
ok(test_getter_setter($t,'end',50));

my $exon = Bio::EnsEMBL::Exon->new();
$exon->start(10);
$exon->end(20);
$exon->strand(1);
$exon->phase(0);
$exon->end_phase( -1 );

$t->start_Exon($exon);
ok($t);

$t->end_Exon($exon);
ok($t);


#
# Tests for the translation adaptor
##################################

my $ta = $db->get_TranslationAdaptor();
my $ids = $ta->list_dbIDs();
ok (@{$ids});

my $stable_ids = $ta->list_stable_ids();
ok (@{$stable_ids});


my $tra = $db->get_TranscriptAdaptor();

my $transcript = $tra->fetch_by_stable_id('ENST00000201961');

#
# test fetch_by_Transcript
#
my $translation = $ta->fetch_by_Transcript($transcript);

ok($translation && $translation->stable_id eq 'ENSP00000201961');
ok($translation && $translation->start_Exon->stable_id eq 'ENSE00000661216');
ok($translation && $translation->end_Exon->stable_id eq 'ENSE00000661212');


#
# test fetch_by_dbID
#
$translation = $ta->fetch_by_dbID(21734);
ok($translation && $translation->stable_id() eq 'ENSP00000201961');


#
# test fetch_by_stable_id
#
$translation = $ta->fetch_by_stable_id('ENSP00000201961');
ok($translation && $translation->dbID() == 21734);

#
# test fetch_by_external_name
#
($translation) = @{$ta->fetch_all_by_external_name('CAC33959')};
ok($translation && $translation->dbID() == 21716);

#
# test get_all_ProteinFeatures
#

my @protein_features = @{$translation->get_all_ProteinFeatures()};
debug("Got " . scalar(@protein_features) ." protein features.");
ok(@protein_features == 3);


#
# test get_all_DomainFeatures
#
my @domain_features = @{$translation->get_all_DomainFeatures()};
debug("Got " . scalar(@domain_features) . " domain features.");
ok(@domain_features == 3);

ok($translation->display_id eq $translation->stable_id);



#
# test length() and seq()
#
my $seq = $translation->seq();
debug("Seq = $seq");
ok($seq);

debug("Lenth = " . $translation->length());
ok(length($seq) == $translation->length());

