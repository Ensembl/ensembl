
use lib 't';

BEGIN { $| = 1;  
	use Test;
	plan tests => 12;
}

use MultiTestDB;
use Bio::EnsEMBL::DnaDnaAlignFeature;
use Bio::EnsEMBL::SeqFeature;
use Bio::EnsEMBL::RawContig;


my($CHR, $START, $END) =  ('20', 30_263_615, 30_275_000);


#
# 1 Test DnaDnaAlignFeature compiles
#
ok(1);



my $multi_db = MultiTestDB->new;
my $db = $multi_db->get_DBAdaptor('core');

my $contig = new Bio::EnsEMBL::RawContig;
$contig->seq('ACTGACTG');
$contig->name('bogus contig');

my @feats;
my $fp = new Bio::EnsEMBL::FeaturePair;
$fp->start(5);
$fp->end  (7);
$fp->strand(1);
$fp->score(10);
$fp->contig($contig);
$fp->hstart(105);
$fp->hend    (107);
$fp->hstrand (1);
$fp->hseqname('dummy-hid');

push(@feats,$fp);


$fp = new Bio::EnsEMBL::FeaturePair;
$fp->start(10);
$fp->end  (14);
$fp->strand(1);
$fp->score(10);
$fp->contig($contig);
$fp->seqname(1);

$fp->hstart  (106);
$fp->hend    (110);
$fp->hstrand (1);
$fp->hseqname('dummy-hid');
push(@feats,$fp);

#
#
# 2 Test DnaDnaAlignFeature::new(-features)
#
$dnaf = Bio::EnsEMBL::DnaDnaAlignFeature->new( -features => \@feats );
ok($dnaf);

#
# 3 Test DnaDnaAlignFeature::seqname
#
ok($dnaf->seqname eq 'bogus contig');

#
# 4 Test DnaDnaAlignFeature::hseqname
#
ok($dnaf->hseqname eq 'dummy-hid');


#
# 5 Test DnaDnaAlignFeature::cigar_string
#
ok($dnaf->cigar_string =~ '3M2I5M');

#
# 6 Test DnaDnaAlignFeature::start
#
ok($dnaf->start == 5);

#
# 7 Test DnaDnaAlignFeature::end
#
ok($dnaf->end == 14);

#
# 8 Test DnaDnaAlignFeature::ungapped_features
#
ok( scalar($dnaf->ungapped_features) == 2);

my $slice = $db->get_SliceAdaptor->fetch_by_chr_start_end($CHR,$START,$END);

#
# 9 Test retrieval from database
#
my $features = $slice->get_all_DnaAlignFeatures;

ok(scalar @$features);

#
# 10 Test transformation to raw contig
#
$features = &transform2rawcontig($features);
ok( scalar @$features );

#
# 11 Test transformation back to slice
#
$features = &transform2slice($features, $slice); 
ok( scalar @$features );

#
# 12 Test transformation onto negative strand slice
#
$features = &transform2slice($features, $slice->invert);
ok(scalar @$features );


########################################################
sub transform2rawcontig {
  my $fs = shift;

  my @out = ();


  foreach my $f (@$fs) {
    my @feats = $f->transform;
    unless(@feats) {
      return [];
    }
    push @out, @feats;
  }
  return \@out;
}

sub transform2slice {
  my ($fs, $slice) = @_;

  my @out = ();

  foreach my $f (@$fs) {
    $f = $f->transform($slice);
    unless($f) {
      return [];
    }
    push @out, $f;
  }
  return \@out;
}
