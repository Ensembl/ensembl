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
use Test::More;
use Test::Warnings qw(warning);
use Bio::EnsEMBL::Test::MultiTestDB;
use Bio::EnsEMBL::DnaPepAlignFeature;

use Bio::EnsEMBL::Test::TestUtils;

# switch on the debug prints

our $verbose = 0;

my($CHR, $START, $END) =  ('20', 30_363_615, 30_475_000);
my $CTG_BOUNDARY       =  62877;

#
# 1 Test DnaPepAlignFeature compiles
#
ok(1);

my $multi_db = Bio::EnsEMBL::Test::MultiTestDB->new;
my $db = $multi_db->get_DBAdaptor('core');



my $chr_slice = $db->get_SliceAdaptor->fetch_by_region('chromosome',
                                                       $CHR,$START,$END);

my $ctg_slice = $db->get_SliceAdaptor->fetch_by_region('contig',
                                                       'AL359765.6.1.13780');

my @feats;


push @feats, Bio::EnsEMBL::FeaturePair->new
  (-start  => 5,
   -end    => 7,
   -score  => 10,
   -strand => 1,
   -slice  => $ctg_slice,
   -hstart => 105,
   -hend   => 105,
   -hstrand => 1,
   -hseqname => 'dummy-hid');

push @feats, Bio::EnsEMBL::FeaturePair->new
  (-start  => 11,
   -end    => 16,
   -score  => 10,
   -strand => 1,
   -slice  => $ctg_slice,
   -hstart => 106,
   -hend    => 107,
   -hstrand => 1,
   -hseqname => 'dummy-hid');

#
#
# Test DnaPepAlignFeature::new(-features)
#
my $dnaf;
warning { $dnaf = Bio::EnsEMBL::DnaPepAlignFeature->new( -features => \@feats ); };
ok(ref($dnaf) && $dnaf->isa('Bio::EnsEMBL::DnaPepAlignFeature'));

#
# Test DnaPepAlignFeature::hseqname
#
ok($dnaf->hseqname eq 'dummy-hid');


#
# Test DnaPepAlignFeature::cigar_string
#
ok($dnaf->cigar_string =~ '3M3I6M');

#
# Test DnaPepAlignFeature::reverse_complement
#
my $strand = $dnaf->strand;
my $hstrand = $dnaf->hstrand;
$dnaf->reverse_complement;
ok($dnaf->cigar_string =~ '6M3I3M');
ok(($strand*-1) == $dnaf->strand);
ok(($hstrand*-1) == $dnaf->hstrand); 



#
# Test DnaPepAlignFeature::start
#
ok($dnaf->start == 5);

#
# Test DnaPepAlignFeature::end
#
ok($dnaf->end == 16);

#
# Test DnaPepAlignFeature::ungapped_features
#
ok( scalar($dnaf->ungapped_features) == 2);


#
# 12 Test retrieval from database
#
my $features = $chr_slice->get_all_ProteinAlignFeatures;

ok(scalar @$features);

#
# Test transformation to raw contig
#
my $f = $features->[0];
my @fs = $f->transform('contig');
ok( scalar @fs );

#
# Test transformation back to slice
#
($f) = @fs;
$f = $f->transfer($chr_slice);
ok($f);

#
#  Test transformation onto negative strand slice
#
$f = $f->transfer($chr_slice->invert);
ok($f);


##
## 16-21 create a dnaalign feature on a slice across a contig boundary
##       and convert to raw contig coordinates
##       (+ve strand, +ve hitstrand)
##
#@feats = ();
#$fp = new Bio::EnsEMBL::FeaturePair;
#$fp->start($CTG_BOUNDARY - 5);
#$fp->end  ($CTG_BOUNDARY );
#$fp->strand(1);
#$fp->score(10);
#$fp->contig($slice);
#$fp->hstart(104);
#$fp->hend  (105);
#$fp->hstrand (1);
#$fp->hseqname('dummy-hid');
#push(@feats,$fp);

#$fp = new Bio::EnsEMBL::FeaturePair;
#$fp->start($CTG_BOUNDARY + 4);
#$fp->end  ($CTG_BOUNDARY + 9);
#$fp->strand(1);
#$fp->score(10);
#$fp->contig($slice);
#$fp->hstart  (106);
#$fp->hend    (107);
#$fp->hstrand (1);
#$fp->hseqname('dummy-hid');
#push(@feats,$fp);

#$fp = new Bio::EnsEMBL::FeaturePair;
#$fp->start($CTG_BOUNDARY + 10);
#$fp->end  ($CTG_BOUNDARY + 12);
#$fp->strand(1);
#$fp->score(10);
#$fp->contig($slice);
#$fp->hstart  (110);
#$fp->hend    (110);
#$fp->hstrand (1);
#$fp->hseqname('dummy-hid');
#push(@feats,$fp);

#$dnaf = Bio::EnsEMBL::DnaPepAlignFeature->new( -features => \@feats );
#ok($dnaf);
#ok($dnaf->cigar_string eq '6M3I6M6D3M');
#ok($dnaf->validate || 1); #validate doesn't return true but throws on fail

#@dnafs = $dnaf->transform;
#ok(scalar(@dnafs) == 2);
#ok($dnafs[0]->validate || 1); 
#ok($dnafs[1]->validate || 1);




##
## 22-27 create a dnaalign feature on a slice across a contig boundary
##       and convert to raw contig coordinates
##       (-ve strand, +ve hitstrand)
##
#@feats = ();

#$fp = new Bio::EnsEMBL::FeaturePair;
#$fp->start($CTG_BOUNDARY + 8);
#$fp->end  ($CTG_BOUNDARY + 10);
#$fp->strand(-1);
#$fp->score(10);
#$fp->contig($slice);
#$fp->hstart  (100);
#$fp->hend    (100);
#$fp->hstrand (1);
#$fp->hseqname('dummy-hid');
#push(@feats,$fp);

#$fp = new Bio::EnsEMBL::FeaturePair;
#$fp->start($CTG_BOUNDARY - 1);
#$fp->end  ($CTG_BOUNDARY + 4);
#$fp->strand(-1);
#$fp->score(10);
#$fp->contig($slice);
#$fp->hstart(101);
#$fp->hend    (102);
#$fp->hstrand (1);
#$fp->hseqname('dummy-hid');
#push(@feats,$fp);

#$fp = new Bio::EnsEMBL::FeaturePair;
#$fp->start($CTG_BOUNDARY - 4);
#$fp->end  ($CTG_BOUNDARY - 2);
#$fp->strand(-1);
#$fp->score(10);
#$fp->contig($slice);
#$fp->seqname(1);
#$fp->hstart  (105);
#$fp->hend    (105);
#$fp->hstrand (1);
#$fp->hseqname('dummy-hid');
#push(@feats,$fp);


#$dnaf = Bio::EnsEMBL::DnaPepAlignFeature->new( -features => \@feats );
#ok($dnaf);
#ok($dnaf->cigar_string eq '3M3I6M6D3M');
#ok($dnaf->validate || 1); #validate doesn't return true but throws on fail

#@dnafs = $dnaf->transform;
#ok(scalar(@dnafs) == 2);

#debug( "Feature 0 dump" );
#while( my ($k, $v) = each %{$dnafs[0]} ) {
#  debug( "  ->".$k." = ".$v );
#}

#ok($dnafs[0]->validate || 1); 

#debug( "Feature 1 dump" );
#while( my ($k, $v) = each %{$dnafs[1]} ) {
#  debug( "  ->".$k." = ".$v );
#}
#ok($dnafs[1]->validate || 1);




##
##
## Do the same tests again on the negative strand slice
##
##
#$CTG_BOUNDARY = $slice->length - $CTG_BOUNDARY + 1;
#$slice = $slice->invert;



##
## 28-33 create a dnaalign feature on a slice across a contig boundary
##       and convert to raw contig coordinates
##       (+ve strand, +ve hitstrand)
##
#@feats = ();
#$fp = new Bio::EnsEMBL::FeaturePair;
#$fp->start($CTG_BOUNDARY - 2);
#$fp->end  ($CTG_BOUNDARY);
#$fp->strand(1);
#$fp->score(10);
#$fp->contig($slice);
#$fp->hstart(105);
#$fp->hend  (105);
#$fp->hstrand (1);
#$fp->hseqname('dummy-hid');
#push(@feats,$fp);

#$fp = new Bio::EnsEMBL::FeaturePair;
#$fp->start($CTG_BOUNDARY + 4);
#$fp->end  ($CTG_BOUNDARY + 9);
#$fp->strand(1);
#$fp->score(10);
#$fp->contig($slice);
#$fp->hstart  (106);
#$fp->hend    (107);
#$fp->hstrand (1);
#$fp->hseqname('dummy-hid');
#push(@feats,$fp);

#$fp = new Bio::EnsEMBL::FeaturePair;
#$fp->start($CTG_BOUNDARY + 10);
#$fp->end  ($CTG_BOUNDARY + 12);
#$fp->strand(1);
#$fp->score(10);
#$fp->contig($slice);
#$fp->hstart  (110);
#$fp->hend    (110);
#$fp->hstrand (1);
#$fp->hseqname('dummy-hid');
#push(@feats,$fp);

#$dnaf = Bio::EnsEMBL::DnaPepAlignFeature->new( -features => \@feats );
#ok($dnaf);
#ok($dnaf->cigar_string eq '3M3I6M6D3M');
#ok($dnaf->validate || 1); #validate doesn't return true but throws on fail

#@dnafs = $dnaf->transform;
#ok(scalar(@dnafs) == 2);
#ok($dnafs[0]->validate || 1); 
#ok($dnafs[1]->validate || 1);



##
## 34-39 create a dnaalign feature on a slice across a contig boundary
##       and convert to raw contig coordinates
##       (-ve strand, +ve hitstrand)
##
#@feats = ();

#$fp = new Bio::EnsEMBL::FeaturePair;
#$fp->start($CTG_BOUNDARY + 8);
#$fp->end  ($CTG_BOUNDARY + 10);
#$fp->strand(-1);
#$fp->score(10);
#$fp->contig($slice);
#$fp->hstart  (100);
#$fp->hend    (100);
#$fp->hstrand (1);
#$fp->hseqname('dummy-hid');
#push(@feats,$fp);

#$fp = new Bio::EnsEMBL::FeaturePair;
#$fp->start($CTG_BOUNDARY - 1);
#$fp->end  ($CTG_BOUNDARY + 4);
#$fp->strand(-1);
#$fp->score(10);
#$fp->contig($slice);
#$fp->hstart(101);
#$fp->hend    (102);
#$fp->hstrand (1);
#$fp->hseqname('dummy-hid');
#push(@feats,$fp);

#$fp = new Bio::EnsEMBL::FeaturePair;
#$fp->start($CTG_BOUNDARY - 4);
#$fp->end  ($CTG_BOUNDARY - 2);
#$fp->strand(-1);
#$fp->score(10);
#$fp->contig($slice);
#$fp->seqname(1);
#$fp->hstart  (105);
#$fp->hend    (105);
#$fp->hstrand (1);
#$fp->hseqname('dummy-hid');
#push(@feats,$fp);


#$dnaf = Bio::EnsEMBL::DnaPepAlignFeature->new( -features => \@feats );
#ok($dnaf);
#ok($dnaf->cigar_string eq '3M3I6M6D3M');
#ok($dnaf->validate || 1); #validate doesn't return true but throws on fail

#@dnafs = $dnaf->transform;
#ok(scalar(@dnafs) == 2);
#ok($dnafs[0]->validate || 1); 
#ok($dnafs[1]->validate || 1);

done_testing();

