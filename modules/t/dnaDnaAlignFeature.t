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
use Bio::EnsEMBL::DnaDnaAlignFeature;
use Bio::EnsEMBL::FeaturePair;


my($CHR, $START, $END) =  ('20', 30_363_615, 30_475_000);
my $CTG_BOUNDARY       =  62877;

#
# 1 Test DnaDnaAlignFeature compiles
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
  (-START => 5,
   -END   => 7,
   -STRAND => 1,
   -SCORE => 10,
   -SLICE => $ctg_slice,
   -HSTART => 105,
   -HEND   => 107,
   -HSTRAND => 1,
   -HSEQNAME => 'dummy-hid');

push @feats, Bio::EnsEMBL::FeaturePair->new
  (-start   => 10,
   -end     => 14,
   -strand  => 1,
   -score   => 10,
   -slice   => $ctg_slice,
   -hstart  => 108,
   -hend    => 112,
   -hstrand => 1,
   -hseqname => 'dummy-hid');

#
#
# Test DnaDnaAlignFeature::new(-features)
#
my $dnaf;
warning { $dnaf = Bio::EnsEMBL::DnaDnaAlignFeature->new( -features => \@feats ); };
ok(ref($dnaf) && $dnaf->isa('Bio::EnsEMBL::DnaDnaAlignFeature'));

#
# Test DnaDnaAlignFeature::hseqname
#
ok($dnaf->hseqname eq 'dummy-hid');

#
# 5 Test DnaDnaAlignFeature::cigar_string
#
ok($dnaf->cigar_string eq '3M2I5M');

#
# 6-8 Test DnaDnaAlignFeature::reverse_complement
#
my $strand = $dnaf->strand;
my $hstrand = $dnaf->hstrand;
$dnaf->reverse_complement;
ok($dnaf->cigar_string eq '5M2I3M');
ok(($strand*-1) == $dnaf->strand);
ok(($hstrand*-1) == $dnaf->hstrand);



#
# 9 Test DnaDnaAlignFeature::start
#
ok($dnaf->start == 5);

#
# 10 Test DnaDnaAlignFeature::end
#
ok($dnaf->end == 14);

#
# 11 Test DnaDnaAlignFeature::ungapped_features
#
ok( scalar($dnaf->ungapped_features) == 2);

#
# Test alignment strings
#
my $alignment_strings = $dnaf->alignment_strings('NO_HSEQ');
is($alignment_strings->[0], 'AAATTAAGGG', 'Retrieved query sequence');
is($alignment_strings->[1], '', 'No target sequence');
$alignment_strings = $dnaf->alignment_strings('FIX_SEQ', 'NO_HSEQ');
is($alignment_strings->[0], 'AAATTAAGGG', 'Retrieved fixed query sequence');

#
# Test restrict_between_positions
#
my $new_daf;
warning { $new_daf = $dnaf->restrict_between_positions($dnaf->start + 2, $dnaf->end, 'SEQ'); };
is($new_daf->start, $dnaf->start + 2, 'New daf start');
is($new_daf->end, $dnaf->end, 'End has not changed');
warning { $new_daf = $dnaf->restrict_between_positions(1, $dnaf->end - 2, 'SEQ'); };
is($new_daf->start, $dnaf->start, 'No start change if beyond the start');
is($new_daf->end, $dnaf->end - 2, 'New daf end');


#
# 12 Test retrieval from database
#
my $features = $chr_slice->get_all_DnaAlignFeatures;

ok(scalar @$features);

#
# 13 Test transformation to raw contig
#
my $f = $features->[0];
my @fs = $f->transform('contig');
ok( scalar @fs );

#
# 14 Test transformation back to slice
#
($f) = @fs;
$f = $f->transfer($chr_slice);
ok($f);

#
# 15 Test transformation onto negative strand slice
#
$f = $f->transfer($chr_slice->invert);
ok($f);



done_testing();
