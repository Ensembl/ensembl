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
use Test::Warnings;

use Bio::EnsEMBL::Test::TestUtils;

use Bio::EnsEMBL::Feature;
use Bio::EnsEMBL::Slice;
use Bio::EnsEMBL::Analysis;
use Bio::EnsEMBL::CoordSystem;
use Bio::EnsEMBL::Test::MultiTestDB;

our $verbose= 0; #turn on or off debug statements


my $multi_db = Bio::EnsEMBL::Test::MultiTestDB->new;
my $db = $multi_db->get_DBAdaptor('core');

my $coord_system = Bio::EnsEMBL::CoordSystem->new
  (-NAME    => 'chromosome',
   -VERSION => 'NCBI34',
   -DBID    => 123,
   -RANK => 1);

my $analysis = Bio::EnsEMBL::Analysis->new(-LOGIC_NAME => 'test');

my $slice = Bio::EnsEMBL::Slice->new(-COORD_SYSTEM    => $coord_system,
                                     -SEQ_REGION_NAME => 'X',
                                     -SEQ_REGION_LENGTH => 15e6,
                                     -START           => 1_000_000,
                                     -END             => 2_000_000);

#
# Test new and getters
#
my $start  = 10;
my $end    = 100;
my $strand = -1;
my $feature = Bio::EnsEMBL::Feature->new(-START => 10,
                                         -END   => 100,
                                         -STRAND => -1,
                                         -ANALYSIS => $analysis,
                                         -SLICE => $slice);


ok($feature && $feature->isa('Bio::EnsEMBL::Feature'));

ok($feature->start == $start);
ok($feature->end == $end);
ok($feature->strand == $strand);
ok($feature->analysis == $analysis);
ok($feature->slice == $slice);

#
# Test setters
#
$analysis = Bio::EnsEMBL::Analysis->new(-LOGIC_NAME => 'new analysis');
$slice = Bio::EnsEMBL::Slice->new(-COORD_SYSTEM    => $coord_system,
                                  -SEQ_REGION_NAME => 'Y',
                                  -SEQ_REGION_LENGTH => 15e6,
                                  -START           => 1_000_000,
                                  -END             => 2_000_000);
ok(&test_getter_setter($feature, 'start', 1000));
ok(&test_getter_setter($feature, 'end', 2000));
ok(&test_getter_setter($feature, 'strand', 1));
ok(&test_getter_setter($feature, 'analysis', $analysis));
ok(&test_getter_setter($feature, 'slice', $slice));


#
# Test length
#

ok($feature->length == ($end-$start+1));


#
# Test move
#

$feature->move(300, 500, 1);

ok($feature->start == 300);
ok($feature->end == 500);
ok($feature->strand == 1);


#################
# Test transform
#################

$slice = $db->get_SliceAdaptor->fetch_by_region('chromosome',
                                                 '20',
                                                 30_249_935,
                                                 31_000_000);

$feature->slice($slice);

#
# Test Transform chromosome -> contig
#

$feature = $feature->transform('contig');

debug("\nchromosome -> contig");
debug("start  = " . $feature->start());
debug("end    = " . $feature->end());
debug("strand = " . $feature->strand());
debug("seq_region = " . $feature->slice->seq_region_name());

ok($feature->start() == 400);
ok($feature->end() == 600);
ok($feature->strand() == 1);
ok($feature->slice->seq_region_name() eq 'AL359765.6.1.13780');
ok($feature->slice->coord_system->name() eq 'contig');


#
# Test Transform contig -> supercontig
#

$feature = $feature->transform('supercontig');

debug("\ncontig -> supercontig");
debug("start  = " . $feature->start());
debug("end    = " . $feature->end());
debug("strand = " . $feature->strand());
debug("seq_region = " . $feature->slice->seq_region_name());

ok($feature->start()  == 658269);
ok($feature->end()    == 658469);
ok($feature->strand() == 1);
ok($feature->slice->seq_region_name() eq 'NT_028392');
ok($feature->slice->coord_system->name() eq 'supercontig');


#
# Test Transform supercontig -> contig
#

$feature = $feature->transform('contig');
debug("\nsupercontig -> contig");
debug("start  = " . $feature->start());
debug("end    = " . $feature->end());
debug("strand = " . $feature->strand());
debug("seq_region = " . $feature->slice->seq_region_name());

ok($feature->start() == 400);
ok($feature->end() == 600);
ok($feature->strand() == 1);
ok($feature->slice->seq_region_name() eq 'AL359765.6.1.13780');
ok($feature->slice->coord_system->name() eq 'contig');




#
# Test Transform contig -> clone
#

$feature = $feature->transform('clone');
debug("\ncontig -> clone");
debug("start  = " . $feature->start());
debug("end    = " . $feature->end());
debug("strand = " . $feature->strand());
debug("seq_region = " . $feature->slice->seq_region_name());

ok($feature->start() == 400);
ok($feature->end() == 600);
ok($feature->strand() == 1);
ok($feature->slice->seq_region_name() eq 'AL359765.6');
ok($feature->slice->coord_system->name() eq 'clone');



#
# Test transform clone -> toplevel
#

$feature = $feature->transform('toplevel');

debug("\nclone -> toplevel");
debug("start  = " . $feature->start());
debug("end    = " . $feature->end());
debug("strand = " . $feature->strand());
debug("seq_region = " . $feature->slice->seq_region_name());


ok($feature->start() == 300 + $slice->start() - 1);
ok($feature->end()   == 500 + $slice->start() - 1);
ok($feature->strand() == 1);
ok($feature->slice->coord_system->name() eq 'chromosome');
ok($feature->slice->seq_region_name() eq '20');

#put back to clone
$feature = $feature->transform('clone');


#
# Test transform to into gap
#
$slice = $db->get_SliceAdaptor->fetch_by_region('chromosome',
                                                 '20',
                                                 1_000_000,
                                                 2_000_000);
$feature->slice($slice);

ok(!defined($feature->transform('contig')));



###############
# Test transfer
###############
$slice = $db->get_SliceAdaptor->fetch_by_region('chromosome',
                                                 '20',
                                                 30_249_935,
                                                 31_000_000);

$feature->slice($slice);


#
# Transfer to expanded inverted chr slice
#

#get larger slice on other strand
$slice = $slice->invert()->expand(100,100);

$feature = $feature->transfer($slice);
debug("\ntransfer to inverted chromosomal slice");
debug("start  = " . $feature->start());
debug("end    = " . $feature->end());
debug("strand = " . $feature->strand());
debug("seq_region = " . $feature->slice->seq_region_name());

ok($feature->start()  == 749567);
ok($feature->end()    == 749767);
ok($feature->strand() == -1);


#
# Transfer to contig slice
#

$slice = $db->get_SliceAdaptor->fetch_by_region('contig',
                                                'AL359765.6.1.13780',
                                                '30', '3000');

$feature = $feature->transfer($slice);
debug("\ntransfer to contig slice");
debug("start  = " . $feature->start());
debug("end    = " . $feature->end());
debug("strand = " . $feature->strand());
debug("seq_region = " . $feature->slice->seq_region_name());


ok($feature->start()  == 471);
ok($feature->end()    == 671);
ok($feature->strand() == 1);


#
# Transfer to supercontig slice
#
$slice  = $db->get_SliceAdaptor->fetch_by_region('supercontig',
                                                 'NT_028392');
debug( "----------------------------" );

$feature = $feature->transfer($slice);
debug("\ntransfer to supercontig slice");
debug("start  = " . $feature->start());
debug("end    = " . $feature->end());
debug("strand = " . $feature->strand());
debug("seq_region = " . $feature->slice->seq_region_name());

ok($feature->start()  == 658369);
ok($feature->end()    == 658569);
ok($feature->strand() == 1);


#
# Project over a contig boundary
#

$feature->move(671600,671800);

debug( "Before projection to contig" );
debug("start  = " . $feature->start());
debug("end    = " . $feature->end());
debug("strand = " . $feature->strand());
debug("seq_region = " . $feature->slice->seq_region_name());

my @projection = @{$feature->project('contig')};

debug( "After project to contig" );
foreach my $segment (@projection) {
  ($start, $end, $slice) = @$segment;
  debug("[$start-$end] -> " . $slice->seq_region_name
        . ' ' . $slice->start . '-' . $slice->end . '('.$slice->strand().')');
}
debug( "-----------------------" );

ok(@projection == 2);

my ($seg1, $seg2) = @projection;

ok($seg1 && $seg1->[0] == 1);
ok($seg1 && $seg1->[1]   == 50);
ok($seg1 && $seg1->[2]->start == 13731);
ok($seg1 && $seg1->[2]->end   == 13780);
ok($seg1 && $seg1->[2]->seq_region_name eq 'AL359765.6.1.13780');

ok($seg2 && $seg2->[0] == 51);
ok($seg2 && $seg2->[1] == 201);
ok($seg2 && $seg2->[2]->start == 101);
ok($seg2 && $seg2->[2]->end   == 251);
ok($seg2 && $seg2->[2]->seq_region_name eq 'AL031658.11.1.162976');

debug('forward strand slice feature projection');
foreach my $segment (@projection) {
  ($start, $end, $slice) = @$segment;
  debug("[$start-$end] -> " . $slice->seq_region_name
        . ' ' . $slice->start . '-' . $slice->end . '('.$slice->strand().')');
}

#
# Try a projection using an inverted slice
#

$slice = $feature->slice()->invert();
$feature->slice($slice);
$feature->move(4_391_950, 4_392_150);

@projection = @{$feature->project('contig')};

ok(@projection == 2);

($seg1, $seg2) = @projection;

ok($seg1 && $seg1->[0] == 1 );
ok($seg1 && $seg1->[1]   ==  8 );
ok($seg1 && $seg1->[2]->start == 101);
ok($seg1 && $seg1->[2]->end   == 108);
ok($seg1 && $seg1->[2]->seq_region_name eq 'AL031658.11.1.162976');

ok($seg2 && $seg2->[0] == 9);
ok($seg2 && $seg2->[1] == 201);
ok($seg2 && $seg2->[2]->start == 13588);
ok($seg2 && $seg2->[2]->end   == 13780);
ok($seg2 && $seg2->[2]->seq_region_name eq 'AL359765.6.1.13780');

debug('negative strand slice feature projection');
foreach my $segment (@projection) {
  ($start, $end, $slice) = @$segment;
  debug("[$start-$end] -> " . $slice->seq_region_name
        . ' ' . $slice->start . '-' . $slice->end . '('.$slice->strand().')');
}

#
# Try project to same coord system
#


@projection = @{$feature->project('supercontig')};

ok(@projection == 1);

($seg1) = @projection;
$slice = $feature->slice();
ok($seg1 && $seg1->[0] == 1 );
ok($seg1 && $seg1->[1]   == $feature->length() );
ok($seg1 && $seg1->[2]->start == $slice->end - $feature->end() + 1);
ok($seg1 && $seg1->[2]->end   == $slice->end - $feature->start() + 1);
ok($seg1 && $seg1->[2]->strand == $slice->strand * $feature->strand());
ok($seg1 && $seg1->[2]->seq_region_name eq 'NT_028392');


debug('projection to same coord system');
foreach my $segment (@projection) {
  ($start, $end, $slice) = @$segment;
  debug("[$start-$end] -> " . $slice->seq_region_name
        . ' ' . $slice->start . '-' . $slice->end . '('.$slice->strand().')');
}

ok($feature->display_id eq '');


#
# test seq function call
#

$slice = $db->get_SliceAdaptor->fetch_by_region('chromosome',
                                                 '20',
                                                 30_249_935,
                                                 31_000_000);
$feature = new Bio::EnsEMBL::Feature( -slice => $slice,
				    -start => 100,
				    -end => 200,
				    -strand => 1,
				    -analysis => $analysis );

my $seq = $feature->seq();
debug( "Feature ".$feature->start()."-".$feature->end() );
debug( "Feature sequence is: ".$seq );

ok( length( $seq ) == 101 );
ok( $seq =~ /^[GCATgcat]+$/ );

$feature->move( 100,200,-1 );
my $seq2 = $feature->seq();
debug( "Reverse Feature seq: ".$seq2 );

$seq2 = reverse( $seq2 );
$seq2 =~ tr/acgtACGT/tgcaTGCA/;

ok( $seq eq $seq2 );

#
#


my $slice_adaptor = $db->get_SliceAdaptor();
$slice = $slice_adaptor->fetch_by_region('chromosome', '20', 1e6, 10e6);

$feature->slice($slice);

###
# test seq_region functions with forward strand slice
#
ok($feature->seq_region_start   == $feature->start + $slice->start - 1);
ok($feature->seq_region_end     == $feature->end   + $slice->start - 1);
ok($feature->seq_region_strand  == $feature->strand);
ok($feature->seq_region_name    eq $slice->seq_region_name());


#####
# test seq_region functions with reverse strand slice
#
$slice = $slice->invert();
$feature->slice($slice);
ok($feature->seq_region_start   == $slice->end() - $feature->end() + 1);
ok($feature->seq_region_end     == $slice->end() - $feature->start() + 1);
ok($feature->seq_region_strand  == $feature->strand() * -1);
ok($feature->seq_region_name    eq $slice->seq_region_name());

##### 
# test seq_region functions with no slice
#
$feature->slice(undef);
ok(!defined($feature->seq_region_start));
ok(!defined($feature->seq_region_end));
ok(!defined($feature->seq_region_strand));
ok(!defined($feature->seq_region_name));


#
# test overlaps function
#

my $chr_slice = $db->get_SliceAdaptor->fetch_by_region('chromosome',
						       '20',
						       30_249_935,
						       31_000_000);
my $all_chr_slice =  $db->get_SliceAdaptor->fetch_by_region('chromosome', 20);

my $hap_chr_slice = $db->get_SliceAdaptor->fetch_by_region('chromosome',
                   '20_HAP1',
                   30_249_935,
                   31_000_000);

# my $sctg_slice  = $db->get_SliceAdaptor->fetch_by_region('supercontig',
# 							 'NT_028392');

# my $ctg_slice = $db->get_SliceAdaptor->fetch_by_region('contig',
# 						       'AL359765.6.1.13780',
# 						       '30', '3000');


my $f1 = new Bio::EnsEMBL::Feature( -start => 1,
				    -end => 10,
				    -strand => 1,
				    -slice => $chr_slice,
				    -analysis => $analysis
				  );

my $f2 = new Bio::EnsEMBL::Feature( -start => 10,
				    -end => 20,
				    -strand => -1,
				    -slice => $chr_slice,
				    -analysis => $analysis
				  );

#
# simple same coord system overlap
#

ok( $f1->overlaps( $f2 ), 'Overlaps: feature 1 <-> feature 2 overlap (1 ends @ 10. 2 starts @ 10)');
ok( $f1->overlaps_local( $f2 ), 'Overlaps Local: feature 1 <-> feature 2 overlap (1 ends @ 10. 2 starts @ 10)');

$f2 = new Bio::EnsEMBL::Feature( -start => 11,
				 -end => 20,
				 -strand => -1,
				 -slice => $chr_slice,
				 -analysis => $analysis
			       );

ok( ! $f1->overlaps( $f2 ), 'Overlaps: feature 1 <-> new feature 2 do not overlap (1 ends @ 10. 2 starts @ 11)');
ok( ! $f1->overlaps_local( $f2 ), 'Overlaps Local: feature 1 <-> feature 2 overlap (1 ends @ 10. 2 starts @ 10)');

#
# Testing when seq_region overlaps do not return true but local ones do
#
my $f1_absolute = Bio::EnsEMBL::Feature->new(-start => 30_249_935, -end => 30_249_945, -strand => 1, -slice => $all_chr_slice);
ok($f1->overlaps($f1_absolute), 'Overlaps: feature 1 overlaps as both it and absolute feature are on the same chromsomal coordinate');
ok(!$f1->overlaps_local($f1_absolute), 'Overlaps Local: feature 1 does not overlap absolute feature; their starts and ends are different');

# 
# Testing alternative seq region overlap
#
my $hap_feature_1 = Bio::EnsEMBL::Feature->new(-start => 1, -end => 10, -strand => 1, -slice => $hap_chr_slice);
warns_like 
  { ok( ! $hap_feature_1->overlaps( $f1 ), 'Overlaps: feature 1 <-> new feature 2 do not overlap (different seq regions)') }
  qr/different seq regions/,
  'Error when using different seq regions as expected';

warns_like 
  { ok( ! $hap_feature_1->overlaps_local( $f1 ), 'Overlaps Local: feature 1 <-> feature 2 overlap (different seq regions)') }
  qr/different seq regions/,
  'Error when using different seq regions as expected';

#
# test get_all_alt_locations
#
$chr_slice = $slice_adaptor->fetch_by_region('chromosome', 'Y', 10000000,
                                             10500000);

$f1 = new Bio::EnsEMBL::Feature(-start => 10, -end => 20, -strand => 1,
                                -slice => $chr_slice);

my @alt = @{$f1->get_all_alt_locations()};

ok(@alt == 1);
ok($alt[0]->slice->seq_region_name eq '20' &&
   $alt[0]->start == 30300009 && $alt[0]->end == 30300019);

print_locations($f1);


print_locations($alt[0]);

@alt = @{$alt[0]->get_all_alt_locations()};
debug( "There are ".scalar( @alt )." alternative locations" );

ok($alt[0]->slice->seq_region_name eq 'Y' &&
   $alt[0]->start == 10000009 && $alt[0]->end == 10000019);


$chr_slice = $db->get_SliceAdaptor->fetch_by_region('chromosome', 
                                                    '20_HAP1');

$f1 = new Bio::EnsEMBL::Feature(-start => 10, -end => 20, -strand => 1,
                                -slice => $chr_slice);

@alt = @{$f1->get_all_alt_locations()};
ok( scalar( @alt ) == 1 );

debug( "Alt locations from a HAP thing ".scalar( @alt ) );

$chr_slice = $db->get_SliceAdaptor->fetch_by_region('chromosome', 
                                                    '20');

$f1 = new Bio::EnsEMBL::Feature(-start => 10, -end => 20, -strand => 1,
                                -slice => $chr_slice);

@alt = @{$f1->get_all_alt_locations()};

debug( "Alt locations from non HAP thing ".scalar( @alt ) );

ok( scalar( @alt ) == 0 );

#
# other coord system overlaps
#

# $f2 = new Bio::EnsEMBL::Feature( -start => 1,
# 				 -end => 1000_000,
# 				 -strand => -1,
# 				 -slice => $sctg_slice,
# 				 -analysis => $analysis
# 			       );

# debug( "Other coord overlaps ".( $f1->overlaps( $f2 )));
# ok( $f1->overlaps( $f2 ));


# #
# # Just seqname features
# #


# $f1->seqname( "minigenewise" );
# $f1->slice( undef );

# $f2 = new Bio::EnsEMBL::Feature( -start => 10,
# 				 -end => 20,
# 				 -strand => -1,
# 				 -seqname => "minigenewise",
# 				 -analysis => $analysis
# 			       );


# debug( "seqname f1<->f2 overlap = ".($f1->overlaps( $f2 )));
# ok( ($f1->overlaps( $f2 )));


# get_overlapping_Genes

my $gene_adaptor = $db->get_GeneAdaptor;
my $gene = $gene_adaptor->fetch_by_stable_id('ENSG00000171456');
my $o_genes = $gene->get_overlapping_Genes;
is($o_genes->[0]->stable_id,'ENSG00000171456','Check overlap with a feature that overlaps nothing');
$chr_slice = $o_genes->[0]->feature_Slice;

my $f3 = new Bio::EnsEMBL::Feature( 
  -start => 1,
  -end => 100,
  -strand => 1,
  -slice => $chr_slice,
  -analysis => $analysis,
  -adaptor => $db->get_GeneAdaptor,
); 

$o_genes = $f3->get_overlapping_Genes;
debug('Got some overlapping things:'.scalar @$o_genes);
foreach my $g (@$o_genes) {
  debug('Overlapping: '.$g->stable_id);
}
is($o_genes->[0]->stable_id,'ENSG00000171456','Check overlap with a fake feature on top of an existing gene');
# Note, doesn't find itself as it isn't in the test database.
sub print_locations {
  my $f = shift;

  debug("* ". $f->slice->seq_region_name(). " ".
        $f->seq_region_start. "-". $f->seq_region_end(). "\n");

  foreach my $af (@{$f->get_all_alt_locations()}) {
    debug("  ". $af->slice->seq_region_name(). " ".
        $af->seq_region_start. "-". $af->seq_region_end(). "\n");
  }
}

done_testing();
