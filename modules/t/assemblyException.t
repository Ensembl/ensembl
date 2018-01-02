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
use Test::Warnings;

use Bio::EnsEMBL::Test::MultiTestDB;

use Bio::EnsEMBL::DnaDnaAlignFeature;
use Bio::EnsEMBL::Test::TestUtils;

our $verbose = 0;

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new();
ok(1);

my $db = $multi->get_DBAdaptor('patch');
my $sfa = $db->get_SimpleFeatureAdaptor();
my $aexc_adaptor = $db->get_AssemblyExceptionFeatureAdaptor();
my $slice_adaptor = $db->get_SliceAdaptor();

my $slice = $slice_adaptor->fetch_by_region('chromosome', 'Y', 1, 400000);
my $feats = $sfa->fetch_all_by_Slice($slice);
is( @$feats, 94, "Returned 94 features" );

print_features($feats);


$slice = $slice_adaptor->fetch_by_region('chromosome', 'Y',
                                            1, 200000);
my $org_slice = $slice_adaptor->fetch_by_region('chromosome', 'X',
                                            1, 200000);


$feats = $sfa->fetch_all_by_Slice($slice);

is( @$feats, 24 , "Fetched 24 features");

print_features($feats);

$multi->hide( "core", "simple_feature" );
$multi->save( "core", "meta_coord" );

for my $f ( @$feats ) {
  $f->dbID( undef );
  $f->adaptor( undef );
  $sfa->store( $f );

  $f->dbID( undef );
  $f->adaptor( undef );
  $f->slice( $org_slice );
  $sfa->store( $f );
    
}


$slice = $slice_adaptor->fetch_by_region('chromosome', 'Y',
					 1, 200000);
$feats = $sfa->fetch_all_by_Slice( $slice );

debug( "After storing retrieval" );
print_features($feats);
is(@$feats, 24, "Fetched 24 features");



#
# sequence getting tests
#

my $hap_slice = $slice_adaptor->fetch_by_region('chromosome', 'Y',
                                             10001, 400000);

$org_slice = $slice_adaptor->fetch_by_region('chromosome', 'X',
                                             60001, 400000);

my ( $fhs, $bhs, $fos, $bos );

debug( "Front hap seq: ".($fhs = $hap_slice->subseq( 99_991, 100_000 )));
debug( "Back hap seq: ".($bhs = $hap_slice->subseq( 400_001, 400_010 )));
debug( "Front org seq: ".( $fos = $org_slice->subseq( 99_991, 100_000 )));
debug( "Back org seq: ".( $bos = $org_slice->subseq( 400_001, 400_010 )));

is( $fhs, $fos, "Front sequences for haplotype and reference are the same" );
is( $bhs, $bos, "Back sequences for haplotype and reference are the same" );


$slice = $slice_adaptor->fetch_by_region('chromosome', 'Y',
					 299_998, 300_002);

is( $slice->seq(), "ACAGA", "Fetched haplotype subslice");

$multi->restore();


#try projecting a hap slice to the contig coordinate system
debug("Org slice projection");
my $projection = $org_slice->project('contig');
print_projection($projection);

debug("Hap slice projection");
$projection = $hap_slice->project('contig');
print_projection($projection);


sub print_features {
  return if(!$verbose);
  my $features = shift;
  foreach my $f (@$features) {
    if(defined($f)) {
      my $seqname = $f->slice->seq_region_name();
      my $analysis = $f->analysis->logic_name();
      debug($seqname . ' ' . $f->start().'-'.$f->end().'('.$f->strand().
            ') ['. $f->dbID.'] '.$f->score() ." ($analysis) ");
    } else {
      debug('UNDEF');
    }
  }
}



sub print_projection {
  my $proj = shift;
  foreach my $seg (@$proj) {
    my ($start, $end, $seq_reg) = ($seg->[0],$seg->[1],$seg->[2]->seq_region_name());
    debug("[$start-$end] $seq_reg");
  }
}

done_testing();
