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

use Bio::EnsEMBL::FeaturePair;
use Bio::EnsEMBL::Slice;
use Bio::EnsEMBL::Analysis;
use Bio::EnsEMBL::CoordSystem;
use Bio::EnsEMBL::Test::MultiTestDB;

our $verbose= 0; #turn on or off debug statements

my $coord_system = Bio::EnsEMBL::CoordSystem->new
  (-NAME    => 'chromosome',
   -VERSION => 'NCBI34',
   -DBID    => 123,
   -RANK    => 1);

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
my $hstart = 1;
my $hend   = 90;
my $hstrand = 1;
my $hseqname = 'RF1231';
my $percent_id = 90.8;
my $p_value = '1.52';
my $score   = 50;
my $species = 'Homo_sapiens';
my $hspecies = 'Mus_musculus';

my $fp = Bio::EnsEMBL::FeaturePair->new(-START    => $start,
                                    -END      => $end,
                                    -STRAND   => $strand,
                                    -ANALYSIS => $analysis,
                                    -SLICE    => $slice,
                                    -HSTART   => $hstart,
                                    -HSTRAND  => $hstrand,
                                    -HEND     => $hend,
                                    -HSEQNAME => $hseqname,
                                    -PERCENT_ID => $percent_id,
                                    -P_VALUE  => $p_value,
                                    -SCORE    => $score,
                                    -SPECIES  => $species,
                                    -HSPECIES => $hspecies);


ok($fp && $fp->isa('Bio::EnsEMBL::FeaturePair'));

ok($fp->start == $start);
ok($fp->end == $end);
ok($fp->strand == $strand);
ok($fp->analysis == $analysis);
ok($fp->slice == $slice);

ok($fp->hstart == $hstart);
ok($fp->hend   == $hend);
ok($fp->hseqname eq $hseqname);
ok($fp->percent_id == $percent_id);
ok($fp->p_value == $p_value);
ok($fp->score == $score);
ok($fp->species eq $species);
ok($fp->hspecies eq $hspecies);


#
# Test setters
#
ok(test_getter_setter($fp, 'hstart', 1002));
ok(test_getter_setter($fp, 'hend',   2002));
ok(test_getter_setter($fp, 'hseqname', 'test'));
ok(test_getter_setter($fp, 'percent_id', 123.2));
ok(test_getter_setter($fp, 'p_value', 23.0));
ok(test_getter_setter($fp, 'score', 123));
ok(test_getter_setter($fp, 'species', 'Rattus_norvegicus'));
ok(test_getter_setter($fp, 'hspecies', 'Danio_rerio'));


ok($fp->display_id eq $fp->hseqname());

#
# Test invert
#

$fp->invert();

ok($fp->start   == $hstart);
ok($fp->end     == $hend);
ok($fp->strand  == $hstrand);
ok($fp->hstart  == $start);
ok($fp->hend    == $end);
ok($fp->hstrand == $strand);

ok($fp->seqname eq $hseqname);
ok($fp->hseqname eq $slice->name());
ok(!defined($fp->slice()));

ok($fp->species eq $hspecies);
ok($fp->hspecies eq $species);

done_testing();
