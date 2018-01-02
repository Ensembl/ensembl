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

use Bio::EnsEMBL::ProteinFeature;
use Bio::EnsEMBL::Analysis;

our $verbose = 0; #turn on or off debug statements


#
# Test new and getters
#
my $start  = 10;
my $end    = 100;
my $hstart = 1;
my $hend   = 90;
my $hstrand = 1;
my $hseqname = 'RF1231';
my $percent_id = 90.8;
my $p_value = '1.52';
my $score   = 50;
my $species = 'Homo_sapiens';
my $hspecies = 'Mus_musculus';
my $hdes = "Hit description";

my $idesc = 'interpro description';
my $ilabel = 'interpro label';
my $interpro_ac = 'interpro accession';
my $translation_id = 1234;

my $analysis_db = 'test_db';

my $analysis = Bio::EnsEMBL::Analysis->new(-LOGIC_NAME => 'test', -DB => $analysis_db);


my $f = Bio::EnsEMBL::ProteinFeature->new
  (-START       => $start,
   -END         => $end,
   -ANALYSIS    => $analysis,
   -HSTART      => $hstart,
   -HEND        => $hend,
   -HSEQNAME    => $hseqname,
   -PERCENT_ID  => $percent_id,
   -P_VALUE     => $p_value,
   -SCORE       => $score,
   -SPECIES     => $species,
   -HSPECIES    => $hspecies,
   -HDESCRIPTION=> $hdes,
   -IDESC       => $idesc,
   -ILABEL      => $ilabel,
   -INTERPRO_AC => $interpro_ac,
   -TRANSLATION_ID     => $translation_id);



ok($f && $f->isa('Bio::EnsEMBL::ProteinFeature'));

ok($f->start == $start);
ok($f->end == $end);
ok($f->analysis == $analysis);
ok($f->translation_id == $translation_id);

ok($f->hstart == $hstart);
ok($f->hend   == $hend);
ok($f->hseqname eq $hseqname);
ok($f->percent_id == $percent_id);
ok($f->p_value == $p_value);
ok($f->score == $score);
ok($f->species eq $species);
ok($f->hspecies eq $hspecies);
ok($f->hdescription eq $hdes);

ok($f->idesc eq $idesc);
ok($f->ilabel eq $ilabel);
ok($f->interpro_ac eq $interpro_ac);

# check that the strand is 0
ok($f->strand == 0);

# Test summary_as_hash
my $summary = $f->summary_as_hash();
ok($summary->{'type'} = $analysis_db);
ok($summary->{'id'} = $hseqname);
ok($summary->{'start'} = $start);
ok($summary->{'end'} = $end);
ok($summary->{'interpro'} = $interpro_ac);
ok($summary->{'description'} = $idesc);
ok($summary->{'hit_start'} = $hstart);
ok($summary->{'hit_end'} = $hend);

#
# Test setters
#
ok(test_getter_setter($f, 'idesc', 'interpro desc1'));
ok(test_getter_setter($f, 'ilabel', 'interpro label1'));
ok(test_getter_setter($f, 'interpro_ac',   'interpro ac1'));


done_testing();

