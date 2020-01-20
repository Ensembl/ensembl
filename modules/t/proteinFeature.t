# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2020] EMBL-European Bioinformatics Institute
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
use Bio::EnsEMBL::Test::MultiTestDB;
use Bio::EnsEMBL::BaseAlignFeature;
use Bio::EnsEMBL::FeaturePair;
use Bio::EnsEMBL::ProteinFeature;
use Bio::EnsEMBL::Analysis;
use Data::Dumper;
use Test::Exception;

our $verbose = 0; #turn on or off debug statements

my $multi_db = Bio::EnsEMBL::Test::MultiTestDB->new;
my $db = $multi_db->get_DBAdaptor('core');
my $pfa = $db->get_ProteinFeatureAdaptor();



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


#
# Test alignment length and alignment string

$analysis = Bio::EnsEMBL::Analysis->new(-LOGIC_NAME => 'gifts_import', -DB => $analysis_db);

my $seq_start  = 1;
my $seq_end    = 111;
my $hit_start = 0;
my $hit_end   = 0;
my $hit_name = "P20366";
my $cigar_string = "MD:Z:0G105";
my $align_type = "mdtag";
my $transl_id = 21739;

#Create ProteinFeature
my $pf = Bio::EnsEMBL::ProteinFeature->new
  (-START       => $seq_start,
   -END         => $seq_end,
   -ANALYSIS    => $analysis,
   -HSTART      => $hit_start,
   -HEND        => $hit_end,
   -HSEQNAME    => $hit_name,
   -CIGAR_STRING => $cigar_string,
   -ALIGN_TYPE => $align_type,
   -TRANSLATION_ID     => $transl_id,
   -ADAPTOR           => $pfa
   );

ok($pf && $pf->isa('Bio::EnsEMBL::ProteinFeature'));

# Test summary_as_hash
my $pf_summary = $pf->summary_as_hash();

ok($pf_summary->{'type'} eq $analysis_db);
ok($pf_summary->{'id'} eq $hit_name);
ok($pf_summary->{'start'} eq $seq_start);
ok($pf_summary->{'end'} = $seq_end);
ok($pf_summary->{'interpro'} eq '');
ok($pf_summary->{'description'} eq '');
ok($pf_summary->{'hit_start'} eq $hit_start);
ok($pf_summary->{'hit_end'} eq $hit_end);
ok($pf_summary->{'cigar_string'} eq $cigar_string);
ok($pf_summary->{'align_type'} eq $align_type);
ok($pf_summary->{'hseqname'} eq $hit_name);
ok($pf_summary->{'translation_id'} eq $transl_id);


# Test alignment strings
my $alignment_string = $pf->alignment_strings();
is($alignment_string->[0], "DNSSLSGEERLKCKLGKSFLLEKSLGKGMLIHCSLGVSMGKGKPPSPLTLTSFPPFCDLAKSAFHVVLTTTGVKLTMIPYSRSRLMSSEDLAEIPQLQKLSIPHGF", "Got query string right");
is($alignment_string->[1], "GNSSLSGEERLKCKLGKSFLLEKSLGKGMLIHCSLGVSMGKGKPPSPLTLTSFPPFCDLAKSAFHVVLTTTGVKLTMIPYSRSRLMSSEDLAEIPQLQKLSIPHGF", "Got target string right");


my $mdz_alignment_length_test = {
  "MD:Z:14" => 14,
  "MD:Z:323C10" => 334,
  "MD:Z:3C10" => 14,
  "MD:Z:0C12T0" => 14,
  "MD:Z:96^RHKTDSFVGLMGKRALNS0V14"=> 111,
  "MD:Z:50^EIGVLAKAFIDQGKLIPDDVMTRLALHELKNLTQYSWLLD137" => 187,
  "MD:Z:422^PLNGSGQLKMPSHCLSAQMLAPPPPGLPRLAL2A0T0K1A0^TTSEGGATSPTSPSY0S1P0D1S1A0^NRSFVGLGPRDPAGIYQAQS1Y0^LG" => 439,
  "MD:Z:108^GITRKERPPLDVDEMLERFKTEAQ930" => 1038,
  "MD:Z:0^MP0P0T1D1F0Q0Q0P0T1D0N0^DDSYLGELRA0S0K142" => 157,
  "MD:Z:212S43" => 256,
  "MD:Z:35S0R0Y0I0V0M0G0H0I1^HKRRQLPTA0L1Q0V0L0R0G0R0L0R0P0G0D0G0L0L0R0S0S0S0S0Y0V0K1F0N0R0K0R0^EGQIQGAIHTQCI" => 75,
  "MD:Z:0G105" => 106,
  "MD:Z:35^VIVALE31^GRPLIQPRRKKAYQLEHTFQGLLGKRSLFTE10 " => 76
  
};

while(my ($mdz_string, $length) = each %$mdz_alignment_length_test){
  $pf->{'_alignment_length'} = undef;
  $pf->cigar_string($mdz_string);
  is($pf->align_type, "mdtag", "Got the right align_type");
  is($pf->cigar_string, $mdz_string);
  is( $pf->alignment_length(), $length, "Got back correct alignment length $length for mdz string ". $pf->cigar_string);
}

my $mdz_alignment_test = {
  "MD:Z:96^RHKTDSFVGLMGKRALNS0V14" => # mdz string ENSP00000289576 P20366
    [ 
      "MKILVALAVFFLVSTQLFAEEIGANDDLNYWSDWYDSDQIKEELPEPFEHLLQRIARRPKPQQFFGLMGKRDADSSIEKQVALLKALYGHGQISHKMAYERSAMQNYERRR",  # input seq
      "MKILVALAVFFLVSTQLFAEEIGANDDLNYWSDWYDSDQIKEELPEPFEHLLQRIARRPKPQQFFGLMGKRDADSSIEKQVALLKALYGHGQISHK------------------MAYERSAMQNYERRR", # target seq
      "MKILVALAVFFLVSTQLFAEEIGANDDLNYWSDWYDSDQIKEELPEPFEHLLQRIARRPKPQQFFGLMGKRDADSSIEKQVALLKALYGHGQISHKRHKTDSFVGLMGKRALNSVAYERSAMQNYERRR" # query seq
    ],
  "MD:Z:35S0R0Y0I0V0M0G0H0I1^HKRRQLPTA0L1Q0V0L0R0G0R0L0R0P0G0D0G0L0L0R0S0S0S0S0Y0V0K1F0N0R0K0R0^EGQIQGAIHTQCI" =>     #ENSP00000334741  Q2M2W7
    [
      "MNRLYLTPDGFFFRVHMLALDSSSCNKPCPEFKPGIETDLNDAAYVLYTTVCNVGATARAVGRPAFFWERWETMT",
      "MNRLYLTPDGFFFRVHMLALDSSSCNKPCPEFKPGIETDLNDAAY---------VLYTTVCNVGATARAVGRPAFFWERWETMT-------------",
      "MNRLYLTPDGFFFRVHMLALDSSSCNKPCPEFKPGSRYIVMGHIYHKRRQLPTALLQVLRGRLRPGDGLLRSSSSYVKRFNRKREGQIQGAIHTQCI"
    ],
  "MD:Z:0G105" => #ENSP00000374857 P0DOY3
    [
      "XQPKAAPSVTLFPPSSEELQANKATLVCLISDFYPGAVTVAWKADSSPVKAGVETTTPSKQSNNKYAASSYLSLTPEQWKSHKSYSCQVTHEGSTVEKTVAPTECS",
      "XQPKAAPSVTLFPPSSEELQANKATLVCLISDFYPGAVTVAWKADSSPVKAGVETTTPSKQSNNKYAASSYLSLTPEQWKSHKSYSCQVTHEGSTVEKTVAPTECS",
      "GQPKAAPSVTLFPPSSEELQANKATLVCLISDFYPGAVTVAWKADSSPVKAGVETTTPSKQSNNKYAASSYLSLTPEQWKSHKSYSCQVTHEGSTVEKTVAPTECS"
    ],
  "MD:Z:35^VIVALE31^GRPLIQPRRKKAYQLEHTFQGLLGKRSLFTE10" => #ENSP00000340461 Q86UU9
    [
      "MLPCLALLLLMELSVCTVAGDGGEEQTLSTEAETWEGAGPSIQLQLQEVKTGKASQFFGLMGKRVGGREDEAQGSE",
      "MLPCLALLLLMELSVCTVAGDGGEEQTLSTEAETW------EGAGPSIQLQLQEVKTGKASQFFGLMGKRVG-------------------------------GREDEAQGSE",
      "MLPCLALLLLMELSVCTVAGDGGEEQTLSTEAETWVIVALEEGAGPSIQLQLQEVKTGKASQFFGLMGKRVGGRPLIQPRRKKAYQLEHTFQGLLGKRSLFTEGREDEAQGSE"
    ]
};

while(my ($mdz_string, $expected_results) = each %$mdz_alignment_test){
  my $alignment_strs = $pf->_mdz_alignment_string($expected_results->[0],$mdz_string);
  my $qc_status = qc_check_sequence($expected_results->[0], $$alignment_strs[0]);
  ok($qc_status == 1);
  is($expected_results->[1],$alignment_strs->[0], "Got the right target seq");
  is($expected_results->[2],$alignment_strs->[1], "Got the right query seq");
}


#dummy test
dies_ok { $pf->_mdz_alignment_string("dummy","MD:Z:") } '_mdz_alignment_string() dies ok with dummy sequence';
dies_ok { $pf->transform() } 'transform() dies ok with no features';


sub qc_check_sequence{
  my ($input_seq, $target_seq) = @_;
  
  $target_seq =~ s/-//g;
  return 0 if $input_seq ne $target_seq;
  
  return 1;

}


done_testing();

