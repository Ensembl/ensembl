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
use Bio::EnsEMBL::BaseAlignFeature;
use Bio::EnsEMBL::FeaturePair;
use Bio::EnsEMBL::ProteinFeature;
use Bio::EnsEMBL::Analysis;


my $multi_db = Bio::EnsEMBL::Test::MultiTestDB->new;
my $db = $multi_db->get_DBAdaptor('core');
my $pfa = $db->get_ProteinFeatureAdaptor();

my $analysis_db = 'test_db';

my $analysis = Bio::EnsEMBL::Analysis->new(-LOGIC_NAME => 'gifts_import', -DB => $analysis_db);

my $seq_start  = 1;
my $seq_end    = 111;
my $hit_start = 0;
my $hit_end   = 0;
my $hit_name = "P20366";
my $cigar_string = "MD:Z:96^RHKTDSFVGLMGKRALNS0V14";
my $align_type = "mdtag";
my $translation_id = 21741;

my $paf = Bio::EnsEMBL::ProteinFeature->new
  (-START       => $seq_start,
   -END         => $seq_end,
   -ANALYSIS    => $analysis,
   -HSTART      => $hit_start,
   -HEND        => $hit_end,
   -HSEQNAME    => $hit_name,
   -CIGAR_STRING => $cigar_string,
   -ALIGN_TYPE => $align_type,
   -TRANSLATION_ID     => $translation_id,
   -ADAPTOR           => $pfa
   );


ok($paf && $paf->isa('Bio::EnsEMBL::ProteinFeature'));



# Test alignment strings

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
  $paf->{'_alignment_length'} = undef;
  $paf->cigar_string($mdz_string);
  ok($paf->cigar_string eq $mdz_string);
  ok( $paf->alignment_length() == $length, "Got back correct alignment length $length for mdz string ". $paf->cigar_string);
}


my $mdz_alignment_test = {
  "MD:Z:96^RHKTDSFVGLMGKRALNS0V14" => # mdz string
    [ 
      "MKILVALAVFFLVSTQLFAEEIGANDDLNYWSDWYDSDQIKEELPEPFEHLLQRIARRPKPQQFFGLMGKRDADSSIEKQVALLKALYGHGQISHKMAYERSAMQNYERRR",  # input seq
      "MKILVALAVFFLVSTQLFAEEIGANDDLNYWSDWYDSDQIKEELPEPFEHLLQRIARRPKPQQFFGLMGKRDADSSIEKQVALLKALYGHGQISHK------------------MAYERSAMQNYERRR", # target seq
      "MKILVALAVFFLVSTQLFAEEIGANDDLNYWSDWYDSDQIKEELPEPFEHLLQRIARRPKPQQFFGLMGKRDADSSIEKQVALLKALYGHGQISHKRHKTDSFVGLMGKRALNSVAYERSAMQNYERRR" # query seq
    ],
  "MD:Z:35S0R0Y0I0V0M0G0H0I1^HKRRQLPTA0L1Q0V0L0R0G0R0L0R0P0G0D0G0L0L0R0S0S0S0S0Y0V0K1F0N0R0K0R0^EGQIQGAIHTQCI" => 
    [
      "MNRLYLTPDGFFFRVHMLALDSSSCNKPCPEFKPGIETDLNDAAYVLYTTVCNVGATARAVGRPAFFWERWETMT",
      "MNRLYLTPDGFFFRVHMLALDSSSCNKPCPEFKPGIETDLNDAAY---------VLYTTVCNVGATARAVGRPAFFWERWETMT-------------",
      "MNRLYLTPDGFFFRVHMLALDSSSCNKPCPEFKPGSRYIVMGHIYHKRRQLPTALLQVLRGRLRPGDGLLRSSSSYVKRFNRKREGQIQGAIHTQCI"
    ],
  "MD:Z:0G105" => 
    [
      "XQPKAAPSVTLFPPSSEELQANKATLVCLISDFYPGAVTVAWKADSSPVKAGVETTTPSKQSNNKYAASSYLSLTPEQWKSHKSYSCQVTHEGSTVEKTVAPTECS",
      "XQPKAAPSVTLFPPSSEELQANKATLVCLISDFYPGAVTVAWKADSSPVKAGVETTTPSKQSNNKYAASSYLSLTPEQWKSHKSYSCQVTHEGSTVEKTVAPTECS",
      "GQPKAAPSVTLFPPSSEELQANKATLVCLISDFYPGAVTVAWKADSSPVKAGVETTTPSKQSNNKYAASSYLSLTPEQWKSHKSYSCQVTHEGSTVEKTVAPTECS"
    ],
  "MD:Z:35^VIVALE31^GRPLIQPRRKKAYQLEHTFQGLLGKRSLFTE10" => 
    [
      "MLPCLALLLLMELSVCTVAGDGGEEQTLSTEAETWEGAGPSIQLQLQEVKTGKASQFFGLMGKRVGGREDEAQGSE",
      "MLPCLALLLLMELSVCTVAGDGGEEQTLSTEAETW------EGAGPSIQLQLQEVKTGKASQFFGLMGKRVG-------------------------------GREDEAQGSE",
      "MLPCLALLLLMELSVCTVAGDGGEEQTLSTEAETWVIVALEEGAGPSIQLQLQEVKTGKASQFFGLMGKRVGGRPLIQPRRKKAYQLEHTFQGLLGKRSLFTEGREDEAQGSE"
    
    ]
};

while(my ($mdz_string, $expected_results) = each %$mdz_alignment_test){
  my $alignment_strs = $paf->_mdz_alignment_string($$expected_results[0],$mdz_string);
  my $qc_status = qc_check_sequence($$expected_results[0], $$alignment_strs[0]);
  ok($qc_status == 1);
  ok($$expected_results[1] eq $$alignment_strs[0], "Got the right target seq");
  ok($$expected_results[2] eq $$alignment_strs[1], "Got the right query seq");

}

sub qc_check_sequence{
  my ($input_seq, $target_seq) = @_;
  
  $target_seq =~ s/-//g;
  return 0 if $input_seq ne $target_seq;
  
  return 1;

}

done_testing();
