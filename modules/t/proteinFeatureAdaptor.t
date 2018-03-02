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

use Bio::EnsEMBL::Test::TestUtils;

use Test::More;
use Test::Warnings;

use Bio::EnsEMBL::Test::MultiTestDB;
use Test::Exception;

our $verbose = 0;

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new;

# get a core DBAdaptor
my $dba = $multi->get_DBAdaptor("core");

#
# Test get_ProteinFeatureAdaptor works
#
my $pfa = $dba->get_ProteinFeatureAdaptor();

ok($pfa && ref($pfa) && $pfa->isa('Bio::EnsEMBL::DBSQL::ProteinFeatureAdaptor'));

my $pfs = $pfa->fetch_all_by_translation_id(21724);

print_features($pfs);

ok(@$pfs == 15);

sub print_features {
  my $features = shift;
  foreach my $f (@$features) {
	if (defined($f)) {
	  debug($f->start . '-' . $f->end . ' -> ' . $f->hseqname . ':' . $f->hstart . '-' . $f->hend);
	}
  }
}

# test adding and retrieving a new feature
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
my $interpro_ac = 'interpro accession';

my $analysis_db = 'test_db';

my $analysis = Bio::EnsEMBL::Analysis->new(-LOGIC_NAME => 'test', -DB => $analysis_db);
$multi->save('core', 'protein_feature', 'meta_coord');


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
   -INTERPRO_AC => $interpro_ac);
   

my $summary = $f->summary_as_hash();
is($summary->{'type'}, $analysis_db);
is($summary->{'id'}, $hseqname);
is($summary->{'start'}, $start);
is($summary->{'end'}, $end);
is($summary->{'interpro'}, $interpro_ac);
is($summary->{'description'}, $idesc);
is($summary->{'hit_start'}, $hstart);
is($summary->{'hit_end'}, $hend);   
   
my $dbID = $pfa->store($f,21724);

ok($dbID, "New object created has dbID");

#fetch protein feature object by dbID and compare the core attributes
#deep comparison is not possible because all the attributes are not stored in the protein_feature table
my $f_by_dbID = $pfa->fetch_by_dbID($dbID);

is($f_by_dbID->translation_id, 21724);
is($f_by_dbID->start, $start);
is($f_by_dbID->end, $end);
is($f_by_dbID->dbID, $dbID);
is($f_by_dbID->analysis->logic_name, 'test');
is($f_by_dbID->hstart, $hstart);
is($f_by_dbID->hend, $hend);
is($f_by_dbID->hseqname, $hseqname);
is($f_by_dbID->hdescription, $hdes);
is($f_by_dbID->score, $score);
is($f_by_dbID->percent_id, $percent_id);
is($f_by_dbID->p_value, $p_value);


$pfs = $pfa->fetch_all_by_translation_id(21724);

ok(@$pfs == 16);

my @pfs = grep{$_->hdescription() eq $hdes} @$pfs;

ok(scalar @pfs > 0);

$multi->restore('core', 'protein_feature');

$pfs = $pfa->fetch_all();
is(@$pfs, 157, "Retrieved all protein features");

$pfs = $pfa->fetch_all_by_logic_name('pfscan');
is(@$pfs, 156, "Retrieved pfscan features");


my $test_logic_name = 'gifts_import';
my $test_align_type = 'mdtag';
my $test_cigar_string = 'MD:Z:35^VIVALE31^GRPLIQPRRKKAYQLEHTFQGLLGKRSLFTE10';
my $test_dbID = 242847;
my $test_hitname = 'Q86UU9';
my $transl_id = 21739;

# fetch_all_by_translation_id
ok($pfa && $pfa->isa('Bio::EnsEMBL::DBSQL::ProteinFeatureAdaptor'));
my $trl_features = $pfa->fetch_all_by_translation_id($transl_id);

my @trl_gift_features = grep { $_->analysis()->logic_name eq $test_logic_name } @$trl_features;
my $trl_feature = shift @trl_gift_features;

ok($trl_feature && $trl_feature->isa('Bio::EnsEMBL::ProteinFeature'));
ok($trl_feature->analysis()->logic_name eq $test_logic_name, 'Got the right logic name ' .$test_logic_name);

ok($trl_feature->cigar_string eq $test_cigar_string, 'Got the right cigar string ' . $test_cigar_string);
ok($trl_feature->align_type eq $test_align_type, 'Got the right align type ' . $test_align_type );

# fetch_all_by_translation_id and logic_name
$trl_features = $pfa->fetch_all_by_translation_id($transl_id, "gifts_import");

@trl_gift_features = grep { $_->analysis()->logic_name eq $test_logic_name } @$trl_features;
$trl_feature = shift @trl_gift_features;

ok($trl_feature && $trl_feature->isa('Bio::EnsEMBL::ProteinFeature'));
ok($trl_feature->analysis()->logic_name eq $test_logic_name, 'Got the right logic name ' .$test_logic_name);

ok($trl_feature->cigar_string eq $test_cigar_string, 'Got the right cigar string ' . $test_cigar_string);
ok($trl_feature->align_type eq $test_align_type, 'Got the right align type ' . $test_align_type );


# fetch_all_by_logic_name
my $logic_name_features = $pfa->fetch_all_by_logic_name('gifts_import');
my $logic_name_feature = shift @$logic_name_features;
ok($logic_name_feature->cigar_string eq $test_cigar_string);
ok($logic_name_feature->align_type eq $test_align_type);
ok($logic_name_feature->analysis()->logic_name eq $test_logic_name);

# fetch_all_by_dbID_list
my $dbid_features = $pfa->fetch_all_by_dbID_list ([$test_dbID]);
my $dbid_feature = shift @$dbid_features;
ok($dbid_feature->cigar_string eq $test_cigar_string);
ok($dbid_feature->align_type eq $test_align_type);
ok($dbid_feature->analysis()->logic_name eq $test_logic_name);

$dbid_feature = $pfa->fetch_by_dbID($test_dbID);
ok($dbid_feature->cigar_string eq $test_cigar_string);
ok($dbid_feature->align_type eq $test_align_type);
ok($dbid_feature->analysis()->logic_name eq $test_logic_name);


# fetch_all_by_hit_name. Test with uniprot accession
my $hitname_features = $pfa->fetch_all_by_hit_name($test_hitname, $test_logic_name);
my $hitname_feature = shift @$hitname_features;
ok($hitname_feature->cigar_string eq $test_cigar_string);
ok($hitname_feature->align_type eq $test_align_type);
ok($hitname_feature->analysis()->logic_name eq $test_logic_name, "Got the right logicname " . $test_logic_name );
ok($hitname_feature->translation_id eq $transl_id, "Got the right transl_id " . $transl_id );

# fetch_all_by_uniprot_acc. Test with uniprot accession
my $hitname_uniprot_features = $pfa->fetch_all_by_uniprot_acc($test_hitname);
my $hitname_uniprot_feature = shift @$hitname_uniprot_features;
ok($hitname_uniprot_feature->cigar_string eq $test_cigar_string);
ok($hitname_uniprot_feature->align_type eq $test_align_type);
ok($hitname_uniprot_feature->analysis()->logic_name eq $test_logic_name, "Got the right logicname " . $test_logic_name );
ok($hitname_uniprot_feature->translation_id eq $transl_id, "Got the right transl_id " . $transl_id );

# Test inherited methods from BaseAlignFeatureAdaptor
dies_ok { $pfa->fetch_all_by_Slice_and_hcoverage() } 'fetch_all_by_Slice_and_hcoverage() dies ok with no features';

dies_ok { $pfa->fetch_all_by_Slice_and_external_db() } 'fetch_all_by_Slice_and_hcoverage() dies ok with no features';

dies_ok { $pfa->fetch_all_by_Slice() } 'fetch_all_by_Slice() dies ok with no features';

dies_ok { $pfa->fetch_all_by_Slice_and_score() } 'fetch_all_by_Slice_and_score() dies ok with no features';

dies_ok { $pfa->fetch_all_by_Slice_constraint() } 'fetch_all_by_Slice_constraint() dies ok with no features';

dies_ok { $pfa->fetch_all_by_stable_id_list() } 'fetch_all_by_stable_id_list() dies ok with no features';

dies_ok { $pfa->count_by_Slice_constraint() } 'count_by_Slice_constraint() dies ok with no features';

dies_ok { $pfa->remove_by_Slice() } 'remove_by_Slice() dies ok with no features';

dies_ok { $pfa->get_seq_region_id_internal() } 'get_seq_region_id_internal() dies ok with no features';

dies_ok { $pfa->get_seq_region_id_external() } 'get_seq_region_id_external() dies ok with no features';





done_testing();
