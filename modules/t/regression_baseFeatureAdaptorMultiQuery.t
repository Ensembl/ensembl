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
use Bio::EnsEMBL::Test::MultiTestDB;
use DBI qw/:sql_types/;

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new();
my $db = $multi->get_DBAdaptor( "core" );

my $sa = $db->get_SliceAdaptor();
my $ga = $db->get_GeneAdaptor();

my $name  = '20_HAP1';
my $start = undef;
my $end = undef;
my $strand = 1;

my $slice = $sa->fetch_by_region('toplevel', $name, $start, $end, $strand);
$ga->bind_param_generic_fetch('protein_coding', SQL_VARCHAR);
my $genes = $ga->fetch_all_by_Slice_constraint($slice, 'g.biotype =?');

# Logic here is that we divide our slice into 3; prior HAP (Chr20), HAP & post HAP (Chr20)
# and each one needs the bind parameters set for each query. By the time we
# get into BaseFeatureAdaptor's internals we can record and replay them everytime
# we run the same query on a different Slice

is(scalar(@{$genes}), 14, 'Should have 14 genes when we query slices which need normalisation');

done_testing();
