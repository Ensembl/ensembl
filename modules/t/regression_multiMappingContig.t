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

# TEST CASE TO RECREATE HUMAN CONTIG MAPPING ISSUE
# 
# GRCh38 used the same contig in full in multiple locations on the same chromosome. The
# original mapper code would only give the 1st mapping & then ignores anything else.
# Here we have 3 fake contigs each composed of 4 bases; A, C & T. These map to chr1
# and only contigs A & C are repeated.

use strict;
use warnings;

use Test::More;
use Test::Warnings;
use Test::Differences;
use Bio::EnsEMBL::Test::MultiTestDB;

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new('mapping');
my $dba = $multi->get_DBAdaptor('core');
ok(defined $dba, 'Test database instatiated');

my $seqlevel = $dba->get_CoordSystemAdaptor()->fetch_sequence_level();
my $slice = $dba->get_SliceAdaptor->fetch_by_toplevel_location('chr1:1-20');

my $expected_seq = 'AAAACCCCAAAACCCCTTTT';
is(length($slice->seq), length($expected_seq), "Checking sequence emitted is the right length");
is($slice->seq(), $expected_seq, 'Checking expected sequence emitted');

my $projection = $slice->project($seqlevel->name(), $seqlevel->version());
my $expected_projections = [
  [1, 4, 'contig::A:1:4:1'],
  [5, 8, 'contig::C:1:4:1'],
  [9, 12, 'contig::A:1:4:1'],
  [13, 16, 'contig::C:1:4:1'],
  [17, 20, 'contig::T:1:4:1'],
];
my @converted_projections = map { [$_->[0], $_->[1], $_->[2]->name()] } @{$projection};
eq_or_diff(\@converted_projections, $expected_projections, 'Checking projections has all values');

done_testing();
