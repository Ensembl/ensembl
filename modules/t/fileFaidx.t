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
use Test::Differences;
use Test::Exception;

use Bio::EnsEMBL::Utils::IO qw/slurp/;
use Bio::EnsEMBL::Test::TestUtils;
use File::Temp qw/:seekable tempdir/;

use Bio::EnsEMBL::Utils::IO::FileFaidx;

my $dir = tempdir( CLEANUP => 1 );

# Test a working FASTA file with an known index
{
  my $fasta = <<'FASTA';
>TEST1
ACTGG
ACTGt
A

>TEST2
AAAA
AAAA
AAAA
>TEST3
TTTTTTTTTTTTTTTTTTTT
FASTA

  my ($tmp_fh, $fsa) = to_fsa($fasta);
  $fsa->write_index_to_disk(1);
  my $expected_lookup = {TEST1 => [11,7,5,6,'TEST1'], TEST2=> [12,29,4,5,'TEST2'], TEST3 => [20,51,20,21,'TEST3']};
  eq_or_diff($fsa->lookup(), $expected_lookup, 'Lookup generated correctly');
  
  #Check that the index was written correct
  my $expected_index = <<INDEX;
TEST1\t11\t7\t5\t6
TEST2\t12\t29\t4\t5
TEST3\t20\t51\t20\t21
INDEX
  my $actual_index = slurp("${tmp_fh}.fai");
  eq_or_diff($actual_index, $expected_index, 'Serialised index conforms to the expected format');
  
  # Check a non-generating system reading the above index back in
  
  my $non_generating_fsa = Bio::EnsEMBL::Utils::IO::FileFaidx->new("$tmp_fh");
  $non_generating_fsa->no_generation(1);
  {
    #make sure we properly mess up the code which loads the data from FASTA and the test will still work
    no warnings qw/redefine/;
    local *Bio::EnsEMBL::Utils::IO::FileFaidx::load_faindex_from_fasta = sub {
      die "This cannot happen";
    };
    eq_or_diff($non_generating_fsa->lookup(), $expected_lookup, 'Lookup was read from the file system');
  }

  my $actual_seq = $non_generating_fsa->fetch_seq('TEST1', 6, 6);
  eq_or_diff(${$actual_seq}, 'ACTGTA', 'Fetched sequence as expected');
  
  # Checking that fetches end at the limits of the sequence and we ignore lowercasing
  eq_or_diff(${$non_generating_fsa->fetch_seq('TEST1', 10, 20)}, 'TA', 'Reads off the end do not crash and we ignore case from the file');

  # Check that we can still get lower-cased characters if we need to
  $non_generating_fsa->uppercase_sequence(0);
  eq_or_diff(${$non_generating_fsa->fetch_seq('TEST1', 10, 11)}, 'tA', 'uppercase_sequence() is off so we can get lower-cased values back');
}

# 1st broken index. In-record whitespace
{
  my $fasta = <<'FASTA';
>TEST
ACTG

ACTG
A
FASTA
  my ($fh, $fsa) = to_fsa($fasta);
  throws_ok { $fsa->lookup() } qr/Line 3 is blank/, 'Inline whitespace detected';
}

# 2nd broken index. Mismatched length records
{
  my $fasta = <<'FASTA';
>Test
ACTGTTCG
ATTG
ATCGTCCC
ATTG
FASTA
  my ($fh, $fsa) = to_fsa($fasta);
  throws_ok { $fsa->lookup() } qr/Line 3 is different to the detected record length 8/, 'Mismatched record length detected';
}

# 3rd broken index. Leading whitespace
{
  my $fasta = <<'FASTA';


>Test
ACTGTTCG
FASTA
  my ($fh, $fsa) = to_fsa($fasta);
  #we warn twice; once for line 1 and once for line 2
  warns_like { $fsa->lookup() } qr/Found whitespace at line 1.+Consider trimming.+line 2/s, 'Leading FASTA whitespace detcted';
}

# 4th broken index. Not there
{
  throws_ok { my $fsa = Bio::EnsEMBL::Utils::IO::FileFaidx->new("/a/random/location/yeah") } qr/No file found at/, 'No file means exception';
}

sub to_fsa {
  my ($string) = @_;
  my $tmp = File::Temp->new(UNLINK => 1, SUFFIX => '.fa', DIR => $dir);
  print $tmp $string;
  $tmp->seek(0, SEEK_SET); #rewind to start
  my $fsa = Bio::EnsEMBL::Utils::IO::FileFaidx->new("$tmp");
  return ($tmp, $fsa);
}

done_testing();
