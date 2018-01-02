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
use Test::Exception;

use Bio::EnsEMBL::Test::TestUtils;
use Bio::EnsEMBL::Test::MultiTestDB;

use File::Spec;
use File::Temp qw/:seekable tempdir/;
use Bio::EnsEMBL::Utils::IO::FASTASerializer;
use Bio::EnsEMBL::DBSQL::FastaSequenceAdaptor;
use Bio::EnsEMBL::Utils::IO::FileFaidx;
use Bio::EnsEMBL::Utils::IO qw/work_with_file/;

# "Globals"
my $multi_db = Bio::EnsEMBL::Test::MultiTestDB->new;
my $dba = $multi_db->get_DBAdaptor('core');

my $CHR           = '20';
my $START         = 30_220_000;
my $END           = 31_200_000;
my $STRAND        = 1;

my $DIR = tempdir( CLEANUP => 1 );
my $FASTA_FILE = File::Spec->catfile($DIR, 'chromosome.fa');
my $SA = $dba->get_SliceAdaptor();

# Write sequence to a file and index
{
  work_with_file($FASTA_FILE, 'w', sub {
    my ($fh) = @_;
    my $fasta = Bio::EnsEMBL::Utils::IO::FASTASerializer->new($fh, sub {
      my ($slice) = @_;
      return $slice->seq_region_name();
    });
    my $slice = $SA->fetch_by_region('chromosome', $CHR);
    $fasta->print_Seq($slice);
  });
}

# Now index and compare
{
  my $fidx = Bio::EnsEMBL::Utils::IO::FileFaidx->new($FASTA_FILE);
  my $slice = $SA->fetch_by_region('chromosome', $CHR, $START, $END, $STRAND);
  my $db_seq = $slice->seq();

  # Setup switchable adaptor and test
  my $power = 20;
  my $fsa = Bio::EnsEMBL::DBSQL::FastaSequenceAdaptor->new($fidx, $power); #request 1MB cache chunks
  $dba->switch_adaptor('sequence', $fsa, sub {
    my $fasta_seq = $slice->seq();
    is($db_seq, $fasta_seq, 'Checking both adaptors return the same sequence');
    
    my $sa = $dba->get_SequenceAdaptor();
    ok($sa->isa('Bio::EnsEMBL::DBSQL::FastaSequenceAdaptor'), 'Checking replacement worked as expected');
    dies_ok { $dba->get_SequenceAdaptor()->store() } 'Call should die due to FASTASequenceAdaptor being in place';

    # Drill into the object and look at the cache
    my ($cache_key) = keys %{$sa->{seq_cache}};
    my $expected_length = 1 << $power;
    is(length(${$sa->{seq_cache}->{$cache_key}}), $expected_length, 'Checking length of cached sequence');
  });
}

done_testing();
