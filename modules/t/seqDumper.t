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

use IO::String;
use Bio::EnsEMBL::Slice;
use Bio::EnsEMBL::Utils::SeqDumper;
use Bio::EnsEMBL::Test::TestUtils;
use Bio::EnsEMBL::Test::MultiTestDB;

our $verbose = 0;

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new();

my $db = $multi->get_DBAdaptor( "core" );


my $seq_dumper = Bio::EnsEMBL::Utils::SeqDumper->new();

ok(ref($seq_dumper) && $seq_dumper->isa('Bio::EnsEMBL::Utils::SeqDumper'));

my $file;

if($verbose) {
  $file = undef;
} else {
  $file = '/dev/null';
}


#do not dump snps they are not in core db
$seq_dumper->disable_feature_type('variation');

my $slice_adaptor = $db->get_SliceAdaptor();


my $slice =
  $slice_adaptor->fetch_by_region('contig', 'AL031658.11.1.162976');

$seq_dumper->dump($slice, 'EMBL', $file);
ok(1);

$seq_dumper->dump($slice, 'GENBANK', $file);
ok(1);

$seq_dumper->dump($slice, 'FASTA', $file);
ok(1);




$slice =
  $slice_adaptor->fetch_by_region('chromosome', '20', 30_500_000, 30_600_000);

$seq_dumper->dump($slice, 'EMBL', $file);
ok(1);

$seq_dumper->dump($slice, 'GENBANK', $file);
ok(1);

$seq_dumper->dump($slice, 'FASTA', $file);
ok(1);

my $index_fh = sub {
  my ($fh, $substr) = @_;
  $fh->setpos(0);
  my @lines;
  while(my $line = <$fh>) {
    chomp $line;
    push(@lines, $line) if index($line, $substr) == 0;
  }
  return \@lines;
};

my $index_count_fh = sub {
  my ($fh, $substr) = @_;
  return scalar(@{$index_fh->($fh, $substr)});
};

{
  my $sd = Bio::EnsEMBL::Utils::SeqDumper->new();
  $sd->{feature_types}->{$_} = 0 for keys %{$sd->{feature_types}};
  
  {
    my $fh = IO::String->new();
    $sd->dump_embl($slice, $fh);
    my $lines = $index_fh->($fh, 'SQ ');
    is(scalar(@{$lines}), 1, 'Expect only 1 EMBL SQ line describing a sequence');
    is($lines->[0], 'SQ   Sequence     100001 BP;      24986 A;      24316 C;      24224 G;      26475 T;          0 other;', 'Formatting of SQ as expected');
  }

  # check if transl_table is included
  $sd = Bio::EnsEMBL::Utils::SeqDumper->new();
  $sd->{feature_types}->{$_} = 0 for keys %{$sd->{feature_types}};
  $sd->{feature_types}->{'gene'} = 1;
  
  {
    my $mt_slice = $slice_adaptor->fetch_by_region('chromosome', 'MT_NC_001807', 10060, 10405);
    my $fh = IO::String->new();
    $sd->dump_embl($mt_slice, $fh);
    my $lines = $index_fh->($fh, 'FT ');
    like( $lines->[9], qr/FT\s+\/transl_table=2/,  "Expected transl_table line at FT  CDS");
  }

  {
    my $fh = IO::String->new();
    $sd->dump_genbank($slice, $fh);
    my $lines = $index_fh->($fh, 'BASE COUNT');
    is(@{$lines}, 1, 'Expect only 1 Genbank BASE COUNT line describing a sequence');
    is($lines->[0], 'BASE COUNT       24986 a      24316 c      24224 g      26475 t          0 n', 'Formatting of BASE COUNT as expected');
  }
  
}

done_testing();
