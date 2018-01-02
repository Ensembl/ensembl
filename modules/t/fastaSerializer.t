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

use IO::String;
use Bio::EnsEMBL::Utils::IO::FASTASerializer;
use Bio::EnsEMBL::Test::TestUtils qw/warns_like/;
use Bio::EnsEMBL::Slice;
use Bio::EnsEMBL::CoordSystem;
use Bio::Seq;

#
# TEST - Slice creation from adaptor
#
my $coord_system = Bio::EnsEMBL::CoordSystem->new(
  -NAME => 'chromosome', -RANK => 1
);

# instantiate slice
#SEQ COORD_SYSTEM SEQ_REGION_NAME SEQ_REGION_LENGTH
#                      START END STRAND ADAPTOR EMPTY

my $slice = Bio::EnsEMBL::Slice->new( 
    -SEQ_REGION_NAME => "top_banana",
    -COORD_SYSTEM => $coord_system,
    -STRAND => 1,
    -START => 110,
    -END => 199,
    -SEQ_REGION_LENGTH => 90,
    -SEQ => "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAATGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACGCGCGCGCGGGA",
);

# ensure Serializer produces output identical to the well-used SeqDumper. 
my $fh_SeqDumper = IO::String->new();

Bio::EnsEMBL::Utils::SeqDumper->dump_fasta( $slice, $fh_SeqDumper);

my $fh_Serializer = IO::String->new();
my $serializer = Bio::EnsEMBL::Utils::IO::FASTASerializer->new($fh_Serializer);
$serializer->print_Seq($slice);

my $SeqDumper_output = ${$fh_SeqDumper->string_ref()};
my $Serializer_output = ${$fh_Serializer->string_ref()};

$fh_SeqDumper->close;
$fh_Serializer->close;

#print STDERR $Serializer_output."\n";

is ($SeqDumper_output,$Serializer_output,"Outputs should match from both serializers");

# Test custom header capabilities

my $custom_header = sub {
    my $slice = shift; 
    return "It's a FASTA header";
};

$fh_Serializer = IO::String->new();
$serializer = Bio::EnsEMBL::Utils::IO::FASTASerializer->new($fh_Serializer, $custom_header);
$serializer->print_Seq($slice);
$Serializer_output = ${$fh_Serializer->string_ref()};

note $Serializer_output;

is ($Serializer_output,">It's a FASTA header\nAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAATGAAAAAAAAAAAAAAAAAAAAAAAAAA\nAAAAAAAAAAAAAAAAACGCGCGCGCGGGA\n", "Serializer custom header should override correctly.");
$fh_Serializer->close;

{
  my $seq = 'A'x120;
  my $s = Bio::EnsEMBL::Slice->new(-SEQ_REGION_NAME => 'a', -COORD_SYSTEM => $coord_system, -SEQ => $seq, -SEQ_REGION_LENGTH => 120, -START => 1, -END => 120);
  my $header = sub { return 'a'; };
  my $io = IO::String->new();
  my $ser = Bio::EnsEMBL::Utils::IO::FASTASerializer->new($io, $header, 600, 20); # chunk is smaller than line width
  my $expected = <<'FASTA';
>a
AAAAAAAAAAAAAAAAAAAA
AAAAAAAAAAAAAAAAAAAA
AAAAAAAAAAAAAAAAAAAA
AAAAAAAAAAAAAAAAAAAA
AAAAAAAAAAAAAAAAAAAA
AAAAAAAAAAAAAAAAAAAA
FASTA
  $ser->print_Seq($s);
  is(${$io->string_ref()}, $expected, 'Testing round number serialisation');
}

{
  my $seq = 'A'x21;
  my $s = Bio::EnsEMBL::Slice->new(-SEQ_REGION_NAME => 'a', -COORD_SYSTEM => $coord_system, -SEQ => $seq, -SEQ_REGION_LENGTH => 21, -START => 1, -END => 21);
  my $header = sub { return 'a'; };
  my $io = IO::String->new();
  my $ser = Bio::EnsEMBL::Utils::IO::FASTASerializer->new($io, $header, 1, 20);
  my $expected = <<'FASTA';
>a
AAAAAAAAAAAAAAAAAAAA
A
FASTA
  $ser->print_Seq($s);
  is(${$io->string_ref()}, $expected, 'Testing odd (as in strange) number line length serialisation');
}

{
  my $seq = Bio::Seq->new(-SEQ => 'M', -DISPLAY_ID => 'A');
  my $io = IO::String->new();
  my $ser = Bio::EnsEMBL::Utils::IO::FASTASerializer->new($io);
  my $expected = <<'FASTA';
>A
M
FASTA
  $ser->print_Seq($seq);
  is(${$io->string_ref()}, $expected, 'Testing single base serialisation');
  my $header = sub { return 'A'; };
  $io = IO::String->new();
  $ser = Bio::EnsEMBL::Utils::IO::FASTASerializer->new($io, $header, 1, '100000');
  $ser->print_Seq($seq);
  is(${$io->string_ref()}, $expected, 'Test with crazy line length');
  $io->close;
}
throws_ok( sub {
  my $io = IO::String->new;
  my $ser = Bio::EnsEMBL::Utils::IO::FASTASerializer->new($io, 'Hello', 1, '1000000000000');
}, qr/Maximum line width/, 'Test really crazy line length');

throws_ok( sub { Bio::EnsEMBL::Utils::IO::FASTASerializer->new(undef, 'Hello', 1, '+INF'); }, qr/Maximum line width/,'Stupid line length input');


done_testing();
