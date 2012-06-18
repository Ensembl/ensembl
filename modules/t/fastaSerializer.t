use strict;
use warnings;

use Test::More;

use IO::String;
use Bio::EnsEMBL::Utils::IO::FASTASerializer;
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
  my $ser = Bio::EnsEMBL::Utils::IO::FASTASerializer->new($io, $header, 6, 20);
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
}

done_testing();