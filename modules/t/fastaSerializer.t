use strict;
use warnings;

use Test::More;
use File::Temp qw/tempfile/;

use Bio::EnsEMBL::Utils::IO::FASTASerializer;
use Bio::EnsEMBL::Slice;
use Bio::EnsEMBL::CoordSystem;

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
my $fh_SeqDumper = tempfile();

Bio::EnsEMBL::Utils::SeqDumper->dump_fasta( $slice, $fh_SeqDumper);

my $fh_Serializer = tempfile();
my $serializer = Bio::EnsEMBL::Utils::IO::FASTASerializer->new($fh_Serializer);
$serializer->print_Seq($slice);

$fh_SeqDumper->seek(0,0);
$fh_Serializer->seek(0,0);

local $/ = undef;

my $SeqDumper_output = <$fh_SeqDumper>;
my $Serializer_output = <$fh_Serializer>;

$fh_SeqDumper->close;
$fh_Serializer->close;

#print STDERR $Serializer_output."\n";

is ($SeqDumper_output,$Serializer_output,"Outputs should match from both serializers");

# Test custom header capabilities

my $custom_header = sub {
    my $slice = shift; 
    return "It's a FASTA header";
};

$fh_Serializer = tempfile();
$serializer = Bio::EnsEMBL::Utils::IO::FASTASerializer->new($fh_Serializer, $custom_header);
$serializer->print_Seq($slice);
$fh_Serializer->seek(0,0);
$Serializer_output = <$fh_Serializer>;

diag $Serializer_output;

is ($Serializer_output,">It's a FASTA header\nAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAATGAAAAAAAAAAAAAAAAAAAAAAAAAA\nAAAAAAAAAAAAAAAAACGCGCGCGCGGGA\n", "Serializer custom header should override correctly.");
$fh_Serializer->close;

done_testing();