package XrefParser::CCDSParser;

use strict;
use File::Basename;

use XrefParser::BaseParser;

use vars qw(@ISA);
@ISA = qw(XrefParser::BaseParser);

# Parse file of CCDS records and assign direct xrefs
# All assumed to be linked to transcripts

sub run {

  my ($self, $file, $source_id, $species_id) = @_;

  open(CCDS,"<".$file) || die "Could not open $file\n";

  my $count = 0;

  while (<CCDS>) {

    # TODO extract stable ID, CCDS identifier
    my ($stable_id, $ccds) = split;

    my $xref_id = $self->add_xref($ccds, 1, $ccds, "", $source_id, $species_id);
    $self->add_direct_xref($xref_id, $stable_id, "transcript", "");
    $count++;

  }

  print "Parsed $count CCDS identifiers from $file\n";

}


sub new {

  my $self = {};
  bless $self, "XrefParser::CCDSParser";
  return $self;

}

1;
