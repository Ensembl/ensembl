package XrefParser::CeleraTranscriptParser;

use strict;

use XrefParser::CeleraParser;

use vars qw(@ISA);
@ISA = qw(XrefParser::CeleraParser);

# See CeleraParser for details

sub get_sequence_type() {

  return 'dna';

}


sub new {

  my $self = {};
  bless $self, "XrefParser::CeleraTranscriptParser";
  return $self;

}

1;
