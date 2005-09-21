package XrefParser::CeleraProteinParser;

use strict;

use XrefParser::CeleraParser;

use vars qw(@ISA);
@ISA = qw(XrefParser::CeleraParser);

# See CeleraParser for details

sub get_sequence_type() {

  return 'peptide';

}


sub new {

  my $self = {};
  bless $self, "XrefParser::CeleraProteinParser";
  return $self;

}

1;
