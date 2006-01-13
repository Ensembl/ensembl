
package XrefParser::JGI_ProteinParser;

use strict;

use XrefParser::JGI_Parser;

use vars qw(@ISA);
@ISA = qw(XrefParser::JGI_Parser);

# See JGI_Parser for details
sub get_sequence_type() {
  return 'peptide';
}


sub new {
  my $self = {};
  bless $self, "XrefParser::JGI_ProteinParser";
  return $self;
}

1;
