
package XrefParser::JGI_ProteinParser;

use strict;

use base qw( XrefParser::JGI_Parser );

# See JGI_Parser for details

sub get_sequence_type {
  return 'peptide';
}

1;
