package XrefParser::CeleraProteinParser;

use strict;

use base qw( XrefParser::CeleraParser );

# See CeleraParser for details

sub get_sequence_type
{
    return 'peptide';
}

1;
