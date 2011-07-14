package XrefParser::CeleraTranscriptParser;

use strict;

use base qw( XrefParser::CeleraParser );

# See CeleraParser for details

sub get_sequence_type
{
    return 'dna';
}

1;
