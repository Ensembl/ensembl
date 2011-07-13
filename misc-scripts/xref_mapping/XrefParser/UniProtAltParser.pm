# See UniProtParser.pm for docs.

package XrefParser::UniProtAltParser;

use strict;
use POSIX qw(strftime);
use File::Basename;

use base qw( XrefParser::BaseParser );
use base qw( XrefParser::UniProtParser );
use vars qw(@ISA);
@ISA = qw(XrefParser::UniProtParser);
my $verbose;



sub get_name {
  my $self = shift;
  my $acc  = shift;
  my $label = shift;

  return $label;
}


1;
