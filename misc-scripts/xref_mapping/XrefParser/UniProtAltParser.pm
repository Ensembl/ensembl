=head1 LICENSE

Copyright [1999-2014] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=cut

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
