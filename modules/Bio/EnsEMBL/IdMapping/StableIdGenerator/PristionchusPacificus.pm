=head1 LICENSE

  Copyright 2010, The WormBase Consortium. All rights reserved.

  Redistribution and use in source and binary forms, with or without modification, are
  permitted provided that the following conditions are met:

   1. Redistributions of source code must retain the above copyright notice, this list of
      conditions and the following disclaimer.

   2. Redistributions in binary form must reproduce the above copyright notice, this list
      of conditions and the following disclaimer in the documentation and/or other materials
      provided with the distribution.

  THIS SOFTWARE IS PROVIDED BY THE WORMBASE CONSORTIUM ``AS IS'' AND ANY EXPRESS OR IMPLIED
  WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
  FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> OR
  CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
  ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
  NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
  ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

  The views and conclusions contained in the software and documentation are those of the
  authors and should not be interpreted as representing official policies, either expressed
  or implied, of The WormBase Consortium.

=head1 CONTACT

  Please email comments or questions to the public WormBase
  help list at <help@wormbase.org>.

=cut

package Bio::EnsEMBL::IdMapping::StableIdGenerator::PristionchusPacificus;

# Class to generate WormBase conform Pristionchus IDs to be injected into the mapper

use strict;
use warnings;

use base qw(Bio::EnsEMBL::IdMapping::StableIdGenerator::EnsemblGeneric);

# PPAxxxxx
# ups the id by one
sub increment_stable_id {
  my ( $self, $lastId ) = @_;

  throw("invalid stable ID: $lastId.") unless ($lastId=~/PPA/);

  $lastId =~ /^PPA(\d+)/;

  my $number = $1+1;
  my $stable_id = sprintf("PPA%05d",$number);

  return $stable_id;
}

# just in case it is actually used somewhere
sub is_valid {
  my ( $self, $stableId ) = @_;

  if ( defined($stableId) && $stableId =~ /^PPA\d+/)
     {return 1} else {return undef}
}

1;
