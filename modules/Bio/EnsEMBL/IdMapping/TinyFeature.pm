package Bio::EnsEMBL::IdMapping::TinyFeature;

=head1 NAME


=head1 SYNOPSIS


=head1 DESCRIPTION


=head1 METHODS


=head1 LICENCE

This code is distributed under an Apache style licence. Please see
http://www.ensembl.org/info/about/code_licence.html for details.

=head1 AUTHOR

Patrick Meidl <meidl@ebi.ac.uk>, Ensembl core API team

=head1 CONTACT

Please post comments/questions to the Ensembl development list
<ensembl-dev@ebi.ac.uk>

=cut


# internal data structure (array indices):
#
#  0  dbID
#  1  stable_id
#
# other instance variables differ by subclass implementation, so look there.


use strict;
use warnings;
no warnings 'uninitialized';

use Bio::EnsEMBL::Utils::Exception qw(throw warning);


sub new_fast {
  my $class = shift;
  my $array_ref = shift;
  return bless $array_ref, $class;
}


sub id {
  my $self = shift;
  $self->[0] = shift if (@_);
  return $self->[0];
}


sub stable_id {
  my $self = shift;
  $self->[1] = shift if (@_);
  return $self->[1];
}


1;

