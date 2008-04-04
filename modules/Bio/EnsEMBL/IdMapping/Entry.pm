package Bio::EnsEMBL::IdMapping::Entry;

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

use strict;
use warnings;
no warnings 'uninitialized';

use Bio::EnsEMBL::Utils::Exception qw(throw warning);


sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;

  my $self = [];
  bless ($self, $class);

  return $self;
}


sub new_fast {
  my $class = shift;
  my $array_ref = shift;
  return bless $array_ref, $class;
}


sub source {
  my $self = shift;
  $self->[0] = shift if (@_);
  return $self->[0];
}


sub target {
  my $self = shift;
  $self->[1] = shift if (@_);
  return $self->[1];
}


sub score {
  my $self = shift;
  $self->[2] = shift if (@_);
  return $self->[2];
}


sub to_string {
  my $self = shift;
  return sprintf('%-10s%-10s%-5.6f', $self->source, $self->target, $self->score);
}


sub to_string_compat {
  my $self = shift;
  return join(" ", $self->source, $self->target, $self->score);
}


1;

