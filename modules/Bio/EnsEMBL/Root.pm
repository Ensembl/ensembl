
=head1 DESCRIPTION

This package provides a bare-bones root class for Ensembl objects to
inheriet from. In particular it allows us to put some magic in here
to make sure we can manage 1.0 <-> 0.7 transitions for bioperl

=cut



package Bio::EnsEMBL::Root;

use vars qw(@ISA);
use strict;
use Bio::Root::RootI;

@ISA = qw(Bio::Root::RootI);


sub new{
  my($caller,@args) = @_;
  
  my $self = {};
  bless $self, $caller;

  return $self;
}



sub DESTROY{

  my ($self) = @_;


  #print STDERR "destroying $self in Bio::EnsEMBL::Root\n";
  

}

1;
