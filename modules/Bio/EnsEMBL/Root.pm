package Bio::EnsEMBL::Root;

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
