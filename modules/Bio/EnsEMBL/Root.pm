package Bio::EnsEMBL::Root;

use Bio::Root::Root;

@ISA = qw(Bio::Root::Root);


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
