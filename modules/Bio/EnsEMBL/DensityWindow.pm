package Bio::EnsEMBL::DensityWindow;

use strict;
use vars qw(@ISA);
use Bio::Root::RootI;
@ISA = qw(Bio::Root::RootI);



sub new {
    my ($class) = @_;

    my $self = {};
    bless $self,$class;

    return $self;
}




=head2 start

 Title   : start
 Usage   : $self->start($newval)
 Function: 
 Returns : value of start
 Args    : newvalue (optional)


=cut

sub start{
   my $self = shift;
   if( @_ ) {
      my $value = shift;
      $self->{'start'} = $value;
    }
    return $self->{'start'};

}

=head2 end

 Title   : end
 Usage   : $self->end($newval)
 Function: 
 Returns : value of end
 Args    : newvalue (optional)


=cut

sub end{
   my $self = shift;
   if( @_ ) {
      my $value = shift;
      $self->{'end'} = $value;
    }
    return $self->{'end'};

}

=head2 value

 Title   : value
 Usage   : $self->value($newval)
 Function: 
 Returns : value of end
 Args    : newvalue (optional)


=cut

sub value{
   my $self = shift;
   if( @_ ) {
      my $value = shift;
      $self->{'value'} = $value;
    }
    return $self->{'value'};

}
