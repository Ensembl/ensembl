# Let the code begin...
package Bio::EnsEMBL::Temp;
use vars qw(@ISA);
use strict;

use Bio::Root::Object;

@ISA = qw(Bio::Root::Object);
# new() is inherited from Bio::Root::Object

# _initialize is where the heavy stuff will happen when new is called

sub _initialize {
  my($self,@args) = @_;

  my $make = $self->SUPER::_initialize(@args);
# set stuff in self from @args
  return $make; # success - we hope!
}

=head2 gene

 Title   : gene
 Usage   : $obj->gene($newval)
 Function: stores the gene of the deleted object
 Returns : value of gene
 Args    : newvalue (optional)


=cut

sub gene{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'gene'} = $value;
    }
    return $obj->{'gene'};

}

=head2 transcript

 Title   : transcript
 Usage   : $obj->transcript($newval)
 Function: stores the transcript of the deleted object
 Returns : value of transcript
 Args    : newvalue (optional)


=cut

sub transcript{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'transcript'} = $value;
    }
    return $obj->{'transcript'};

}

