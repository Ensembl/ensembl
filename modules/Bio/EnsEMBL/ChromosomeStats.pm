package Bio::EnsEMBL::ChromosomeStats;

use strict;
use vars qw(@ISA);
use Bio::EnsEMBL::Root;
@ISA = qw(Bio::EnsEMBL::Root);



sub new {
    my ($class) = @_;

    my $self = {};
    bless $self,$class;

    return $self;
}




=head2 chromosome_id

 Title   : chromosome_id
 Usage   : $self->chromosome_id($newval)
 Function: 
 Returns : chromosome_id
 Args    : newvalue (optional)


=cut

sub chromosome_id{
   my $self = shift;
   if( @_ ) {
      my $value = shift;
      $self->{'chromosome_id'} = $value;
    }
    return $self->{'chromosome_id'};

}

=head2 id

 Title   : id
 Usage   : $self->id($newval)
 Function: 
 Returns : id
 Args    : newvalue (optional)


=cut

sub id{
   my $self = shift;
   if( @_ ) {
      my $value = shift;
      $self->{'id'} = $value;
    }
    return $self->{'id'};

}

=head2 name

 Title   : name
 Usage   : $self->name($newval)
 Function: 
 Returns : name
 Args    : newvalue (optional)


=cut

sub name{
   my $self = shift;
   if( @_ ) {
      my $value = shift;
      $self->{'name'} = $value;
    }
    return $self->{'name'};

}

=head2 species

 Title   : species
 Usage   : $self->species($newval)
 Function: 
 Returns : species
 Args    : newvalue (optional)


=cut

sub species{
   my $self = shift;
   if( @_ ) {
      my $value = shift;
      $self->{'species'} = $value;
    }
    return $self->{'species'};

}





=head2 known_genes

 Title   : known_genes
 Usage   : $self->known_genes($newval)
 Function: 
 Returns : known_genes
 Args    : newvalue (optional)


=cut

sub known_genes{
   my $self = shift;
   if( @_ ) {
      my $value = shift;
      $self->{'known_genes'} = $value;
    }
    return $self->{'known_genes'};

}

=head2 unknown_genes

 Title   : unknown_genes
 Usage   : $self->unknown_genes($newval)
 Function: 
 Returns : unknown_genes
 Args    : newvalue (optional)


=cut

sub unknown_genes{
   my $self = shift;
   if( @_ ) {
      my $value = shift;
      $self->{'unknown_genes'} = $value;
    }
    return $self->{'unknown_genes'};

}

=head2 snps

 Title   : snps
 Usage   : $self->snps($newval)
 Function: 
 Returns : snps
 Args    : newvalue (optional)


=cut

sub snps{
   my $self = shift;
   if( @_ ) {
      my $value = shift;
      $self->{'snps'} = $value;
    }
    return $self->{'snps'};

}

=head2 length

 Title   : length
 Usage   : $self->length($newval)
 Function: 
 Returns : length
 Args    : newvalue (optional)


=cut

sub length{
   my $self = shift;
   if( @_ ) {
      my $value = shift;
      $self->{'length'} = $value;
    }
    return $self->{'length'};

}

