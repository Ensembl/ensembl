#
# EnsEMBL module for Bio::EnsEMBL::DBOLD::MetaContainer
#
# Cared for by Arne Stabenau
#
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

  Bio::EnsEMBL::DBOLD::MeatContainer - Encapsulates all access to database meta information

=head1 SYNOPSIS


=head1 DESCRIPTION


=head1 CONTACT



=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are 
usually preceded with a _

=cut


# Let the code begin...

package Bio::EnsEMBL::DBOLD::MetaContainer;

use vars qw(@ISA);
use strict;

# Object preamble - inheriets from Bio::Root::Object


use Bio::EnsEMBL::DBOLD::BaseAdaptor;
use Bio::EnsEMBL::DBOLD::DBAdaptor;

@ISA = qw(Bio::EnsEMBL::DBOLD::BaseAdaptor);

# new() is inherited from Bio::EnsEMBL::DBOLD::BaseAdaptor


=head2 list_value_by_key

 Title   : list_value_by_key
 Usage   : 
 Function: gets a value for a key. Can be anything
 Example : $metaContainer->list_value_by_key
 Returns : a list of values
 Args    : a key string

=cut


sub list_value_by_key {
  my ($self,$key) = @_;
  my @result;
  
  my $sth = $self->prepare( "select meta_value from meta where meta_key = ? order by meta_id" );
  $sth->execute( $key );
  while( my $arrRef = $sth->fetchrow_arrayref() ) {
    push( @result, $arrRef->[0] );
  }
  
  return @result;
}


sub store_key_value {
  my ( $self, $key, $value ) = @_;
  # store it
  my $sth = $self->prepare( "insert into meta( meta_key, meta_value) values( ?, ? )" );
  my $res = $sth->execute( $key, $value );
  return;
}


sub create_tables {
  my $self = shift;
  my $sth = $self->prepare( "drop table if exists meta" );
  $sth->execute();

  $sth = $self->prepare( "
     CREATE TABLE meta (
        meta_id INT not null auto_increment,
        meta_key varchar( 40 ) not null,
        meta_value varchar( 255 ) not null,

        PRIMARY KEY( meta_id ),
        KEY meta_key_index ( meta_key ),
        KEY meta_value_index ( meta_value ))
  " );
  $sth->execute();
}



# add well known meta info get-functions below

sub get_Species {
  my $self = shift;
  my $sth = $self->prepare( "select meta_value from meta where meta_key = 'species.common_name'" );
  $sth->execute;
  my $common_name;
  if( my $arrRef = $sth->fetchrow_arrayref() ) {
    $common_name = $arrRef->[0];
  } else {
    return undef;
  }
  my @classification = $self->list_value_by_key( 'species.classification' );
  if( scalar(@classification) == 0 ) {
    return undef;
  }
  my $species = new Bio::Species;
  $species->common_name( $common_name );
  $species->classification( @classification );

  return $species;
}


1;

