#
# Ensembl module for Bio::EnsEMBL::AffyProbe
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::AffyProbe - a module to represent the data for an affymetrix oligo spot.
The same oligo as part of the same Set can be found on many Arrays.


=head1 SYNOPSIS

use Bio::EnsEMBL::AffyProbe;

$feature = Bio::EnsEMBL::AffyProbe->new 
    
      -probeset => 'some setname'
)

=head1 DESCRIPTION

Affyprobe represent an oligo (probe) on one or more Affymetrix Arrays. 

=head1 CONTACT

Post comments/questions to the ensembl development list: ensembl-dev@ebi.ac.uk

=head1 METHODS

=cut


use strict;
use warnings;

package Bio::EnsEMBL::AffyProbe;

use Bio::EnsEMBL::Utils::Argument qw( rearrange ) ;
use Bio::EnsEMBL::Utils::Exception qw( throw warning );

use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Storable);

sub new {
  my $caller = shift;

  #allow constructor to be called as class or object method
  my $class = ref($caller) || $caller;

  my $self = $class->SUPER::new(@_);

  my($arrays, $probenames, $probeset, $arraynames, $arrayname, $probename, $array ) =
      rearrange(['ARRAYS', 'NAMES', 'PROBESET', 'ARRAYNAMES', 'ARRAYNAME', 'NAME', 'ARRAY' ],
		@_);

  my $i;
  if( $probenames && ref( $probenames ) eq "ARRAY") {
      $i = scalar( @$probenames );
  } elsif( defined $probename ) {
    if( defined $arrayname ) {
      $self->add_arrayname_probename( $arrayname, $probename );
    } elsif( defined $array ) {
      $self->add_Array_probename( $array, $probename );
    } else {
      throw( "Provide array or arrayname with probename" );
    }
  } else {
      throw( "You need to provide probenames to generat a probe" );
  }

  # check the parameters and fill internal data structures
  if( $arrays &&  ref( $arrays ) eq "ARRAY" ) {
      if( scalar( @$arrays ) != $i ) {
	  throw( "Need a probename for each array" );
      }
      for( $i = 0 ; $i < scalar( @$probenames ); $i++ ) {
	  $self->add_Array_probename( $arrays->[$i], $probenames->[$i] );
      }
  } elsif( $arraynames && ref( $arraynames )) {
      if( scalar( @$arraynames ) != $i ) {
	  throw( "Need a probename for each array" );
      }
      for( $i = 0 ; $i < scalar( @$probenames ); $i++ ) {
	  $self->add_arrayname_probename( $arraynames->[$i], $probenames->[$i] );
      }
  } elsif( ! defined $probename ) {
      throw( "Need to provide arrays or names for a probe" );
  }
  
  $self->{'probeset'} = $probeset;

  return $self;
}


sub add_Array_probename {
    my $self = shift;
    my ( $array, $probename ) = @_;
    $self->{'arrays'} ||= {};
    $self->{'arrays'}->{$array->name()} = $array;
    $self->{'probenames'}->{$array->name()} = $probename;
}

sub get_all_AffyFeatures {
    my $self = shift;
    if( $self->adaptor() && $self->dbID() ) {
	
    } elsif( defined $self->{'features'} ) {
	return $self->{'features'};
    } else {
	warning( "Need to have attached features or database connection for that" );
	return [];
    }    
}

sub get_all_AffyArrays {
    my $self = shift;
    if( defined $self->{'arrays'} ) {
	return [values %{$self->{'arrays'}}];
    } elsif( $self->adaptor() && $self->dbID()) { 
	# retrieve them by name
	my $array_adaptor = $self->adaptor()->db->get_AffyArrayAdaptor();
    } else {
	warning( "Need database connection to make Arrays from names" );
	return [];
    }
}


sub get_all_complete_names {
    my $self = shift;
    my @result = ();

    while( my ( $arrayname, $probename ) = each %{$self->{'probenames'}} ) {
	push( @result , ($arrayname.":".$self->probeset().":".$probename));
    }
    return \@result;
}

sub get_complete_name {
    my $self = shift;
    my $arrayname = shift;

    my $probename = $self->{'probenames'}->{$arrayname};
    if(  ! defined $probename ) {
	throw( "Unknown array name" );
    } 
    $probename = $arrayname.":".$self->probeset().":".$probename;
    return $probename;
}

sub get_all_probenames {
    my $self = shift;
    return [ values %{$self->{'probenames'}} ]
}

sub get_probename {
    my $self = shift;
    my $arrayname = shift;
    my $probename = $self->{'probenames'}->{$arrayname};
    if( ! defined $probename ) {
	throw( "Unknown array name" );
    }
    return $probename;
}

sub add_arrayname_probename {
    my $self = shift;
    my ( $arrayname, $probename ) = @_;
    $self->{'probenames'}->{$arrayname} = $probename;
}



# if the probe is set, take it from there
sub probeset {
    my $self = shift;
    $self->{'probeset'} = shift if( @_ );
    return $self->{'probeset'};
}

1;







