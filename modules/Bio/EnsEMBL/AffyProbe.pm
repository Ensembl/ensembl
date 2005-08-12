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

AffyProbes can be part of more than one Array, although not part of more than
one probeset. Onn each different Array the Probe has a slightly different name
( a number combination). A complete probename consists of the Arrayname, the setname
and the name of the probe. This probes can have a number of different names and complete
Names depending on the Array they are placed on.

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
use Bio::EnsEMBL::Storable;
use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Storable);


=head2 new

  Arg [ARRAYS] [ARRAYNAMES] :
               Probes can be generated either with a list of AffyArrays or with just the names.
               The latter is, to make the creation process easier, although you will need
               fully database persistent Arrays to store the probe. Either parameter will require
               a listref of the data.
  Arg [ARRAY] [ARRAYNAME] :
               Instead of many you can just provide one name/array.  
  Arg [NAME] [NAMES] :
               Either a single probename or a list can be provided. If a list (ref) the order
               has to be the same as in the ARRAYS/ARRAYNAMES list.
  Arg [PROBESET] :
               Each probe is part of exactly one probeset which has to be provided here.
  Example    : none
  Description: Constructor for an array probe
  Returntype : Bio::EnsEMBL::AffyProbe
  Exceptions : throws if each probe does not have a name
  Caller     : general
  Status     : Medium Risk
             : may be replaced with none affy specific methods

=cut


sub new {
  my $caller = shift;

  #allow constructor to be called as class or object method
  my $class = ref($caller) || $caller;

  my $self = $class->SUPER::new(@_);

  my(
    $arrays, $probenames, $probeset, $arraynames, $arrayname, $name, $array ) =
    rearrange(
      ['ARRAYS', 'PROBENAMES', 'PROBESET', 'ARRAYNAMES', 'ARRAYNAME', 'NAME', 'ARRAY' ], 
      @_
    );

  my $i;
  
  if( $probenames && ref( $probenames ) eq "ARRAY") {
      $i = scalar( @$probenames );
  } elsif( defined $name ) {
    if( defined $arrayname ) {
      $self->add_arrayname_probename( $arrayname, $name );
    } elsif( defined $array ) {
      $self->add_Array_probename( $array, $name );
    } else {
      throw( "Provide array or arrayname with probename" );
    }
  } else {
      throw( "You need to provide probenames to generate a probe" );
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
  } elsif( ! defined $name ) {
      throw( "Need to provide arrays or names for a probe" );
  }
  
#  $self->dbID( $probe_id ) if $probe_id;
#  $self->adaptor( $adaptor ) if $adaptor;
  $self->{'probeset'} = $probeset;

  return $self;
}



=head2 add_Array_probename

  Arg [1]    : Bio::EnsEMBL::AffyArray $array
  Arg [2]    : string $probename
  Example    : none
  Description: The pair of array and probename are added to this probe. This allows incremental
               generation of the probe.
  Returntype : none
  Exceptions : none
  Caller     : general, constructor, AffyProbeAdaptor->_obj_from_sth()
  Status     : Medium Risk
             : may be replaced with none affy specific methods

=cut


sub add_Array_probename {
    my $self = shift;
    my ( $array, $probename ) = @_;
    $self->{'arrays'} ||= {};
    $self->{'arrays'}->{$array->name()} = $array;
    $self->{'probenames'}->{$array->name()} = $probename;
}


=head2 get_all_AffyFeatures

  Args       : none
  Example    : none
  Description: All features for this probe. The probe needs to be database persistent for
               this to work.
  Returntype : listref Bio::EnsEMBL:AffyFeature
  Exceptions : none
  Caller     : general
  Status     : Medium Risk
             : may be replaced with none affy specific methods

=cut

sub get_all_AffyFeatures {
  my $self = shift;
  if( $self->adaptor() && $self->dbID() ) {
    return $self->adaptor()->db()->get_AffyFeatureAdaptor()->fetch_all_by_AffyProbe( $self );
  } else {
    warning( "Need to have attached features or database connection for that" );
    return [];
  }    
}


=head2 get_all_AffyArrays

  Args       : none
  Example    : none
  Description: Returns all arrays that are set in this probe. If there are none, and 
               the probe is persistent it tries to retrieve them from the database.
               This shouldnt really happen, if you got this from the database, its already
               properly filled. (and consequently is not implemented yet)
  Returntype : listref Bio::EnsEMBL::AffyArray
  Exceptions : none
  Caller     : general
  Status     : Medium Risk
             : may be replaced with none affy specific methods

=cut

sub get_all_AffyArrays {
    my $self = shift;
    if( defined $self->{'arrays'} ) {
	return [values %{$self->{'arrays'}}];
    } elsif( $self->adaptor() && $self->dbID()) { 
	# retrieve them by name
	my $array_adaptor = $self->adaptor()->db->get_AffyArrayAdaptor();
	warn( "Not yet implemented" );
	return [];
    } else {
	warning( "Need database connection to make Arrays from names" );
	return [];
    }
}



=head2 get_all_complete_names

  Args       : none
  Example    : none
  Description: Retreives all complete probenames from this probe. This concats the 
               Arrayname, the probesetname and the probename.
  Returntype : listref of strings
  Exceptions : none
  Caller     : general
  Status     : Medium Risk
             : may be replaced with none affy specific methods

=cut


sub get_all_complete_names {
    my $self = shift;
    my @result = ();

    while( my ( $arrayname, $probename ) = each %{$self->{'probenames'}} ) {
	push( @result , ($arrayname.":".$self->probeset().":".$probename));
    }
    return \@result;
}


=head2 get_complete_name

  Arg [1]    : string $arrayname
  Example    : none
  Description: For a given arrayname generate the completed probename and returns
               it.
  Returntype : string
  Exceptions : throws if the arrayname is not known in this probe
  Caller     : general
  Status     : Medium Risk
             : may be replaced with none affy specific methods


=cut

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


=head2 get_all_probenames

  Args       : none
  Example    : none
  Description: Probably useless method to return a list of the short probenames.
               (They really only make sense with the array)
  Returntype : listref of strings
  Exceptions : none
  Caller     : general
  Status     : Medium Risk
             : may be replaced with none affy specific methods

=cut

sub get_all_probenames {
    my $self = shift;
    return [ values %{$self->{'probenames'}} ]
}


=head2 get_probename

  Arg [1]    : string $arrayname
  Example    : none
  Description: For given arrayname return the specific probename
  Returntype : string
  Exceptions : throw if arrayname is not known for this probe
  Caller     : general
  Status     : Medium Risk
             : may be replaced with none affy specific methods

=cut


sub get_probename {
    my $self = shift;
    my $arrayname = shift;
    my $probename = $self->{'probenames'}->{$arrayname};
    if( ! defined $probename ) {
      throw( "Unknown array name" );
    }
    return $probename;
}


=head2 add_arrayname_probename

  Arg [1]    : string $arrayname
  Arg [2]    : string $probename
  Example    : none
  Description: Adds given pair of array and probename to this probe. Useless if this needs to
               be stored in the database, better to use Array objects then.
  Returntype : none
  Exceptions : none
  Caller     : general
  Status     : Medium Risk
             : may be replaced with none affy specific methods

=cut

sub add_arrayname_probename {
    my $self = shift;
    my ( $arrayname, $probename ) = @_;
    $self->{'probenames'}->{$arrayname} = $probename;
}


=head2 probeset

  Arg   [1]  : string $probeset
  Example    : none
  Description: getter setter for this probes probeset attribute
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Medium Risk
             : may be replaced with none affy specific methods

=cut

sub probeset {
    my $self = shift;
    $self->{'probeset'} = shift if( @_ );
    return $self->{'probeset'};
}

1;







