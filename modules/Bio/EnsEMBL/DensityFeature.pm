=head1 LICENSE

  Copyright (c) 1999-2012 The European Bioinformatics Institute and
  Genome Research Limited.  All rights reserved.

  This software is distributed under a modified Apache license.
  For license details, please see

    http://www.ensembl.org/info/about/code_licence.html

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <dev@ensembl.org>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk@ensembl.org>.

=cut

=head1 NAME

Bio::EnsEMBL::DensityFeature - A feature representing a density, or
precentage coverage etc. in a given region.

=head1 SYNOPSIS

  use Bio::EnsEMBL::DensityFeature;

  $feature = Bio::EnsEMBL::DensityFeature->new(
    -seq_region    => $region,
    -start         => 1,
    -end           => 1e6,
    -density_type  => $dt,
    -density_value => 98.5
  );

=head1 DESCRIPTION

A density feature represents a count, density, or percentage coverage,
etc. for a given region.

This module is part of the Ensembl project http://www.ensembl.org

=head1 METHODS

=cut


use strict;
use warnings;

package Bio::EnsEMBL::DensityFeature;

use Bio::EnsEMBL::Feature;
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Utils::Exception qw(throw);
use Bio::EnsEMBL::DensityType;

use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Feature);


=head2 new

  Arg [SEQ_REGION] : the sequence over which the density was calculated.

  Arg [START] : start point on the seq at which density was calulated.

  Arg [END] : end point on the seq at which density was calulated.

  Arg [DENSITY_TYPE] : the type of density calculated.

  Arg [DENSITY_VALUE] : the density.

  Arg [...]  : Named arguments passed to superclass
  Example    : $feature = Bio::EnsEMBL::DensityFeature->new
                            (-seq_region    => $region,
			     -start         => 1,
                             -end           => 1e6,
                             -density_type  => $dt,
			     -density_value => 98.5)

  Description: Creates a new density feature.
  Returntype : Bio::EnsEMBL::DensityFeature
  Exceptions : throw if invalid density value type is provided
  Caller     : general
  Status     : Stable

=cut

sub new {
  my $caller = shift;

  #allow constructor to be called as class or object method
  my $class = ref($caller) || $caller;

  my $self = $class->SUPER::new(@_);

  my($seq_region, $start, $end, $dt, $dv) =
    rearrange(['SEQ_REGION', 'START', 'END', 'DENSITY_TYPE', 'DENSITY_VALUE'],
              @_);

  throw("Density value must be >= 0.") if($dv < 0);

  if(!defined($dt)){
    throw("Density Type is NOT optional.");
  }

  $self->{'density_type'} = $dt;
  $self->{'density_value'} = $dv;

  $self->{'slice'}    = $seq_region;
  $self->{'start'} = $start;
  $self->{'end'}   = $end;
  

  return $self;
}


=head2 new_fast

  Arg [...]  : none
  Example    : $feature = Bio::EnsEMBL::DensityFeature->new_fast();  
  Description: Creates a new density feature.
  Returntype : Bio::EnsEMBL::DensityFeature
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub new_fast{
  my $caller = shift;

  #allow constructor to be called as class or object method
  my $class = ref($caller) || $caller;

  my $self = $class->SUPER::new(@_);

  return $self;
}


=head2 strand

  Arg [1]    : none
  Example    : $strand = $df->strand();
  Description: Getter fot the strand attribute. Density features always have
               strand 0 and this attribute is not settable.
  Returntype : int (always 0)
  Exceptions : warning if an attempt is made to set the strand
  Caller     : general
  Status     : Stable

=cut

sub strand {
  my $self = shift;
  warning("DensityFeature strand is not settable") if(@_);
  return 0;
}



=head2 density_value

  Arg [1]    : (optional) float $density_value
  Example    : $dv = $density_feature->density_value();
  Description: Getter/Setter for the density value of this DensityFeature.
               The density value may be a count, a percentage, or a coverage
               of a feature type in the area defined by this feature.
  Returntype : float
  Exceptions : throw if a negative density value is provided
  Caller     : general
  Status     : Stable

=cut

sub density_value {
  my $self = shift;

  if(@_) {
    my $density_value = shift;
    throw("Density value must be >= 0.") if($density_value < 0);
    $self->{'density_value'} = $density_value;
  }

  return $self->{'density_value'};
}



=head2 analysis

  Arg [1]    : (optional) Bio::EnsEMBL::Analysis $analysis
               New value for the analysis of the attached DensityType
  Example    : print $df->analysis()->logic_name();
  Description: Overridden superclass analysis method, to chain to analysis
               method on attached DensityType.
  Returntype : Bio::EnsEMBL::Analysis
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub analysis {
  my $self = shift;

  my $dt = $self->density_type();

  return undef if(!$dt);

  return $dt->analysis(@_);
}



=head2 density_type

  Arg [1]    : string $newval (optional) 
               The new value to set the density_value_type attribute to
  Example    : $density_value_type = $obj->density_value_type()
  Description: Getter/Setter for the density_value_type attribute
  Returntype : Bio::EnsEMBL::DensityType
  Exceptions : if object passed is not of type DensityType
  Caller     : general
  Status     : Stable

=cut

sub density_type{
  my $self = shift;
  if(@_) {
    my $type = shift;
    if( !ref $type || !$type->isa("Bio::EnsEMBL::DensityType")){
      throw("object passed must be an ensembl DensityType ". 
	    "not a [".ref($type)."]");
    }
    else{
      $self->{'density_type'}=$type;
    }
  }
  return $self->{'density_type'};
}


###BG########

=head2 scaledvalue

 Title   : scaledvalue
 Usage   : $obj->scaledvalue($newval)
 Function:
 Returns : scalar - object's scaled value
 Args    : newvalue (optional)
 Status  : Stable

=cut

sub scaledvalue{
  my $obj = shift;
  if( @_ ) {
    my $scaledvalue = shift;
    $obj->{'scaledvalue'} = $scaledvalue;
  }
  return $obj->{'scaledvalue'};
}



=head2 url

 Title   : url
 Usage   : $obj->url($newval)
 Function:
 Returns : String containing this object's url
 Args    : newvalue (optional)
 Status  : Stable


=cut

sub url{
   my $obj = shift;
   if( @_ ) {
      my $url = shift;
      $obj->{'url'} = $url;
    }
    return $obj->{'url'};

}


1;



