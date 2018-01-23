=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=cut


=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=cut

=head1 NAME

Bio::EnsEMBL::DensityFeatureSet -
A feature representing a set of density features

=head1 SYNOPSIS

  use Bio::EnsEMBL::DensityFeatureSet;

  my $densitySet = Bio::EnsEMBL::DensityFeatureSet->new(
    -bin_array = \@out,
    -stretch   = 1,
  );

=head1 DESCRIPTION

A density feature set is a wrap around a array of density features with
additional information about the collective density feature set, such as
max_min_values and scale factors etc. a given region.

This module is part of the Ensembl project http://www.ensembl.org

=head1 METHODS

=cut


package Bio::EnsEMBL::DensityFeatureSet;

use strict;
use warnings;

use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Utils::Exception qw(throw);

=head2 new

  Description: Creates a new density feature set.
  Returntype : Bio::EnsEMBL::DensityFeatureSet
  Exceptions : throw if invalid density value type is provided
  Caller     : general
  Status     : Stable

=cut

sub new {
  my $class = shift;

  my $max_value = undef;
  my $min_value = undef;

  my($dfeats, $stretch, $scale_to_fit) =
      rearrange(['FEATURES', 'STRETCH', 'SCALE_TO_FIT'], @_);
  foreach (@$dfeats){
	  my $value = $_->density_value;
	  $max_value = $value if (!defined($max_value) || $value > $max_value); 
	  $min_value = $value if (!defined($min_value) || $value < $min_value);
  }

  return bless {'bin_array'    => $dfeats,
  	            'stretch'      => $stretch,
                'scale_to_fit' => $scale_to_fit,
                'min_value'    => $min_value,
                'max_value'    => $max_value}, $class;
}


=head2 stretch

 Title   : stretch
 Usage   : $obj->stretch($newval)
 Function: gets/sets a boolean for whether we should stretch the data over the
 range (i.e. from min to max rather than absolute numbers).
 Returns : value of _stretch
 Args    : newvalue (optional)
 Status  : Stable

=cut

sub stretch{
   my $self = shift;
   $self->{'stretch'} = shift if(@_);
   return $self->{'stretch'};
}


=head2 scale_to_fit

 Title   : scale_to_fit
 Usage   : $obj->scale_to_fit($newval)
 Function: gets/sets the number that the BinValues are to be scaled against -
 i.e. the greatest BinValue->value will be scaled to this number, and the rest
 scaled in proportion.
 Returns : scale_to_fit value
 Args    : newvalue (optional)
 Status  : Stable


=cut

sub scale_to_fit{
   my $self = shift;
   $self->{'scale_to_fit'} = shift if (@_);
   return $self->{'scale_to_fit'};

}

=head2 colour

 Title   : colour
 Usage   : $obj->colour($newval)
 Function: 
 Returns : value of colour
 Args    : newvalue (optional)
 Status  : Stable


=cut


sub colour{
   my $self = shift;
   $self->{'color'} = shift if(@_);
   return $self->{'color'};

}

=head2 label

 Title   : label
 Usage   : $obj->label($newval)
 Function: 
 Returns : String containing label
 Args    : newvalue (optional)
 Status  : Stable


=cut

sub label{
   my $self = shift;
   $self->{'label'} = shift if (@_);
   return $self->{'label'};

}


=head2 label2

 Title   : label2
 Usage   : $obj->label2($newval)
 Function: 
 Returns : String containing label2
 Args    : newvalue (optional)
 Status  : Stable


=cut

sub label2{
    my $self = shift;
   $self->{'label2'} = shift if (@_);
   return $self->{'label2'};
}



=head2 get_all_binvalues

  Arg [1]    : none
  Example    : @binvalues = @{$dfs->get_all_binvalues};
  Description: Scales all of the contained DensityFeatures by $scalefactor
               and returns them.
  Returntype : reference to a list of DensityFeatures
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_all_binvalues{
  my $self = shift;
  my $max_value = $self->max_value();
  my $min_value = $self->min_value();

 	return [] if(!@{$self->{'bin_array'}});

  my $width = $self->scale_to_fit();
  return [] unless defined($width);
  # throw("Cannot scale values - scale_to_fit has not been set");

  if ($self->stretch && ($max_value-$min_value) ){
    foreach my $bv (@{ $self->{'bin_array'}}){
      my $scaledval = (($bv->density_value - $min_value) /
                       ($max_value-$min_value) )* $width;
      $bv->scaledvalue($scaledval);
    }
  } elsif($max_value) {
    foreach my $bv (@{ $self->{'bin_array'}}){
      my $scaledval = ($bv->density_value / $max_value) * $width;
      $bv->scaledvalue($scaledval);
    }
  } else {
    foreach my $bv (@{ $self->{'bin_array'}}){
      $bv->scaledvalue(0);
    }
  }

  return $self->{'bin_array'};
}


=head2 max_value

  Arg [1]    : none
  Example    : my $max = $dfs->max_value();
  Description: Returns the maximum density feature value from the density
               feature set
  Returntype : int
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub max_value{ $_[0]->{'max_value'};}


=head2 min_value

  Arg [1]    : none
  Example    : my $min = $dfs->min_value();
  Description: Returns the minimum density feature value from the density
               feature set.
  Returntype : int
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub min_value{ $_[0]->{'min_value'};}



=head2 size

  Arg [1]    : none
  Example    : my $num_features = $dfs->size();
  Description: Returns the number of density features in this density feature
               set.
  Returntype : int
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub size {
    my $self = shift;
    return scalar @{$self->{'bin_array'}};
}

1;



