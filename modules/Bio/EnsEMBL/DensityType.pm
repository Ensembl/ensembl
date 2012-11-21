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

Bio::EnsEMBL::DensityType - A type representing a density, or percentage
coverage etc. in a given region.

=head1 SYNOPSIS

  use Bio::EnsEMBL::DensityType;

  $type = Bio::EnsEMBL::DensityType->new(
    -analysis   => $analysis,
    -blocksize  => 1000000,
    -value_type => $type
  );

=head1 DESCRIPTION

A density type represents a density, or percentage coverage etc. in a
given region.

=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::DensityType;

use Bio::EnsEMBL::Storable;
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Utils::Exception qw(throw);

use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Storable);

=head2 new

  Arg [..]   :  Takes a set of named arguments
  Example    : $dt = new Bio::EnsEMBL::DensityType::DensityType(
				-analysis  =>  $analysis,
				-blocksize =>  1e6,
				-value_type => 'sum')

  Description: Creates a new Density Type object
  Returntype : Bio::EnsEMBL::DensityType
  Exceptions : blocksize > 0, 
               valuetype must be 'sum' or 'ratio',
               valid analysis object must be passed
  Caller     : general
  Status     : Stable

=cut


sub new {
  my $class = shift;

  my $self = $class->SUPER::new(@_);

  my ($analysis, $block_size, $value_type, $region_features) = 
    rearrange(['ANALYSIS','BLOCK_SIZE','VALUE_TYPE','REGION_FEATURES'],@_);

  if($analysis) {
    if(!ref($analysis) || !$analysis->isa('Bio::EnsEMBL::Analysis')) {
      throw('-ANALYSIS argument must be a Bio::EnsEMBL::Analysis not '.
            $analysis);
    }
  }

  if($value_type ne "sum" and $value_type ne "ratio"){
    throw('-VALUE_TYPE argument must be "ratio" or "sum" not *'.
	  $value_type."*");
  }

  $block_size ||= 0;
  $region_features ||= 0;

  if($block_size && $region_features){
    throw('Set either -BLOCK_SIZE or -REGION_FEATURES, not both'); 
  }

  if( $block_size <0 or $region_features < 0 ) {
    throw( 'No negative values for -BLOCK_SIZE or -REGION_FEATURES' );
  }


  $self->{'analysis'} = $analysis;
  $self->{'block_size'} = $block_size;
  $self->{'value_type'} = $value_type;
  $self->{'region_features'} = $region_features;

  return $self;
}

=head2 analysis

  Arg [1]    : Bio::EnsEMBL::Analysis
  Description: get/set for attribute analysis
  Returntype : Bio::EnsEMBL::Analysis
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub analysis{
  my $self = shift;

  if(@_) {
    my $a = shift;
    if(defined($a) && (!ref($a) || !$a->isa('Bio::EnsEMBL::Analysis'))) {
      throw("Argument must be undef or a Bio::EnsEMBL::Analysis object.");
    }
    $self->{'analysis'} = $a;
  }
  return $self->{'analysis'};
}

=head2 value_type

  Arg [1]    : string $value_type
  Description: gettter/setter for the type 
  Returntype : float
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub value_type{
  my $self = shift;
  $self->{'value_type'} = shift if(@_);
  return $self->{'value_type'};
}


=head2 block_size

  Arg [1]    : int
  Description: getter/setter for attribute block_size
  Returntype : int
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub block_size{
  my $self = shift;
  $self->{'block_size'} = shift if(@_);
  return $self->{'block_size'};
}


=head2 region_features

  Arg [1]    : int $region_features
  Example    : The number of features per seq_region inside this density_type..
  Description: get/set for attribute region_features
  Returntype : string
  Exceptions : none
  Caller     : general

=cut

sub region_features {
   my $self = shift;
  $self->{'region_features'} = shift if( @_ );
  return $self->{'region_features'};
}


1;
