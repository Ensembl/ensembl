#
# Ensembl module for Bio::EnsEMBL::DensityType
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::DensityType - A type representing a density, or percentage
coverage etc. in a given region.

=head1 SYNOPSIS

use Bio::EnsEMBL::DensityType;

$type = Bio::EnsEMBL::DensityType->new(-analysis => $analysis,
				       -blocksize => 1000000,
				       -value_type => $type);

=head1 DESCRIPTION

A density type represents a density, or percentage
coverage etc. in a given region.

				       
This module is part of the Ensembl project http://www.ensembl.org

=head1 CONTACT

Post comments/questions to the ensembl development list: ensembl-dev@ebi.ac.uk

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

  my ($analysis, $block_size, $value_type) = 
    rearrange(['ANALYSIS','BLOCK_SIZE','VALUE_TYPE'],@_);

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

  if($block_size <=0 ){
    throw('-BLOCK_SIZE must be greater than 0'); 
  }

  $self->{'analysis'} = $analysis;
  $self->{'block_size'} = $block_size;
  $self->{'value_type'} = $value_type;

  return $self;
}

=head2 analysis

  Arg [1]    : Bio::EnsEMBL::Analysi
  Example    : none
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
  Example    : none
  Description: gettter/setter for the 
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
  Example    : none
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


1;
