#
# Ensembl module for Bio::EnsEMBL::DensityType
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::DensityType - A type representing a density, or precentage
coverage etc. in a given region.

=head1 SYNOPSIS

use Bio::EnsEMBL::DensityType;

$type = Bio::EnsEMBL::DensityType->new(-analysis => $analysis,
				       -blocksize => 1000000,
				       -vlaue_type => $type);

=head1 DESCRIPTION

A density type represents a density, or precentage
coverage etc. in a given region.

				       
This module is part of the Ensembl project http://www.ensembl.org

=head1 CONTACT

Post comments/questions to the ensembl development list: ensembl-dev@ebi.ac.uk

=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::DensityType;
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
               valuetype must be 'sum' or ratio',
               valid analysis object must be passed
  Caller     : general

=cut


sub new {
  my $caller = shift;

  my $class = ref($caller) || $caller;

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

  return bless({'analysis'   => $analysis,
		'block_size' => $block_size,
		'value_type' => $value_type});
}

=head2 analysis

  Arg [1]    : Bio::EnsEMBL::Analysi
  Example    : none
  Description: get/set for attribute analysis
  Returntype : Bio::EnsEMBL::Analysis
  Exceptions : none
  Caller     : general

=cut
  
sub analysis{
  my ($self,$arg) = @_;
  
  if (defined($arg)) {
    $self->{analysis} = $arg;
  }


  return $self->{analysis};
}

=head2 value_type

  Arg [1]    : string $analysis
  Example    : none
  Description: get/set for attribute analysis
  Returntype : Bio::EnsE
  Exceptions : none
  Caller     : general

=cut

sub value_type{
  my ($self,$arg) = @_;
  
  if (defined($arg)) {
    $self->{'value_type'} = $arg;
  }


  return $self->{'value_type'};
}

=head2 block_size

  Arg [1]    : string
  Example    : none
  Description: get/set for attribute block_size
  Returntype : string
  Exceptions : none
  Caller     : general

=cut

sub block_size{
  my ($self,$arg) = @_;
  
  if (defined($arg)) {
    $self->{'block_size'} = $arg;
  }


  return $self->{'block_size'};
}
  
