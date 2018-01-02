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

Bio::EnsEMBL::AlignStrainSlice - Represents the slice of the genome aligned with certain strains (applying the variations/indels)

=head1 SYNOPSIS

  $sa = $db->get_SliceAdaptor;

  $slice =
    $sa->fetch_by_region( 'chromosome', 'X', 1_000_000, 2_000_000 );

  $strainSlice1 = $slice->get_by_Strain($strain_name1);
  $strainSlice2 = $slice->get_by_Strain($strain_name2);

  my @strainSlices;
  push @strainSlices, $strainSlice1;
  push @strainSlices, $strainSlice2;

  $alignSlice = Bio::EnsEMBL::AlignStrainSlice->new(
    -SLICE   => $slice,
    -STRAINS => \@strainSlices
  );

  # Get coordinates of variation in alignSlice
  my $alleleFeatures = $strainSlice1->get_all_AlleleFeature_Slice();

  foreach my $af ( @{$alleleFeatures} ) {
    my $new_feature = $alignSlice->alignFeature( $af, $strainSlice1 );
    print( "Coordinates of the feature in AlignSlice are: ",
      $new_feature->start, "-", $new_feature->end, "\n" );
  }

=head1 DESCRIPTION

A AlignStrainSlice object represents a region of a genome align for
certain strains.  It can be used to align certain strains to a reference
slice.

=head1 METHODS

=cut

package Bio::EnsEMBL::AlignStrainSlice;
use strict;

use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Mapper;
use Bio::EnsEMBL::Mapper::RangeRegistry;
use Bio::EnsEMBL::Utils::Exception qw(throw deprecate warning);

=head2 new

    Arg[1]      : Bio::EnsEMBL::Slice $Slice
    Arg[2]      : listref of Bio::EnsEMBL::StrainSlice $strainSlice
    Example     : push @strainSlices, $strainSlice1;
                  push @strainSlices, $strainSlice2;
                  .....
                  push @strainSlices, $strainSliceN;
                  $alignStrainSlice = Bio::EnsEMBL::AlignStrainSlice->new(-SLICE => $slice,
									  -STRAIN => \@strainSlices);
    Description : Creates a new Bio::EnsEMBL::AlignStrainSlice object that will contain a mapper between
                  the Slice object, plus all the indels from the different Strains
    ReturnType  : Bio::EnsEMBL::AlignStrainSlice
    Exceptions  : none
    Caller      : general

=cut

sub new{
    my $caller = shift;
    my $class = ref($caller) || $caller;

    my ($slice, $strainSlices) = rearrange([qw(SLICE STRAINS)],@_);

    #check that both StrainSlice and Slice are identical (must have been defined in the same slice)
    foreach my $strainSlice (@{$strainSlices}){
	if (($strainSlice->start != $slice->start) || ($strainSlice->end != $slice->end) || ($strainSlice->seq_region_name ne $slice->seq_region_name)){
	    warning("Not possible to create Align object from different Slices");
	    return [];
	}
    }

    return bless{'slice' => $slice,
		 'strains' => $strainSlices}, $class;
}

=head2 alignFeature

    Arg[1]      : Bio::EnsEMBL::Feature $feature
    Arg[2]      : Bio::EnsEMBL::StrainSlice $strainSlice
    Example     : $new_feature = $alignSlice->alignFeature($feature, $strainSlice);
    Description : Creates a new Bio::EnsEMBL::Feature object that aligned to 
                  the AlignStrainSlice object.
    ReturnType  : Bio::EnsEMBL::Feature
    Exceptions  : none
    Caller      : general

=cut

sub alignFeature{
    my $self = shift;
    my $feature = shift;

    #check that the object is a Feature
    if (!ref($feature) || !$feature->isa('Bio::EnsEMBL::Feature')){	
	throw("Bio::EnsEMBL::Feature object expected");
    }
    #and align it to the AlignStrainSlice object
    my $mapper_strain = $self->mapper();

    my @results;
  
    if ($feature->start > $feature->end){
	#this is an Indel, map it with the special method
	@results = $mapper_strain->map_indel('Slice',$feature->start, $feature->end, $feature->strand,'Slice');
	#and modify the coordinates according to the length of the indel
	$results[0]->end($results[0]->start + $feature->length_diff -1);
    }
    else{
	@results = $mapper_strain->map_coordinates('Slice',$feature->start, $feature->end, $feature->strand,'Slice');
     }
    #get need start and end of the new feature, aligned ot AlignStrainSlice
    my @results_ordered = sort {$a->start <=> $b->start} @results;

    my %new_feature = %$feature; #make a shallow copy of the Feature
    $new_feature{'start'}= $results_ordered[0]->start();
    $new_feature{'end'} = $results_ordered[-1]->end();  #get last element of the array, the end of the slice

    return bless \%new_feature, ref($feature);
    
}


#getter for the mapper between the Slice and the different StrainSlice objects
sub mapper{
    my $self = shift;
    
    if (!defined $self->{'mapper'}){
	#get the alleleFeatures in all the strains
	if (!defined $self->{'indels'}){
	    #when the list of indels is not defined, get them
	    $self->{'indels'} = $self->_get_indels();
	}
	my $indels = $self->{'indels'}; #gaps in reference slice
	my $mapper = Bio::EnsEMBL::Mapper->new('Slice', 'AlignStrainSlice');
	my $start_slice = 1;
	my $end_slice;
	my $start_align = 1;
	my $end_align;
	my $length_indel = 0;
	my $length_acum_indel = 0;
	foreach my $indel (@{$indels}){
	    $end_slice = $indel->[0] - 1;
	    $end_align = $indel->[0] - 1 + $length_acum_indel; #we must consider length previous indels

	    $length_indel = $indel->[1] - $indel->[0] + 1;
	    

	    $mapper->add_map_coordinates('Slice',$start_slice,$end_slice,1,'AlignStrainSlice',$start_align,$end_align);
	    
	    $mapper->add_indel_coordinates('Slice',$end_slice + 1,$end_slice,1,'AlignStrainSlice',$end_align + 1,$end_align + $length_indel);
	    $start_slice = $end_slice + 1;
	    $start_align = $indel->[1] + 1 + $length_acum_indel; #we must consider legnth previous indels
	    
	    $length_acum_indel += $length_indel;
	}
	if ($start_slice <= $self->length){
	    $mapper->add_map_coordinates('Slice',$start_slice,$self->length,1,'AlignStrainSlice',$start_align,$start_align + $self->length - $start_slice)
	}
	$self->{'mapper'} = $mapper;
	
    }
    return $self->{'mapper'};
}

#returns the length of the AlignSlice: length of the Slice plus the gaps
sub length{
    my $self = shift;
    my $length;
    if (!defined $self->{'indels'}){
	#when the list of indels is not defined, get them
	$self->{'indels'} = $self->_get_indels();	
    }
    $length = $self->{'slice'}->length;
    map {$length += ($_->[1] - $_->[0] + 1)} @{$self->{'indels'}};
    return $length;
}

=head2 strains

  Args       : None
  Description: Returns list with all strains used to
               define this AlignStrainSlice object
  Returntype : listref of Bio::EnsEMBL::StrainSlice objects
  Exceptions : none
  Caller     : general

=cut

sub strains{
    my $self = shift;

    return $self->{'strains'};
}

=head2 Slice

  Args       : None
  Description: Returns slice where the AlignStrainSlice
               is defined
  Returntype : Bio::EnsEMBL::Slice object
  Exceptions : none
  Caller     : general

=cut

sub Slice{
    my $self = shift;
    return $self->{'slice'};
}
#method to retrieve, in order, a list with all the indels in the different strains
sub _get_indels{
    my $self = shift;
    
    #go throuh all the strains getting ONLY the indels (length_diff <> 0)
    my @indels;
    foreach my $strainSlice (@{$self->strains}){
	my $differences = $strainSlice->get_all_AlleleFeatures_Slice(); #need to check there are differences....
	foreach my $af (@{$differences}){
	    #if length is 0, but is a -, it is still a gap in the strain
	    if (($af->length_diff != 0) || ($af->length_diff == 0 && $af->allele_string =~ /-/)){
		push @indels, $af;
	    }
	}
    }
    #need to overlap the gaps using the RangeRegistry module
    my $range_registry = Bio::EnsEMBL::Mapper::RangeRegistry->new();
    foreach my $indel (@indels){
	#in the reference and the strain there is a gap
	$range_registry->check_and_register(1,$indel->start,$indel->start) if ($indel->length_diff == 0);
	#deletion in reference slice
	$range_registry->check_and_register(1,$indel->start, $indel->end ) if ($indel->length_diff < 0);
	#insertion in reference slice
	$range_registry->check_and_register(1,$indel->start,$indel->start + $indel->length_diff - 1) if ($indel->length_diff > 0);
    }
    #and return all the gap coordinates....
    return $range_registry->get_ranges(1);
}

=head2 get_all_Slices

  Args       : none
  Description: This Slice is made of several Bio::EnsEMBL::StrainSlices
               sequence. This method returns these StrainSlices (or part of
               them) with the original coordinates 
  Returntype : listref of Bio::EnsEMBL::StrainSlice objects
  Exceptions : end should be at least as big as start
  Caller     : general

=cut

sub get_all_Slices {
  my $self = shift;

  my @strains;
  #add the reference strain
  my $dbVar = $self->Slice->adaptor->db->get_db_adaptor('variation');
  unless($dbVar) {
	warning("Variation database must be attached to core database to " .
		"retrieve variation information" );
	return '';
    }
  my $indAdaptor = $dbVar->get_IndividualAdaptor();
  my $ref_name =  $indAdaptor->get_reference_strain_name;
  my $ref_strain = Bio::EnsEMBL::StrainSlice->new(
					  -START   => $self->Slice->{'start'},
					  -END     => $self->Slice->{'end'},
					  -STRAND  => $self->Slice->{'strand'},
					  -ADAPTOR => $self->Slice->adaptor(),
					  -SEQ    => $self->Slice->{'seq'},
					  -SEQ_REGION_NAME => $self->Slice->{'seq_region_name'},
					  -SEQ_REGION_LENGTH => $self->Slice->{'seq_region_length'},
					  -COORD_SYSTEM    => $self->Slice->{'coord_system'},
					  -STRAIN_NAME     => $ref_name,
					   );
  #this is a fake reference alisce, should not contain any alleleFeature
  undef $ref_strain->{'alleleFeatures'};
  
  push @strains, @{$self->strains};
  my $new_feature;
  my $indel;
  my $aligned_features;
  my $indels = (); #reference to a hash containing indels in the different strains
  #we need to realign all Features in the different Slices and add '-' in the reference Slice
  foreach my $strain (@{$self->strains}){
      foreach my $af (@{$strain->get_all_AlleleFeatures_Slice()}){
	  $new_feature = $self->alignFeature($af); #align feature in AlignSlice coordinates
	  push @{$aligned_features},$new_feature if($new_feature->seq_region_start <= $strain->end); #some features might map outside slice
	  if ($af->start != $af->end){ #an indel, need to add to the reference, and realign in the strain	     
              #make a shallow copy of the indel - clear it first!
          $indel = undef;
	      %{$indel} = %{$new_feature};
	      bless $indel, ref($new_feature);
	      $indel->allele_string('-');
	      push @{$indels},$indel; #and include in the list of potential indels
	  }
      }
      next if (!defined $aligned_features);
      undef $strain->{'alleleFeatures'}; #remove all features before adding new aligned features
      push @{$strain->{'alleleFeatures'}}, @{$aligned_features};
      undef $aligned_features;
  }  
  push @strains, $ref_strain;
  #need to add indels in the different strains, if not present
  if (defined $indels){
      foreach my $strain (@strains){
	  #inlcude the indels in the StrainSlice object
	  push @{$strain->{'alignIndels'}},@{$indels};
      }
  }
  return \@strains;
}

1;
