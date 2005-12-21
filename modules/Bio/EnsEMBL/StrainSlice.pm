#
# Ensembl module for Bio::EnsEMBL::StrainSlice
#
#
# Copyright Team Ensembl
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::StrainSlice - SubClass of the Slice. Represents the slice of the genome for a certain strain (applying the variations)

=head1 SYNOPSIS

   $sa = $db->get_SliceAdaptor;

   $slice = $sa->fetch_by_region('chromosome', 'X', 1_000_000, 2_000_000);

   $strainSlice = $slice->get_by_Strain($strain_name);

   #get the sequence from the Strain Slice
   my $seq = $strainSlice->seq();
   print $seq;




=head1 DESCRIPTION

A StrainSlice object represents a region of a genome for a certain strain.  It can be used to retrieve
sequence or features from a strain.

=head1 CONTACT

This modules is part of the Ensembl project http://www.ensembl.org

Questions can be posted to the ensembl-dev mailing list:
ensembl-dev@ebi.ac.uk

=head1 METHODS

=cut

package Bio::EnsEMBL::StrainSlice;
use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp);
use Bio::EnsEMBL::Slice;
use Bio::EnsEMBL::Mapper;
use Bio::EnsEMBL::Utils::Exception qw(throw deprecate warning);

use Data::Dumper;
@ISA = qw(Bio::EnsEMBL::Slice);


=head2 new

    Arg[1]      : Bio::EnsEMBL::Slice $slice
    Arg[2]      : string $strain_name
    Example     : $strainSlice = Bio::EnsEMBL::StrainSlice->new(-.... => ,
							      -strain_name => $strain_name);
    Description : Creates a new Bio::EnsEMBL::StrainSlice object that will contain a shallow copy of the
                  Slice object, plus additional information such as the Strain this Slice refers to
                  and listref of Bio::EnsEMBL::Variation::VariationFeatures of differences with the
                  reference sequence
    ReturnType  : Bio::EnsEMBL::StrainSlice
    Exceptions  : none
    Caller      : general

=cut

sub new{
    my $caller = shift;
    my $class = ref($caller) || $caller;

    my ($strain_name) = rearrange(['STRAIN_NAME'],@_);

    my $self = $class->SUPER::new(@_);

    $self->{'strain_name'} = $strain_name;

    if(!$self->adaptor()) {
	warning('Cannot get new StrainSlice features without attached adaptor');
	return '';
    }
    my $variation_db = $self->adaptor->db->get_db_adaptor('variation');

    unless($variation_db) {
	warning("Variation database must be attached to core database to " .
		"retrieve variation information" );
	return '';
    }
    
    my $vf_adaptor = $variation_db->get_VariationFeatureAdaptor;
    
    if( $vf_adaptor ) {
	#get the Population for the given strain
	my $pop_adaptor = $variation_db->get_PopulationAdaptor;

	if ($pop_adaptor){
	    my $population = $pop_adaptor->fetch_by_name($self->{'strain_name'});
	    #check that the population returned is a strain
	    if ((defined $population) && ($population->is_strain)){
		my $variation_features = $vf_adaptor->fetch_all_by_Slice_Population($self,$population);
		$self->{'variationFeatures'} = $variation_features;
		return $self;
	    }
	    else{ 
		warning("Strain not in the database");
		return '';
	    }
	}
	else{
	    warning("Not possible to retrieve PopulationAdaptor from the variation database");
	    return '';
	}
    } else {
	warning("Not possible to retrieve VariationFeatureAdaptor from variation database");
	return '';
    }
}

=head2 seq

  Arg [1]    : none
  Example    : print "SEQUENCE = ", $strainSlice->seq();
  Description: Returns the sequence of the region represented by this
               slice formatted as a string in the strain.
  Returntype : string
  Exceptions : none
  Caller     : general

=cut

sub seq {
  my $self = shift;

  # special case for in-between (insert) coordinates
  return '' if($self->start() == $self->end() + 1);

  return $self->{'seq'} if($self->{'seq'});

  if($self->adaptor()) {
    my $seqAdaptor = $self->adaptor()->db()->get_SequenceAdaptor();
    my $reference_sequence = $seqAdaptor->fetch_by_Slice_start_end_strand($self,1,undef,1); #get the reference sequence for that slice
    #apply all differences to the reference sequence

    # sort edits in reverse order to remove complication of
    # adjusting downstream edits
    my @variation_features_ordered = sort {$b->start() <=> $a->start()} @{$self->{'variationFeatures'}} if (defined $self->{'variationFeatures'});

    foreach my $vf (@variation_features_ordered){
	$vf->apply_edit($reference_sequence); #change, in the reference sequence, the vf
    }
    return substr(${$reference_sequence},0,1) if ($self->length == 1); 
    return ${$reference_sequence}; #returns the reference sequence, applying the variationFeatures
  }

  # no attached sequence, and no db, so just return Ns
  return 'N' x $self->length();
}

=head2 get_all_differences_Slice

    Args        : nonre
    Example     : my $differences = $strainSlice->get_all_differences_Slice()
    Description : Gets all differences between the StrainSlice object and the Slice is defined
    ReturnType  : listref of Bio::EnsEMBL::Variation::VariationFeature
    Exceptions  : none
    Caller      : general

=cut

sub get_all_differences_Slice{
    my $self = shift;

    return $self->{'variationFeatures'};
}

=head2 get_all_differences_StrainSlice

    Arg[1]      : Bio::EnsEMBL::StrainSlice $ss
    Example     : my $differences = $strainSlice->get_all_differences_StrainSlice($ss)
    Description : Gets differences between 2 StrainSlice objects
    ReturnType  : listref of Bio::EnsEMBL::Variation::VariationFeature
    Exceptions  : thrown on bad argument
    Caller      : general

=cut

sub get_all_differences_StrainSlice{
    my $self = shift;
    my $strainSlice = shift;

    if (!ref($strainSlice) || !$strainSlice->isa('Bio::EnsEMBL::StrainSlice')){
	throw('Bio::EnsEMBL::StrainSlice arg expected');
    }
    if ( @{$self->{'variationFeatures'}} == 0 && @{$strainSlice->{'variationFeatures'}} == 0){
	return undef; #there are no differences in any of the Strains
	
    }
    my $differences; #differences between strains
    if (@{$strainSlice->{'variationFeatures'}} == 0){
	#need to create a copy of VariationFeature
	foreach my $difference (@{$self->{'variationFeatures'}}){
	    my %vf = %$difference;
	    push @{$differences},bless \%vf,ref($difference);
	}
    }
    elsif (@{$self->{'variationFeatures'}} == 0){
	#need to create a copy of VariationFeature, but changing the allele by the allele in the reference
	foreach my $difference (@{$strainSlice->{'variationFeatures'}}){
	    push @{$differences}, $strainSlice->_convert_difference($difference);
	}
    }
    else{
	#both strains have differences
	#create a hash with the differences in the self strain slice
	my %variation_features_self = map {$_->start => $_} @{$self->{'variationFeatures'}};
	foreach my $difference (@{$strainSlice->{'variationFeatures'}}){
	    #there is no difference in the other strain slice, convert the allele
	    if (!defined $variation_features_self{$difference->start}){
		push @{$differences},$strainSlice->_convert_difference($difference);
	    }		
	    else{
		#if it is defined and have the same allele, delete from the hash
		if ($variation_features_self{$difference->start}->allele_string eq $difference->allele_string){
		    delete $variation_features_self{$difference->start};
		}
	    }
	}	
	#and copy the differences that in the self
	foreach my $difference (values %variation_features_self){
	    my %vf = %$difference;
	    push @{$differences},bless \%vf,ref($difference);
	}

    }
    #need to map differences to the self 
    my $mapper = $self->mapper(); #now that we have the differences, map them in the StrainSlice
#    print Dumper($mapper);
    my @results;
    foreach my $difference (@{$differences}){	
	@results = $mapper->map_coordinates('Slice',$difference->start,$difference->end,$difference->strand,'Slice');
	#we can have 3 possibilities:
	#the difference is an insertion and when mapping returns the boundaries of the insertion in the StrainSlice
	if (@results == 2){
	    #the first position in the result is the beginning of the insertion
	    if($results[0]->start < $results[1]->start){
		$difference->start($results[0]->end+1);
		$difference->end($results[1]->start-1);
	    }
	    else{
		$difference->start($results[1]->end+1);
		$difference->end($results[0]->start-1);
	    }
	    $difference->strand($results[0]->strand);
	}
	else{
	    #it can be either a SNP or a deletion, and we have the coordinates in the result, etither a Bio::EnsEMBL::Mapper::Coordinate
	    # or a Bio::EnsEMBL::Mapper::IndelCoordinate
#	    print "Difference: ",$difference->start, "-", $difference->end,"strand ",$difference->strand,"\n";
	    $difference->start($results[0]->start);
	    $difference->end($results[0]->end);
	    $difference->strand($results[0]->strand);
	}
    }

    return $differences;
}

#for a given VariationFeatures, converts the allele into the reference allele and returns a new list with
#the converted VariationFeatures
sub _convert_difference{
    my $self = shift;
    my $difference = shift;
    my %new_vf = %$difference; #make a copy of the variationFeature
    #and change the allele with the one from the reference Slice
    $new_vf{'allele_string'} = $self->SUPER::subseq($difference->start,$difference->end,$difference->strand);
    return bless \%new_vf,ref($difference);
}

=head2 sub_Slice

  Arg   1    : int $start
  Arg   2    : int $end
  Arge [3]   : int $strand
  Example    : none
  Description: Makes another StrainSlice that covers only part of this slice
               with the appropriate differences to the reference Slice
               If a slice is requested which lies outside of the boundaries
               of this function will return undef.  This means that
               behaviour will be consistant whether or not the slice is
               attached to the database (i.e. if there is attached sequence
               to the slice).  Alternatively the expand() method or the
               SliceAdaptor::fetch_by_region method can be used instead.
  Returntype : Bio::EnsEMBL::StrainSlice or undef if arguments are wrong
  Exceptions : thrown when trying to get the subSlice in the middle of a
               insertion
  Caller     : general

=cut

sub sub_Slice {
  my ( $self, $start, $end, $strand ) = @_;
  my $mapper = $self->mapper();
  #finally map from the Slice to the Strain
  my @results = $mapper->map_coordinates('StrainSlice',$start,$end,$strand,'StrainSlice');
  my $new_start;
  my $new_end;
  my $new_strand;
  my $new_seq;
 
  #Get need start and end for the subSlice of the StrainSlice
  my @results_ordered = sort {$a->start <=> $b->start} @results;
  $new_start = $results_ordered[0]->start();
  $new_strand = $results_ordered[0]->strand() if (ref($results_ordered[0]) eq 'Bio::EnsEMBL::Mapper::Coordinate');
  $new_strand = $results_ordered[-1]->strand() if (ref($results_ordered[-1]) eq 'Bio::EnsEMBL::Mapper::Coordinate');
  $new_end = $results_ordered[-1]->end();  #get last element of the array, the end of the slice

  my $subSlice = $self->SUPER::sub_Slice($new_start,$new_end,$new_strand);
  $subSlice->{'strain_name'} = $self->{'strain_name'};

  my $new_variations; #reference to an array that will contain the variationFeatures in the new subSlice
  #update the VariationFeatures in the sub_Slice of the Strain
  my $vf_start;
  my $vf_end;
  my $offset = $subSlice->start - $self->start;

  foreach my $variationFeature (@{$self->{'variationFeatures'}}){
      #calculate the new position of the variation_feature in the subSlice
      $vf_start = $variationFeature->start - $offset;
      $vf_end = $variationFeature->end - $offset;
      if ($vf_start  >= 1 and $vf_end <= $subSlice->length){
	  #copy the variationFeature
	  my %new_vf;
	  %new_vf = %$variationFeature;
	  #and shift to the new coordinates
	  $new_vf{'start'} = $vf_start;
	  $new_vf{'end'} = $vf_end;	
	  my $test = bless \%new_vf, ref($variationFeature);
	  push @{$new_variations}, $test;
      }
  }
  $subSlice->{'variationFeatures'} = $new_variations;
  return $subSlice;

}

=head2 subseq

  Arg  [1]   : int $startBasePair
               relative to start of slice, which is 1.
  Arg  [2]   : int $endBasePair
               relative to start of slice.
  Arg  [3]   : (optional) int $strand
               The strand of the slice to obtain sequence from. Default
               value is 1.
  Description: returns string of dna sequence
  Returntype : txt
  Exceptions : end should be at least as big as start
               strand must be set
  Caller     : general

=cut

sub subseq {
  my ( $self, $start, $end, $strand ) = @_;

  if ( $end+1 < $start ) {
    throw("End coord + 1 is less then start coord");
  }

  # handle 'between' case for insertions
  return '' if( $start == $end + 1);

  $strand = 1 unless(defined $strand);

  if ( $strand != -1 && $strand != 1 ) {
    throw("Invalid strand [$strand] in call to Slice::subseq.");
  }

  my $subseq;
  my $seq;
  if($self->adaptor){
      $seq = $self->seq;
      reverse_comp(\$seq) if ($strand == -1);
      $subseq = substr($seq,$start-1,$end);
  } 
  else {
      ## check for gap at the beginning and pad it with Ns
      if ($start < 1) {
	  $subseq = "N" x (1 - $start);
	  $start = 1;
      }
      $subseq .= substr ($self->seq(), $start-1, $end - $start + 1);
      ## check for gap at the end and pad it with Ns
    if ($end > $self->length()) {
	$subseq .= "N" x ($end - $self->length());
    }
      reverse_comp(\$subseq) if($strand == -1);
  }
  return $subseq;
  
}


sub mapper{
    my $self = shift;
   
    if (@_) {
	delete $self->{'mapper'};
    }
    if(!defined $self->{'mapper'}){
	#create the mapper between the Slice and StrainSlice
	my $mapper = Bio::EnsEMBL::Mapper->new('Slice','StrainSlice');
	#align with Slice
	#get all the VariationFeatures in the strain Slice, from start to end in the Slice
	my @variation_features_ordered = sort {$a->start() <=> $b->start()} @{$self->{'variationFeatures'}} if (defined $self->{'variationFeatures'});
	
	my $start_slice = 1;
	my $end_slice;
	my $start_strain = 1;
	my $end_strain;
	my $length_allele;
	foreach my $variation_feature (@variation_features_ordered){
	    #we have a insertion/deletion: marks the beginning of new slice move coordinates	    
	    if ($variation_feature->length_diff != 0){
		$length_allele = $variation_feature->length + $variation_feature->length_diff();	    
		$end_slice = $variation_feature->start() - 1;

		if ($end_slice >= $start_slice){
		    $end_strain = $end_slice - $start_slice + $start_strain;
		    #add the sequence that maps
		    $mapper->add_map_coordinates('Slice',$start_slice,$end_slice,1,'StrainSlice',$start_strain,$end_strain);
		    #add the indel
		    $mapper->add_indel_coordinates('Slice',$end_slice+1,$end_slice + $variation_feature->length,1,'StrainSlice',$end_strain+1,$end_strain + $length_allele);
		    $start_strain = $end_strain + $length_allele + 1;
		}
		else{
		    #add the indel
		    $mapper->add_indel_coordinates('Slice',$end_slice+1,$end_slice + $variation_feature->length,1,'StrainSlice',$end_strain+1,$end_strain + $length_allele);
		    $start_strain += $length_allele;
		}
		$start_slice = $end_slice + $variation_feature->length+ 1;
	    }
	}
	if ($start_slice <= $self->length){
	    $mapper->add_map_coordinates('Slice',$start_slice,$self->length,1,'StrainSlice',$start_strain,$start_strain + $self->length - $start_slice);
	}
	$self->{'mapper'} = $mapper;
    }
    return $self->{'mapper'};
}
1;
