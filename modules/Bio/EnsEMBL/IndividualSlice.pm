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

Bio::EnsEMBL::IndividualSlice - SubClass of the Slice. Represents the
slice of the genome for a certain individual (applying the alleles for
this individual)

=head1 SYNOPSIS

  $sa = $db->get_SliceAdaptor;

  $slice =
    $sa->fetch_by_region( 'chromosome', 'X', 1_000_000, 2_000_000 );

  $individualSlice = $slice->get_by_Individual($individual_name);

  # Get the sequence from the Individual Slice: will contain IUPAC codes
  # for SNPs and Ensembl ambiguity codes for indels
  my $seq = $individualSlice->seq();
  print $seq;

  # Get a subSlice of the Strain
  my $subSlice_individual =
    $individualSlice->sub_Slice( 5_000, 8_000, 1 )

  # Compare two different individuals in the same Slice
  my $sliceIndividual2 = $slice->get_by_Individual($individual_name2);
  my $differences =
    $individualSlice->get_all_differences_IndividualSlice(
    $sliceIndividual2);

  foreach my $af ( @{$differences} ) {
    print
      "There is a difference between $individual_name "
      . "and $individual_name2 at ",
      $af->start, "-", $af->end,
      " with allele ", $af->allele_string(), "\n";
  }

=head1 DESCRIPTION

A IndividualSlice object represents a region of a genome for a certain
individual.  It can be used to retrieve sequence or features from a
individual.

=head1 METHODS

=cut

package Bio::EnsEMBL::IndividualSlice;
use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp);
use Bio::EnsEMBL::Slice;
use Bio::EnsEMBL::Mapper;
use Bio::EnsEMBL::Utils::Exception qw(throw deprecate warning);

@ISA = qw(Bio::EnsEMBL::Slice);


=head2 new

    Arg [1..N]  : List of named arguments
                  Bio::EnsEMBL::CoordSystem COORD_SYSTEM
                  string SEQ_REGION_NAME,
                  int    START,
                  int    END,
                  string VERSION (optional, defaults to '')
                  int    STRAND, (optional, defaults to 1)
                  Bio::EnsEMBL::DBSQL::SliceAdaptor ADAPTOR (optional)
    Arg[N+1]      : string $individual_name
    Example     : $individualSlice = Bio::EnsEMBL::IndividualSlice->new(-coord_system => $cs,
									-start => 1,
									-end => 10000,
									-strand => 1,
									-seq_region_name => 'X',
									-seq_region_length => 12e6,
									-individual_name => $individual_name);
    Description : Creates a new Bio::EnsEMBL::IndividualSlice object that will contain a shallow copy of the
                  Slice object, plus additional information such as the individual this Slice refers to
                  and listref of Bio::EnsEMBL::Variation::AlleleFeatures of differences with the
                  reference sequence
    ReturnType  : Bio::EnsEMBL::IndividualSlice
    Exceptions  : none
    Caller      : general

=cut

sub new{
    my $caller = shift;
    my $class = ref($caller) || $caller;

    #create the IndividualSlice object as the Slice, plus the individual attribute
    my ($individual_name, $sample_id) = rearrange(['INDIVIDUAL', 'SAMPLE_ID'],@_);

    my $self = $class->SUPER::new(@_);

    $self->{'individual_name'} = $individual_name;
    $self->{'sample_id'} = $sample_id;

    return $self;

}

=head2 individual_name

    Arg [1]     : (optional) string $individual_name
    Example     : my $individual_name = $individualSlice->individual_name();
    Description : Getter/Setter for the name of the individual in the slice
    ReturnType  : string
    Exceptions  : none
    Caller      : general

=cut

sub individual_name{
   my $self = shift;
   if (@_){
       $self->{'individual_name'} = shift @_;
   }
   return $self->{'individual_name'};
}

=head2 seq

  Arg [1]    : none
  Example    : print "SEQUENCE = ", $strainSlice->seq();
  Description: Returns the sequence of the region represented by this
               StrainSlice formatted as a string.
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
    my @allele_features_ordered;
    @allele_features_ordered = sort {$b->start() <=> $a->start() || $b->end() <=> $a->end()} @{$self->{'alleleFeatures'}} if (defined $self->{'alleleFeatures'});

    foreach my $af (@allele_features_ordered){
	$af->apply_edit($reference_sequence); #change, in the reference sequence, the af
    }
#    return substr(${$reference_sequence},0,1) if ($self->length == 1); 
    return ${$reference_sequence}; #returns the reference sequence, applying the alleleFeatures
  }

  # no attached sequence, and no db, so just return Ns
  return 'N' x $self->length();
}

=head2 get_all_differences_Slice

    Args        : none
    Example     : my $differences = $individualSlice->get_all_differences_Slice()
    Description : Gets all differences between the IndividualSlice object and the Slice is defined
    ReturnType  : listref of Bio::EnsEMBL::Variation::AlleleFeature
    Exceptions  : none
    Caller      : general

=cut

sub get_all_differences_Slice{
    my $self = shift;
    my $differences; #reference to the array with the differences between Slice and StrainSlice
    my $ref_allele;
    foreach my $difference (@{$self->{'alleleFeatures'}}){
	if ($difference->length_diff == 0){
	    #the difference is a SNP, check if it is the same as the reference allele
	    $ref_allele = $self->SUPER::subseq($difference->start,$difference->end,$difference->strand);
	    $ref_allele = '-' if ($ref_allele eq '');
	    if ($ref_allele ne $difference->allele_string){
		#when the alleleFeature is different from the reference allele, add to the differences list
		push @{$differences},$difference;
	    }
	}
	else{
	    push @{$differences},$difference;
	}
    }

    return $differences;

}

=head2 get_all_differences_IndividualSlice

    Arg[1]      : Bio::EnsEMBL::IndividualSlice $is
    Example     : my $differences = $individualSlice->get_all_differences_IndividualSlice($individualslice)
    Description : Gets differences between 2 IndividualSlice objects
    ReturnType  : listref of Bio::EnsEMBL::Variation::AlleleFeature
    Exceptions  : thrown on bad argument
    Caller      : general

=cut

sub get_all_differences_IndividualSlice{
    my $self = shift;
    my $individualSlice = shift;

    if (!ref($individualSlice) || !$individualSlice->isa('Bio::EnsEMBL::IndividualSlice')){
	throw('Bio::EnsEMBL::IndividualSlice arg expected');
    }
    if ( @{$self->{'alleleFeatures'}} == 0 && @{$individualSlice->{'alleleFeatures'}} == 0){
	return undef; #there are no differences in any of the Individuals
	
    }
    my $differences; #differences between individuals
    if (@{$individualSlice->{'alleleFeatures'}} == 0){
	#need to create a copy of alleleFeature for the first Individual
	foreach my $difference (@{$self->{'alleleFeatures'}}){
	    my %vf = %$difference;
	    push @{$differences},bless \%vf,ref($difference);
	}
    }
    elsif (@{$self->{'alleleFeatures'}} == 0){
	#need to create a copy of AlleleFeature, but changing the allele by the allele in the reference sequence
	foreach my $difference (@{$individualSlice->{'alleleFeatures'}}){
	    push @{$differences}, $individualSlice->_convert_difference($difference);
	}
    }
    else{
	#both individuals have differences
	#create a hash with the differences in the first slice
	my %allele_features_self = map {$_->start.'-'.$_->end => $_} @{$self->{'alleleFeatures'}};
	foreach my $difference (@{$individualSlice->{'alleleFeatures'}}){
	    #there is no difference in the other individual slice, convert the allele
	  if (!defined $allele_features_self{$difference->start.'-'.$difference->end}){
	      push @{$differences},$individualSlice->_convert_difference($difference);
	    }
	    else{
		#if it is defined and have the same allele, delete from the hash since it is not a difference
		#between the individuals
		if ($allele_features_self{$difference->start.'-'.$difference->end}->allele_string eq $difference->allele_string){
		  delete $allele_features_self{$difference->start.'-'.$difference->end};
		}
	    }
	}	
	#and finally, make a shallow copy of the differences in the first individual
	foreach my $difference (values %allele_features_self){
	  my %vf = %$difference;
	  push @{$differences},bless \%vf,ref($difference);
	}
	
    }
    #need to map differences to the first individual, self, since the coordinates are in the Slice coordinate system
    my $mapper = $self->mapper(); #now that we have the differences, map them in the IndividualSlice
    my @results;
    foreach my $difference (@{$differences}){
      @results = $mapper->map_coordinates('Slice',$difference->start,$difference->end,$difference->strand,'Slice');
      #we can have 3 possibilities:
      #the difference is an insertion and when mapping returns the boundaries of the insertion in the IndividualSlice
      if (@results == 2){
	#the first position in the result is the beginning of the insertion
	if($results[0]->start < $results[1]->start){
	  $difference->start($results[0]->end+1);
	  $difference->end($results[1]->start-1);
	}
	else{
	  #it is the second position the beginning of the insertion
	  $difference->start($results[1]->end+1);
	  $difference->end($results[0]->start-1);
	}
	$difference->strand($results[0]->strand);
      }
      else{
	#it can be either a SNP or a deletion, and we have the coordinates in the result, etither a Bio::EnsEMBL::Mapper::Coordinate
	# or a Bio::EnsEMBL::Mapper::IndelCoordinate
	$difference->start($results[0]->start);
	$difference->end($results[0]->end);
	$difference->strand($results[0]->strand);
      }
    }
    
    return $differences;
}

#for a given AlleleFeature, converts the allele into the reference allele and returns
#the converted AlleleFeature

sub _convert_difference{
    my $self = shift;
    my $difference = shift;
    my %new_af = %$difference; #make a copy of the alleleFeature
    #and change the allele with the one from the reference Slice
    $new_af{'allele_string'} = $self->SUPER::subseq($difference->start,$difference->end,$difference->strand);
    return bless \%new_af,ref($difference);
}

=head2 mapper

  Args       : none
  Description: Getter for the mapper between the between the IndividualSlice and the Slice it refers to. 
               It is done automatically when necessary to create subSlice or to get the differences between individuals
  Returntype : Bio::EnsEMBL::Mapper
  Exceptions : none
  Caller     : Internal function

=cut

sub mapper{
    my $self = shift;
   
    if (@_) {
	#allow to create again the mapper
	delete $self->{'mapper'};
    }
    if(!defined $self->{'mapper'}){
	#create the mapper between the Slice and StrainSlice
	my $mapper = Bio::EnsEMBL::Mapper->new('Slice','IndividualSlice');
	#align with Slice
	#get all the VariationFeatures in the Individual Slice, from start to end in the Slice
  my @allele_features_ordered;
	@allele_features_ordered = sort {$a->start() <=> $b->start() || $b->end() <=> $a->end()} @{$self->{'alleleFeatures'}} if (defined $self->{'alleleFeatures'});
	
	my $start_slice = 1;
	my $end_slice;
	my $start_individual = 1;
	my $end_individual;
	my $length_allele;
	my $total_length_diff = 0;
	#we will walk from left to right in the slice object, updating the start and end individual every time
	#there is a new alleleFeature in the Individual
	foreach my $allele_feature (@allele_features_ordered){
	    #we have a insertion/deletion: marks the beginning of new slice move coordinates
	    if ($allele_feature->length_diff != 0){
	        $total_length_diff += $allele_feature->length_diff;
		$length_allele = $allele_feature->length + $allele_feature->length_diff(); #length of the allele in the Individual
		$end_slice = $allele_feature->start() - 1; #set the end of the slice before the alleleFeature
		if ($end_slice >= $start_slice){
		    #normal cases (not with gaps)
		    $end_individual = $end_slice - $start_slice + $start_individual; #set the end of the individual from the beginning plus the offset
		    #add the sequence that maps
		    $mapper->add_map_coordinates('Slice',$start_slice,$end_slice,1,'IndividualSlice',$start_individual,$end_individual);
		    #and add the indel
		    $mapper->add_indel_coordinates('Slice',$end_slice+1,$end_slice + $allele_feature->length,1,'IndividualSlice',$end_individual+1,$end_individual + $length_allele);
		    $start_individual = $end_individual + $length_allele + 1; #set the beginning of the individual after the allele
		}
		else{
		    #add the indel
		    $mapper->add_indel_coordinates('Slice',$end_slice+1,$end_slice + $allele_feature->length,1,'IndividualSlice',$end_individual+1,$end_individual + $length_allele);
		    $start_individual += $length_allele;
		}
		$start_slice = $end_slice + $allele_feature->length+ 1; #set the beginning of the slice after the variation feature
	    }
	}
	if ($start_slice <= $self->length){
	    #if we haven't reached the end of the IndividualSlice, add the final map coordinates between the individual and the slice
	    $mapper->add_map_coordinates('Slice',$start_slice,$self->length,1,'IndividualSlice',$start_individual,$start_individual + $self->length - $start_slice);
	}

	$mapper->add_map_coordinates('Slice', -$self->start+1, 0,1, 'IndividualSlice', -$self->start +1,0) if ($self->start > 0); #before individualSlice
	$mapper->add_map_coordinates('Slice', $self->length + 1,$self->seq_region_length - ($self->length +1),1, 'IndividualSlice', $self->length + 1 + $total_length_diff,$self->seq_region_length + $total_length_diff - ($self->length +1) ) if ($self->length <= $self->seq_region_length); #after strainSlice
	$self->{'mapper'} = $mapper;
    }
    return $self->{'mapper'};
}

=head2 sub_Slice

  Arg   1    : int $start
  Arg   2    : int $end
  Arge [3]   : int $strand
  Example    : none
  Description: Makes another IndividualSlice that covers only part of this IndividualSlice
               with the appropriate differences to the reference Slice
               If a slice is requested which lies outside of the boundaries
               of this function will return undef.  This means that
               behaviour will be consistant whether or not the slice is
               attached to the database (i.e. if there is attached sequence
               to the slice).  Alternatively the expand() method or the
               SliceAdaptor::fetch_by_region method can be used instead.
  Returntype : Bio::EnsEMBL::IndividualSlice or undef if arguments are wrong
  Exceptions : thrown when trying to get the subSlice in the middle of a
               insertion
  Caller     : general

=cut

sub sub_Slice {
  my ( $self, $start, $end, $strand ) = @_;
  my $mapper = $self->mapper();
  #map from the Individual to the Slice to get the sub_Slice, and then, apply the differences in the subSlice
  my @results = $mapper->map_coordinates('IndividualSlice',$start,$end,$strand,'IndividualSlice');
  my $new_start;
  my $new_end;
  my $new_strand;
  my $new_seq;
  #Get need start and end for the subSlice of the IndividualSlice
  my @results_ordered = sort {$a->start <=> $b->start} grep {ref($_) eq 'Bio::EnsEMBL::Mapper::Coordinate'} @results;
  $new_start = $results_ordered[0]->start();
  $new_strand = $results_ordered[0]->strand() if (ref($results_ordered[0]) eq 'Bio::EnsEMBL::Mapper::Coordinate');
#  $new_strand = $results_ordered[-1]->strand() if (ref($results_ordered[-1]) eq 'Bio::EnsEMBL::Mapper::Coordinate');
  $new_end = $results_ordered[-1]->end();  #get last element of the array, the end of the slice

  my $subSlice = $self->SUPER::sub_Slice($new_start,$new_end,$new_strand);
  $subSlice->{'individual_name'} = $self->{'individual_name'};

  my $new_alleles; #reference to an array that will contain the variationFeatures in the new subSlice
  #update the VariationFeatures in the sub_Slice of the Individual
  my %af;
  my $new_allele_feature;
  foreach my $alleleFeature (@{$self->{'alleleFeatures'}}){      
      $new_allele_feature = $alleleFeature->transfer($subSlice);
      #only transfer the coordinates to the SubSlice that are within the boundaries
      if ($new_allele_feature->start >= 1 && $new_allele_feature->end <= $subSlice->length){
	  push @{$new_alleles}, $new_allele_feature;
      }
  }
  $subSlice->{'alleleFeatures'} = $new_alleles;
  return $subSlice;

}

=head2 subseq

  Arg  [1]   : int $startBasePair
               relative to start of slice, which is 1.
  Arg  [2]   : int $endBasePair
               relative to start of slice.
  Arg  [3]   : (optional) int $strand
               The strand of the individual slice to obtain sequence from. Default
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
    throw("End coord + 1 is less than start coord");
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
    my $seqAdaptor = $self->adaptor()->db()->get_SequenceAdaptor();
    $subseq = ${$seqAdaptor->fetch_by_Slice_start_end_strand($self,$start,$end,$strand)}; #get the reference sequence for that slice
    #apply all differences to the reference sequence
    # sort edits in reverse order to remove complication of
    # adjusting downstream edits
    my @allele_features_ordered;
    @allele_features_ordered = sort {$b->start() <=> $a->start() || $b->end() <=> $a->end()} @{$self->{'alleleFeatures'}} if (defined $self->{'alleleFeatures'});
    my $af_start;
    my $af_end;
    foreach my $af (@allele_features_ordered){
	if (($af->start - $start +1 > 0) && ($end - $af->end > 0)){
	    #save the current start and end of the alleleFeature before changing for apply_edit
	    $af_start = $af->start;
	    $af_end = $af->end;
	    #apply the difference if the feature is in the new slice
	    $af->start($af->start - $start +1);
	    $af->end($af->end - $start +1);
	    $af->apply_edit(\$subseq); #change, in the reference sequence, the af
	    #restore the initial values of alleleFeature start and end
	    $af->start($af_start);
	    $af->end($af_end);
	    
	}
    }
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

=head2 get_all_Transcripts

  Args       : None
  Example    : @transcripts = @{$individualslice->get_all_Transcripts)};
  Description: Gets all transcripts which overlap this Individual Slice.  If you want to
               specify a particular analysis or type, then you are better off
               using get_all_Genes or get_all_Genes_by_type and iterating
               through the transcripts of each gene.
  Returntype : reference to a list of Bio::EnsEMBL::Transcripts
  Exceptions : none
  Caller     : general

=cut

sub get_all_Transcripts {
  my $self = shift;

  my $transcripts = $self->SUPER::get_all_Transcripts(1);
  $self->map_to_Individual($transcripts);

  return $transcripts;
}


=head2 get_all_Exons

  Arg [1]    : (optional) string $dbtype
               The dbtype of exons to obtain.  This assumes that the db has
               been added to the DBAdaptor under this name (using the
               DBConnection::add_db_adaptor method).
  Example    : @exons = @{$individualSlice->get_all_Exons};
  Description: Gets all exons which overlap this IndividualSlice.  Note that these exons
               will not be associated with any transcripts, so this may not
               be terribly useful.
  Returntype : reference to a list of Bio::EnsEMBL::Exons
  Exceptions : none
  Caller     : general

=cut

sub get_all_Exons {
  my $self = shift;
  my $dbtype = shift;

  my $exons = $self->SUPER::get_all_Exons($dbtype);
  $self->map_to_Individual($exons); #map the exons to the Individual

  return $exons;
}

=head2 get_all_Genes

  Arg [1]    : (optional) string $logic_name
               The name of the analysis used to generate the genes to retrieve
  Arg [2]    : (optional) string $dbtype
               The dbtype of genes to obtain.  This assumes that the db has
               been added to the DBAdaptor under this name (using the
               DBConnection::add_db_adaptor method).
  Example    : @genes = @{$individualSlice->get_all_Genes};
  Description: Retrieves all genes that overlap this slice.
  Returntype : listref of Bio::EnsEMBL::Genes
  Exceptions : none
  Caller     : none

=cut

sub get_all_Genes{
  my ($self, $logic_name, $dbtype) = @_;

  my $genes = $self->SUPER::get_all_Genes($logic_name, $dbtype, 1);

  $self->map_to_Individual($genes);

  foreach my $gene (@{$genes}){
      $self->map_to_Individual($gene->get_all_Exons); #map the Exons to the Individual
      $self->map_to_Individual($gene->get_all_Transcripts); #map the Transcripts to the Individual
  }

  return $genes;
}

=head2 map_to_Individual

    Arg[1]      : ref $features
    Example     : $individualSlice->map_to_Individual($exons);
    Description : Gets the features from the Slice and maps it in the IndividualSlice, using the mapper
                  between Slice and IndividualSlice
    ReturnType  : None
    Exceptions  : None
    Caller      : general
    
=cut

sub map_to_Individual{
    my $self = shift;
    my $features = shift;

    my $mapper = $self->mapper();
    my (@results, @results_ordered, $new_start, $new_end, $new_strand);
    #foreach of the transcripts, map them to the IndividualSlice and replace the Slice with the IndividualSlice
    foreach my $feature (@{$features}){
	$feature->slice($self); #replace the IndividualSlice as the Slice for this feature (the Slice plus the AlleleFeatures)
	#map from the Slice to the Individual Slice
	my @results = $mapper->map_coordinates('Slice',$feature->start,$feature->end,$feature->strand,'Slice');
	#from the results, order them but filter out those that are not coordinates
	@results_ordered = sort {$a->start <=> $b->start} grep {ref($_) eq 'Bio::EnsEMBL::Mapper::Coordinate'} @results;
	$new_start = $results_ordered[0]->start();
	$new_strand = $results_ordered[0]->strand();
	$new_end = $results_ordered[-1]->end();  #get last element of the array, the end of the slice
	$feature->start($new_start);  #update new coordinates
	$feature->end($new_end);
	$feature->strand($new_strand);
  }
}

sub alleleFeatures{
    my $self = shift;
    return $self->{'alleleFeatures'};
}

sub add_AlleleFeature{
    my $self = shift;

    if (@_){
	if(!ref($_[0]) || !$_[0]->isa('Bio::EnsEMBL::Variation::AlleleFeature')) {
	    throw("Bio::EnsEMBL::Variation::AlleleFeature argument expected");
	}
	#add the alleleFeature to the individualSlice
	push @{$self->{'alleleFeatures'}},shift;
    }
}
1;
