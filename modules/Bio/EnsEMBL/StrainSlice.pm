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

Bio::EnsEMBL::StrainSlice - SubClass of the Slice. Represents the slice
of the genome for a certain strain (applying the variations)

=head1 SYNOPSIS

  $sa = $db->get_SliceAdaptor;

  $slice =
    $sa->fetch_by_region( 'chromosome', 'X', 1_000_000, 2_000_000 );

  $strainSlice = $slice->get_by_strain($strain_name);

  # get the sequence from the Strain Slice
  my $seq = $strainSlice->seq();
  print $seq;

  # get allele features between this StrainSlice and the reference
  my $afs = $strainSlice->get_all_AlleleFeatures_Slice();
  foreach my $af ( @{$afs} ) {
    print "AlleleFeature in position ", $af->start, "-", $af->end,
      " in strain with allele ", $af->allele_string, "\n";
  }

  # compare a strain against another strain
  my $strainSlice_2 = $slice->get_by_strain($strain_name_2);
  my $differences =
    $strainSlice->get_all_differences_StrainSlice($strainSlice_2);

  foreach my $difference ( @{$differences} ) {
    print "Difference in position ", $difference->start, "-",
      $difference->end(),           " in strain with allele ",
      $difference->allele_string(), "\n";
  }

=head1 DESCRIPTION

A StrainSlice object represents a region of a genome for a certain
strain.  It can be used to retrieve sequence or features from a strain.

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

@ISA = qw(Bio::EnsEMBL::Slice);


=head2 new

    Arg[1]      : Bio::EnsEMBL::Slice $slice
    Arg[2]      : string $strain_name
    Example     : $strainSlice = Bio::EnsEMBL::StrainSlice->new(-.... => ,
							      -strain_name => $strain_name);
    Description : Creates a new Bio::EnsEMBL::StrainSlice object that will contain a shallow copy of the
                  Slice object, plus additional information such as the Strain this Slice refers to
                  and listref of Bio::EnsEMBL::Variation::AlleleFeatures of differences with the
                  reference sequence
    ReturnType  : Bio::EnsEMBL::StrainSlice
    Exceptions  : none
    Caller      : general

=cut

sub new{
  my $caller = shift;
  deprecate("new is deprecated and will be removed in e88. Please use Bio::EnsEMBL::Variation::StrainSlice::new instead");
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

  my $af_adaptor = $variation_db->get_AlleleFeatureAdaptor;

  if( $af_adaptor ) {
    #get the Sample for the given strain
    my $sample_adaptor = $variation_db->get_SampleAdaptor;

    if ($sample_adaptor){
      my $sample = shift @{$sample_adaptor->fetch_all_by_name($self->{'strain_name'})}; #the name should be unique for a strain
      #check that the sample returned isin the database

      if (defined $sample){
        my $allele_features = $af_adaptor->fetch_all_by_Slice($self, $sample);
        #warning("No strain genotype data available for Slice ".$self->name." and Strain ".$sample->name) if ! defined $allele_features->[0];
        my $vf_ids = {}; #hash containing the relation vf_id->af
        $self->{'_strain'} = $sample;		
        map {defined $_->{'_variation_feature_id'} ? $vf_ids->{$_->{'_variation_feature_id'}} = $_ : ''} @{$allele_features};
        #		my $new_allele_features = $self->_filter_af_by_coverage($allele_features);
        #		$self->{'alleleFeatures'} = $new_allele_features;
        $self->{'alleleFeatures'} = $allele_features || [];
        $self->{'_vf_ids'} = $vf_ids;
        return $self;
      }
      else{ 
        warning("Strain ($self->{strain_name}) not in the database");
        return $self;
      }
    }
    else{
      warning("Not possible to retrieve SampleAdaptor from the variation database");
      return '';
    }
  } else {
    warning("Not possible to retrieve VariationFeatureAdaptor from variation database");
    return '';
  }
}

=head2 _filter_af_by_coverage

    Arg [1]     : listref to Bio::EnsEMBL::Variation::AlleleFeatures  $allele_features
    Example     : my $new_list_allele_features = $strainSlice->_filter_af_by_coverage($allele_features);
    Description : For a list of allele features, gets a new list where they are filter depending on coverage
    ReturnType  : listref of Bio::EnsEMBL::Variation::AlleleFeature
    Exceptions  : none
    Caller      : internal function

=cut

sub _filter_af_by_coverage{
  my $self = shift;
  deprecate("_filter_af_by_coverage is deprecated and will be removed in e88. Please use Bio::EnsEMBL::Variation::StrainSlice::_filter_af_by_coverage instead");
  my $allele_features = shift;
  
  my $variation_db = $self->adaptor->db->get_db_adaptor('variation');
  
  unless($variation_db) {
    warning("Variation database must be attached to core database to " .
    "retrieve variation information" );
    return '';
  }  
  
  return $allele_features unless $self->sample->has_coverage;
  
  my $rc_adaptor = $variation_db->get_ReadCoverageAdaptor();
  #this is ugly, but ReadCoverage is always defined in the positive strand
  
  ### EK : - it looks like the arguments to fetch_all_by_Slice_Sample_depth have changed
  ###  passing 1 will only get you the coverage of level 1
  ###  by omitting the parameter we take into account all coverage regions 
  #    my $rcs = $rc_adaptor->fetch_all_by_Slice_Sample_depth($self,$self->{'_strain'},1);
  my $rcs = $rc_adaptor->fetch_all_by_Slice_Sample_depth($self,$self->{'_strain'});
  my $new_af;
  foreach my $af (@{$allele_features}){
    foreach my $rc (@{$rcs}){
      if ($af->start <= $rc->end and $af->start >= $rc->start){
        push @{$new_af}, $af;
        last;
      }
    }
  }
  
  return $new_af;
}


=head2 strain_name

    Arg [1]     : (optional) string $strain_name
    Example     : my $strain_name = $strainSlice->strain_name();
    Description : Getter/Setter for the name of the strain
    ReturnType  : string
    Exceptions  : none
    Caller      : general

=cut

sub strain_name{
   my $self = shift;
  deprecate("strain_name is deprecated and will be removed in e88. Please use Bio::EnsEMBL::Variation::StrainSlice::strain_name instead");
   if (@_){
       $self->{'strain_name'} = shift @_;
   }
   return $self->{'strain_name'};
}


=head2 sample

    Example     : my $sample = $strainSlice->sample();
    Description : Getter for the Sample object associated
    ReturnType  : Bio::EnsEMBL::Variation::Sample
    Exceptions  : none
    Caller      : general

=cut

sub sample {
  deprecate("sample is deprecated and will be removed in e88. Please use Bio::EnsEMBL::Variation::StrainSlice::sample instead");
  return $_[0]->{_strain};
}


=head2 display_Slice_name

    Args        : none
    Example     : my $strain_name = $strainSlice->display_Slice_name();
    Description : Getter for the name of the strain
    ReturnType  : string
    Exceptions  : none
    Caller      : webteam

=cut

sub display_Slice_name{
    my $self = shift;
  deprecate("display_Slice_name is deprecated and will be removed in e88. Please use Bio::EnsEMBL::Variation::StrainSlice::display_Slice_name instead");
    return $self->strain_name;
}

=head2 seq

  Arg [1]    : int $with_coverage (optional)
  Example    : print "SEQUENCE = ", $strainSlice->seq();
  Description: Returns the sequence of the region represented by this
               slice formatted as a string in the strain. If flag with_coverage
               is set to 1, returns sequence if there is coverage in the region
  Returntype : string
  Exceptions : none
  Caller     : general

=cut

sub seq {
  my $self = shift;
  deprecate("seq is deprecated and will be removed in e88. Please use Bio::EnsEMBL::Variation::StrainSlice::seq instead");
  my $with_coverage = shift;

  $with_coverage ||= 0;

  # special case for in-between (insert) coordinates
  return '' if($self->start() == $self->end() + 1);

  return $self->{'seq'} if($self->{'seq'});

  if($self->adaptor()) {
    my $seqAdaptor = $self->adaptor()->db()->get_SequenceAdaptor();
    my $reference_sequence = $seqAdaptor->fetch_by_Slice_start_end_strand($self,1,undef,1); #get the reference sequence for that slice
    
    # default to lowercase sequence
    # this gets overwritten with uppercase sequence in covered regions
    # remains lowercase in strains with implicit coverage
    $$reference_sequence = lc($$reference_sequence) if $with_coverage;
	
    #apply all differences to the reference sequence
    #first, in case there are any indels, create the new sequence (containing the '-' bases)
    # sort edits in reverse order to remove complication of
    # adjusting downstream edits
    my @indels_ordered;
    @indels_ordered = sort {$b->start() <=> $a->start()} @{$self->{'alignIndels'}} if (defined $self->{'alignIndels'});

    foreach my $vf (@indels_ordered){
      $vf->apply_edit($reference_sequence); #change, in the reference sequence, the vf
    }
	
    #need to find coverage information if diffe
    # sort edits in reverse order to remove complication of
    # adjusting downstream edits
    my @variation_features_ordered;
    @variation_features_ordered = sort {$b->start() <=> $a->start()} @{$self->{'alleleFeatures'}} if (defined $self->{'alleleFeatures'});

    foreach my $vf (@variation_features_ordered){
      $vf->apply_edit($reference_sequence); #change, in the reference sequence, the vf
    }

    #need to find coverage information if different from reference
    my $sampleAdaptor = $self->adaptor->db->get_db_adaptor('variation')->get_SampleAdaptor;
    my $ref_strain = $sampleAdaptor->get_reference_strain_name;
    $self->_add_coverage_information($reference_sequence) if ($with_coverage == 1 && $self->strain_name ne $ref_strain);
    return substr(${$reference_sequence},0,1) if ($self->length == 1); 
    return substr(${$reference_sequence},0,$self->expanded_length); #returns the reference sequence, applying the variationFeatures. Remove additional bases added due to indels
  }

  # no attached sequence, and no db, so just return Ns
  return 'N' x $self->length();
}

sub expanded_length {
	my $self = shift;
  deprecate("expanded_length is deprecated and will be removed in e88. Please use Bio::EnsEMBL::Variation::StrainSlice::expanded_length instead");
	my $length = $self->SUPER::length();
	
	foreach my $af(@{$self->{'alleleFeatures'}}) {
		$length += $af->length_diff() if $af->length_diff > 0;
	}
	
	return $length;
}



sub _add_coverage_information{
  my $self = shift;
  deprecate("_add_coverage_information is deprecated and will be removed in e88. Please use Bio::EnsEMBL::Variation::StrainSlice::_add_coverage_information instead");
  my $reference_sequence = shift;
  
  my $variation_db = $self->adaptor->db->get_db_adaptor('variation');
  
  unless($variation_db) {
    warning("Variation database must be attached to core database to " .
    "retrieve variation information" );
    return '';
  }
  
  # only fetch RC data if sample is flagged as having coverage
  return unless $self->sample->has_coverage;
  
  my $rc_adaptor = $variation_db->get_ReadCoverageAdaptor();
  ### EK : - it looks like the arguments to fetch_all_by_Slice_Sample_depth have changed
  ###  passing 1 will only get you the coverage of level 1
  ###  by omitting the parameter we take into account all coverage regions 
  #    my $rcs = $rc_adaptor->fetch_all_by_Slice_Sample_depth($self,$self->{'_strain'},1);
  my $rcs = $rc_adaptor->fetch_all_by_Slice_Sample_depth($self,$self->{'_strain'});
  my $rcs_sorted;
  @{$rcs_sorted} = sort {$a->start <=> $b->start} @{$rcs} if ($self->strand == -1);
  $rcs = $rcs_sorted if ($self->strand == -1);
  my $start = 1;
	
	
	# wm2 - new code to mask sequence, instead starts with masked string
	# and unmasks seq where there is read coverage
	
	# get all length-changing vars
  my @indels_ordered;
	@indels_ordered = sort {$a->start() <=> $b->start()} @{$self->{'alignIndels'}} if (defined $self->{'alignIndels'});
	
	my $masked_seq = '~' x length($$reference_sequence);
	
	foreach my $rc(@{$rcs}) {
		my ($start, $end) = ($rc->start, $rc->end);
		
		# adjust region for indels
		foreach my $indel(@indels_ordered) {
			next if $rc->start > $end;
			
			# if within RC region, only need adjust the end
			$start += $indel->length_diff unless $indel->start > $start;
			$end += $indel->length_diff;
		}
		
		# adjust coords for seq boundaries
		$start = 1 if $start < 1;
		$end = CORE::length($masked_seq) if $end > CORE::length($masked_seq);
		
		# now unmask the sequence using $$reference_sequence
		substr($masked_seq, $start - 1, $end - $start + 1) = uc(substr($$reference_sequence, $start - 1, $end - $start + 1));
	}
	
	# copy the masked sequence to the reference sequence
	$$reference_sequence = $masked_seq;
}


=head2 get_AlleleFeature

    Arg[1]      : Bio::EnsEMBL::Variation::VariationFeature $vf
    Example     : my $af = $strainSlice->get_AlleleFeature($vf);
    Description : Returns the AlleleFeature object associated with the VariationFeature (if any)
    ReturnType  : Bio::EnsEMBL::Variation::AlleleFeature
    Exceptions  : none
    Caller      : general

=cut

sub get_AlleleFeature{
    my $self = shift;
  deprecate("get_all_AlleleFeature is deprecated and will be removed in e88. Please use Bio::EnsEMBL::Variation::StrainSlice::get_all_AlleleFeature instead");
    my $vf = shift;
    
    my $af;
    #look at the hash containing the relation vf_id->alleleFeature, if present, return object, otherwise, undef
    $af = $self->{'_vf_ids'}->{$vf->dbID} if (defined $self->{'_vf_ids'}->{$vf->dbID});
    return $af;
}


=head2 get_all_AlleleFeatures_Slice

    Arg[1]      : int $with_coverage (optional)
    Example     : my $af = $strainSlice->get_all_AlleleFeatures_Slice()
    Description : Gets all AlleleFeatures between the StrainSlice object and the Slice is defined.
                  If argument $with_coverage set to 1, returns only AF if they have coverage information
    ReturnType  : listref of Bio::EnsEMBL::Variation::AlleleFeature
    Exceptions  : none
    Caller      : general

=cut

sub get_all_AlleleFeatures_Slice{
  my $self = shift;
  deprecate("get_all_AlleleFeatures_Slice is deprecated and will be removed in e88. Please use Bio::EnsEMBL::Variation::StrainSlice::get_all_AlleleFeatures_Slice instead");
  my $with_coverage = shift;

  my $variation_db = $self->adaptor->db->get_db_adaptor('variation');

  unless($variation_db) {
    warning("Variation database must be attached to core database to " .
    "retrieve variation information" );
    return '';
  }
  my $sampleAdaptor = $variation_db->get_SampleAdaptor();
  my $ref_name =  $sampleAdaptor->get_reference_strain_name;
  return [] if ($self->strain_name eq $ref_name);
  $with_coverage ||= 0; #by default, get all AlleleFeatures
  if ($with_coverage == 1){
    my $new_allele_features = $self->_filter_af_by_coverage($self->{'alleleFeatures'});
    return $new_allele_features || [];
  }

  return $self->{'alleleFeatures'} || [];
}

=head2 get_all_differences_StrainSlice

    Arg[1]      : Bio::EnsEMBL::StrainSlice $ss
    Example     : my $differences = $strainSlice->get_all_differences_StrainSlice($ss)
    Description : Gets differences between 2 StrainSlice objects
    ReturnType  : listref of Bio::EnsEMBL::Variation::AlleleFeature
    Exceptions  : thrown on bad argument
    Caller      : general

=cut

sub get_all_differences_StrainSlice{
    my $self = shift;
  deprecate("get_all_differences_StrainSlice is deprecated and will be removed in e88. Please use Bio::EnsEMBL::Variation::StrainSlice::get_all_differences_StrainSlice instead");

    my $strainSlice = shift;

    if (!ref($strainSlice) || !$strainSlice->isa('Bio::EnsEMBL::StrainSlice')){
	throw('Bio::EnsEMBL::StrainSlice arg expected');
    }
    if ( @{$self->{'alleleFeatures'}} == 0 && @{$strainSlice->{'alleleFeatures'}} == 0){
	return undef; #there are no differences in any of the Strains
	
    }
    my $differences; #differences between strains
    if (@{$strainSlice->{'alleleFeatures'}} == 0){
	#need to create a copy of VariationFeature
	foreach my $difference (@{$self->{'alleleFeatures'}}){
	    my %vf = %$difference;
	    push @{$differences},bless \%vf,ref($difference);
	}
    }
    elsif (@{$self->{'alleleFeatures'}} == 0){
	#need to create a copy of VariationFeature, but changing the allele by the allele in the reference
	foreach my $difference (@{$strainSlice->{'alleleFeatures'}}){
	    push @{$differences}, $strainSlice->_convert_difference($difference);
	}
    }
    else{
	#both strains have differences
	#create a hash with the differences in the self strain slice
	my %variation_features_self = map {$_->start => $_} @{$self->{'alleleFeatures'}};
	foreach my $difference (@{$strainSlice->{'alleleFeatures'}}){
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
  deprecate("_convert_difference is deprecated and will be removed in e88. Please use Bio::EnsEMBL::Variation::StrainSlice::_convert_difference instead");
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
  deprecate("sub_Slice is deprecated and will be removed in e88. Please use Bio::EnsEMBL::Variation::StrainSlice::sub_Slice instead");
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

  foreach my $variationFeature (@{$self->{'alleleFeatures'}}){
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
  $subSlice->{'alleleFeatures'} = $new_variations;
  return $subSlice;

}

=head2 ref_subseq

  Arg  [1]   : int $startBasePair
               relative to start of slice, which is 1.
  Arg  [2]   : int $endBasePair
               relative to start of slice.
  Arg  [3]   : (optional) int $strand
               The strand of the slice to obtain sequence from. Default
               value is 1.
  Description: returns string of dna from reference sequence
  Returntype : txt
  Exceptions : end should be at least as big as start
               strand must be set
  Caller     : general

=cut

sub ref_subseq{
  my $self = shift;
  deprecate("ref_subseq is deprecated and will be removed in e88. Please use Bio::EnsEMBL::Variation::StrainSlice::ref_subseq instead");
  my $start = shift;
  my $end = shift;
  my $strand = shift;
  # special case for in-between (insert) coordinates
  return '' if($start == $end + 1);

  my $subseq;
  if($self->adaptor){
    my $seqAdaptor = $self->adaptor->db->get_SequenceAdaptor();
    $subseq = ${$seqAdaptor->fetch_by_Slice_start_end_strand
      ( $self, $start,
        $end, $strand )};
  } else {
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
  deprecate("subseq is deprecated and will be removed in e88. Please use Bio::EnsEMBL::Variation::StrainSlice::subseq instead");

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


      $seq = $self->seq;
      reverse_comp(\$seq) if ($strand == -1);
      $subseq = substr($seq,$start-1,$end - $start + 1);
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
  deprecate("mapper is deprecated and will be removed in e88. Please use Bio::EnsEMBL::Variation::StrainSlice::mapper instead");
  
    if (@_) {
	delete $self->{'mapper'};
    }
    if(!defined $self->{'mapper'}){
	#create the mapper between the Slice and StrainSlice
	my $mapper = Bio::EnsEMBL::Mapper->new('Slice','StrainSlice');
	#align with Slice
	#get all the VariationFeatures in the strain Slice, from start to end in the Slice
  my @variation_features_ordered;
	@variation_features_ordered = sort {$a->start() <=> $b->start()} @{$self->{'alleleFeatures'}} if (defined $self->{'alleleFeatures'});
	
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


=head2 get_all_VariationFeatures

    Arg[1]     : int $with_coverage (optional)
    Description :returns all alleleFeatures features on this slice. 
    ReturnType : listref of Bio::EnsEMBL::Variation::AlleleFeature
    Exceptions : none
    Caller     : contigview, snpview

=cut

sub get_all_VariationFeatures {
  my $self = shift;
  deprecate("get_all_VariationFeatures is deprecated and will be removed in e88. Please use Bio::EnsEMBL::Variation::StrainSlice::get_all_VariationFeatures instead");
  my $with_coverage = shift;
  $with_coverage ||= 0;
  return $self->get_all_AlleleFeatures_Slice($with_coverage);
}

=head2 get_original_seq_region_position

  Arg  [1]   : int $position
               relative to start of slice, which is 1.
  Description: Placeholder method - this method has no explicit use beyond
			   providiing compatibility with AlignSlice. To map positions
			   between the StrainSlice and the reference slice, use the
			   mapper and its methods.
  Returntype : ($strainSlice, $seq_region_position), an array where the first
               element is a Bio::EnsEMBL::StrainSlice and the second one is the
               requested seq_region_position.
  Exceptions : none
  Caller     : general

=cut

sub get_original_seq_region_position {
    my $self = shift;
  deprecate("get_original_seq_region_position is deprecated and will be removed in e88. Please use Bio::EnsEMBL::Variation::StrainSlice::get_original_seq_region_position instead");
    my $position = shift;
    #coordinates in the AlignSlice and Slice are the same, so far will return the same Slice
    #and coordinate
    return ($self,$position);
}


=head2 remove_indels

    Args        : none
    Example     : $strainSlice->remove_indels();
    Description : Removes insertions and deletions from the allele features
				  of this object
    ReturnType  : none
    Exceptions  : none
    Caller      : webteam

=cut

sub remove_indels {
	my $self = shift;
  deprecate("remove_indels is deprecated and will be removed in e88. Please use Bio::EnsEMBL::Variation::StrainSlice::remove_indels instead");
	my @new_afs = grep { $_->variation->var_class ne 'in-del' } @{$self->{'alleleFeatures'}};
	
	$self->{'alleleFeatures'} = \@new_afs;
}

1;
