#
# Ensembl module for Bio::EnsEMBL::StrainSlice
#
#
# Copyright Team Ensembl
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::StrainSlice - SubClass of the Slice. Represents the slice of the genome 
for a certain strain (applying the alleles for a this strain)

=head1 SYNOPSIS

   $sa = $db->get_SliceAdaptor;

   $slice = $sa->fetch_by_region('chromosome', 'X', 1_000_000, 2_000_000);

   $strainSlice = $slice->get_by_Strain($strain_name);

   #get the sequence from the Strain Slice
   my $seq = $strainSlice->seq();
   print $seq;

   #get a subSlice of the Strain
   my $subSlice_strain = $strainSlice->sub_Slice(5_000,8_000,1)

   #compare two different strains in the same Slice
   my $sliceStrain2 = $slice->get_by_Strain($strain_name2);
   my $differences = $sliceStrain->get_all_differences_StrainSlice($strain_name2);
   foreach my $af (@{$differences}){
      print "There is a difference between $strain_name and $strain_name2 at ", $af->start,"-",$af->end,
            " with allele ", $af->allele_string(),"\n";
   }

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

    Arg [1..N]  : List of named arguments
                  Bio::EnsEMBL::CoordSystem COORD_SYSTEM
                  string SEQ_REGION_NAME,
                  int    START,
                  int    END,
                  string VERSION (optional, defaults to '')
                  int    STRAND, (optional, defaults to 1)
                  Bio::EnsEMBL::DBSQL::SliceAdaptor ADAPTOR (optional)
    Arg[N+1]      : string $strain_name
    Example     : $strainSlice = Bio::EnsEMBL::StrainSlice->new(-coord_system => $cs,
								-start => 1,
								-end => 10000,
								-strand => 1,
								-seq_region_name => 'X',
								-seq_region_length => 12e6,
								-strain_name => $strain_name);
    Description : Creates a new Bio::EnsEMBL::StrainSlice object that will contain a shallow copy of the
                  Slice object, plus additional information such as the train_name this Slice refers to
                  and listref of Bio::EnsEMBL::Variation::VariationFeatures of differences with the
                  reference sequence
    ReturnType  : Bio::EnsEMBL::StrainSlice
    Exceptions  : none
    Caller      : general

=cut

sub new{
    my $caller = shift;
    my $class = ref($caller) || $caller;

    #create the StrainSlice object as the Slice, plust the strain_name attribute
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
    #get the AlleleFeatures in the Strain
    my $af_adaptor = $variation_db->get_AlleleFeatureAdaptor;
    
    if( $af_adaptor ) {
	#get the Population for the given strain
	my $pop_adaptor = $variation_db->get_PopulationAdaptor;

	if ($pop_adaptor){
	    my $population = $pop_adaptor->fetch_by_name($self->{'strain_name'});
	    #check that the population returned is a strain
	    if ((defined $population) && ($population->is_strain)){
		#get all the VariationFeatures in the $population and the Slice given
		my $allele_features = $af_adaptor->fetch_all_by_Slice_Population($self,$population);
		$self->{'alleleFeatures'} = $allele_features;
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
	warning("Not possible to retrieve AlleleFeatureAdaptor from variation database");
	return '';
    }
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
    my @allele_features_ordered = sort {$b->start() <=> $a->start()} @{$self->{'alleleFeatures'}} if (defined $self->{'alleleFeatures'});

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

    Args        : nonre
    Example     : my $differences = $strainSlice->get_all_differences_Slice()
    Description : Gets all differences between the StrainSlice object and the Slice is defined
    ReturnType  : listref of Bio::EnsEMBL::Variation::AlleleFeature
    Exceptions  : none
    Caller      : general

=cut

sub get_all_differences_Slice{
    my $self = shift;

    return $self->{'alleleFeatures'};
}

=head2 get_all_differences_StrainSlice

    Arg[1]      : Bio::EnsEMBL::StrainSlice $ss
    Example     : my $differences = $strainSlice->get_all_differences_StrainSlice($strainslice)
    Description : Gets differences between 2 StrainSlice objects
    ReturnType  : listref of Bio::EnsEMBL::Variation::AlleleFeature
    Exceptions  : thrown on bad argument
    Caller      : general

=cut

sub get_all_differences_StrainSlice{
    my $self = shift;
    my $strainSlice = shift;

    if (!ref($strainSlice) || !$strainSlice->isa('Bio::EnsEMBL::StrainSlice')){
	throw('Bio::EnsEMBL::StrainSlice arg expected');
    }
    if ( @{$self->{'alleleFeatures'}} == 0 && @{$strainSlice->{'alleleFeatures'}} == 0){
	return undef; #there are no differences in any of the Strains
	
    }
    my $differences; #differences between strains
    if (@{$strainSlice->{'alleleFeatures'}} == 0){
	#need to create a copy of alleleFeature for the first Strain
	foreach my $difference (@{$self->{'alleleFeatures'}}){
	    my %vf = %$difference;
	    push @{$differences},bless \%vf,ref($difference);
	}
    }
    elsif (@{$self->{'alleleFeatures'}} == 0){
	#need to create a copy of AlleleFeature, but changing the allele by the allele in the reference sequence
	foreach my $difference (@{$strainSlice->{'alleleFeatures'}}){
	    push @{$differences}, $strainSlice->_convert_difference($difference);
	}
    }
    else{
	#both strains have differences
	#create a hash with the differences in the first slice
	my %allele_features_self = map {$_->start => $_} @{$self->{'alleleFeatures'}};
	foreach my $difference (@{$strainSlice->{'alleleFeatures'}}){
	    #there is no difference in the other strain slice, convert the allele
	    if (!defined $allele_features_self{$difference->start}){
		push @{$differences},$strainSlice->_convert_difference($difference);
	    }		
	    else{
		#if it is defined and have the same allele, delete from the hash since it is not a difference
		#between the strains
		if ($allele_features_self{$difference->start}->allele_string eq $difference->allele_string){
		    delete $allele_features_self{$difference->start};
		}
	    }
	}	
	#and finally, make a shallow copy of the differences in the first strain
	foreach my $difference (values %allele_features_self){
	    my %vf = %$difference;
	    push @{$differences},bless \%vf,ref($difference);
	}

    }
    #need to map differences to the first strain, self, since the coordinates are in the Slice coordinate system
    my $mapper = $self->mapper(); #now that we have the differences, map them in the StrainSlice
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

=head2 sub_Slice

  Arg   1    : int $start
  Arg   2    : int $end
  Arge [3]   : int $strand
  Example    : none
  Description: Makes another StrainSlice that covers only part of this StrainSlice
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
  #map from the Strain to the Slice to get the sub_Slice, and then, apply the differences in the subSlice
  my @results = $mapper->map_coordinates('StrainSlice',$start,$end,$strand,'StrainSlice');
  my $new_start;
  my $new_end;
  my $new_strand;
  my $new_seq;
  #Get need start and end for the subSlice of the StrainSlice
  my @results_ordered = sort {$a->start <=> $b->start} grep {ref($_) eq 'Bio::EnsEMBL::Mapper::Coordinate'} @results;
  $new_start = $results_ordered[0]->start();
  $new_strand = $results_ordered[0]->strand() if (ref($results_ordered[0]) eq 'Bio::EnsEMBL::Mapper::Coordinate');
#  $new_strand = $results_ordered[-1]->strand() if (ref($results_ordered[-1]) eq 'Bio::EnsEMBL::Mapper::Coordinate');
  $new_end = $results_ordered[-1]->end();  #get last element of the array, the end of the slice

  my $subSlice = $self->SUPER::sub_Slice($new_start,$new_end,$new_strand);
  $subSlice->{'strain_name'} = $self->{'strain_name'};

  my $new_alleles; #reference to an array that will contain the variationFeatures in the new subSlice
  #update the VariationFeatures in the sub_Slice of the Strain
  my @results;
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
               The strand of the strain slice to obtain sequence from. Default
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
      $subseq = substr($seq,$start-1,$end - $start +1);
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

=head2 mapper

  Args       : none
  Description: Getter for the mapper between the between the StrainSlice and the Slice it refers to. 
               It is done automatically when necessary to create subSlice or to get the differences between strains
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
	my $mapper = Bio::EnsEMBL::Mapper->new('Slice','StrainSlice');
	#align with Slice
	#get all the VariationFeatures in the strain Slice, from start to end in the Slice
	my @allele_features_ordered = sort {$a->start() <=> $b->start()} @{$self->{'alleleFeatures'}} if (defined $self->{'alleleFeatures'});
	
	my $start_slice = 1;
	my $end_slice;
	my $start_strain = 1;
	my $end_strain;
	my $length_allele;
	#we will walk from left to right in the slice object, updating the start and end strain every time
	#there is a new alleleFeature in the Strain
	foreach my $allele_feature (@allele_features_ordered){
	    #we have a insertion/deletion: marks the beginning of new slice move coordinates	    
	    if ($allele_feature->length_diff != 0){
		$length_allele = $allele_feature->length + $allele_feature->length_diff(); #length of the allele in the Strain		    
		$end_slice = $allele_feature->start() - 1; #set the end of the slice before the alleleFeature
		if ($end_slice >= $start_slice){
		    #normal cases (not with gaps)
		    $end_strain = $end_slice - $start_slice + $start_strain; #set the end of the strain from the beginning plus the offset
		    #add the sequence that maps
		    $mapper->add_map_coordinates('Slice',$start_slice,$end_slice,1,'StrainSlice',$start_strain,$end_strain);
		    #and add the indel
		    $mapper->add_indel_coordinates('Slice',$end_slice+1,$end_slice + $allele_feature->length,1,'StrainSlice',$end_strain+1,$end_strain + $length_allele);
		    $start_strain = $end_strain + $length_allele + 1; #set the beginning of the strain after the allele
		}
		else{
		    #add the indel
		    $mapper->add_indel_coordinates('Slice',$end_slice+1,$end_slice + $allele_feature->length,1,'StrainSlice',$end_strain+1,$end_strain + $length_allele);
		    $start_strain += $length_allele;
		}
		$start_slice = $end_slice + $allele_feature->length+ 1; #set the beginning of the slice after the variation feature
	    }
	}
	if ($start_slice <= $self->length){
	    #if we haven't reached the end of the StrainSlice, add the final map coordinates between the strain and the slice
	    $mapper->add_map_coordinates('Slice',$start_slice,$self->length,1,'StrainSlice',$start_strain,$start_strain + $self->length - $start_slice);
	}
	$self->{'mapper'} = $mapper;
    }
    return $self->{'mapper'};
}


=head2 get_all_Transcripts

  Args       : None
  Example    : @transcripts = @{$slice->get_all_Transcripts)_};
  Description: Gets all transcripts which overlap this Strain Slice.  If you want to
               specify a particular analysis or type, then you are better off
               using get_all_Genes or get_all_Genes_by_type and iterating
               through the transcripts of each gene.
  Returntype : reference to a list of Bio::EnsEMBL::Transcripts
  Exceptions : none
  Caller     : general

=cut

sub get_all_Transcripts {
  my $self = shift;

  if(!$self->adaptor()) {
    warning('Cannot get Transcripts without attached adaptor');
    return [];
  }

  my $ta = $self->adaptor()->db()->get_TranscriptAdaptor();
  #get first all Transcript in the Slice
  my $transcripts = $ta->fetch_all_by_Slice($self, 1);
  #get the mapper between the Slice and the StrainSlice
  my $mapper = $self->mapper();
  my (@results, @results_ordered, $new_start, $new_end, $new_strand);
  #foreach of the transcripts, map them to the StrainSlice and replace the Slice with the StrainSlice
  foreach my $transcript (@{$transcripts}){
      $transcript->{'slice'} = $self; #add the StrainSlice as the Slice for this Transcript
      #map from the Slice to the Strain Slice
      my @results = $mapper->map_coordinates('Slice',$transcript->start,$transcript->end,$transcript->strand,'Slice');
      #from the results, order them but filter out those that are not coordinates
      @results_ordered = sort {$a->start <=> $b->start} grep {ref($_) eq 'Bio::EnsEMBL::Mapper::Coordinate'} @results;
      $new_start = $results_ordered[0]->start();
      $new_strand = $results_ordered[0]->strand();
      $new_end = $results_ordered[-1]->end();  #get last element of the array, the end of the slice
      $transcript->start($new_start);  #update the new coordinates for the transcript
      $transcript->end($new_end);
      $transcript->strand($new_strand);
      
  }
  return $transcripts;
}


=head2 get_all_Exons

  Arg [1]    : none
  Example    : @exons = @{$slice->get_all_Exons};
  Description: Gets all exons which overlap this StrainSlice.  Note that these exons
               will not be associated with any transcripts, so this may not
               be terribly useful.
  Returntype : reference to a list of Bio::EnsEMBL::Exons
  Exceptions : none
  Caller     : general

=cut

sub get_all_Exons {
  my $self = shift;

  if(!$self->adaptor()) {
    warning('Cannot get Exons without attached adaptor');
    return [];
  }
  my $exons = $self->adaptor->db->get_ExonAdaptor->fetch_all_by_Slice($self);
  #get the mapper between the Slice and the StrainSlice
  my $mapper = $self->mapper();
  my (@results, @results_ordered, $new_start, $new_end, $new_strand);
  #foreach of the transcripts, map them to the StrainSlice and replace the Slice with the StrainSlice
  foreach my $exon (@{$exons}){
      $exon->{'slice'} = $self; #add the StrainSlice as the Slice for this Exon
      #map from the Slice to the Strain Slice
      my @results = $mapper->map_coordinates('Slice',$exon->start,$exon->end,$exon->strand,'Slice');
      #from the results, order them but filter out those that are not coordinates
      @results_ordered = sort {$a->start <=> $b->start} grep {ref($_) eq 'Bio::EnsEMBL::Mapper::Coordinate'} @results;
      $new_start = $results_ordered[0]->start();
      $new_strand = $results_ordered[0]->strand();
      $new_end = $results_ordered[-1]->end();  #get last element of the array, the end of the slice
      $exon->start($new_start);  #update new coordinates for the Exon
      $exon->end($new_end);
      $exon->strand($new_strand);
  }
  return $exons;
}

1;
