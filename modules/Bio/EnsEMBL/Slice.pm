#
# Ensembl module for Bio::EnsEMBL::Slice
#
#
# Copyright Team Ensembl
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::Slice - Arbitary Slice of a genome

=head1 SYNOPSIS

   $sa = $db->get_SliceAdaptor;

   $slice = $sa->fetch_by_region('chromosome', 'X', 1_000_000, 2_000_000);

   #get some attributes of the slice
   my $seqname = $slice->seq_region_name();
   my $start = $slice->start();
   my $end   = $slice->end();
   print "$seqname, $start-$end\n";

   #get the sequence from the slice
   my $seq = $slice->seq();
   print $seq;

   #get some features from the slice
   foreach $gene ( @{$slice->get_all_Genes} ) {
      # do something with a gene
   }

   foreach my $feature ( @{$slice->get_all_DnaAlignFeatures}) {
     # do something with dna-dna alignments
   }


=head1 DESCRIPTION

A Slice object represents a region of a genome.  It can be used to retrieve
sequence or features from an area of interest.

=head1 CONTACT

This modules is part of the Ensembl project http://www.ensembl.org

Questions can be posted to the ensembl-dev mailing list:
ensembl-dev@ebi.ac.uk

=head1 METHODS

=cut

package Bio::EnsEMBL::Slice;
use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Root; #included for backwards compatibility
use Bio::EnsEMBL::Tile; #included for backwards compatibility
use Bio::EnsEMBL::Chromosome; #included for backwards compatibility
use Bio::EnsEMBL::RawContig; #included for backwards compatibility
use Bio::PrimarySeqI;# included only for backwards compatibility

use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Utils::Exception qw(throw deprecate warning);

use Bio::EnsEMBL::RepeatMaskedSlice;
use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp);

#inheritance to Bio::EnsEMBL::Root will eventually be removed
@ISA = qw(Bio::EnsEMBL::Root Bio::PrimarySeqI);


=head2 new

  Arg [...]  : List of named arguments
               string COORD_SYSTEM
               string SEQ_REGION_NAME,
               int    START,
               int    END,
               string VERSION (optional, defaults to '')
               int    STRAND, (optional, defaults to 1)
               Bio::EnsEMBL::DBSQL::SliceAdaptor ADAPTOR (optional)
  Example    : $slice = Bio::EnsEMBL::Slice->new(-coord_system => $cs,
                                                 -start => 1,
                                                 -end => 10000,
                                                 -strand => 1,
                                                 -seq_region_name => 'X',
                                                 -seq_region_length => 12e6,
                                                 -adaptor => $slice_adaptor);
  Description: Creates a new slice object.  A slice represents a region
               of sequence in a particular coordinate system.  Slices can be
               used to retrieve sequence and features from an area of
               interest in a genome.

               Coordinates start at 1 and are inclusive.  Negative
               coordinates or coordinates exceeding the length of the
               seq_region are permitted.  Start must be less than or equal.
               to end regardless of the strand.

               Slice objects are immutable. Once instantiated their attributes
               (with the exception of the adaptor) may not be altered.  To
               change the attributes a new slice must be created.
  Returntype : Bio::EnsEMBL::Slice
  Exceptions : none
  Caller     : general, Bio::EnsEMBL::SliceAdaptor

=cut

sub new {
  my $caller = shift;

  #new can be called as a class or object method
  my $class = ref($caller) || $caller;

  my ($seq, $coord_system, $seq_region_name, $seq_region_length,
      $start, $end, $strand, $adaptor, $empty) =
        rearrange([qw(SEQ COORD_SYSTEM SEQ_REGION_NAME SEQ_REGION_LENGTH
                      START END STRAND ADAPTOR EMPTY)], @_);

  #empty is only for backwards compatibility
  if($empty) {
    deprecate("Creation of empty slices is no longer needed" .
              "and is deprecated");
    return bless({'empty' => 1, 'adaptor' => $adaptor}, $class);
  }

  $seq_region_name  || throw('SEQ_REGION_NAME argument is required');
  defined($start)   || throw('START argument is required');
  defined($end)     || throw('END argument is required');
  ($start <= $end)  || throw('start must be less than or equal to end');

  if(!defined($seq_region_length)) {
    warning("SEQ_REGION_LENGTH argument not provided.\nUsing end of Slice ".
            "which may not be correct.");
    $seq_region_length = $end;
  } else {
    throw('SEQ_REGION_LENGTH must be > 0') if($seq_region_length < 1);
  }

  if($seq && length($seq) != ($end - $start + 1)){
      throw('SEQ must be the same length as the defined LENGTH not '.
            length($seq).' compared to '.($end - $start + 1));
  }
  if(!$coord_system || !ref($coord_system) ||
     !$coord_system->isa('Bio::EnsEMBL::CoordSystem')) {
    throw('COORD_SYSTEM argument must be a Bio::EnsEMBL::CoordSystem');
  } elsif($coord_system->is_top_level()) {
    throw('Cannot create slice on toplevel CoordSystem.');
  }

  #strand defaults to 1 if not defined
  $strand ||= 1;

  if($strand != 1 && $strand != -1) {
    throw('STRAND argument must be -1 or 1');
  }

  if(defined($adaptor)) {
    if(!ref($adaptor) || !$adaptor->isa('Bio::EnsEMBL::DBSQL::SliceAdaptor')) {
      throw('ADAPTOR argument must be a Bio::EnsEMBL::DBSQL::SliceAdaptor');
    }
  }

  return bless {'coord_system'      => $coord_system,
                'seq'               => $seq,
                'seq_region_name'   => $seq_region_name,
                'seq_region_length' => $seq_region_length,
                'start'             => $start,
                'end'               => $end,
                'strand'            => $strand,
                'adaptor'           => $adaptor}, $class;
}



=head2 adaptor

  Arg [1]    : (optional) Bio::EnsEMBL::DBSQL::SliceAdaptor $adaptor
  Example    : $adaptor = $slice->adaptor();
  Description: Getter/Setter for the slice object adaptor used
               by this slice for database interaction.
  Returntype : Bio::EnsEMBL::DBSQL::SliceAdaptor
  Exceptions : none
  Caller     : general

=cut

sub adaptor{
   my $self = shift;

   if(@_) {
     my $ad = shift;
     if(defined($ad)) {
       if(!ref($ad) || !$ad->isa('Bio::EnsEMBL::DBSQL::SliceAdaptor')) {
         throw('Argument must be a Bio::EnsEMBL::DBSQL::SliceAdaptor');
       }
     }
     $self->{'adaptor'} = $ad;
   }

   return $self->{'adaptor'};
}



=head2 seq_region_name

  Arg [1]    : none
  Example    : $seq_region = $slice->seq_region_name();
  Description: Returns the name of the seq_region that this slice is on. For
               example if this slice is in chromosomal coordinates the
               seq_region_name might be 'X' or '10'.

               This function was formerly named chr_name, but since slices can
               now be on coordinate systems other than chromosomal it has been
               changed.
  Returntype : string
  Exceptions : none
  Caller     : general

=cut

sub seq_region_name {
  my $self = shift;
  return $self->{'seq_region_name'};
}



=head2 seq_region_length

  Arg [1]    : none
  Example    : $seq_region_length = $slice->seq_region_length();
  Description: Returns the length of the entire seq_region that this slice is
               on. For example if this slice is on a chromosome this will be
               the length (in basepairs) of the entire chromosome.
  Returntype : int
  Exceptions : none
  Caller     : general

=cut

sub seq_region_length {
  my $self = shift;
  return $self->{'seq_region_length'};
}


=head2 coord_system

  Arg [1]    : none
  Example    : print $slice->coord_system->name();
  Description: Returns the coordinate system that this slice is on.
  Returntype : Bio::EnsEMBL::CoordSystem
  Exceptions : none
  Caller     : general

=cut

sub coord_system {
  my $self = shift;
  return $self->{'coord_system'};
}



=head2 start

  Arg [1]    : none
  Example    : $start = $slice->start();
  Description: Returns the start position of this slice relative to the
               start of the sequence region that it was created on.
               Coordinates are inclusive and start at 1.  Negative coordinates
               or coordinates exceeding the length of the sequence region are
               permitted.  Start is always less than or equal to end 
               regardless of the orientation of the slice.
  Returntype : int
  Exceptions : none
  Caller     : general

=cut

sub start {
  my $self = shift;
  return $self->{'start'};
}



=head2 end

  Arg [1]    : none
  Example    : $end = $slice->end();
  Description: Returns the end position of this slice relative to the
               start of the sequence region that it was created on.
               Coordinates are inclusive and start at 1.  Negative coordinates
               or coordinates exceeding the length of the sequence region are
               permitted.  End is always greater than or equal to start
               regardless of the orientation of the slice.
  Returntype : int
  Exceptions : none
  Caller     : general

=cut

sub end {
  my $self = shift;
  return $self->{'end'};
}



=head2 strand

  Arg [1]    : none
  Example    : $strand = $slice->strand();
  Description: Returns the orientation of this slice on the seq_region it has
               been created on
  Returntype : int (either 1 or -1)
  Exceptions : none
  Caller     : general, invert

=cut

sub strand{
  my $self = shift;
  return $self->{'strand'};
}





=head2 name

  Arg [1]    : none
  Example    : my $results = $cache{$slice->name()};
  Description: Returns the name of this slice. The name is formatted as a colon
               delimited string with the following attributes:
               coord_system:version:seq_region_name:start:end:strand

               Slices with the same name are equivalent and thus the name can
               act as a hash key.
  Returntype : string
  Exceptions : none
  Caller     : general

=cut

sub name {
  my $self = shift;

  return join(':',
              $self->{'coord_system'}->name(),
              $self->{'coord_system'}->version(),
              $self->{'seq_region_name'},
              $self->{'start'},
              $self->{'end'},
              $self->{'strand'});
}



=head2 length

  Arg [1]    : none
  Example    : $length = $slice->length();
  Description: Returns the length of this slice in basepairs
  Returntype : int
  Exceptions : none
  Caller     : general

=cut

sub length {
  my ($self) = @_;

  return $self->{'end'} - $self->{'start'} + 1;
}



=head2 invert

  Arg [1]    : none
  Example    : $inverted_slice = $slice->invert;
  Description: Creates a copy of this slice on the opposite strand and
               returns it.
  Returntype : Bio::EnsEMBL::Slice
  Exceptions : none
  Caller     : general

=cut

sub invert {
  my $self = shift;

  #make a shallow copy of the slice via a hash copy,
  my %s = %$self;
  if(!$self->adaptor){
    warning("You can't invert a Slice which isn't attached to a database");
    return;
  }
  #flip the strand,
  $s{'strand'} = $self->{'strand'} * -1;
  #bless and return the copy
  return  bless \%s, ref $self;
}



=head2 seq

  Arg [1]    : none
  Example    : print "SEQUENCE = ", $slice->seq();
  Description: Returns the sequence of the region represented by this
               slice formatted as a string.
  Returntype : string
  Exceptions : none
  Caller     : general

=cut

sub seq {
  my $self = shift;

  if(!$self->{'seq'}){
    my $seqAdaptor = $self->adaptor->db->get_SequenceAdaptor;
    return ${$seqAdaptor->fetch_by_Slice_start_end_strand( $self, 1,undef, 1 )};
   }

  return $self->{'seq'};
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
  
  if ( $end < $start ) {
    throw("End coord is less then start coord");
  }

  $strand = 1 unless(defined $strand);

  if ( $strand != -1 && $strand != 1 ) {
    throw("Invalid strand [$strand] in call to Slice::subseq.");
  }
  my $subseq;
  if($self->adaptor){
    my $seqAdaptor = $self->adaptor->db->get_SequenceAdaptor();
    $subseq = ${$seqAdaptor->fetch_by_Slice_start_end_strand
      ( $self, $start,
        $end, $strand )};
  } else {
    $subseq = substr ($self->seq(), $start-1, $end - $start + 1);
    reverse_comp(\$subseq) if($strand == -1);
  }
  return $subseq;
}



=head2 get_base_count

  Arg [1]    : none
  Example    : $c_count = $slice->get_base_count->{'c'};
  Description: Retrieves a hashref containing the counts of each bases in the
               sequence spanned by this slice.  The format of the hash is :
               { 'a' => num,
                 'c' => num,
                 't' => num,
                 'g' => num,
                 'n' => num,
                 '%gc' => num }

               All bases which are not in the set [A,a,C,c,T,t,G,g] are
               included in the 'n' count.  The 'n' count could therefore be
               inclusive of ambiguity codes such as 'y'.
               The %gc is the ratio of GC to AT content as in:
               total(GC)/total(ACTG) * 100
               This function is conservative in its memory usage and scales to
               work for entire chromosomes.
  Returntype : hashref
  Exceptions : none
  Caller     : general

=cut

sub get_base_count {
  my $self = shift;

  my $a = 0; my $c = 0; my $t = 0; my $g = 0;

  my $start = 1;
  my $end;
  my $RANGE = 100_000;
  my $len = $self->length;
  my $seq;

  while($start <= $len) {
    $end = $start + $RANGE - 1;
    
    $end = $len if($end > $len);
    $seq = $self->subseq($start, $end);
    $a += $seq =~ tr/Aa/Aa/;
    $c += $seq =~ tr/Cc/Cc/;
    $t += $seq =~ tr/Tt/Tt/;
    $g += $seq =~ tr/Gg/Gg/;

    $start = $end + 1;
  }

  my $gc_content = 0;
  if($a || $g || $c || $t) {  #avoid divide by 0
    $gc_content = sprintf( "%1.2f", (($g + $c)/($a + $g + $t + $c)) * 100);
  }

  return {'a' => $a,
	  'c' => $c,
	  't' => $t,
	  'g' => $g,
	  'n' => $len - $a - $c - $t - $g,
	  '%gc' => $gc_content};
}



=head2 project

  Arg [1]    : string $name
               The name of the coordinate system to project this slice onto
  Arg [2]    : string $version
               The version of the coordinate system (such as 'NCBI34') to
               project this slice onto
  Example    :
    my $clone_projection = $slice->project('clone');

    foreach my $segment (@$clone_projection) {
      my ($start, $end, $clone) = @$segment;
      print $slice->seq_region_name, ':', $start, '-', $end , ' -> ',
            $clone->seq_region_name, ':', $clone->start, '-', $clone->end,
            $clone->strand, "\n";
    }
  Description: Returns the results of 'projecting' this slice onto another
               coordinate system.  Projecting to a coordinate system that
               the slice is assembled from is analagous to retrieving a tiling
               path.  This method may also be used to 'project up' to a higher
               level coordinate system, however.

               This method returns a listref of triplets [start,end,slice]
               which represents the projection.  The start and end defined the
               region of this slice which is made up of the third value of
               the triplet: a slice in the requested coordinate system.
  Returntype : list reference of [$start,$end,$slice] triplets
  Exceptions : none
  Caller     : general

=cut

sub project {
  my $self = shift;
  my $cs_name = shift;
  my $cs_version = shift;

  throw('Coord_system name argument is required') if(!$cs_name);

  #obtain a mapper between this coordinate system and the requested one
  if(!$self->adaptor){
    warning("You can't project a slice which isn't attached to a ".
            "database on to another coordinate system");
    return;
  }
  my $slice_adaptor = $self->adaptor();

  if(!$slice_adaptor) {
    warning("Cannot project slice without an attached SliceAdaptor");
    return [];
  }

  my $db = $slice_adaptor->db();
  my $csa = $db->get_CoordSystemAdaptor();
  my $cs = $csa->fetch_by_name($cs_name, $cs_version);
  my $slice_cs = $self->coord_system();

  if(!$cs) {
    throw("Cannot project to unknown coordinate system " .
          "[$cs_name $cs_version]");
  }

  #no mapping is needed if the requested coord system is the one we are in
  #but we do need to check if some of the slice is outside of defined regions
  if($slice_cs->equals($cs)) {
    my $entire_slice = $slice_adaptor->fetch_by_region($cs->name(),
                                                   $self->{'seq_region_name'},
                                                   undef, undef, undef,
                                                   $cs->version());

    my $entire_len = $entire_slice->length();

    #if the slice has negative coordinates or coordinates exceeding the 
    #exceeding length of the sequence region we want to shrink the slice to
    #the defined region

    if($self->{'start'} > $entire_len || $self->{'end'} < 1) {
      #none of this slice is in a defined region
      return [];
    }

    my $right_contract = 0;
    my $left_contract  = 0;
    if($self->{'end'} > $entire_len) {
      $right_contract = $entire_len - $self->{'end'};
    }
    if($self->{'start'} < 1) {
      $left_contract = $self->{'start'} - 1;
    }

    my $new_slice;
    if($left_contract || $right_contract) {
      $new_slice = $self->expand($left_contract, $right_contract);
    } else {
      $new_slice = $self;
    }

    return [bless [1-$left_contract, $self->length()+$right_contract, $new_slice], "Bio::EnsEMBL::ProjectionSegment" ];
  }

  my @projection;
  my $current_start = 1;

  #decompose this slice into its symlinked components.
  #this allows us to handle haplotypes and PARs
  my $normal_slice_proj = 
    $slice_adaptor->fetch_normalized_slice_projection($self);
  foreach my $segment (@$normal_slice_proj) {
    my $normal_slice = $segment->[2];

    $slice_cs = $normal_slice->coord_system();

    my $asma = $db->get_AssemblyMapperAdaptor();
    my $asm_mapper = $asma->fetch_by_CoordSystems($slice_cs, $cs);

    # perform the mapping between this slice and the requested system
    my @coords = $asm_mapper->map($normal_slice->seq_region_name(),
				  $normal_slice->start(),
				  $normal_slice->end(),
				  $normal_slice->strand(),
				  $slice_cs);

    #construct a projection from the mapping results and return it

    foreach my $coord (@coords) {
      my $coord_start  = $coord->start();
      my $coord_end    = $coord->end();
      my $length       = $coord_end - $coord_start + 1;

      #skip gaps
      if($coord->isa('Bio::EnsEMBL::Mapper::Coordinate')) {
        my $coord_cs     = $coord->coord_system();

        #create slices for the mapped-to coord system
        my $slice = $slice_adaptor->fetch_by_region($coord_cs->name(),
                                                    $coord->id(),
                                                    $coord_start,
                                                    $coord_end,
                                                    $coord->strand(),
                                                    $coord_cs->version());

        my $current_end = $current_start + $length - 1;
	
        push @projection, bless([$current_start, $current_end, $slice], "Bio::EnsEMBL::ProjectionSegment");
      }

      $current_start += $length;
    }
  }

  return \@projection;
}



=head2 expand

  Arg [1]    : (optional) int $five_prime_expand
               The number of basepairs to shift this slices five_prime
               coordinate by.  Positive values make the slice larger,
               negative make the slice smaller.
               coordinate left.
               Default = 0.
  Arg [2]    : (optional) int $three_prime_expand
               The number of basepairs to shift this slices three_prime
               coordinate by. Positive values make the slice larger,
               negative make the slice smaller.
               Default = 0.
  Example    : my $expanded_slice      = $slice->expand( 1000, 1000);
               my $contracted_slice    = $slice->expand(-1000,-1000);
               my $shifted_right_slice = $slice->expand(-1000, 1000);
               my $shifted_left_slice  = $slice->expand( 1000,-1000);
  Description: Returns a slice which is a resized copy of this slice.  The
               start and end are moved outwards from the center of the slice
               if positive values are provided and moved inwards if negative
               values are provided. This slice remains unchanged.  A slice
               may not be contracted below 1bp but may grow to be arbitrarily
               large.
  Returntype : Bio::EnsEMBL::Slice
  Exceptions : warning if an attempt is made to contract the slice below 1bp
  Caller     : general

=cut

sub expand {
  my $self = shift;
  my $five_prime_shift = shift || 0;
  my $three_prime_shift = shift || 0;
  
  if($self->{'seq'}){
    warning("You can't expand a slice which already has an attached ".
            "sequence ");
    return;
  }

  my $new_start;
  my $new_end;

  if($self->{'strand'} == 1) {
    $new_start = $self->{'start'} - $five_prime_shift;
    $new_end   = $self->{'end'} + $three_prime_shift;
  } else {
    $new_end = $self->{'end'} + $five_prime_shift;
    $new_start = $self->{'start'} - $three_prime_shift;
  }

  if($new_start > $new_end) {
    throw('Slice start cannot be greater than slice end');
  }

  #fastest way to copy a slice is to do a shallow hash copy
  my %new_slice = %$self;
  $new_slice{'start'} = $new_start;
  $new_slice{'end'}   = $new_end;

  return bless \%new_slice, ref($self);
}



=head2 get_seq_region_id

  Arg [1]    : none
  Example    : my $seq_region_id = $slice->get_seq_region_id();
  Description: Gets the internal identifier of the seq_region that this slice
               is on. Note that this function will not work correctly if this
               slice does not have an attached adaptor. Also note that it may
               be better to go through the SliceAdaptor::get_seq_region_id 
               method if you are working with multiple databases since is 
               possible to work with slices from databases with different
               internal seq_region identifiers.
  Returntype : int
  Exceptions : warning if slice is not associated with a SliceAdaptor
  Caller     : assembly loading scripts, general

=cut

sub get_seq_region_id {
  my ($self) = @_;
  
  if($self->adaptor) {
    return $self->adaptor->get_seq_region_id($self);
  } else {
    warning('Cannot retrieve seq_region_id without attached adaptor.');
    return 0;                                                            
  }
}


=head2 get_all_attribute_types

  Arg [1]    : none
  Example    : my @attribute_types = $slice->get_attribute_types();
  Description: Gets a list of the types of attributes associated with the
               seq_region this slice is on.
  Returntype : list of strings
  Exceptions : warning if slice does not have attached adaptor
  Caller     : general

=cut

sub get_attribute_types {
  my $self = shift;

  if(!$self->adaptor()) {
    warning('Cannot get attributes without an adaptor.');
    return ();
  }
  my %attrib_hash = %{$self->adaptor->get_seq_region_attribs($self)};
  return keys %attrib_hash;
}


=head2 get_attribute

  Arg [1]    : string $attrib_type
               The code of the attribute type to retrieve values for.
  Example    : ($htg_phase) = $slice->get_attribute('htg_phase');
               @synonyms    = $slice->get_attribute('synonyms');
  Description: Gets a list of values for a given attribute of this
               slice''s seq_region.
  Returntype : list of strings
  Exceptions : warning is slice does not have attached adaptor
               throw if argument is not provided
  Caller     : general

=cut

sub get_attribute {
  my $self = shift;
  my $attrib_code = shift;

  throw('Attrib code argument is required.') if(!$attrib_code);

  if(!$self->adaptor()) {
    warning('Cannot get attributes without an adaptor.');
    return ();
  }

  my %attribs = %{$self->adaptor->get_seq_region_attribs($self)};

  return @{$attribs{$attrib_code} || []};
}


=head2 get_all_PredictionTranscripts

  Arg [1]    : (optional) string $logic_name
               The name of the analysis used to generate the prediction
               transcripts obtained.
  Example    : @transcripts = @{$slice->get_all_PredictionTranscripts};
  Description: Retrieves the list of prediction transcripts which overlap
               this slice with logic_name $logic_name.  If logic_name is 
               not defined then all prediction transcripts are retrieved.
  Returntype : listref of Bio::EnsEMBL::PredictionTranscript
  Exceptions : none
  Caller     : none

=cut

sub get_all_PredictionTranscripts {
   my ($self,$logic_name) = @_;

   my $pta = $self->adaptor()->db()->get_PredictionTranscriptAdaptor();

   return $pta->fetch_all_by_Slice($self, $logic_name);
}



=head2 get_all_DnaAlignFeatures

  Arg [1]    : (optional) string $logic_name
               The name of the analysis performed on the dna align features
               to obtain.
  Arg [2]    : (optional) float $score
               The mimimum score of the features to retrieve
  Arg [3]    : (optional) string $dbtype
               The name of an attached database to retrieve the features from
               instead.  e.g. 'est'.
  Example    : @dna_dna_align_feats = @{$slice->get_all_DnaAlignFeatures};
  Description: Retrieves the DnaDnaAlignFeatures which overlap this slice with
               logic name $logic_name and with score above $score.  If 
               $logic_name is not defined features of all logic names are 
               retrieved.  If $score is not defined features of all scores are
               retrieved.
  Returntype : listref of Bio::EnsEMBL::DnaDnaAlignFeatures
  Exceptions : none
  Caller     : general

=cut

sub get_all_DnaAlignFeatures {
   my ($self, $logic_name, $score, $dbtype) = @_;

   my $db;

   if($dbtype) {
     $db = $self->adaptor->db->get_db_adaptor($dbtype);
     if(!$db) {
       warning("Don't have db $dbtype returning empty list\n");
       return [];
     }
   } else {
     $db = $self->adaptor->db;
   }

   my $dafa = $db->get_DnaAlignFeatureAdaptor();

   return $dafa->fetch_all_by_Slice_and_score($self,$score, $logic_name);
}



=head2 get_all_ProteinAlignFeatures

  Arg [1]    : (optional) string $logic_name
               The name of the analysis performed on the protein align features
               to obtain.
  Arg [2]    : (optional) float $score
               The mimimum score of the features to retrieve
  Arg [3]    : (optional) string $dbtype
               The name of an attached datbase to retrieve features from 
               instead.
  Example    : @dna_pep_align_feats = @{$slice->get_all_ProteinAlignFeatures};
  Description: Retrieves the DnaPepAlignFeatures which overlap this slice with
               logic name $logic_name and with score above $score.  If 
               $logic_name is not defined features of all logic names are 
               retrieved.  If $score is not defined features of all scores are
               retrieved.
  Returntype : listref of Bio::EnsEMBL::DnaPepAlignFeatures
  Exceptions : none
  Caller     : general

=cut

sub get_all_ProteinAlignFeatures {
  my ($self, $logic_name, $score, $dbtype) = @_;

  my $db;

  if($dbtype) {
    $db = $self->adaptor->db->get_db_adaptor($dbtype);
    if(!$db) {
      warning("Don't have db $dbtype returning empty list\n");
      return [];
    }
  } else {
    $db = $self->adaptor->db;
  }

  my $pafa = $db->get_ProteinAlignFeatureAdaptor();
  return $pafa->fetch_all_by_Slice_and_score($self, $score, $logic_name);
}



=head2 get_all_SimilarityFeatures

  Arg [1]    : (optional) string $logic_name
               the name of the analysis performed on the features to retrieve
  Arg [2]    : (optional) float $score
               the lower bound of the score of the features to be retrieved
  Example    : @feats = @{$slice->get_all_SimilarityFeatures};
  Description: Retrieves all dna_align_features and protein_align_features
               with analysis named $logic_name and with score above $score.
               It is probably faster to use get_all_ProteinAlignFeatures or
               get_all_DnaAlignFeatures if a sepcific feature type is desired.
               If $logic_name is not defined features of all logic names are 
               retrieved.  If $score is not defined features of all scores are
               retrieved.
  Returntype : listref of Bio::EnsEMBL::BaseAlignFeatures
  Exceptions : none
  Caller     : general

=cut

sub get_all_SimilarityFeatures {
  my ($self, $logic_name, $score) = @_;

  my @out = ();

  push @out, @{$self->get_all_ProteinAlignFeatures($logic_name, $score) };
  push @out, @{$self->get_all_DnaAlignFeatures($logic_name, $score) };

  return \@out;
}



=head2 get_all_SimpleFeatures

  Arg [1]    : (optional) string $logic_name
               The name of the analysis performed on the simple features
               to obtain.
  Arg [2]    : (optional) float $score
               The mimimum score of the features to retrieve
  Example    : @simple_feats = @{$slice->get_all_SimpleFeatures};
  Description: Retrieves the SimpleFeatures which overlap this slice with
               logic name $logic_name and with score above $score.  If 
               $logic_name is not defined features of all logic names are 
               retrieved.  If $score is not defined features of all scores are
               retrieved.
  Returntype : listref of Bio::EnsEMBL::SimpleFeatures
  Exceptions : none
  Caller     : general

=cut

sub get_all_SimpleFeatures {
  my ($self, $logic_name, $score) = @_;

  my $sfa = $self->adaptor()->db()->get_SimpleFeatureAdaptor();

  return $sfa->fetch_all_by_Slice_and_score($self, $score, $logic_name);
}



=head2 get_all_RepeatFeatures

  Arg [1]    : (optional) string $logic_name
               The name of the analysis performed on the repeat features
               to obtain.
  Example    : @repeat_feats = @{$slice->get_all_RepeatFeatures}
  Description: Retrieves the RepeatFeatures which overlap  with
               logic name $logic_name and with score above $score.  If 
               $logic_name is not defined features of all logic names are 
               retrieved.
  Returntype : listref of Bio::EnsEMBL::RepeatFeatures
  Exceptions : none
  Caller     : general

=cut

sub get_all_RepeatFeatures {
   my ($self, $logic_name) = @_;

   my $rpfa = $self->adaptor()->db()->get_RepeatFeatureAdaptor();

   return $rpfa->fetch_all_by_Slice($self, $logic_name);
}



=head2 get_all_SNPs

  Args      : none
  Function  : returns all SNPs on this slice. This function will only work
              correctly if the SNP database or the lite database has been
              attached to the core database.  This can been done through
              a call to DBAdaptor::add_db_adaptor.
  Returntype: listref of Bio::EnsEMBL::External::Variation
  Exceptions: none
  Caller    : contigview, snpview

=cut

sub get_all_SNPs {
  my $self = shift;

  my $snpa = $self->adaptor()->db()->get_SNPAdaptor();
  if( $snpa ) {
    return $snpa->fetch_all_by_Slice($self);
  } else {
    return [];
  }
}



=head2 get_all_Genes

  Arg [1]    : (optional) string $logic_name
               The name of the analysis used to generate the genes to retrieve 
  Arg [2]    : (optional) boolean $empty_flag 
  Example    : @genes = @{$slice->get_all_Genes};
  Description: Retrieves all genes that overlap this slice.  The empty flag is 
               used by the web code and is used to retrieve light weight genes
               that only have a start, end and strand (only works if lite db
               is available).  If the lite database has been attached to the
               core database this method will use the lite database (and 
               genes will not be as full featured).
  Returntype : listref of Bio::EnsEMBL::Genes
  Exceptions : none
  Caller     : none

=cut

sub get_all_Genes{
   my ($self, $logic_name, $empty_flag) = @_;

   #caching is performed on a per slice basis in the GeneAdaptor
   return $self->adaptor->db->get_GeneAdaptor->fetch_all_by_Slice($self,
								  $logic_name,
								  $empty_flag);
}



=head2 get_all_Genes_by_type


  Arg [1]    : string $type 
  Arg [2]    : (optional) string $logic_name
  Arg [3]    : (optional) boolean $empty_flag
  Example    : @genes = @{$slice->get_all_Genes_by_type($type, 
							'ensembl')};
  Description: Retrieves genes that overlap this slice of type $type.
               This is primarily used by the genebuilding code when several
               types of genes are used.

               The logic name is the analysis of the genes that are retrieved.
               If not provided all genes will be retrieved instead.

               The empty flag indicates light weight genes that only have a
               start, end and strand should be used (only works if lite db is
               available). If the lite database has
               been attached to the core database this method will use the
               lite database (and genes will not be as full featured).
  Returntype : listref of Bio::EnsEMBL::Genes
  Exceptions : none
  Caller     : genebuilder

=cut

sub get_all_Genes_by_type{
  my ($self, $type, $logic_name, $empty_flag) = @_;

  my @out = grep { $_->type eq $type } @{ $self->get_all_Genes($logic_name,
							       $empty_flag)};

  return \@out;
}


=head2 get_all_Transcripts

  Arg [1]    : none
  Example    : @transcripts = @{$slice->get_all_Transcripts};
  Description: Gets all transcripts which overlap this slice.  If you want to
               specify a particular analysis or type, then you are better off
               using get_all_Genes or get_all_Genes_by_type and iterating
               through the transcripts of each gene.
  Returntype : reference to a list of Bio::EnsEMBL::Transcripts
  Exceptions : none
  Caller     : general

=cut

sub get_all_Transcripts {
  my $self = shift;
  return $self->adaptor->db->get_TranscriptAdaptor->fetch_all_by_Slice($self);
}


=head2 get_all_Exons

  Arg [1]    : none
  Example    : @exons = @{$slice->get_all_Exons};
  Description: Gets all exons which overlap this slice.  Note that these exons
               will not be associated with any transcripts, so this may not
               be terribly useful.
  Returntype : reference to a list of Bio::EnsEMBL::Exons
  Exceptions : none
  Caller     : general

=cut

sub get_all_Exons {
  my $self = shift;
  return $self->adaptor->db->get_ExonAdaptor->fetch_all_by_Slice($self);
}



=head2 get_all_QtlFeatures

  Args       : none
  Example    : none
  Description: returns overlapping QtlFeatures
  Returntype : listref Bio::EnsEMBL::Map::QtlFeature
  Exceptions : none
  Caller     : general

=cut

sub get_all_QtlFeatures {
  my $self = shift;

  my $qfAdaptor;
  if( $self->adaptor()) {
    $qfAdaptor = $self->adaptor()->db()->get_QtlFeatureAdaptor();
  } else {
    return [];
  }

  return $qfAdaptor->fetch_all_by_Slice_constraint( $self );
}




=head2 get_all_KaryotypeBands

  Arg [1]    : none
  Example    : @kary_bands = @{$slice->get_all_KaryotypeBands};
  Description: Retrieves the karyotype bands which this slice overlaps.
  Returntype : listref oif Bio::EnsEMBL::KaryotypeBands
  Exceptions : none
  Caller     : general, contigview

=cut

sub get_all_KaryotypeBands {
  my ($self) = @_;

  my $kadp = $self->adaptor->db->get_KaryotypeBandAdaptor();
  return $kadp->fetch_all_by_Slice($self);
}




=head2 get_repeatmasked_seq

  Arg [1]    : listref of strings $logic_names (optional)
  Arg [2]    : int $soft_masking_enable (optional)
  Example    : $rm_slice = $slice->get_repeatmasked_seq()
               $softrm_slice->get_repeatmasked_seq(['RepeatMask'],1)
  Description: Returns Bio::EnsEMBL::Slice that can be used to create repeat
               masked sequence instead of the regular sequence.
               Sequence returned by this new slice will have repeat regions
               hardmasked by default (sequence replaced by N) or
               or soft-masked when arg[2] = 1 (sequence in lowercase)
               Will only work with database connection to get repeat features.
  Returntype : Bio::EnsEMBL::RepeatMaskedSlice
  Exceptions : none
  Caller     : general

=cut

sub get_repeatmasked_seq {
    my ($self,$logic_names,$soft_mask) = @_;

    #creat a slice that 
    return Bio::EnsEMBL::RepeatMaskedSlice->new
      (-START   => $self->{'start'},
       -END     => $self->{'end'},
       -STRAND  => $self->{'strand'},
       -ADAPTOR => $self->{'adaptor'},
       -SEQ     => $self->{'seq'},
       -SEQ_REGION_NAME => $self->{'seq_region_name'},
       -SEQ_REGION_LENGTH => $self->{'seq_region_length'},
       -COORD_SYSTEM    => $self->{'coord_system'},
       -REPEAT_MASK     => $logic_names,
       -SOFT_MASK       => $soft_mask);
}



=head2 _mask_features

  Arg [1]    : reference to a string $dnaref
  Arg [2]    : array_ref $repeats
               reference to a list Bio::EnsEMBL::RepeatFeature
               give the list of coordinates to replace with N or with 
               lower case
  Arg [3]    : int $soft_masking_enable (optional)
  Example    : none
  Description: replaces string positions described in the RepeatFeatures
               with Ns (default setting), or with the lower case equivalent 
               (soft masking).  The reference to a dna string which is passed
               is changed in place.
  Returntype : none
  Exceptions : none
  Caller     : seq

=cut

sub _mask_features {
  my ($self,$dnaref,$repeats,$soft_mask) = @_;

  # explicit CORE::length call, to avoid any confusion with the Slice
  # length method
  my $dnalen = CORE::length($$dnaref);

 REP:foreach my $f (@{$repeats}) {
    my $start  = $f->start;
    my $end    = $f->end;
    my $length = ($end - $start) + 1;

    # check if we get repeat completely outside of expected slice range
    if ($end < 1 || $start > $dnalen) {
	  warning("Unexpected: Repeat completely outside slice coordinates.");
      next REP;
    }

    # repeat partly outside slice range, so correct
    # the repeat start and length to the slice size if needed
    if ($start < 1) { 
      $start = 1;
      $length = ($end - $start) + 1;
    }

    # repeat partly outside slice range, so correct
    # the repeat end and length to the slice size if needed
    if ($end > $dnalen) {
      $end = $dnalen;
      $length = ($end - $start) + 1;
    }

    $start--;

    my $padstr;

    if ($soft_mask) {
      $padstr = lc substr ($$dnaref,$start,$length);
    } else {
      $padstr = 'N' x $length;
    }
    substr ($$dnaref,$start,$length) = $padstr;
  }
}


=head2 get_all_SearchFeatures

  Arg [1]    : scalar $ticket_ids
  Example    : $slice->get_all_SearchFeatures('BLA_KpUwwWi5gY');
  Description: Retreives all search features for stored blast
               results for the ticket that overlap this slice
  Returntype : listref of Bio::EnsEMBL::SeqFeatures
  Exceptions : none
  Caller     : general (webby!)

=cut

sub get_all_SearchFeatures {
    my $self = shift;
    my $ticket = shift;
    local $_;

    unless($ticket) {
      throw("ticket argument is required");
    }

    my $sfa = $self->adaptor()->db()->get_db_adaptor('blast');

    my $offset = $self->chr_start-1;

    my $features = $sfa ? $sfa->get_all_SearchFeatures($ticket, $self->chr_name, $self->chr_start, $self->chr_end) : [];

    foreach( @$features ) { 
      $_->start( $_->start-$offset );
      $_->end(   $_->end-$offset );
    };
    return $features;
}



=head2 get_all_MiscFeatures

  Arg [1]    : string $set (optional)
  Example    : $slice->get_all_MiscFeatures('cloneset');
  Description: Retreives all misc features which overlap this slice. If
               a set code is provided only features which are members of
               the requested set are returned.
  Returntype : listref of Bio::EnsEMBL::MiscFeatures
  Exceptions : none
  Caller     : general

=cut

sub get_all_MiscFeatures {
  my $self = shift;
  my $misc_set = shift;

  my $adaptor = $self->adaptor();

  if(!$adaptor) {
    warning('Cannot retrieve features without attached adaptor.');
    return [];
  }

  my $mfa = $adaptor->db->get_MiscFeatureAdaptor();

  if($misc_set) {
    return $mfa->fetch_all_by_Slice_and_set_code($self,$misc_set);
  }

  return $mfa->fetch_all_by_Slice($self);
}




=head2 get_all_MarkerFeatures

  Arg [1]    : (optional) string logic_name
               The logic name of the marker features to retrieve 
  Arg [2]    : (optional) int $priority 
               Lower (exclusive) priority bound of the markers to retrieve
  Arg [3]    : (optional) int $map_weight 
               Upper (exclusive) priority bound of the markers to retrieve
  Example    : my @markers = @{$slice->get_all_MarkerFeatures(undef,50, 2)};
  Description: Retrieves all markers which lie on this slice fulfilling the 
               specified map_weight and priority parameters (if supplied).
  Returntype : reference to a list of Bio::EnsEMBL::MarkerFeatures
  Exceptions : none
  Caller     : contigview, general

=cut

sub get_all_MarkerFeatures {
  my ($self, $logic_name, $priority, $map_weight) = @_;

  my $ma = $self->adaptor->db->get_MarkerFeatureAdaptor;

  my $feats = $ma->fetch_all_by_Slice_and_priority($self, 
					      $priority, 
					      $map_weight, 
					      $logic_name);
  return $feats;
}



=head2 get_all_compara_DnaAlignFeatures

  Arg [1]    : string $qy_species
               The name of the species to retrieve similarity features from
  Arg [2]    : string $qy_assembly
               The name of the assembly to retrieve similarity features from
  Arg [3]    : string $type
               The type of the alignment to retrieve similarity features from
  Example    : $fs = $slc->get_all_compara_DnaAlignFeatures('Mus musculus',
							    'MGSC3',
							    'WGA');
  Description: Retrieves a list of DNA-DNA Alignments to the species specified
               by the $qy_species argument.
               The compara database must be attached to the core database
               for this call to work correctly.  As well the compara database
               must have the core dbadaptors for both this species, and the
               query species added to function correctly.
  Returntype : reference to a list of Bio::EnsEMBL::DnaDnaAlignFeatures
  Exceptions : warning if compara database is not available
  Caller     : contigview

=cut

sub get_all_compara_DnaAlignFeatures {
  my ($self, $qy_species, $qy_assembly, $alignment_type) = @_;

  unless($qy_species && $qy_assembly && $alignment_type) {
    throw("Query species and assembly and alignmemt type arguments are required");
  }

  my $compara_db = $self->adaptor->db->get_db_adaptor('compara');

  unless($compara_db) {
    warning("Compara database must be attached to core database to " .
		"retrieve compara information");
    return [];
  }

  my $dafa = $compara_db->get_DnaAlignFeatureAdaptor;
  return $dafa->fetch_all_by_Slice($self, $qy_species, $qy_assembly, $alignment_type);
}

sub get_all_compara_Syntenies {
  my ($self, $qy_species ) = @_;

  unless($qy_species) {
    throw("Query species and assembly arguments are required");
  }

  my $compara_db = $self->adaptor->db->get_db_adaptor('compara');

  unless($compara_db) {
    warning("Compara database must be attached to core database to " .
		"retrieve compara information");
    return [];
  }

  my $sa = $compara_db->get_SyntenyAdaptor;
  $sa->setSpecies("XX",$self->adaptor->db->get_MetaContainer->get_Species->binomial, $qy_species );
  return $sa->get_synteny_for_chromosome( $self->seq_region_name,$self->start, $self->end );
}


=head2 get_all_Haplotypes

  Arg [1]    : (optional) boolean $lite_flag
               if true lightweight haplotype objects are used
  Example    : @haplotypes = $slice->get_all_Haplotypes;
  Description: Retrieves all of the haplotypes on this slice.  Only works
               if the haplotype adaptor has been attached to the core adaptor
               via $dba->add_db_adaptor('haplotype', $hdba); 
  Returntype : listref of Bio::EnsEMBL::External::Haplotype::Haplotypes
  Exceptions : warning is Haplotype database is not available
  Caller     : contigview, general

=cut

sub get_all_Haplotypes {
  my($self, $lite_flag) = @_;

  my $haplo_db = $self->adaptor->db->get_db_adaptor('haplotype');

  unless($haplo_db) {
    warning("Haplotype database must be attached to core database to " .
		"retrieve haplotype information" );
    return [];
  }

  my $haplo_adaptor = $haplo_db->get_HaplotypeAdaptor;

  my $haplotypes = $haplo_adaptor->fetch_all_by_Slice($self, $lite_flag);

  return $haplotypes;
}



=head2 get_all_DASFeatures

  Arg [1]    : none
  Example    : $features = $slice->get_all_DASFeatures;
  Description: Retreives a hash reference to a hash of DAS feature
               sets, keyed by the DNS, NOTE the values of this hash
               are an anonymous array containing:
                (1) a pointer to an array of features;
                (2) a pointer to the DAS stylesheet
  Returntype : hashref of Bio::SeqFeatures
  Exceptions : ?
  Caller     : webcode

=cut

sub get_all_DASFeatures{
   my ($self,@args) = @_;
  
  my %genomic_features =
      map { ( $_->_dsn => [ $_->fetch_all_by_Slice($self) ]  ) }
         $self->adaptor()->db()->_each_DASFeatureFactory;
  return \%genomic_features;

}


=head2 get_all_ExternalFeatures

  Arg [1]    : (optional) string $track_name
               If specified only features from ExternalFeatureAdaptors with 
               the track name $track_name are retrieved.  
               If not set, all features from every ExternalFeatureAdaptor are 
               retrieved.
  Example    : @x_features = @{$slice->get_all_ExternalFeatures}
  Description: Retrieves features on this slice from external feature adaptors 
  Returntype : listref of Bio::SeqFeatureI implementing objects in slice 
               coordinates 
  Exceptions : none
  Caller     : general

=cut

sub get_all_ExternalFeatures {
   my ($self, $track_name) = @_;

   my $features = [];

   my $xfa_hash = $self->adaptor->db->get_ExternalFeatureAdaptors;
   my @xf_adaptors = ();

   if($track_name) {
     #use a specific adaptor
     push @xf_adaptors, $xfa_hash->{$track_name};
   } else {
     #use all of the adaptors
     push @xf_adaptors, values %$xfa_hash;
   }


   foreach my $xfa (@xf_adaptors) {
     push @$features, @{$xfa->fetch_all_by_Slice($self)};
   }

   return $features;
}





# GENERIC FEATURES (See DBAdaptor.pm)

=head2 get_generic_features

  Arg [1]    : (optional) List of names of generic feature types to return.
               If no feature names are given, all generic features are
               returned.
  Example    : my %features = %{$slice->get_generic_features()};
  Description: Gets generic features via the generic feature adaptors that
               have been added via DBAdaptor->add_GenericFeatureAdaptor (if 
               any)
  Returntype : Hash of named features.
  Exceptions : none
  Caller     : none

=cut

sub get_generic_features() {

  my ($self, @names) = @_;

  my $db = $self->adaptor()->db();

  my %features = ();   # this will hold the results

  # get the adaptors for each feature
  my %adaptors = %{$db->get_GenericFeatureAdaptors(@names)};

  foreach my $adaptor_name (keys(%adaptors)) {
		
    my $adaptor_obj = $adaptors{$adaptor_name};
    # get the features and add them to the hash
    my $features_ref = $adaptor_obj->fetch_all_by_Slice($self);
    # add each feature to the hash to be returned
    foreach my $feature (@$features_ref) {
      $features{$adaptor_name} = $feature;
    }
  }

  return \%features;

}



# sub DEPRECATED METHODS #
###############################################################################

=head1 DEPRECATED METHODS

=cut



=head2 get_all_supercontig_Slices

  Description: DEPRECATED use get_tiling_path("NTcontig") instead

=cut


sub get_all_supercontig_Slices {
  my $self = shift;

  deprecate("Use get_tiling_path('NTcontig') instead");

  my $result = [];

  if( $self->adaptor() ) {
    my $superctg_names = 
      $self->adaptor()->list_overlapping_supercontigs( $self );

    for my $name ( @$superctg_names ) {
      my $slice;
      $slice = $self->adaptor()->fetch_by_supercontig_name( $name );
      $slice->name( $name );
      push( @$result, $slice );
    }
  } else {
    warning( "Slice needs to be attached to a database to get supercontigs" );
  }

  return $result;
}


=head2 get_all_Genes_by_source

  Description: DEPRECATED use get_all_Genes instead

=cut

sub get_all_Genes_by_source {
  deprecate("Use get_all_Genes() instead");
  return get_all_Genes(@_);
}


=head2 id

  Description: DEPRECATED use name() instead

=cut

sub id {
  deprecate("Use name() instead");
  name(@_);
}



=head2 display_id

  Description: DEPRECATED use name() instead

=cut


sub display_id {
  deprecate("Use name() instead");
  name(@_);
}


=head2 desc

  Description: DEPRECATED use name() instead

=cut

sub desc{
  deprecate("Use name() instead");
}

=head2 moltype

  Description: DEPRECATED do not use

=cut

sub moltype {
  deprecate("Do not use this method");
  return 'dna';
}

=head2 alphabet

  Description: DEPRECATED do not use

=cut

sub alphabet {
  deprecate("Do not use this method");
  return 'dna';
}

=head2 accession_number

  Description: DEPRECATED use name() instead

=cut

sub accession_number {
  deprecate("Use name() instead");
  name(@_);
}

=head2 get_Chromosome

  Description: DEPRECATED use this instead:
               $slice_adp->fetch_by_region('chromosome',
                                           $slice->seq_region_name)

=cut

sub get_Chromosome {
  my $self = shift @_;

  deprecate("Use SliceAdaptor::fetch_by_region('chromosome'," .
            '$slice->seq_region_name) instead');

  my $csa = $self->adaptor->db->get_CoordSystemAdaptor();
  my ($top_cs) = @{$csa->fetch_all()};

  my $chr = $self->adaptor->fetch_by_region($top_cs->name(),
                                            $self->seq_region_name(),
                                           undef,undef,undef,
                                            $top_cs->version());

  return bless($chr, 'Bio::EnsEMBL::Chromosome');
}



=head2 chr_name

  Description: DEPRECATED use seq_region_name() instead

=cut

sub chr_name{
  deprecate("Use seq_region_name() instead");
  seq_region_name(@_);
}



=head2 chr_start

  Description: DEPRECATED use start() instead

=cut

sub chr_start{
  deprecate('Use start() instead');
  start(@_);
}



=head2 chr_end

  Description: DEPRECATED use end() instead
  Returntype : int
  Exceptions : none
  Caller     : SliceAdaptor, general

=cut

sub chr_end{
  deprecate('Use end() instead');
  end(@_);
}


=head2 assembly_type

  Description: DEPRECATED use version instead

=cut

sub assembly_type{
  my $self = shift;
  deprecate('Use $slice->coord_system()->version() instead.');
  return $self->coord_system->version();
}


=head2 get_tiling_path

  Description: DEPRECATED use project instead

=cut

sub get_tiling_path {
  my $self = shift;

  deprecate('Use $slice->project("seqlevel") instead.');

  my $csa = $self->adaptor()->db()->get_CoordSystemAdaptor();

  #assume that they want the tiling path to the sequence coord system
  #this might not work well if this isn't a chromosomal slice
  my $projection = $self->project('seqlevel');

  my @tiling_path;

  foreach my $segment (@$projection) {

    my ($slice_start, $slice_end, $contig) = @$segment;
    my $contig_ori   = $contig->strand();
    my $contig_start = $contig->start();
    my $contig_end   = $contig->end();

    #the old get_tiling_path always gave back entire contigs in the forward
    #strand
    $contig = $contig->adaptor->fetch_by_region('seqlevel',
                                                $contig->seq_region_name());

    #bless contigs into contigs
    $contig = bless($contig, 'Bio::EnsEMBL::RawContig');


    push @tiling_path, Bio::EnsEMBL::Tile->new_fast($self,
                                                    $slice_start,
                                                    $slice_end,
                                                    $contig,
                                                    $contig_start,
                                                    $contig_end,
                                                    $contig_ori);
  }

  return \@tiling_path;
}


=head2 dbID

  Description: DEPRECATED use SliceAdaptor::get_seq_region_id instead

=cut

sub dbID {
  my $self = shift;
  deprecate('Use SliceAdaptor::get_seq_region_id instead.');
  if(!$self->adaptor) {
    warning('Cannot retrieve seq_region_id without attached adaptor.');
    return 0;
  }
  return $self->adaptor->get_seq_region_id($self);
}


=head2 get_all_MapFrags

  Description: DEPRECATED use get_all_MiscFeatures instead

=cut

sub get_all_MapFrags {
  my $self = shift;
  deprecate('Use get_all_MiscFeatures instead');
  return $self->get_all_MiscFeatures(@_);
}

=head2 has_MapSet

  Description: Not sure how this function is used, but it probably
               is not necessary.

=cut

sub has_MapSet {
  my( $self, $mapset_name ) = @_;
  deprecate('Not sure what this function was meant to do.');
  my $mfs = $self->get_all_MiscFeatures($mapset_name);
  return (@$mfs > 0);
}



1;
