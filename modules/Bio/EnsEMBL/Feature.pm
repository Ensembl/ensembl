# EnsEMBL module for Bio::EnsEMBL::Feature
#
# Copyright (c) 2003 EnsEMBL
#


=head1 NAME

Bio::EnsEMBL::Feature - Ensembl specific sequence feature.

=head1 SYNOPSIS

    my $feat = new Bio::EnsEMBL::Feature(-start   => 100,
                                       -end     => 220,
                                       -strand  => -1,
                                       -slice   => $slice
                                       -analysis => $analysis
                                      );

    my $start  = $feat->start;
    my $end    = $feat->end;
    my $strand = $feat->strand;

    #move the feature to the chromosomal coordinate system
    $feature = $feature->transform('chromosome');

    #move the feature to a different slice (possibly on another coord system)
    $feature = $feature->transfer($new_slice);

    #project the feature onto another coordinate system possibly across
    #boundaries:
    @projection = @{$feature->project('contig')};

    #change the start, end, and strand of the feature in place
    $feature->move($new_start, $new_end, $new_strand);

=head1 DESCRIPTION

This is the Base feature class from which all EnsEMBL features inherit.  It
provides a bare minimum functionality that all features require.  It basically
describes a location on a sequence of in an arbitrary coordinate system.

=head1 CONTACT

Post questions to the EnsEMBL development list: ensembl-dev@ebi.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


package Bio::EnsEMBL::Feature;

use Bio::EnsEMBL::Root;
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Utils::Exception qw(throw);

#This only inherits from root for backwards compatibility.
#In the future this inheritence will be removed
@ISA = qw(Bio::EnsEMBL::Root);


=head2 new

  Arg [-SLICE]: Bio::EnsEMBL::SLice - Represents the sequence that this
                feature is on. The coordinates of the created feature are
                relative to the start of the slice.
  Arg [-START]: The start coordinate of this feature relative to the start
                of the slice it is sitting on.  Coordinates start at 1 and
                are inclusive.
  Arg [-END]  : The end coordinate of this feature relative to the start of
                the slice it is sitting on.  Coordinates start at 1 and are
                inclusive.
  Arg [-STRAND]: The orientation of this feature.  Valid values are 1,-1,0.
  Example    : $feature = Bio::EnsEMBL::Feature->new(-start    => 1, 
                                                     -end      => 100,
                                                     -strand   => 1,
                                                     -slice    => $slice,
                                                     -analysis => $analysis);
  Description: Constructs a new Bio::EnsEMBL::Feature.  Generally subclasses
               of this method are instantiated, rather than this class itself.
  Returntype : Bio::EnsEMBL::Feature
  Exceptions : Thrown on invalid -SLICE, -ANALYSIS, -STRAND arguments
  Caller     : general, subclass constructors

=cut


sub new {
  my $caller = shift;

  my $class = ref($caller) || $caller;

  my($start, $end, $strand, $slice, $analysis) =
    rearrange(['START','END','STRAND','SLICE','ANALYSIS'], @_);

  if(defined($slice)) {
    if(!ref($slice) || !$slice->isa('Bio::EnsEMBL::Slice')) {
      throw('-SLICE argument must be a Bio::EnsEMBL::Slice');
    }
  }

  if(defined($analysis)) {
    if(!ref($analysis) || !$analysis->isa('Bio::EnsEMBL::Analysis')) {
      throw('-ANALYSIS argument must be a Bio::EnsEMBL::Analysis');
    }
  }

  if(defined($strand)) {
    if(!($strand == 1) && !($strand == -1) && !($strand == 0)) {
      throw('-STRAND argument must be 1, -1, or 0');
    }
  }

  if(defined($start) && defined($end)) {
    if($end < $start) {
      throw('-START argument must be less than -END argument');
    }
  }

  return bless({'start'    => $start,
                'end'      => $end,
                'strand'   => $strand,
                'slice'    => $slice,
                'analysis' => $analysis}, $class);
}



=head2 start

  Arg [1]    : (optional) int $start
               The start of this feature relative to the start of the slice
               that it is on.
  Example    : $start = $feat->start()
  Description: Getter/Setter for the start of this feature relative to the 
               start of the slice it is on.  Note that negative values, or
               values exceeding the length of the slice are permitted.
               Start must be less than or equal to the end regardless of the 
               strand. Coordinate values start at 1 and are inclusive.
  Returntype : int
  Exceptions : none
  Caller     : general

=cut

sub start {
  my $self = shift;
  $self->{'start'} = shift if(@_);
  return $self->{'start'};
}




=head2 end

  Arg [1]    : (optional) int $end
  Example    : $end = $feat->end();
  Description: Getter/Setter for the end of this feature relative to the
               start of the slice that it is on.  Note that negative values,
               of values exceeding the length of the slice are permitted.  End
               must be greater than or equal to start regardless of the strand.
               Coordinate values start at 1 and are inclusive.
  Returntype : int
  Exceptions : none
  Caller     : general

=cut

sub end {
  my $self = shift;
  $self->{'end'} = shift if(@_);
  return $self->{'end'};
}




=head2 strand

  Arg [1]    : (optional) int $strand
  Example    : $feat->strand(-1);
  Description: Getter/Setter for the strand of this feature relative to the
               slice it is on.  0 is an unknown or non-applicable strand.  
               -1 is the reverse (negative) strand and 1 is the forward 
               (positive) strand.  No other values are permitted.
  Returntype : int
  Exceptions : thrown if an invalid strand argument is passed
  Caller     : general

=cut

sub strand {
  my $self = shift;

  if(@_) {
    my $strand = shift || 0;
    if(defined($strand) && $strand != 0 && $strand != 1 && $strand != -1) {
      throw('strand argument must be 0, -1 or 1');
    }

    $self->{'strand'} = $strand;
  }

  return $self->{'strand'};
}



=head2 move

  Arg [1]    : int start
  Arg [2]    : int end
  Arg [3]    : (optional) int strand
  Example    : None
  Description: Sets the start, end and strand in one call rather than in 
               3 seperate calls to the start(), end() and strand() methods.
               This is for convenience and for speed when this needs to be
               done within a tight loop.
  Returntype : none
  Exceptions : Thrown is invalid arguments are provided
  Caller     : general

=cut

sub move {
  my $self = shift;

  throw('start and end arguments are required') if(@_ < 2);

  my $start  = shift;
  my $end    = shift;
  my $strand = shift;

  if(defined($start) && defined($end) && $end < $start) {
    throw('start must be less than or equal to end');
  }
  if(defined($strand) && $strand != 0 && $strand != -1 && $strand != 1) {
    throw('strand must be 0, -1 or 1');
  }

  $self->{'start'} = $start;
  $self->{'end'} = $end;
  $self->{'strand'} = $strand if(defined($strand));
}



=head2 length

  Arg [1]    : none
  Example    : $length = $feat->length();
  Description: Returns the length of this feature
  Returntype : int
  Exceptions : none
  Caller     : general

=cut

sub length {
  my $self = shift;
  return $self->{'start'} - $self->{'end'} + 1;
}




=head2 analysis

  Arg [1]    : (optional) Bio::EnsEMBL::Analysis $analysis
  Example    : $feature->analysis(new Bio::EnsEMBL::Analysis(...))
  Description: Getter/Setter for the analysis that is associated with 
               this feature.  The analysis describes how this feature 
               was derived.
  Returntype : Bio::EnsEMBL::Analysis
  Exceptions : thrown if an invalid argument is passed
  Caller     : general

=cut

sub analysis {
  my $self = shift;

  if(@_) {
    my $an = shift;
    if(defined($an) && (!ref($an) || !$an->isa('Bio::EnsEMBL::Analysis'))) {
      throw('analysis argument must be a Bio::EnsEMBL::Analysis');
    }
    $self->{'analysis'} = $an;
  }

  return $self->{'analysis'};
}



=head2 slice

  Arg [1]    : (optional) Bio::EnsEMBL::Slice $slice
  Example    : $seqname = $feature->slice()->name();
  Description: Getter/Setter for the Slice that is associated with this 
               feature.  The slice represents the underlying sequence that this
               feature is on.  Note that this method call is analagous to the
               old SeqFeature methods contig(), entire_seq(), attach_seq(),
               etc.
  Returntype : Bio::EnsEMBL::Slice
  Exceptions : thrown if an invalid argument is passed
  Caller     : general

=cut

sub slice {
  my $self = shift;

  if(@_) {
    my $sl = shift;
    if(defined($sl) && (!ref($sl) || !$sl->isa('Bio::EnsEMBL::Slice'))) {
      throw('slice argument must be a Bio::EnsEMBL::Slice');
    }

    $self->{'slice'} = $sl;
  }

  return $self->{'slice'};
}


=head2 transform

  Arg [1]    : string $coord_system
               The coord system to transform this feature to.
  Arg [2]    : string $version (optional)
               The version of the coord system to transform this feature to.
  Example    : $feature = $feature->transform('contig');
               next if(!defined($feature));
  Description: Returns a copy of this feature, but converted to a different
               coordinate system. The converted feature will be placed on a
               slice which spans an entire sequence region of the new
               coordinate system.  If The requested coordinate system is the
               same coordinate system (but may be ona different slice if the
               original slice did not span the entire seq_region).

               If a feature spans a boundary in the new coordinate system, 
               undef is returned instead.

               For example, transforming an exon in contig coordinates to one 
               in chromosomal coodinates will place the exon on a slice of an 
               entire chromosome.
  Returntype : Bio::EnsEMBL::Feature (or undef)
  Exceptions : thrown if an invalid coordinate system is provided
               thrown if this feature is not already on a slice
  Caller     : general

=cut

sub transform {
  my $self = shift;
  my $cs_name = shift;
  my $cs_version = shift;

  ### for backwards compatibility it might be a good idea to transform to
  ### contig coords if no args are provided and to chain this method to
  ### transfer() if a slice is provided

  my $slice = $self->{'slice'};

  if(!$slice) {
    throw('Feature is not associated with a slice and may not be transformed');
  }

  my $db = $self->adaptor->db();
  my $cs = $db->get_CoordSystemAdaptor->fetch_by_name($cs_name, $cs_version);
  my $current_cs = $slice->coord_system();

  #convert this features coords to absolute coords (i.e. relative to the start
  #of the seq_region, not to the slice)
  my $slice_adaptor = $db->get_SliceAdaptor;
  $slice = $slice_adaptor->fetch_by_region($current_cs->name(),
                                           $slice->seq_region(),
                                           undef, #start
                                           undef, #end
                                           1, #strand
                                           $current_cs->version);

  #this also copies the feature so we no longer have to, and can alter it
  #in place
  my $feature = $self->transfer($slice);

  #if the current coord system is the same as the new one, we are already done
  return $feature if($cs->equals($current_cs));

  #otherwise convert this features coordinates to the other coord system
  my $asma = $db->get_AssemblyMapperAdaptor();
  my $asm_mapper = $asma->fetch_by_CoordSystems($current_cs, $cs);

  #remember the coords are already relative to start of seq_region
  my @coords = $asm_mapper->map($slice->seq_region_name(),
                                $feature->{'start'},
                                $feature->{'end'},
                                $feature->{'strand'},
                                $current_cs);

  #we don't deal with gaps or with multiple seq_region mappings
  if(@coords != 1 || $coords[0]->isa('Bio::EnsEMBL::Mapper::Gap')) {
    return undef;
  }
  my $coord = $coords[0];

  $slice = $slice_adaptor->fetch_by_region($cs->name(),
                                           $coord->id(), #seq_region
                                           undef,        #start
                                           undef,        #end
                                           1,            #strand
                                           $cs->version());

  #shift the feature and place it on the appropriate slice
  $feature->{'start'}  = $coord->start();
  $feature->{'end'}    = $coord->end();
  $feature->{'strand'} = $coord->strand();
  $feature->{'slice'}  = $slice;

  return $feature;
}


sub transfer {
  

}




sub project {
  my $self = shift;

}



1;
