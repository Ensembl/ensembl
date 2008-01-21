# Ensembl module for Bio::EnsEMBL::Feature
#
# Copyright (c) 2003 Ensembl
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

This is the Base feature class from which all Ensembl features inherit.  It
provides a bare minimum functionality that all features require.  It basically
describes a location on a sequence in an arbitrary coordinate system.

=head1 CONTACT

Post questions to the Ensembl development list: ensembl-dev@ebi.ac.uk

=head1 METHODS

=cut


use strict;
use warnings;

package Bio::EnsEMBL::Feature;

use Bio::EnsEMBL::Storable;
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Utils::Exception qw(throw deprecate warning);
use Bio::EnsEMBL::Slice;
use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Storable);


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
  Arg [-SEQNAME] : A seqname to be used instead of the default name of the 
                of the slice.  Useful for features that do not have an 
                attached slice such as protein features.
  Arg [-dbID]   : (optional) internal database id
  Arg [-ADAPTOR]: (optional) Bio::EnsEMBL::DBSQL::BaseAdaptor
  Example    : $feature = Bio::EnsEMBL::Feature->new(-start    => 1, 
                                                     -end      => 100,
                                                     -strand   => 1,
                                                     -slice    => $slice,
                                                     -analysis => $analysis);
  Description: Constructs a new Bio::EnsEMBL::Feature.  Generally subclasses
               of this method are instantiated, rather than this class itself.
  Returntype : Bio::EnsEMBL::Feature
  Exceptions : Thrown on invalid -SLICE, -ANALYSIS, -STRAND ,-ADAPTOR arguments
  Caller     : general, subclass constructors
  Status     : Stable

=cut


sub new {
  my $caller = shift;

  my $class = ref($caller) || $caller;

  my($start, $end, $strand, $slice, $analysis,$seqname, $dbID, $adaptor) =
    rearrange(['START','END','STRAND','SLICE','ANALYSIS', 'SEQNAME',
               'DBID', 'ADAPTOR'], @_);

  if($slice) {
    if(!ref($slice) || !$slice->isa('Bio::EnsEMBL::Slice')) {
      throw('-SLICE argument must be a Bio::EnsEMBL::Slice not '.$slice);
    }
  }

  if($analysis) {
    if(!ref($analysis) || !$analysis->isa('Bio::EnsEMBL::Analysis')) {
      throw('-ANALYSIS argument must be a Bio::EnsEMBL::Analysis not '.
            $analysis);
    }
  }

  if(defined($strand)) {
    if(!($strand == 1) && !($strand == -1) && !($strand == 0)) {
      throw('-STRAND argument must be 1, -1, or 0');
    }
  }

  if(defined($start) && defined($end)) {
    if($end+1 < $start) {
      throw('Start must be less than or equal to end+1');
    }
  }

  return bless({'start'    => $start,
                'end'      => $end,
                'strand'   => $strand,
                'slice'    => $slice,
                'analysis' => $analysis,
                'adaptor'  => $adaptor,
                'seqname'  => $seqname,
                'dbID'     => $dbID}, $class);
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
  Status     : Stable

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
  Status     : Stable

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
  Status     : Stable

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
  Status     : Stable

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
  Status     : Stable

=cut

sub length {
  my $self = shift;
  return $self->{'end'} - $self->{'start'} + 1;
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
  Status     : Stable

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
  Status     : Stable

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
               coordinate system. If the requested coordinate system is the
               same coordinate system it is simply placed on a slice which
               spans the entire seq_region (as opposed to the original slice
               which may have only partially covered the seq_region).

               If a feature spans a boundary in the new coordinate system,
               undef is returned instead.

               For example, transforming an exon in contig coordinates to one 
               in chromosomal coodinates will place the exon on a slice of an 
               entire chromosome.
  Returntype : Bio::EnsEMBL::Feature (or undef)
  Exceptions : thrown if an invalid coordinate system is provided
               warning if Feature is not attached to a slice
  Caller     : general, transfer()
  Status     : Stable

=cut

sub transform {
  my $self = shift;
  my $cs_name = shift;
  my $cs_version = shift;
  #
  # For backwards compatibility check if the arguments are old style args
  #
  if(!$cs_name || ref($cs_name)) {
    deprecate('Calling transform without a coord system name is deprecated.');
    return $self->_deprecated_transform($cs_name);
  }

  my $slice = $self->{'slice'};

  if(!$slice) {
    warning("Feature cannot be transformed without attached slice.");
    return undef;
  }

  if(!$slice->adaptor()) {
    warning("Feature cannot be transformed without adaptor on" .
            " attached slice.");
    return undef;
  }

  #use db from slice since this feature may not yet be stored in a database
  my $db = $slice->adaptor->db();
  my $cs = $db->get_CoordSystemAdaptor->fetch_by_name($cs_name, $cs_version);
  my $current_cs = $slice->coord_system();

  if(!$current_cs) {
    warning("Feature cannot be transformed without CoordSystem on " .
            "attached slice.");
    return undef;
  }

  if(!$cs) {
    throw("Cannot transform to unknown coordinate system " .
          "[$cs_name $cs_version]\n");
  }

  # if feature is already in the requested coordinate system, we can just
  # return a copy
  if( $cs->equals( $current_cs ) && $slice->start() == 1 &&
      $slice->strand() == 1 ) {
    my $new_feature;
    %$new_feature = %$self;
    bless $new_feature, ref $self;
    return $new_feature;
  }

  my $projection = $self->project( $cs_name, $cs_version );

  if( @$projection != 1 ) {
    return undef;
  } else {
    my $p_slice = $projection->[0]->[2];
    my $slice_adaptor = $db->get_SliceAdaptor;
    $slice = $slice_adaptor->fetch_by_region($p_slice->coord_system()->name(),
					     $p_slice->seq_region_name(),
					     undef, #start
					     undef, #end
					     1, #strand
					     $p_slice->coord_system()->version);

    my $new_feature;
    %$new_feature = %$self;
    bless $new_feature, ref $self;
    $new_feature->{'start'}  = $p_slice->start();
    $new_feature->{'end'}    = $p_slice->end();
    $new_feature->{'strand'} =
      ($self->{'strand'} == 0) ? 0 : $p_slice->strand();
    $new_feature->{'slice'}  = $slice;
    return $new_feature;
  }
}



=head2 transfer

  Arg [1]    : Bio::EnsEMBL::Slice $slice
               The slice to transfer this feature to
  Example    : $feature = $feature->transfer($slice);
               next if(!defined($feature));
  Description: Returns a copy of this feature which has been shifted onto
               another slice.

               If the new slice is in a different coordinate system the
               feature is transformed first and then placed on the slice.
               If the feature would be split across a coordinate system
               boundary or mapped to a gap undef is returned instead.

               If the feature cannot be placed on the provided slice because
               it maps to an entirely different location, undef is returned
               instead.

  Returntype : Bio::EnsEMBL::Feature (or undef)
  Exceptions : throw on incorrect argument
               throw if feature does not have attached slice
  Caller     : general, transform()
  Status     : Stable

=cut


sub transfer {
  my $self = shift;
  my $slice = shift;

  if(!$slice || !ref($slice) || !$slice->isa('Bio::EnsEMBL::Slice')) {
    throw('Slice argument is required');
  }

  #make a shallow copy of the feature to be transfered
  my $feature;
  %{$feature} = %{$self};
  bless $feature, ref($self);

  my $current_slice = $self->{'slice'};

  if(!$current_slice) {
    warning("Feature cannot be transfered without attached slice.");
    return undef;
  }

  my $cur_cs = $current_slice->coord_system();
  my $dest_cs = $slice->coord_system();

  #if we are not in the same coord system a transformation step is needed first
  if(!$dest_cs->equals($cur_cs)) {
    $feature = $feature->transform($dest_cs->name, $dest_cs->version);
    return undef if(!defined($feature));
    $current_slice = $feature->{'slice'};
  }

  # feature went to entirely different seq_region
  if($current_slice->seq_region_name() ne $slice->seq_region_name()) {
    return undef;
  }

  #if the current feature positions are not relative to the start of the
  #seq region, convert them so they are
  my $cur_slice_start  = $current_slice->start();
  my $cur_slice_strand = $current_slice->strand();
  if($cur_slice_start != 1 || $cur_slice_strand != 1) {
    my $fstart = $feature->{'start'};
    my $fend   = $feature->{'end'};

    if($cur_slice_strand == 1) {
      $feature->{'start'} = $fstart + $cur_slice_start - 1;
      $feature->{'end'}   = $fend   + $cur_slice_start - 1;
    } else {
      my $cur_slice_end = $current_slice->end();
      $feature->{'start'}  = $cur_slice_end - $fend   + 1;
      $feature->{'end'}    = $cur_slice_end - $fstart + 1;
      $feature->{'strand'} *= -1;
    }
  }

  my $fstart = $feature->{'start'};
  my $fend   = $feature->{'end'};

  #convert to destination slice coords
  if($slice->strand == 1) {
    $feature->{'start'} = $fstart - $slice->start() + 1;
    $feature->{'end'}   = $fend   - $slice->start() + 1;
  } else {
    $feature->{'start'} = $slice->end() - $fend   + 1;
    $feature->{'end'}   = $slice->end() - $fstart + 1;
    $feature->{'strand'} *= -1;
  }

  $feature->{'slice'} = $slice;

  return $feature;
}

=head2 project_to_slice

  Arg [1]    : slice to project to


  Example    :
    my $clone_projection = $feature->project_to_slice($slice);

    foreach my $seg (@$clone_projection) {
      my $clone = $seg->to_Slice();
      print "Features current coords ", $seg->from_start, '-',
        $seg->from_end, " project onto clone coords " .
        $clone->seq_region_name, ':', $clone->start, '-', $clone->end,
        $clone->strand, "\n";
    }
  Description: Returns the results of 'projecting' this feature onto another
               slcie . This is useful to see where a feature
               would lie in a coordinate system in which it
               crosses a boundary.

               This method returns a reference to a list of
               Bio::EnsEMBL::ProjectionSegment objects.
               ProjectionSegments are blessed arrays and can also be used as
               triplets [from_start,from_end,to_Slice]. The from_start and
               from_end are the coordinates relative to the feature start.
               For example, if a feature is current 100-200bp on a slice
               then the triplets returned might be:
               [1,50,$slice1],
               [51,101,$slice2]

               The to_Slice is a slice spanning the region on the requested
               coordinate system that this feature projected to.

               If the feature projects entirely into a gap then a reference to
               an empty list is returned.

  Returntype : list reference of Bio::EnsEMBL::ProjectionSegments
               which can also be used as [$start,$end,$slice] triplets
  Exceptions : slice does not have an adaptor
  Caller     : general
  Status     : At Risk

=cut

sub project_to_slice {
  my $self = shift;
  my $to_slice = shift;
  my $slice = $self->{'slice'};

  if(!$slice) {
    warning("Feature cannot be projected without attached slice.");
    return [];
  }


  #get an adaptor from the attached slice because this feature may not yet
  #be stored and may not have its own adaptor
  my $slice_adaptor = $slice->adaptor();

  if(!$slice_adaptor) {
    throw("Cannot project feature because associated slice does not have an " .
          " adaptor");
  }

  my $strand = $self->strand() * $slice->strand();
  #fetch by feature always gives back forward strand slice:
  $slice = $slice_adaptor->fetch_by_Feature($self);
  $slice = $slice->invert if($strand == -1);
  return $slice->project_to_slice($to_slice);
}


=head2 project

  Arg [1]    : string $name
               The name of the coordinate system to project this feature onto
  Arg [2]    : string $version (optional)
               The version of the coordinate system (such as 'NCBI34') to
               project this feature onto
  Example    :
    my $clone_projection = $feature->project('clone');

    foreach my $seg (@$clone_projection) {
      my $clone = $seg->to_Slice();
      print "Features current coords ", $seg->from_start, '-',
        $seg->from_end, " project onto clone coords " .
        $clone->seq_region_name, ':', $clone->start, '-', $clone->end,
        $clone->strand, "\n";
    }
  Description: Returns the results of 'projecting' this feature onto another
               coordinate system.  This is useful to see where a feature
               would lie in a coordinate system in which it
               crosses a boundary.

               This method returns a reference to a list of
               Bio::EnsEMBL::ProjectionSegment objects.
               ProjectionSegments are blessed arrays and can also be used as
               triplets [from_start,from_end,to_Slice]. The from_start and
               from_end are the coordinates relative to the feature start.
               For example, if a feature is current 100-200bp on a slice
               then the triplets returned might be:
               [1,50,$slice1],
               [51,101,$slice2]

               The to_Slice is a slice spanning the region on the requested
               coordinate system that this feature projected to.

               If the feature projects entirely into a gap then a reference to
               an empty list is returned.

  Returntype : list reference of Bio::EnsEMBL::ProjectionSegments
               which can also be used as [$start,$end,$slice] triplets
  Exceptions : slice does not have an adaptor
  Caller     : general
  Status     : Stable

=cut

sub project {
  my $self = shift;
  my $cs_name = shift;
  my $cs_version = shift;

  my $slice = $self->{'slice'};

  if(!$slice) {
    warning("Feature cannot be projected without attached slice.");
    return [];
  }


  #get an adaptor from the attached slice because this feature may not yet
  #be stored and may not have its own adaptor
  my $slice_adaptor = $slice->adaptor();

  if(!$slice_adaptor) {
    throw("Cannot project feature because associated slice does not have an " .
          " adaptor");
  }

  my $strand = $self->strand() * $slice->strand();
  #fetch by feature always gives back forward strand slice:
  $slice = $slice_adaptor->fetch_by_Feature($self);
  $slice = $slice->invert if($strand == -1);
  return $slice->project($cs_name, $cs_version);
}



=head2 seqname

  Arg [1]    : (optional) $seqname
  Example    : $seqname = $feat->seqname();
  Description: Getter/Setter for the name of the sequence that this feature
               is on. Normally you can get away with not setting this value
               and it will default to the name of the slice on which this
               feature is on.  It is useful to set this value on features which
               do not ordinarily sit on features such as ProteinFeatures which
               sit on peptides.
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub seqname {
  my $self = shift;
  
  if(@_) {
    $self->{'seqname'} = shift;
  }

  if(!$self->{'seqname'} && $self->slice()) {
    return $self->slice->name();
  }

  return $self->{'seqname'};
}




=head2 display_id

  Arg [1]    : none
  Example    : print $f->display_id();
  Description: This method returns a string that is considered to be
               the 'display' identifier.  It is overridden by subclasses to
               return an appropriate value for objects of that particular 
               class.  If no appropriate display id is available an empty
               string is returned instead.
  Returntype : string
  Exceptions : none
  Caller     : web drawing code
  Status     : Stable

=cut

sub display_id {
  my $self = shift;
  return '';
}


=head2 feature_Slice

  Args       : none
  Example    : $slice = $feature->feature_Slice()
  Description: This is a convenience method to return a slice that covers the
               Area of this feature. The feature start will be at 1 on it, and
               it will have the length of this feature.
  Returntype : Bio::EnsEMBL::Slice or undef if this feature has no attached
               Slice.
  Exceptions : warning if Feature does not have attached slice.
  Caller     : web drawing code
  Status     : Stable

=cut

sub feature_Slice {
  my $self = shift;

  my $slice = $self->slice();

  if(!$slice) {
    warning('Cannot obtain Feature_Slice for feature without attached slice');
    return undef;
  }

  return Bio::EnsEMBL::Slice->new
    (-seq_region_name   => $slice->seq_region_name,
     -seq_region_length => $slice->seq_region_length,
     -coord_system      => $slice->coord_system,
     -start             => $self->seq_region_start(),
     -end               => $self->seq_region_end(),
     -strand            => $self->seq_region_strand(),
     -adaptor           => $slice->adaptor());

  
}


=head2 seq_region_name

  Arg [1]    : none
  Example    : print $feature->seq_region_name();
  Description: Gets the name of the seq_region which this feature is on.
               Returns undef if this Feature is not on a slice.
  Returntype : string or undef
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub seq_region_name {
  my $self = shift;
  my $slice = $self->{'slice'};

  return ($slice) ? $slice->seq_region_name() : undef;
}


=head2 seq_region_length

  Arg [1]    : none
  Example    : print $feature->seq_region_length();
  Description: Returns the length of the seq_region which this feature is on 
               Returns undef if this Feature is not on a slice.
  Returntype : unsigned int or undef
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut


sub seq_region_length {
  my $self = shift;
  my $slice = $self->{'slice'};

  return ($slice) ? $slice->seq_region_length() : undef;
}


=head2 seq_region_strand

  Arg [1]    : none
  Example    : print $feature->seq_region_strand();
  Description: Returns the strand of the seq_region which this feature is on 
               (i.e. feature_strand * slice_strand)
               Returns undef if this Feature is not on a slice.
  Returntype : 1,0,-1 or undef
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut


sub seq_region_strand {
  my $self = shift;
  my $slice = $self->{'slice'};

  return ($slice) ? $slice->strand() * $self->{'strand'} : undef;
}


=head2 seq_region_start

  Arg [1]    : none
  Example    : print $feature->seq_region_start();
  Description: Convenience method which returns the absolute start of this
               feature on the seq_region, as opposed to the relative (slice) 
               position.

               Returns undef if this feature is not on a slice.
  Returntype : int or undef
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub seq_region_start {
  my $self = shift;
  my $slice = $self->{'slice'};

  return undef if(!$slice);

  if($slice->strand == 1) {
    return undef if(!defined($self->{'start'}));
    return $slice->start() + $self->{'start'} - 1;
  } else {
    return undef if(!defined($self->{'end'}));
    return $slice->end() - $self->{'end'} + 1;
  }
}


=head2 seq_region_end

  Arg [1]    : none
  Example    : print $feature->seq_region_end();
  Description: Convenience method which returns the absolute end of this
               feature on the seq_region, as opposed to the relative (slice)
               position.

               Returns undef if this feature is not on a slice.
  Returntype : int or undef
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub seq_region_end {
  my $self = shift;
  my $slice = $self->{'slice'};

  return undef if(!$slice);

  if($slice->strand == 1) {
    return undef if(!defined($self->{'end'}));
    return $slice->start() + $self->{'end'} - 1;
  } else {
    return undef if(!defined($self->{'start'}));
    return $slice->end() - $self->{'start'} + 1;
  }
}


=head2 coord_system_name

  Arg [1]    : none
  Example    : print $feature->coord_system_name()
  Description: Gets the name of the coord_system which this feature is on.
               Returns undef if this Feature is not on a slice.
  Returntype : string or undef
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub coord_system_name {
  my $self = shift;
  my $slice = $self->{'slice'};
  return ($slice) ? $slice->coord_system_name() : undef;
}


=head2 seq

  Args       : none
  Example    : my $dna_sequence = $simple_feature->seq();
  Description: Returns the dna sequence from the attached slice and 
               attached database that overlaps with this feature.
               Returns undef if there is no slice or no database.
               Returns undef if this feature is unstranded (i.e. strand=0).
  Returntype : undef or string
  Exceptions : warning if this feature is not stranded
  Caller     : general
  Status     : Stable

=cut


sub seq {
  my $self = shift;

  if( ! defined $self->{'slice'} ) {
    return undef;
  }

  if(!$self->strand()) {
    warning("Cannot retrieve sequence for unstranded feature.");
    return undef;
  }

  return $self->{'slice'}->subseq($self->start(), $self->end(),
                                  $self->strand());

}




=head2 get_all_alt_locations

  Arg [1]    : none
  Example    : @features = @{$feature->get_all_alt_locations()};
               foreach $f (@features) {
                 print $f->slice->seq_region_name,' ',$f->start, $f->end,"\n";
               }

  Description: Retrieves shallow copies of this feature in its alternate
               locations.  A feature can be considered to have multiple
               locations when it sits on a alternative structural haplotype
               or when it is on a pseudo autosomal region.  Most features will
               just return a reference to an empty list though.
               The features returned by this method will be on a slice which
               covers the entire alternate region.

               Currently this method does not take into account alternate
               locations on the alternate locations (e.g. a reference
               sequence may have multiple alternate haplotypes.  Asking
               for alternate locations of a feature on one of the alternate
               haplotypes will give you back the reference location, but not
               locations on the other alternate haplotypes).

  Returntype : reference to list of features of the same type of this feature.
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_all_alt_locations {
  my $self = shift;
  my $return_all = shift || 0;

  my $slice = $self->{'slice'} or return [];
  my $sa = $slice->adaptor() or return [];

  # get slice of entire region
  $slice = $sa->fetch_by_seq_region_id($slice->get_seq_region_id);

  my $axfa = $sa->db->get_AssemblyExceptionFeatureAdaptor();
  my $axfs = $axfa->fetch_all_by_Slice($slice);

  my (@haps, @alt);

  foreach my $axf (@$axfs) {
    if(uc($axf->type()) eq 'HAP') {
      push @haps, $axf;
    } elsif(uc($axf->type()) =~ 'PAR') {
      push @alt, $axf;
    } elsif( $axf->type() eq "HAP REF" ) {
      push @haps, $axf if $return_all > 0 ;
      # do nothing when you are on REF
    } else {
      warning("Unknown exception feature type ". $axf->type()."- ignoring.");
    }
  }

  # regions surrounding hap are those of interest, not hap itself
  # convert hap alt. exc. features to regions around haps instead
  foreach my $h (@haps) {
    my $haslice = $h->alternate_slice();
    my $hacs    = $haslice->coord_system();

    if($h->start() > 1 && $haslice->start() > 1) {
      my $aslice = $sa->fetch_by_region($hacs->name(),
                                        $haslice->seq_region_name(),
                                        1,
                                        $haslice->start()-1,
                                        $haslice->strand(),
                                        $hacs->version());

      push @alt, Bio::EnsEMBL::AssemblyExceptionFeature->new
        (-start  => 1,
         -end    => $h->start()-1,
         -alternate_slice => $aslice);
    }

    if($h->end() < $slice->seq_region_length() &&
       $haslice->end < $haslice->seq_region_length()) {
      my $aslice = $sa->fetch_by_region($hacs->name(),
                                        $haslice->seq_region_name(),
                                        $haslice->end()+1,
                                        $haslice->seq_region_length(),
                                        $haslice->strand(),
                                        $hacs->version());

      push @alt, Bio::EnsEMBL::AssemblyExceptionFeature->new
        (-start  => $h->end() + 1,
         -end    => $slice->seq_region_length(),
         -alternate_slice => $aslice);
    }
  }


  # check if exception regions contain our feature

  my @features;

  foreach my $axf (@alt) {
    # ignore other region if feature is not entirely on it
    next if($self->seq_region_start() < $axf->start() ||
            $self->seq_region_end()   > $axf->end());

    # quick shallow copy of the feature
    my $f;
    %$f = %$self;
    bless $f, ref($self);

    my $aslice = $axf->alternate_slice();

    # position feature on entire slice of other region
    $f->{'start'}  = $f->seq_region_start() - $axf->start() + $aslice->start();
    $f->{'end'}    = $f->seq_region_end()   - $axf->start() + $aslice->start();
    $f->{'strand'} *= $aslice->strand();

    $f->{'slice'} = $sa->fetch_by_seq_region_id($aslice->get_seq_region_id());

    push @features, $f;
  }

  return \@features;
}


=head2 overlaps

  Arg [1]    : Bio::EnsEMBL::Feature $f
               The other feature you want to check overlap with this feature
               for.
  Description: This method does a range comparison of this features start and
               end and compares it with another features start and end. It will
               return true if these ranges overlap and the features are on the
               same seq_region.
  Returntype : TRUE if features overlap, FALSE if they don't
  Exceptions : warning if features are on different seq_regions
  Caller     : general
  Status     : Stable

=cut

sub overlaps {
  my $self = shift;
  my $f = shift;

  my $sr1_name = $self->seq_region_name;
  my $sr2_name = $f->seq_region_name;

  if ($sr1_name and $sr2_name and ($sr1_name ne $sr2_name)) {
    warning("Bio::EnsEMBL::Feature->overlaps(): features are on different seq regions.");
    return undef;
  }
  
  return ($self->end >= $f->start and $self->start <= $f->end);
}




##############################################
# Methods included for backwards compatibility
##############################################


# contig
#
# This method is included for backwards compatibility.
# Use slice() instead
#
sub contig {
  deprecate('Use slice() instead');
  slice(@_);
}



# sub_SeqFeature
#
# This method is only for genebuild backwards compatibility.
# Avoid using it if possible
#
sub sub_SeqFeature{
  my ($self) = @_;
  return @{$self->{'_gsf_sub_array'}} if($self->{'_gsf_sub_array'});
}

# add_sub_SeqFeature
#
# This method is only for genebuild backwards compatibility.
# Avoid using it if possible
#
sub add_sub_SeqFeature{
  my ($self,$feat,$expand) = @_;
  my ($p, $f, $l) = caller;
  if( $expand eq 'EXPAND' ) {
    # if this doesn't have start/end set - forget it!
    if( ! $self->start && ! $self->end ) {
      
      $self->start($feat->start());
      $self->end($feat->end());
      $self->strand($feat->strand);
    } else {
      if( $feat->start < $self->start ) {
        $self->start($feat->start);
      }

      if( $feat->end > $self->end ) {
        $self->end($feat->end);
      }
    }
   } else {
     if($self->start > $feat->start || $self->end < $feat->end) {
       throw("$feat is not contained within parent feature, " .
             "and expansion is not valid");
     }
   }

   push(@{$self->{'_gsf_sub_array'}},$feat);
}

# flush_sub_SeqFeature
#
# This method is only for genebuild backwards compatibility.
# Avoid using it isf possible
#
sub flush_sub_SeqFeature {
  my ($self) = @_;
  $self->{'_gsf_sub_array'} = [];
}


sub _deprecated_transform {
  my $self = shift;
  my $arg = shift;

  if(!$arg) {
    warning("Calling transform() with no arguments is deprecated.\n".
          "A coordinate system name argument should be used instead.\n".
          "You probably wanted transform('seqlevel') or transform('contig').");
    return $self->transform('seqlevel');
  }

  if(ref($arg) eq 'Bio::EnsEMBL::Slice') {
    if($arg->{'empty'}) {
      warning("Calling transform with an empty slice is deprecated.\n" .
                "A coordinate system name argument should be used instead.\n".
                "You probably wanted transform('chromosome') or " .
                "transform('toplevel')");
      return $self->transform('toplevel');
    }
    warning("Calling transform with a slice is deprecated.\n" .
              "Use the transfer method instead");
    return $self->transfer($arg);
  }

  warning("Calling transform with a [".ref($arg)."] arg is no longer " .
          "(or never was) supported.  Doing nothing instead.");

  return $self;
}

# id
#
# This method is included for backwards compatibility only.
# Use hseqname or dbID or stable_id instead
#
sub id {
  my $self = shift;
  deprecate("id method is not used - use display_id instead");
  return $self->{'stable_id'} if($self->{'stable_id'});
  return $self->{'hseqname'} if($self->{'hseqname'});
  return $self->{'seqname'}  if($self->{'seqname'});
  return $self->{'dbID'};
}


1;
