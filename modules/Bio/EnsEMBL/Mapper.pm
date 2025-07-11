=head1 LICENSE

See the NOTICE file distributed with this work for additional information
regarding copyright ownership.

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

Bio::EnsEMBL::Mapper

=head1 SYNOPSIS

  $map = Bio::EnsEMBL::Mapper->new( 'rawcontig', 'chromosome' );

  # add a coodinate mapping - supply two pairs or coordinates
  $map->add_map_coordinates(
    $contig_id, $contig_start, $contig_end, $contig_ori,
    $chr_name,  chr_start,     $chr_end
  );

  # map from one coordinate system to another
  my @coordlist =
    $mapper->map_coordinates( 627012, 2, 5, -1, "rawcontig" );

=head1 DESCRIPTION

Generic mapper to provide coordinate transforms between two disjoint
coordinate systems. This mapper is intended to be 'context neutral' - in
that it does not contain any code relating to any particular coordinate
system. This is provided in, for example, Bio::EnsEMBL::AssemblyMapper.

Mappings consist of pairs of 'to-' and 'from-' contigs with coordinates on each.
Orientation is abbreviated to 'ori', 

The contig pair hash is divided into mappings per seq_region, the code 
below makes assumptions about how to filter these results, thus the comparisons 
for some properties are absent in the code but implicit by data structure.

The assembly mapping hash '_pair_last' orders itself by the target seq region and looks like this:

   1 => ARRAY(0x1024c79c0)
      0  Bio::EnsEMBL::Mapper::Pair=HASH(0x1024d6198)
         'from' => Bio::EnsEMBL::Mapper::Unit=HASH(0x1025edf98)
            'end' => 4
            'id' => 4
            'start' => 1
         'ori' => 1
         'to' => Bio::EnsEMBL::Mapper::Unit=HASH(0x1025edf68)
            'end' => 4
            'id' => 1
            'start' => 1
      1  Bio::EnsEMBL::Mapper::Pair=HASH(0x1026c20f0)
         'from' => Bio::EnsEMBL::Mapper::Unit=HASH(0x1025ee3a0)
            'end' => 12
            'id' => 4
            'start' => 9
         'ori' => 1
         'to' => Bio::EnsEMBL::Mapper::Unit=HASH(0x1025ee370)
            'end' => 4
            'id' => 1
            'start' => 1
   2 => ARRAY(0x1025ee460)
      0  Bio::EnsEMBL::Mapper::Pair=HASH(0x1025ee400)
         'from' => Bio::EnsEMBL::Mapper::Unit=HASH(0x1025ee2c8)
            'end' => 8
            'id' => 4
            'start' => 5
         'ori' => 1
         'to' => Bio::EnsEMBL::Mapper::Unit=HASH(0x1025ee2b0)
            'end' => 4
            'id' => 2
            'start' => 1
      1  Bio::EnsEMBL::Mapper::Pair=HASH(0x1025ee658)
         'from' => Bio::EnsEMBL::Mapper::Unit=HASH(0x1025eea48)
            'end' => 16
            'id' => 4
            'start' => 13
         'ori' => 1
         'to' => Bio::EnsEMBL::Mapper::Unit=HASH(0x1025eea18)
            'end' => 4
            'id' => 2
            'start' => 1

The other mapping hash available is the reverse sense, putting the 'from' seq_region as the
sorting key. Here is an excerpt.

0  HASH(0x102690bb8)
   4 => ARRAY(0x1025ee028)
      0  Bio::EnsEMBL::Mapper::Pair=HASH(0x1024d6198)
         'from' => Bio::EnsEMBL::Mapper::Unit=HASH(0x1025edf98)
            'end' => 4
            'id' => 4
            'start' => 1
         'ori' => 1
         'to' => Bio::EnsEMBL::Mapper::Unit=HASH(0x1025edf68)
            'end' => 4
            'id' => 1
            'start' => 1
      1  Bio::EnsEMBL::Mapper::Pair=HASH(0x1025ee400)
         'from' => Bio::EnsEMBL::Mapper::Unit=HASH(0x1025ee2c8)
            'end' => 8
            'id' => 4
            'start' => 5
         'ori' => 1
         'to' => Bio::EnsEMBL::Mapper::Unit=HASH(0x1025ee2b0)
            'end' => 4
            'id' => 2
            'start' => 1


=head1 METHODS

=cut

package Bio::EnsEMBL::Mapper;
use strict;
use integer;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Mapper::Pair;
use Bio::EnsEMBL::Mapper::IndelPair;
use Bio::EnsEMBL::Mapper::Unit;
use Bio::EnsEMBL::Mapper::Coordinate;
use Bio::EnsEMBL::Mapper::IndelCoordinate;
use Bio::EnsEMBL::Mapper::Gap;

use Bio::EnsEMBL::Utils::Exception qw(throw);
use Bio::EnsEMBL::Utils::Interval;
use Bio::EnsEMBL::Utils::Tree::Interval::Immutable;

=head2 new

  Arg [1]    : string $from
               The name of the 'from' coordinate system
  Arg [2]    : string $to
               The name of the 'to' coordinate system
  Arg [3]    : (optional) Bio::EnsEMBL::CoordSystem $from_cs
               The 'from' coordinate system
  Arg [4]    : (optional) Bio::EnsEMBL::CoordSystem $to_cs
  Example    : my $mapper = Bio::EnsEMBL::Mapper->new('FROM', 'TO');
  Description: Constructor.  Creates a new Bio::EnsEMBL::Mapper object.
  Returntype : Bio::EnsEMBL::Mapper
  Exceptions : none
  Caller     : general

=cut

sub new {
  my ( $proto, $from, $to, $from_cs, $to_cs ) = @_;

  if ( !defined($to) || !defined($from) ) {
    throw("Must supply 'to' and 'from' tags");
  }

  my $class = ref($proto) || $proto;

  my $self = bless( { "_pair_$from" => {},
                      "_pair_$to"   => {},
                      "_tree_$from" => {},
                      "_tree_$to"   => {},
                      'pair_count'  => 0,
                      'to'          => $to,
                      'from'        => $from,
                      'to_cs'       => $to_cs,
                      'from_cs'     => $from_cs
                    },
                    $class );

  # do sql to get any componente with muliple assemblys.

  return $self;
}

=head2 flush

  Args       : none
  Example    : none
  Description: removes all cached information out of this mapper
  Returntype : none
  Exceptions : none
  Caller     : AssemblyMapper, ChainedAssemblyMapper

=cut

sub flush {
  my $self = shift;
  my $from = $self->from();
  my $to = $self->to();

  $self->{"_pair_$from"} = {};
  $self->{"_pair_$to"} = {};
  $self->{"_tree_$from"} = {};
  $self->{"_tree_$to"} = {};

  $self->{'pair_count'} = 0;
}



=head2 map_coordinates

    Arg  1      string $id
                id of 'source' sequence
    Arg  2      int $start
                start coordinate of 'source' sequence
    Arg  3      int $end
                end coordinate of 'source' sequence
    Arg  4      int $strand
                raw contig orientation (+/- 1)
    Arg  5      string $type
                nature of transform - gives the type of
                coordinates to be transformed *from*
    Arg  6      boolean (0 or 1) $include_original_region
                option to include original input coordinate region mappings in the result
    Arg  7      int  $cdna_coding_start
                cdna coding start  
    Function    generic map method
    Returntype  if $include_original_region == 0
                  array of mappped Bio::EnsEMBL::Mapper::Coordinate
                  and/or   Bio::EnsEMBL::Mapper::Gap
                if $include_original_region == 1
                  hash of mapped and original Bio::EnsEMBL::Mapper::Coordinate
                  and/or   Bio::EnsEMBL::Mapper::Gap
    Exceptions  none
    Caller      Bio::EnsEMBL::Mapper

=cut

sub map_coordinates {
  my ( $self, $id, $start, $end, $strand, $type, $include_original_region, $cdna_coding_start ) = @_;
  unless (    defined($id)
           && defined($start)
           && defined($end)
           && defined($strand)
           && defined($type) )
  {
    throw("Expecting at least 5 arguments");
  }
  
  $cdna_coding_start = defined $cdna_coding_start ? $cdna_coding_start : 1;
  
  # special case for handling inserts:
  if ( $start == $end + 1 ) {
    return $self->map_insert( $id, $start, $end, $strand, $type, undef, $include_original_region);
  }

  if ( !$self->{'_is_sorted'} ) { $self->_sort() }

  my $hash = $self->{"_pair_$type"};

  my ( $from, $to, $cs );

  if ( $type eq $self->{'to'} ) {
    $from = 'to';
    $to   = 'from';
    $cs   = $self->{'from_cs'};
  } else {
    $from = 'from';
    $to   = 'to';
    $cs   = $self->{'to_cs'};
  }

  unless ( defined $hash ) {
    throw("Type $type is neither to or from coordinate systems");
  }

  my @result;
  my @paired_result;

  if ( !defined $hash->{ uc($id) } ) {
    # one big gap!
    my $gap = Bio::EnsEMBL::Mapper::Gap->new( $start, $end );
    if ($include_original_region) {
      push @paired_result, { 'original' => $gap, 'mapped' => $gap };
      return @paired_result;
    } else {
      push @result, $gap;
      return @result;
    }
  }

  my $last_used_pair;

  # if we don't already have an interval tree built for this id,
  # build one now
  if ( !defined( $self->{"_tree_$type"}->{ uc($id) } )) {
    $self->{"_tree_$type"}->{ uc($id) } = _build_immutable_tree($from,
                                                                $hash->{ uc($id) });
  }
  # query the interval tree (either cached or created new) for overlapping intervals
  my $overlap = $self->{"_tree_$type"}->{ uc($id) }->query($start, $end);

  my $rank              = 0;
  my $orig_start        = $start;
  my $last_target_coord = undef;
  foreach my $i (@{$overlap}) {
    my $pair = $i->data;
    my $self_coord   = $pair->{$from};
    my $target_coord = $pair->{$to};

    #
    # But not the case for haplotypes!! need to test for this case???
    # so removing this till a better solution is found
    #
    #
    #     if($self_coord->{'start'} < $start){
    #       $start = $orig_start;
    #       $rank++;
    #     }

    if ( defined($last_target_coord) && $target_coord->{'id'} ne $last_target_coord ) {
      if ( $self_coord->{'start'} < $start ) {    
      # i.e. the same bit is being mapped to another assembled bit
        $start = $orig_start;
      }
    } else {
      $last_target_coord = $target_coord->{'id'};
    }

    if ( $start < $self_coord->{'start'} ) {
      # gap detected
      my $gap = Bio::EnsEMBL::Mapper::Gap->new( $start, $self_coord->{'start'} - 1, $rank );
      push( @result, $gap );
      push(@paired_result, { 'original' => $gap, 'mapped' => $gap });

      $start = $gap->{'end'} + 1;
    }
    my ( $target_start, $target_end);
    my ( $ori_start, $ori_end);
    
    my $res;
    if ( exists $pair->{'indel'} ) {
      # When next pair is an IndelPair and not a Coordinate, create the
      # new mapping Coordinate, the IndelCoordinate.
      $target_start = $target_coord->{'start'};
      $target_end   = $target_coord->{'end'};
      
      #original coordinates
      $ori_start = 	$self_coord->{'start'};
      $ori_end = 	$self_coord->{'end'};
	  
      #create a Gap object
      my $gap = Bio::EnsEMBL::Mapper::Gap->new( $start,
        ( $self_coord->{'end'} < $end ? $self_coord->{'end'} : $end ) );
      #create the Coordinate object
      my $coord =
        Bio::EnsEMBL::Mapper::Coordinate->new( $target_coord->{'id'},
              $target_start, $target_end, $pair->{'ori'}*$strand, $cs );
      #and finally, the IndelCoordinate object with
      $res = Bio::EnsEMBL::Mapper::IndelCoordinate->new( $gap, $coord );
      
    } else {
      # start is somewhere inside the region
      if ( $pair->{'ori'} == 1 ) {
        $target_start = $target_coord->{'start'} + ( $start - $self_coord->{'start'} );
      } else {
        $target_end = $target_coord->{'end'} - ( $start - $self_coord->{'start'} );
      }
       
      # Either we are enveloping this map or not.  If yes, then end
      # point (self perspective) is determined solely by target.  If
      # not we need to adjust.

      if ( $end > $self_coord->{'end'} ) {
        # enveloped
        if ( $pair->{'ori'} == 1 ) {
          $target_end = $target_coord->{'end'};
        } else {
          $target_start = $target_coord->{'start'};
        }
      } else {
        # need to adjust end
        if ( $pair->{'ori'} == 1 ) {
          $target_end =
            $target_coord->{'start'} +
            ( $end - $self_coord->{'start'} );
        } else {
                  
          $target_start =
            $target_coord->{'end'} - ( $end - $self_coord->{'start'} );
        }
         
      }
	
      $res = Bio::EnsEMBL::Mapper::Coordinate->new( $target_coord->{'id'},
               $target_start, $target_end, $pair->{'ori'}*$strand,$cs, $rank );
               
     #save the ori coordinate info
     $ori_start = ($start - $cdna_coding_start) + 1;
     $ori_end =   $ori_start + ($target_end - $target_start);
     
 
    } ## end else [ if ( exists $pair->{'indel'...})]

    push( @result, $res );
    my $res_ori = Bio::EnsEMBL::Mapper::Coordinate->new( $self_coord->{'id'},
                    $ori_start, $ori_end, $pair->{'ori'}*$strand,$cs, $rank);
    push(@paired_result, { 'original' => $res_ori, 'mapped' => $res });
    

    $last_used_pair = $pair;
    $start          = $self_coord->{'end'} + 1;

  } ## end for ( my $i = $start_idx...)

  if ( !defined $last_used_pair ) {
    my $gap = Bio::EnsEMBL::Mapper::Gap->new( $start, $end );
    push( @result, $gap );
    push(@paired_result, { 'original' => $gap, 'mapped' => $gap });
    

  } elsif ( $last_used_pair->{$from}->{'end'} < $end ) {
    # gap at the end
    my $gap = Bio::EnsEMBL::Mapper::Gap->new(
                                  $last_used_pair->{$from}->{'end'} + 1,
                                  $end, $rank );
    push( @result, $gap );
    
    push(@paired_result, { 'original' => $gap, 'mapped' => $gap });
  }

  if ( $strand == -1 ) {
    @result = reverse(@result);
    @paired_result = reverse(@paired_result);
  }

  if ($include_original_region){
    return @paired_result;
  }else{
    return @result;
  }

} ## end sub map_coordinates



=head2 map_insert

  Arg [1]    : string $id
  Arg [2]    : int $start - start coord. Since this is an insert should always
               be one greater than end.
  Arg [3]    : int $end - end coord. Since this is an insert should always
               be one less than start.
  Arg [4]    : int $strand (0, 1, -1)
  Arg [5]    : string $type - the coordinate system name the coords are from.
  Arg [6]    : boolean $fastmap - if specified, this is being called from
               the fastmap call. The mapping done is not any faster for
               inserts, but the return value is different.
  Example    : 
  Description: This is in internal function which handles the special mapping
               case for inserts (start = end +1).  This function will be called
               automatically by the map function so there is no reason to
               call it directly.
  Returntype : list of Bio::EnsEMBL::Mapper::Coordinate and/or Gap objects
  Exceptions : none
  Caller     : map_coordinates()

=cut

sub map_insert {
  my ($self, $id, $start, $end, $strand, $type, $fastmap, $include_original_region) = @_;

  # swap start/end and map the resultant 2bp coordinate
  ($start, $end) =($end,$start);

  my @coords = $self->map_coordinates($id, $start, $end, $strand, $type, $include_original_region);

  if(@coords == 1) {
    my $c = $coords[0];
    if ($include_original_region) {
      my $m = $c->{'mapped'};
      my $orig = $c->{'original'};
      ($m->{'start'}, $m->{'end'}) = ($m->{'end'}, $m->{'start'});
      ($orig->{'start'}, $orig->{'end'}) = ($orig->{'end'}, $orig->{'start'});
    } else {
      ($c->{'start'}, $c->{'end'}) = ($c->{'end'}, $c->{'start'});
    }
  } else {
    throw("Unexpected: Got ",scalar(@coords)," expected 2.") if(@coords != 2);

    # adjust coordinates, remove gaps
    my ($c1, $c2);

    if($strand == -1) {
      ($c2,$c1) = @coords;
    } else {
      ($c1, $c2) = @coords;
    }
    @coords = ();

    my ($m1, $m2);
    if ($include_original_region) {
      ($m1, $m2) = ($c1->{'mapped'}, $c2->{'mapped'});
    } else {
      ($m1, $m2) = ($c1, $c2);
    }
    if(ref($m1) eq 'Bio::EnsEMBL::Mapper::Coordinate') {
      # insert is after first coord
      if($m1->{'strand'} * $strand == -1) {
        $m1->{'end'}--;
      } else {
        $m1->{'start'}++;
      }
      if ($include_original_region) {
        $c1->{'mapped'} = $m1;
        @coords = ($c1);
      } else {
        @coords = ($m1);
      }
    }
    if(ref($m2) eq 'Bio::EnsEMBL::Mapper::Coordinate') {
      # insert is before second coord
      if($m2->{'strand'} * $strand == -1) {
        $m2->{'start'}++;
      } else {
        $m2->{'end'}--;
      }
      if($strand == -1) {
        if ($include_original_region) {
          $c2->{'mapped'} = $m2;
          unshift @coords, $c2;
        } else {
          unshift @coords, $m2;
        }
      } else {
        if ($include_original_region) {
          $c2->{'mapped'} = $m2;
          push @coords, $c2;
        } else {
          push @coords, $m2;
        }
      }
    }

    # as with @coords == 1 above, if we are in
    # include_original_region mode, then we have
    # to flip back the original region as well
    if ($include_original_region) {
      foreach my $coord (@coords) {
        my $orig = $coord->{'original'};
          ($orig->{'start'}, $orig->{'end'}) =
            ($orig->{'end'}, $orig->{'start'});
      }
    }
  }

  if($fastmap) {
    return undef if(@coords != 1);
    my $c = $include_original_region ? $coords[0]->{'mapped'} : $coords[0];
    return ($c->{'id'}, $c->{'start'}, $c->{'end'},
            $c->{'strand'}, $c->{'coord_system'});
  }

  return @coords;
}






=head2 fastmap

    Arg  1      string $id
                id of 'source' sequence
    Arg  2      int $start
                start coordinate of 'source' sequence
    Arg  3      int $end
                end coordinate of 'source' sequence
    Arg  4      int $strand
                raw contig orientation (+/- 1)
    Arg  5      int $type
                nature of transform - gives the type of
                coordinates to be transformed *from*
    Function    inferior map method. Will only do ungapped unsplit mapping.
                Will return id, start, end strand in a list.
    Returntype  list of results
    Exceptions  none
    Caller      Bio::EnsEMBL::AssemblyMapper

=cut

sub fastmap {
   my ($self, $id, $start, $end, $strand, $type) = @_;

   my ($from, $to, $cs);

   if($end+1 == $start) {
     return $self->map_insert($id, $start, $end, $strand, $type, 1);
   }

   if( ! $self->{'_is_sorted'} ) { $self->_sort() }

   if($type eq $self->{'to'}) {
     $from = 'to';
     $to = 'from';
     $cs = $self->{'from_cs'};
   } else {
     $from = 'from';
     $to = 'to';
     $cs = $self->{'to_cs'};
   }

   my $hash = $self->{"_pair_$type"} or
       throw("Type $type is neither to or from coordinate systems");


   my $pairs = $hash->{uc($id)};

   foreach my $pair (@$pairs) {
     my $self_coord   = $pair->{$from};
     my $target_coord = $pair->{$to};

     # only super easy mapping is done 
     if( $start < $self_coord->{'start'} ||
         $end > $self_coord->{'end'} ) {
       next;
     }

     if( $pair->{'ori'} == 1 ) {
       return ( $target_coord->{'id'},
                $target_coord->{'start'}+$start-$self_coord->{'start'},
                $target_coord->{'start'}+$end-$self_coord->{'start'},
                $strand, $cs );
     } else {
       return ( $target_coord->{'id'},
                $target_coord->{'end'} - ($end - $self_coord->{'start'}),
                $target_coord->{'end'} - ($start - $self_coord->{'start'}),
                -$strand, $cs );
     }
   }

   return ();
}



=head2 add_map_coordinates

    Arg  1      int $id
                id of 'source' sequence
    Arg  2      int $start
                start coordinate of 'source' sequence
    Arg  3      int $end
                end coordinate of 'source' sequence
    Arg  4      int $strand
                relative orientation of source and target (+/- 1)
    Arg  5      int $id
                id of 'target' sequence
    Arg  6      int $start
                start coordinate of 'target' sequence
    Arg  7      int $end
                end coordinate of 'target' sequence
    Function    Stores details of mapping between
                'source' and 'target' regions.
    Returntype  none
    Exceptions  none
    Caller      Bio::EnsEMBL::Mapper

=cut

sub add_map_coordinates {
  my ( $self, $contig_id, $contig_start, $contig_end, $contig_ori,
       $chr_name, $chr_start, $chr_end )
    = @_;

  unless (    defined($contig_id)
           && defined($contig_start)
           && defined($contig_end)
           && defined($contig_ori)
           && defined($chr_name)
           && defined($chr_start)
           && defined($chr_end) )
  {
    throw("7 arguments expected");
  }

  if ( ( $chr_end > $chr_start ) and ( $contig_end - $contig_start ) != ( $chr_end - $chr_start ) ) {
    throw("Cannot deal with mis-lengthed mappings so far");
  }

  my $from = Bio::EnsEMBL::Mapper::Unit->new( $contig_id, $contig_start,
                                              $contig_end );
  my $to =
    Bio::EnsEMBL::Mapper::Unit->new( $chr_name, $chr_start, $chr_end );

  my $pair = Bio::EnsEMBL::Mapper::Pair->new( $from, $to, $contig_ori );

  # place into hash on both ids
  my $map_to   = $self->{'to'};
  my $map_from = $self->{'from'};

  push( @{ $self->{"_pair_$map_to"}->{ uc($chr_name) } },    $pair );
  push( @{ $self->{"_pair_$map_from"}->{ uc($contig_id) } }, $pair );

  # any interval trees for this set of intervals are now invalid
  # so remove them
  $self->{"_tree_$map_to"}->{ uc($contig_id) } = undef;
  $self->{"_tree_$map_from"}->{ uc($contig_id) } = undef;

  $self->{'pair_count'}++;
  $self->{'_is_sorted'} = 0;
} ## end sub add_map_coordinates


=head2 add_indel_coordinates

    Arg  1      int $id
                id of 'source' sequence
    Arg  2      int $start
                start coordinate of 'source' sequence
    Arg  3      int $end
                end coordinate of 'source' sequence
    Arg  4      int $strand
                relative orientation of source and target (+/- 1)
    Arg  5      int $id
                id of 'targe' sequence
    Arg  6      int $start
                start coordinate of 'targe' sequence
    Arg  7      int $end
                end coordinate of 'targe' sequence
    Function    stores details of mapping between two regions:
                'source' and 'target'. Returns 1 if the pair was added, 0 if it
                was already in. Used when adding an indel
    Returntype  int 0,1
    Exceptions  none
    Caller      Bio::EnsEMBL::Mapper

=cut

sub add_indel_coordinates{
  my ($self, $contig_id, $contig_start, $contig_end, 
      $contig_ori, $chr_name, $chr_start, $chr_end) = @_;

  unless(defined($contig_id) && defined($contig_start) && defined($contig_end)
   && defined($contig_ori) && defined($chr_name) && defined($chr_start)
   && defined($chr_end)) {
    throw("7 arguments expected");
  }

  #we need to create the IndelPair object to add to both lists, to and from
  my $from =
    Bio::EnsEMBL::Mapper::Unit->new($contig_id, $contig_start, $contig_end);
  my $to   =
    Bio::EnsEMBL::Mapper::Unit->new($chr_name, $chr_start, $chr_end);

  my $pair = Bio::EnsEMBL::Mapper::IndelPair->new($from, $to, $contig_ori);

  # place into hash on both ids
  my $map_to = $self->{'to'};
  my $map_from = $self->{'from'};

  push( @{$self->{"_pair_$map_to"}->{uc($chr_name)}}, $pair );
  push( @{$self->{"_pair_$map_from"}->{uc($contig_id)}}, $pair );

  # any interval trees for this set of intervals are now invalid
  # so remove them
  $self->{"_tree_$map_to"}->{ uc($chr_name) } = undef;
  $self->{"_tree_$map_from"}->{ uc($contig_id) } = undef;

  $self->{'pair_count'}++;

  $self->{'_is_sorted'} = 0;
  return 1;
}


=head2 map_indel

  Arg [1]    : string $id
  Arg [2]    : int $start - start coord. Since this is an indel should always
               be one greater than end.
  Arg [3]    : int $end - end coord. Since this is an indel should always
               be one less than start.
  Arg [4]    : int $strand (0, 1, -1)
  Arg [5]    : string $type - the coordinate system name the coords are from.
  Example    : @coords = $mapper->map_indel();
  Description: This is in internal function which handles the special mapping
               case for indels (start = end +1). It will be used to map from
               a coordinate system with a gap to another that contains an
               insertion. It will be mainly used by the Variation API.
  Returntype : Bio::EnsEMBL::Mapper::Unit objects
  Exceptions : none
  Caller     : general

=cut

sub map_indel {
  my ( $self, $id, $start, $end, $strand, $type ) = @_;

  # swap start/end and map the resultant 2bp coordinate
  ( $start, $end ) = ( $end, $start );

  if ( !$self->{'_is_sorted'} ) { $self->_sort() }

  my $hash = $self->{"_pair_$type"};

  my ( $from, $to, $cs );

  if ( $type eq $self->{'to'} ) {
    $from = 'to';
    $to   = 'from';
    $cs   = $self->{'from_cs'};
  } else {
    $from = 'from';
    $to   = 'to';
    $cs   = $self->{'to_cs'};
  }

  unless ( defined $hash ) {
    throw("Type $type is neither to or from coordinate systems");
  }
  my @indel_coordinates;

  # if we don't already have an interval tree built for this id,
  # build one now
  if ( !defined $self->{"_tree_$type"}->{ uc($id) } ) {
    $self->{"_tree_$type"}->{ uc($id) } = _build_immutable_tree($from,
                                                                $hash->{ uc($id) });
  }
  # query the interval tree (either cached or created new) for overlapping intervals
  my $overlap = $self->{"_tree_$type"}->{ uc($id) }->query($start, $end);

  foreach my $i (@{$overlap}) {
    my $pair = $i->data;
    my $self_coord   = $pair->{$from};
    my $target_coord = $pair->{$to};

    if ( exists $pair->{'indel'} ) {
      #need to return unit coordinate
      my $to =
        Bio::EnsEMBL::Mapper::Unit->new( $target_coord->{'id'},
                                         $target_coord->{'start'},
                                         $target_coord->{'end'}, );
      push @indel_coordinates, $to;
      last;
    }
  }

  return @indel_coordinates;
} ## end sub map_indel


=head2 add_Mapper

    Arg  1      Bio::EnsEMBL::Mapper $mapper2
    Example     $mapper->add_Mapper($mapper2)
    Function    add all the map coordinates from $mapper to this mapper.
                This object will contain mapping pairs from both the old
                object and $mapper2.
    Returntype  int 0,1
    Exceptions  throw if 'to' and 'from' from both Bio::EnsEMBL::Mappers
                are incompatible
    Caller      $mapper->methodname()

=cut

sub add_Mapper{
  my ($self, $mapper) = @_;

  my $mapper_to = $mapper->{'to'};
  my $mapper_from = $mapper->{'from'};
  if ($mapper_to ne $self->{'to'} or $mapper_from ne $self->{'from'}) {
    throw("Trying to add an incompatible Mapper");
  }

  my $count_a = 0;
  foreach my $seq_name (keys %{$mapper->{"_pair_$mapper_to"}}) {
    push(@{$self->{"_pair_$mapper_to"}->{$seq_name}},
        @{$mapper->{"_pair_$mapper_to"}->{$seq_name}});
    $self->{"_tree_$mapper_to"}->{ uc($seq_name) } = undef;
    $count_a += scalar(@{$mapper->{"_pair_$mapper_to"}->{$seq_name}});
  }
  my $count_b = 0;
  foreach my $seq_name (keys %{$mapper->{"_pair_$mapper_from"}}) {
    push(@{$self->{"_pair_$mapper_from"}->{$seq_name}},
        @{$mapper->{"_pair_$mapper_from"}->{$seq_name}});
    $self->{"_tree_$mapper_from"}->{ uc($seq_name) } = undef;
    $count_b += scalar(@{$mapper->{"_pair_$mapper_from"}->{$seq_name}});
  }

  if ($count_a == $count_b) {
    $self->{'pair_count'} += $count_a;
  } else {
    throw("Trying to add a funny Mapper");
  }

  $self->{'_is_sorted'} = 0;
  return 1;
}



=head2 list_pairs

    Arg  1      int $id
                id of 'source' sequence
    Arg  2      int $start
                start coordinate of 'source' sequence
    Arg  3      int $end
                end coordinate of 'source' sequence
    Arg  4      string $type
                nature of transform - gives the type of
                coordinates to be transformed *from*
    Function    list all pairs of mappings in a region
    Returntype  list of Bio::EnsEMBL::Mapper::Pair
    Exceptions  none
    Caller      Bio::EnsEMBL::Mapper

=cut

sub list_pairs {
  my ( $self, $id, $start, $end, $type ) = @_;

  if ( !$self->{'_is_sorted'} ) { $self->_sort() }

  if ( !defined $type ) {
    throw("Expected 4 arguments");
  }

  if ( $start > $end ) {
    throw(   "Start is greater than end "
           . "for id $id, start $start, end $end\n" );
  }

  my $hash = $self->{"_pair_$type"};

  my ( $from, $to );

  if ( $type eq $self->{'to'} ) {
    $from = 'to';
    $to   = 'from';
  } else {
    $from = 'from';
    $to   = 'to';
  }

  unless ( defined $hash ) {
    throw("Type $type is neither to or from coordinate systems");
  }

  my @list;

  unless ( exists $hash->{ uc($id) } ) {
    return ();
  }

  @list = @{ $hash->{ uc($id) } };

  my @output;
  if ( $start == -1 && $end == -1 ) {
    return @list;
  } else {

    foreach my $p (@list) {

      if ( $p->{$from}->{'end'} < $start ) {
        next;
      }
      if ( $p->{$from}->{'start'} > $end ) {
        last;
      }
      push( @output, $p );
    }
    return @output;
  }
} ## end sub list_pairs


=head2 to

    Arg  1      Bio::EnsEMBL::Mapper::Unit $id
                id of 'source' sequence
    Function    accessor method form the 'source'
                and 'target' in a Mapper::Pair
    Returntype  Bio::EnsEMBL::Mapper::Unit
    Exceptions  none
    Caller      Bio::EnsEMBL::Mapper

=cut

sub to {
  my ( $self, $value ) = @_;

  if ( defined($value) ) {
    $self->{'to'} = $value;
  }

  return $self->{'to'};
}

=head2 from

    Arg  1      Bio::EnsEMBL::Mapper::Unit $id
                id of 'source' sequence
    Function    accessor method form the 'source'
                and 'target' in a Mapper::Pair
    Returntype  Bio::EnsEMBL::Mapper::Unit
    Exceptions  none
    Caller      Bio::EnsEMBL::Mapper

=cut
sub from {
  my ( $self, $value ) = @_;

  if ( defined($value) ) {
    $self->{'from'} = $value;
  }

  return $self->{'from'};
}


#  _dump
#
#     Arg  1      *FileHandle $fh
#     Function    convenience dump function
#                 possibly useful for debugging
#     Returntype  none
#     Exceptions  none
#     Caller      internal
#

sub _dump{
  my ($self,$fh) = @_;

  if( !defined $fh ) {
    $fh = \*STDERR;
  }

  my $from = $self->{'from'};
  my $to = $self->{'to'};

  print $fh "dumping from-hash _pair_$from\n";
  foreach my $id ( keys %{$self->{"_pair_$from"}} ) {
    print $fh "{_pair_$from}->{" . uc($id) . "}:\n";
    foreach my $pair ( @{$self->{"_pair_$from"}->{uc($id)}} ) {
      print $fh "    ",$pair->from->start," ",$pair->from->end,":",$pair->to->start," ",$pair->to->end," ",$pair->to->id,"\n";
    }
    if (defined($self->{"_tree_$from"}->{uc($id)})) {
      print $fh "{_tree_$from}->{" . uc($id) . "} instantiated\n";
    } else {
      print $fh "{_tree_$from}->{" . uc($id) . "} empty\n";
    }
  }
  print $fh "dumping to-hash _pair_$to\n";
  foreach my $id ( keys %{$self->{"_pair_$to"}} ) {
    print $fh "{_pair_$to}->{" . uc($id) . "}:\n";
    foreach my $pair ( @{$self->{"_pair_$to"}->{uc($id)}} ) {
      print $fh "    ",$pair->to->start," ",$pair->to->end,":",$pair->from->start," ",$pair->from->end," ",$pair->from->id,"\n";
    }
    if (defined($self->{"_tree_$to"}->{uc($id)})) {
      print $fh "{_tree_$to}->{" . uc($id) . "} instantiated\n";
    } else {
      print $fh "{_tree_$to}->{" . uc($id) . "} empty\n";
    }
  }
}


# _sort
#
#    Function    sort function so that all
#                mappings are sorted by
#                chromosome start
#    Returntype  none
#    Exceptions  none
#    Caller      internal
#

sub _sort {
  my ($self) = @_;

  my $to   = $self->{'to'};
  my $from = $self->{'from'};

  foreach my $id ( keys %{ $self->{"_pair_$from"} } ) {
    @{ $self->{"_pair_$from"}->{$id} } =
      sort { $a->{'from'}->{'start'} <=> $b->{'from'}->{'start'} }
      @{ $self->{"_pair_$from"}->{$id} };
  }

  foreach my $id ( keys %{ $self->{"_pair_$to"} } ) {
    @{ $self->{"_pair_$to"}->{$id} } =
      sort { $a->{'to'}->{'start'} <=> $b->{'to'}->{'start'} }
      @{ $self->{"_pair_$to"}->{$id} };
  }

  $self->_merge_pairs();
  $self->_is_sorted(1);
}

# this function merges pairs that are adjacent into one
sub _merge_pairs {
  my $self = shift;

  my ( $lr, $lr_from, $del_pair, $next_pair, $current_pair );

  my $map_to = $self->{'to'};
  my $map_from = $self->{'from'};
  $self->{'pair_count'} = 0;

  for my $key ( keys %{$self->{"_pair_$map_to"}} ) {
    # merging pairs means all interval trees for
    # this set of intervals are now invalid, so delete them
    $self->{"_tree_$map_to"}->{ uc($key) } = undef;
    $lr = $self->{"_pair_$map_to"}->{$key}; 
    
    my $i = 0;
    my $next = 1;
    my $length = $#{$lr};
    while( $next <= $length ) {
      $current_pair = $lr->[$i];
      $next_pair = $lr->[$next];
      $del_pair = undef;
      
      if(exists $current_pair->{'indel'} || exists $next_pair->{'indel'}){
        #necessary to modify the merge function to not merge indels
        $next++;
        $i++;
      }
      else{
        # duplicate filter
        if( $current_pair->{'to'}->{'start'} == $next_pair->{'to'}->{'start'} 
           && $current_pair->{'from'}->{'id'} == $next_pair->{'from'}->{'id'} 
           && $current_pair->{'from'}->{'start'} == $next_pair->{'from'}->{'start'}) {
        # Modified in e75 to support GRCh38. Contigs used repeatedly in one seq region were 
        # being pre-emptively deleted as copies. Extra from-start condition above ensures copy
        # deletion is restricted to same-location copies. Even more stringent checks can be made
        # at cost of speed.
            $del_pair = $next_pair;
          } elsif ( ( $current_pair->{'from'}->{'id'} eq $next_pair->{'from'}->{'id'} ) &&
                    ( $next_pair->{'ori'} == $current_pair->{'ori'} ) &&
                    ( $next_pair->{'to'}->{'start'} -1 == $current_pair->{'to'}->{'end'} )) {
          
          if( $current_pair->{'ori'} == 1 ) {
            # check forward strand merge
            if( $next_pair->{'from'}->{'start'} - 1 == $current_pair->{'from'}->{'end'} ) {
              # normal merge with previous element
              $current_pair->{'to'}->{'end'} = $next_pair->{'to'}->{'end'};
              $current_pair->{'from'}->{'end'} = $next_pair->{'from'}->{'end'};
              $del_pair = $next_pair;
            }
          } else {
            # check backward strand merge
            if( $next_pair->{'from'}->{'end'} + 1 == $current_pair->{'from'}->{'start'} ) {
              # yes its a merge
              $current_pair->{'to'}->{'end'} = $next_pair->{'to'}->{'end'};
              $current_pair->{'from'}->{'start'} = $next_pair->{'from'}->{'start'};
              $del_pair = $next_pair;
            }
          }
        }
      
          if( defined $del_pair ) {
            splice( @$lr, $next, 1 );
            $lr_from = $self->{"_pair_$map_from"}->{uc($del_pair->{'from'}->{'id'})};
            for( my $j=0; $j <= $#{$lr_from}; $j++ ) {
              if( $lr_from->[$j] == $del_pair ) {
                splice( @$lr_from, $j, 1 );
                last;
              }
            }
            $length--;
            if( $length < $next ) { last; }
          }
          else {
            $next++;
            $i++;
          }
      }  
      
    }
    $self->{'pair_count'} += scalar( @$lr );
  }
}


# _is_sorted
#
#    Arg  1      int $sorted
#    Function    toggle for whether the (internal)
#                map data are sorted
#    Returntype  int
#    Exceptions  none
#    Caller      internal
#

sub _is_sorted {
   my ($self, $value) = @_;
   $self->{'_is_sorted'} = $value if (defined($value));
   return $self->{'_is_sorted'};
}

# _build_immutable_tree
#
#  Arg 1       string $pair_side - the from or to half of each pair to be
#              the source of intervals
#  Arg 2       listref $pair_list - a list of Bio::EnsEMBL::Mapper::Pair
#  Function    builds a Bio::EnsEMBL::Utils::Tree::Interval::Immutable with
#              intervals corresponding to the chosen side (from or to) of
#              each Pair in $pair_list and a pointer to each Pair
#  Returntype  Bio::EnsEMBL::Utils::Tree::Interval::Immutable
#  Exceptions  none
#  Caller      internal

sub _build_immutable_tree {
  my ($pair_side, $pair_list) = @_;
  # create set of intervals for the tree
  my $from_intervals;

  foreach my $i (@{$pair_list}) {
    my $start = $i->{$pair_side}{start};
    my $end = $i->{$pair_side}{end};

    if ($end < $start) {
        ($end, $start) = ($start, $end);
    }

    push @{$from_intervals}, Bio::EnsEMBL::Utils::Interval->new($start, $end, $i);
  }

  # Create the interval tree defined on the above set of intervals
  #
  # As of release/99, we have more experience with the immutable
  # interval tree, so we will stick with this one. A future
  # refactoring effort may wish to replace this with a dynamically
  # maintained mutable interval tree, rather than simply throwing
  # trees away and rebuilding when the underlying interval set changes
  return Bio::EnsEMBL::Utils::Tree::Interval::Immutable->new($from_intervals);
}

1;
