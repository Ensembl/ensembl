

#
# Ensembl module for Bio::EnsEMBL::Mapper
#
# Written by Ewan Birney <birney@ebi.ac.uk>
#
# Copyright GRL/EBI
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::Mapper

=head1 SYNOPSIS

  $map = Bio::EnsEMBL::Mapper->new('rawcontig', 'chromosome');

  # add a coodinate mapping - supply two pairs or coordinates
  $map->add_map_coordinates(
    $contig_id, $contig_start, $contig_end, $contig_ori,
    $chr_name, chr_start, $chr_end
  );

  # map from one coordinate system to another
  my @coordlist = $mapper->map_coordinates(627012, 2, 5, -1, "rawcontig");

=head1 DESCRIPTION

Generic mapper to provide coordinate transforms between two
disjoint coordinate systems. This mapper is intended to be
'context neutral' - in that it does not contain any code
relating to any particular coordinate system. This is
provided in, for example, Bio::EnsEMBL::AssemblyMapper.

=head1 AUTHOR - Ewan Birney

This module is part of the Ensembl project http://www.ensembl.org

Post general queries to B<ensembl-dev@ebi.ac.uk>

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

package Bio::EnsEMBL::Mapper;
use vars qw(@ISA);
use strict;


use Bio::EnsEMBL::Root;
use Bio::EnsEMBL::Mapper::Pair;
use Bio::EnsEMBL::Mapper::Unit;
use Bio::EnsEMBL::Mapper::Coordinate;
use Bio::EnsEMBL::Mapper::Gap;


@ISA = qw(Bio::EnsEMBL::Root);


sub new {
  my($class,@args) = @_;

  my $from = shift @args;
  my $to   = shift @args;

  my $self = {};
  bless $self,$class;

  if( !defined $to ) {
      $self->throw("Must supply from and to tags");
  }

  $self->{"_pair_$from"} = {};
  $self->{"_pair_$to"} = {};

  $self->to($to);
  $self->from($from);

# set stuff in self from @args
  return $self;
}


=head2 map_coordinates

    Arg  1      int $id
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
    Function    generic map method
    Returntype  array of Bio::EnsEMBL::Mapper::Coordinate
                and/or   Bio::EnsEMBL::Mapper::Gap
    Exceptions  none
    Caller      Bio::EnsEMBL::Mapper

=cut

sub map_coordinates{
   my ($self, $id, $start, $end, $strand, $type) = @_;

   #&eprof_start('map_coordinates');

   unless(defined($id) && defined($start) && defined($end) && 
	  defined($strand) && defined($type) ) {
       $self->throw("Must start,end,strand,id,type as coordinates");
   }


   my $hash = $self->{"_pair_$type"};

   my ($from, $to);

   if($type eq $self->{'to'}) {
     $from = 'to';
     $to = 'from';
   } else {
     $from = 'from';
     $to = 'to';
   }

   unless(defined $hash) {
       $self->throw("Type $type is neither to or from coordinate systems");
   }

   if( $self->{'_is_sorted'} == 0 ) {
       $self->_sort();
   }

   if( !defined $hash->{uc($id)} ) {
       # one big gap!
       my $gap = Bio::EnsEMBL::Mapper::Gap->new($start, $end);
       return $gap;
   }

   my $last_used_pair;
   my @result;

   my ( $start_idx, $end_idx, $mid_idx, $pair, $self_coord );
   my $lr = $hash->{uc($id)};
   
   $start_idx = 0;
   $end_idx = $#$lr;
   
   # binary search the relevant pairs
   # helps if the list is big
   while(( $end_idx - $start_idx ) <= 1 ) {
     $mid_idx = ($start_idx+$end_idx)>>1;
     $pair = $lr->[$mid_idx];
     $self_coord   = $pair->{$from};
     if( $self_coord->{'end'} < $start ) {
       $start_idx = $mid_idx;
     } else {
       $end_idx = $mid_idx;
     }
   }

   my $i = $start_idx;
   for( my $i = $start_idx; $i<=$#$lr; $i++ ) {
     $pair = $lr->[$i];
     my $self_coord   = $pair->{$from};
     my $target_coord = $pair->{$to};

     # if we haven't even reached the start, move on
     if( $self_coord->{'end'} < $start ) {
       next;
     }
     
     # if we have over run, break
     if( $self_coord->{'start'} > $end ) {
       last;
     }

     if( $start < $self_coord->{'start'} ) {
       # gap detected
       my $gap = Bio::EnsEMBL::Mapper::Gap->new($start, 
                                                $self_coord->{'start'}-1);
       
       push(@result,$gap);
       $start = $gap->{'end'}+1;
     }
     
     my ($target_start,$target_end,$target_ori);
     
     # start is somewhere inside the region
     if( $pair->{'ori'} == 1 ) {
       $target_start = 
	   $target_coord->{'start'} + ($start - $self_coord->{'start'});
     } else {
       $target_end = 
	   $target_coord->{'end'} - ($start - $self_coord->{'start'});
     }
     
     # either we are enveloping this map or not. If yes, then end
     # point (self perspective) is determined solely by target. If not
     # we need to adjust

     if( $end > $self_coord->{'end'} ) {
       # enveloped
       if( $pair->{'ori'} == 1 ) {
         $target_end = $target_coord->{'end'};
       } else {
         $target_start = $target_coord->{'start'};
       }
     } else {
       # need to adjust end
       if( $pair->{'ori'} == 1 ) {
         $target_end = 
             $target_coord->{'start'} + ($end - $self_coord->{'start'});
       } else {
         $target_start = 
             $target_coord->{'end'} - ($end - $self_coord->{'start'});
       }
     }

     my $res = Bio::EnsEMBL::Mapper::Coordinate->new($target_coord->{'id'},
						     $target_start,
						     $target_end,
						     $pair->{'ori'} * $strand);
     push(@result,$res);
     
     $last_used_pair = $pair;
     $start = $self_coord->{'end'}+1;
   }

   if( !defined $last_used_pair ) {
       my $gap = Bio::EnsEMBL::Mapper::Gap->new($start, $end);
       push(@result,$gap);

   } elsif( $last_used_pair->{$from}->{'end'} < $end ) {
       # gap at the end
       my $gap = Bio::EnsEMBL::Mapper::Gap->new(
			   $last_used_pair->{$from}->{'end'} + 1,
			   $end);
       push(@result,$gap);
   }

   if ( $strand == -1 ) {
       @result = reverse ( @result);
   }

   #&eprof_end('map_coordinates');

   return @result;
}

=head2 fastmap

    Arg  1      int $id
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

   my ($from, $to);

   if($type eq $self->{'to'}) {
     $from = 'to';
     $to = 'from';
   } else {
     $from = 'from';
     $to = 'to';
   }

   my $hash = $self->{"_pair_$type"} or
       $self->throw("Type $type is neither to or from coordinate systems");

   if( $self->{'_is_sorted'} == 0 ) {
       $self->_sort();
   }

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
		$strand );
     } else {
       return ( $target_coord->{'id'},
		$target_coord->{'end'} - ($end - $self_coord->{'start'}),
		$target_coord->{'end'} - ($start - $self_coord->{'start'}),
		-$strand );
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
                id of 'targe' sequence
    Arg  6      int $start
                start coordinate of 'targe' sequence
    Arg  7      int $end
                end coordinate of 'targe' sequence
    Function    stores details of mapping between two regions:
                'source' and 'target'
    Returntype  none
    Exceptions  none
    Caller      Bio::EnsEMBL::Mapper

=cut

sub add_map_coordinates{
   my ($self, $contig_id, $contig_start, $contig_end, 
       $contig_ori, $chr_name, $chr_start, $chr_end) = @_;

   unless(defined($contig_id) && defined($contig_start) && defined($contig_end)
          && defined($contig_ori) && defined($chr_name) && defined($chr_start)
          && defined($chr_end)) {
       $self->throw("7 arguments expected");
   }

   if( ($contig_end - $contig_start)  != ($chr_end - $chr_start) ) {
       $self->throw("Cannot deal with mis-lengthed mappings so far");
   }

   my $pair = Bio::EnsEMBL::Mapper::Pair->new();

   my $from = Bio::EnsEMBL::Mapper::Unit->new();
   $from->start($contig_start);
   $from->end($contig_end);
   $from->id($contig_id);

   my $to = Bio::EnsEMBL::Mapper::Unit->new();
   $to->start($chr_start);
   $to->end($chr_end);
   $to->id($chr_name);

   $pair->to($to);
   $pair->from($from);

   $pair->ori($contig_ori);

   # place into hash on both ids
   my $map_to = $self->{'to'};
   my $map_from = $self->{'from'};

   $self->{"_pair_$map_to"}->{uc($chr_name)} ||= [];
   push(@{$self->{"_pair_$map_to"}->{uc($chr_name)}},$pair);

   $self->{"_pair_$map_from"}->{uc($contig_id)} ||= [];
   push(@{$self->{"_pair_$map_from"}->{uc($contig_id)}},$pair);

   $self->_is_sorted(0);
}


=head2 list_pairs

    Arg  1      int $id
                id of 'source' sequence
    Arg  2      int $start
                start coordinate of 'source' sequence
    Arg  3      int $end
                end coordinate of 'source' sequence
    Arg  4      int $type
                nature of transform - gives the type of
                coordinates to be transformed *from*
    Function    list all pairs of mappings in a region
    Returntype  list of Bio::EnsEMBL::Mapper::Pair
    Exceptions  none
    Caller      Bio::EnsEMBL::Mapper

=cut

sub list_pairs{
   my ($self, $id, $start, $end, $type) = @_;


   if( !defined $type ) {
       $self->throw("Must start,end,id,type as coordinates");
   }

   if( $start > $end ) {
     $self->throw("Start is greater than end for id $id, start $start, end $end\n");
   }

   # perhaps a little paranoid/excessive
   if( $self->_is_sorted == 0 ) {
       $self->_sort();
   }


   my $hash = $self->{"_pair_$type"};

   my ($from, $to);

   if($type eq $self->{'to'}) {
     $from = 'to';
     $to = 'from';
   } else {
     $from = 'from';
     $to = 'to';
   }
     
   unless(defined $hash) {
       $self->throw("Type $type is neither to or from coordinate systems");
   }

   my @list;

   unless(exists $hash->{uc($id)}) {
     return ();
   }

   @list = @{$hash->{uc($id)}}; 

   my @output;
   if( $start == -1 && $end == -1 ) {
     return @list;
   } else {
     
     foreach my $p ( @list ) {
       
       if( $p->{$from}->{'end'} < $start ) {
	 next;
       }
       if( $p->{$from}->{'start'} > $end ) {
	 last;
       }
       push(@output,$p);
     }
     return @output;
   }
}


=head2 from, to

    Arg  1      Bio::EnsEMBL::Mapper::Unit $id
                id of 'source' sequence
    Function    accessor method form the 'source'
                and 'target' in a Mapper::Pair
    Returntype  Bio::EnsEMBL::Mapper::Unit
    Exceptions  none
    Caller      Bio::EnsEMBL::Mapper

=cut

sub to{
   my ($self,$value) = @_;
   if( defined $value) {
      $self->{'to'} = $value;
    }
    return $self->{'to'};

}

sub from{
   my ($self,$value) = @_;
   if( defined $value) {
      $self->{'from'} = $value;
    }
    return $self->{'from'};

}


=head2 _dump

    Arg  1      *FileHandle $fh
    Function    convenience dump function
                possibly useful for debugging
    Returntype  none
    Exceptions  none
    Caller      internal

=cut

sub _dump{
   my ($self,$fh) = @_;

   if( !defined $fh ) {
       $fh = \*STDERR;
   }

   foreach my $id ( keys %{$self->{'_pair_hash_from'}} ) {
       print $fh "From Hash $id\n";
       foreach my $pair ( @{$self->{'_pair_hash_from'}->{uc($id)}} ) {
	   print $fh "    ",$pair->from->start," ",$pair->from->end,":",$pair->to->start," ",$pair->to->end," ",$pair->to->id,"\n";
       }
   }

}


=head2 _sort

    Function    sort function so that all
                mappings are sorted by
                chromosome start
    Returntype  none
    Exceptions  none
    Caller      internal

=cut

sub _sort{
   my ($self) = @_;

   my $to = $self->{'to'};
   my $from = $self->{'from'};

   foreach my $id ( keys %{$self->{"_pair_$from"}} ) {
       @{$self->{"_pair_$from"}->{$id}} = sort { $a->{'from'}->{'start'} <=> $b->{'from'}->{'start'} } @{$self->{"_pair_$from"}->{$id}};
   }

   foreach my $id ( keys %{$self->{"_pair_$to"}} ) {
       @{$self->{"_pair_$to"}->{$id}} = sort { $a->{'to'}->{'start'} <=> $b->{'to'}->{'start'} } @{$self->{"_pair_$to"}->{$id}};
   }

   $self->_is_sorted(1);

}


=head2 _is_sorted

    Arg  1      int $sorted
    Function    toggle for whether the (internal)
                map data are sorted
    Returntype  int
    Exceptions  none
    Caller      internal

=cut

sub _is_sorted{
   my ($self,$value) = @_;
   if( defined $value) {
      $self->{'_is_sorted'} = $value;
    }
    if (! defined $self->{'_is_sorted'}) {
      $self->{'_is_sorted'} = 0;
    }
    return $self->{'_is_sorted'};

}


1;
