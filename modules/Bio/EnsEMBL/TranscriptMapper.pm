=head1 LICENSE

  Copyright (c) 1999-2012 The European Bioinformatics Institute and
  Genome Research Limited.  All rights reserved.

  This software is distributed under a modified Apache license.
  For license details, please see

    http://www.ensembl.org/info/about/code_licence.html

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <dev@ensembl.org>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk@ensembl.org>.

=cut

=head1 NAME

TranscriptMapper - A utility class used to perform coordinate conversions
between a number of coordinate systems relating to transcripts

=head1 SYNOPSIS

  my $trmapper = Bio::EnsEMBL::TranscriptMapper->new($transcript);

  @coords = $trmapper->cdna2genomic( 123, 554 );

  @coords = $trmapper->genomic2cdna( 141, 500, -1 );

  @coords = $trmapper->genomic2cds( 141, 500, -1 );

  @coords = $trmapper->pep2genomic( 10, 60 );

  @coords = $trmapper->genomic2pep( 123, 400, 1 );

=head1 DESCRIPTION

This is a utility class which can be used to perform coordinate conversions
between a number of coordinate systems relating to transcripts.

=head1 METHODS

=cut

package Bio::EnsEMBL::TranscriptMapper;

use strict;
use warnings;

use Bio::EnsEMBL::Utils::Exception qw(throw);

use Bio::EnsEMBL::Mapper;
use Bio::EnsEMBL::Mapper::Gap;
use Bio::EnsEMBL::Mapper::Coordinate;



=head2 new

  Arg [1]    : Bio::EnsEMBL::Transcript $transcript
               The transcript for which a TranscriptMapper should be created.
  Example    : $trans_mapper = Bio::EnsEMBL::TranscriptMapper->new($transcript)
  Description: Creates a TranscriptMapper object which can be used to perform
               various coordinate transformations relating to transcripts.
               Note that the TranscriptMapper uses the transcript state at the
               time of creation to perform the conversions, and that a new
               TranscriptMapper must be created if the Transcript is altered.
               'Genomic' coordinates are coordinates which are relative to the
               slice that the Transcript is on.
  Returntype : Bio::EnsEMBL::TranscriptMapper
  Exceptions : throws if a transcript is not an argument
  Caller     : Transcript::get_TranscriptMapper
  Status     : Stable

=cut

sub new {
  my $caller = shift;
  my $transcript = shift;

  my $class = ref($caller) || $caller;

  if(!ref($transcript) || !$transcript->isa('Bio::EnsEMBL::Transcript')) {
    throw("Transcript argument is required.");
  }


  my $exons = $transcript->get_all_Exons();
  my $start_phase;
  if(@$exons) {
    $start_phase = $exons->[0]->phase;
  } else {
    $start_phase = -1;
  }

  # Create a cdna <-> genomic mapper and load it with exon coords
  my $mapper = _load_mapper($transcript,$start_phase);

  my $self = bless({'exon_coord_mapper' => $mapper,
                    'start_phase'       => $start_phase,
                    'cdna_coding_start' => $transcript->cdna_coding_start(),
                    'cdna_coding_end'   => $transcript->cdna_coding_end()},
                   $class);
}


=head2 _load_mapper

  Arg [1]    : Bio::EnsEMBL::Transcript $transcript
               The transcript for which a mapper should be created.
  Example    : my $mapper = _load_mapper($transcript);
  Description: loads the mapper
  Returntype : Bio::EnsEMBL::Mapper
  Exceptions : none
  Caller     : Internal
  Status     : Stable

=cut

sub _load_mapper {
  my $transcript = shift;
  my $start_phase = shift;

  my $mapper = Bio::EnsEMBL::Mapper->new( 'cdna', 'genomic');

  my $edits_on = $transcript->edits_enabled();
  my @edits;

  if($edits_on) {
    @edits = @{$transcript->get_all_SeqEdits()};
    @edits = sort {$a->start() <=> $b->start()} @edits;
  }

  my $edit_shift = 0;

  my $cdna_start = undef;

  my $cdna_end = 0;


  foreach my $ex (@{$transcript->get_all_Exons}) {
    my $gen_start = $ex->start();
    my $gen_end   = $ex->end();

    $cdna_start = $cdna_end + 1;
    $cdna_end   = $cdna_start + $ex->length() - 1;

    my $strand = $ex->strand();

    # add deletions and insertions into pairs when SeqEdits turned on
    # ignore mismatches (i.e. treat as matches)
    if($edits_on) {
      while(@edits && $edits[0]->start() + $edit_shift <= $cdna_end) {

        my $edit = shift(@edits);
        my $len_diff = $edit->length_diff();

        if($len_diff) {
          # break pair into two parts, finish first pair just before edit

          my $prev_cdna_end   = $edit->start() + $edit_shift - 1;
          my $prev_cdna_start = $cdna_start;
          my $prev_len        = $prev_cdna_end - $prev_cdna_start + 1;

          my $prev_gen_end;
          my $prev_gen_start;
          if($strand == 1) {
            $prev_gen_start = $gen_start;
            $prev_gen_end   = $gen_start + $prev_len - 1;
          } else {
            $prev_gen_start = $gen_end - $prev_len + 1;
            $prev_gen_end   = $gen_end;
          }

          if($prev_len > 0) { # only create map pair if not boundary case
            $mapper->add_map_coordinates
              ('cdna', $prev_cdna_start, $prev_cdna_end, $strand,
               'genome', $prev_gen_start,$prev_gen_end);
          }

          $cdna_start = $prev_cdna_end  + 1;

          if($strand == 1) {
            $gen_start  = $prev_gen_end + 1;
          } else {
            $gen_end    = $prev_gen_start - 1;
          }

          $cdna_end  += $len_diff;

          if($len_diff > 0) {
            # insert in cdna, shift cdna coords along
            $cdna_start += $len_diff;
          } else {
            # delete in cdna (insert in genomic), shift genomic coords along

            if($strand == 1) {
              $gen_start  -= $len_diff;
            } else {
              $gen_end    += $len_diff;
            }
          }

          $edit_shift += $len_diff;
        }
      }
    }

    my $pair_len = $cdna_end - $cdna_start + 1;

    if($pair_len > 0) {
      $mapper->add_map_coordinates('cdna', $cdna_start, $cdna_end, $strand,
                                   'genome', $gen_start, $gen_end);
    }
  }

  return $mapper;
}


=head2 cdna2genomic

  Arg [1]    : $start
               The start position in cdna coordinates
  Arg [2]    : $end
               The end position in cdna coordinates
  Example    : @cdna_coords = $transcript_mapper->cdna2genomic($start, $end);
  Description: Converts cdna coordinates to genomic coordinates.  The
               return value is a list of coordinates and gaps.
  Returntype : list of Bio::EnsEMBL::Mapper::Coordinate and
               Bio::EnsEMBL::Mapper::Gap objects
  Exceptions : throws if no start or end
  Caller     : general
  Status     : Stable

=cut


sub cdna2genomic {
  my ($self,$start,$end) = @_;

  if( !defined $end ) {
    throw("Must call with start/end");
  }

  my $mapper = $self->{'exon_coord_mapper'};

  return  $mapper->map_coordinates( 'cdna', $start, $end, 1, "cdna" );

}


=head2 genomic2cdna

  Arg [1]    : $start
               The start position in genomic coordinates
  Arg [2]    : $end
               The end position in genomic coordinates
  Arg [3]    : $strand
               The strand of the genomic coordinates (default value 1)
  Example    : @coords = $trans_mapper->genomic2cdna($start, $end, $strnd);
  Description: Converts genomic coordinates to cdna coordinates.  The
               return value is a list of coordinates and gaps.  Gaps
               represent intronic or upstream/downstream regions which do
               not comprise this transcripts cdna.  Coordinate objects
               represent genomic regions which map to exons (utrs included).
  Returntype : list of Bio::EnsEMBL::Mapper::Coordinate and
               Bio::EnsEMBL::Mapper::Gap objects
  Exceptions : throws if start, end or strand not defined
  Caller     : general
  Status     : Stable

=cut

sub genomic2cdna {
  my ($self, $start, $end, $strand) = @_;

  unless(defined $start && defined $end && defined $strand) {
    throw("start, end and strand arguments are required\n");
  }

  my $mapper = $self->{'exon_coord_mapper'};

  return $mapper->map_coordinates("genome", $start, $end, $strand,"genomic");

}


=head2 cds2genomic

  Arg [1]    : int $start
               start position in cds coords
  Arg [2]    : int $end
               end position in cds coords
  Example    : @genomic_coords = $transcript_mapper->cds2genomic(69, 306);
  Description: Converts cds coordinates into genomic coordinates.  The
               coordinates returned are relative to the same slice that the
               transcript used to construct this TranscriptMapper was on.
  Returntype : list of Bio::EnsEMBL::Mapper::Gap and
               Bio::EnsEMBL::Mapper::Coordinate objects
  Exceptions : throws if no end
  Caller     : general
  Status     : at risk

=cut

sub cds2genomic {
  my ( $self, $start, $end ) = @_;

  if ( !( defined($start) && defined($end) ) ) {
    throw("Must call with start and end");
  }

  # Move start end into translate cDNA coordinates now.
  $start = $start +( $self->{'cdna_coding_start'} - 1 ) ;
  $end = $end + ( $self->{'cdna_coding_start'} - 1 );

  return $self->cdna2genomic( $start, $end );
}

=head2 pep2genomic

  Arg [1]    : int $start
               start position in peptide coords
  Arg [2]    : int $end
               end position in peptide coords
  Example    : @genomic_coords = $transcript_mapper->pep2genomic(23, 102);
  Description: Converts peptide coordinates into genomic coordinates.  The
               coordinates returned are relative to the same slice that the
               transcript used to construct this TranscriptMapper was on.
  Returntype : list of Bio::EnsEMBL::Mapper::Gap and
               Bio::EnsEMBL::Mapper::Coordinate objects
  Exceptions : throws if no end
  Caller     : general
  Status     : Stable

=cut

sub pep2genomic {
  my ( $self, $start, $end ) = @_;

  if ( !( defined($start) && defined($end) ) ) {
    throw("Must call with start and end");
  }

  # Take possible N-padding at beginning of CDS into account.
  my $start_phase = $self->{'start_phase'};
  my $shift = ( $start_phase > 0 ) ? $start_phase : 0;

  # Move start end into translate cDNA coordinates now.
  $start = 3*$start - 2 + ( $self->{'cdna_coding_start'} - 1 ) - $shift;
  $end = 3*$end + ( $self->{'cdna_coding_start'} - 1 ) - $shift;

  return $self->cdna2genomic( $start, $end );
}


=head2 genomic2cds

  Arg [1]    : int $start
               The genomic start position
  Arg [2]    : int $end
               The genomic end position
  Arg [3]    : int $strand
               The genomic strand
  Example    : @cds_coords = $trans_mapper->genomic2cds($start, $end, $strand);
  Description: Converts genomic coordinates into CDS coordinates of the
               transcript that was used to create this transcript mapper.
  Returntype : list of Bio::EnsEMBL::Mapper::Coordinate and
               Bio::EnsEMBL::Mapper::Gap objects
  Exceptions : throw if start, end or strand not defined
  Caller     : general
  Status     : Stable

=cut

sub genomic2cds {
 my ($self, $start, $end, $strand) = @_;

 if(!defined($start) || !defined($end) || !defined($strand)) {
   throw("start, end and strand arguments are required");
 }

 if($start > $end + 1) {
   throw("start arg must be less than or equal to end arg + 1");
 }

 my $cdna_cstart = $self->{'cdna_coding_start'};
 my $cdna_cend   = $self->{'cdna_coding_end'};

 #this is a pseudogene if there is no coding region
 if(!defined($cdna_cstart)) {
   #return a gap of the entire requested region, there is no CDS
   return Bio::EnsEMBL::Mapper::Gap->new($start,$end);
 }

 my @coords = $self->genomic2cdna($start, $end, $strand);

 my @out;

 foreach my $coord (@coords) {
   if($coord->isa('Bio::EnsEMBL::Mapper::Gap')) {
     push @out, $coord;
   } else {
     my $start = $coord->start;
     my $end   = $coord->end;

     if($coord->strand == -1 || $end < $cdna_cstart || $start > $cdna_cend) {
       #is all gap - does not map to peptide
       push @out, Bio::EnsEMBL::Mapper::Gap->new($start,$end);
     } else {
       #we know area is at least partially overlapping CDS
	
       my $cds_start = $start - $cdna_cstart + 1;
       my $cds_end   = $end   - $cdna_cstart + 1;

       if($start < $cdna_cstart) {
         #start of coordinates are in the 5prime UTR
         push @out, Bio::EnsEMBL::Mapper::Gap->new($start, $cdna_cstart-1);

         #start is now relative to start of CDS
         $cds_start = 1;
       }
	
       my $end_gap = undef;
       if($end > $cdna_cend) {
         #end of coordinates are in the 3prime UTR
         $end_gap = Bio::EnsEMBL::Mapper::Gap->new($cdna_cend + 1, $end);
         #adjust end to relative to CDS start
         $cds_end = $cdna_cend - $cdna_cstart + 1;
       }

       #start and end are now entirely in CDS and relative to CDS start
       $coord->start($cds_start);
       $coord->end($cds_end);

       push @out, $coord;

       if($end_gap) {
         #push out the region which was in the 3prime utr
         push @out, $end_gap;
       }
     }	
   }
 }

 return @out;

}


=head2 genomic2pep

  Arg [1]    : $start
               The start position in genomic coordinates
  Arg [2]    : $end
               The end position in genomic coordinates
  Arg [3]    : $strand
               The strand of the genomic coordinates
  Example    : @pep_coords = $transcript->genomic2pep($start, $end, $strand);
  Description: Converts genomic coordinates to peptide coordinates.  The
               return value is a list of coordinates and gaps.
  Returntype : list of Bio::EnsEMBL::Mapper::Coordinate and
               Bio::EnsEMBL::Mapper::Gap objects
  Exceptions : throw if start, end or strand not defined
  Caller     : general
  Status     : Stable

=cut

sub genomic2pep {
 my ($self, $start, $end, $strand) = @_;

 unless(defined $start && defined $end && defined $strand) {
   throw("start, end and strand arguments are required");
 }

 my @coords = $self->genomic2cds($start, $end, $strand);

 my @out;

 my $start_phase = $self->{'start_phase'};

 #take into account possible N padding at beginning of CDS
 my $shift = ($start_phase > 0) ? $start_phase : 0;

 foreach my $coord (@coords) {
   if($coord->isa('Bio::EnsEMBL::Mapper::Gap')) {
     push @out, $coord;
   } else {

     #start and end are now entirely in CDS and relative to CDS start

     #convert to peptide coordinates
     my $pep_start = int(($coord->start + $shift + 2) / 3);
     my $pep_end   = int(($coord->end   + $shift + 2) / 3);
     $coord->start($pep_start);
     $coord->end($pep_end);

     push @out, $coord;
   }	
 }

 return @out;
}


1;
