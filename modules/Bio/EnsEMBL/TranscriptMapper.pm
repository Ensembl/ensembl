#
# Ensembl module for TranscriptMapper
#
# Copyright (c) 2004 Ensembl
#
# You may distribute this module under the same terms as perl itself
#

=head1 NAME

TranscriptMapper - A utility class used to perform coordinate conversions
between a number of coordinate systems relating to transcripts

=head1 SYNOPSIS

  my $trmapper = Bio::EnsEMBL::TranscriptMapper->new($transcript);

  @coords = $trmapper->cdna2genomic(123, 554);

  @coords = $trmapper->genomic2cdna(141, 500, -1);

  @coords = $trmapper->genomic2cds(141, 500, -1);

  @coords = $trmapper->pep2genomic(10, 60);

  @coords = $trmapper->genomic2pep(123, 400, 1);


=head1 DESCRIPTION

This is a utility class which can be used to perform coordinate conversions
between a number of coordinate systems relating to transcripts.

=head1 CONTACT

Email questions to the ensembl developer mailing list <ensembl-dev@ebi.ac.uk>

=head1 METHODS

=cut

use strict;
use warnings;

use Bio::EnsEMBL::Utils::Exception qw(throw);

use Bio::EnsEMBL::Mapper;
use Bio::EnsEMBL::Mapper::Gap;
use Bio::EnsEMBL::Mapper::Coordinate;

package Bio::EnsEMBL::TranscriptMapper;



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
  Exceptions : none
  Caller     : Transcript::get_TranscriptMapper

=cut

sub new {
  my $caller = shift;
  my $transcript = shift;

  my $class = ref($caller) || $caller;

  if(!ref($transcript) || !$transcript->isa('Bio::EnsEMBL::Transcript')) {
    throw("Transcript argument is required.");
  }

  #
  # Create a cdna <-> genomic mapper and load it with exon coords
  #

  my $mapper = Bio::EnsEMBL::Mapper->new( 'cdna', 'genomic');

  my $slice = $transcript->slice();

  my $seqname = ($slice) ? $slice->name() : $transcript->seqname();
  $seqname ||= 'genome';

  my $cdna_start = 1;
  foreach my $exon (@{$transcript->get_all_Exons()}) {
    my $cdna_end = $cdna_start + $exon->length() - 1;
    $mapper->add_map_coordinates('cdna', $cdna_start, $cdna_end,
                                 $exon->strand(),
                                 $seqname, $exon->start, $exon->end());
    $cdna_start = $cdna_end + 1;
  }

  my $exons = $transcript->get_all_Exons();
  my $start_phase;
  if(@$exons) {
    $start_phase = $exons->[0]->phase;
  } else {
    $start_phase = -1;
  }

  return bless({'exon_coord_mapper' => $mapper,
                'seqname'          => $seqname,
                'start_phase'       => $start_phase,
                'cdna_coding_start' => $transcript->cdna_coding_start(),
                'cdna_coding_end'   => $transcript->cdna_coding_end()},$class);
}



sub exon_coord_mapper {
  my $self = shift;
  return $self->{'exon_coord_mapper'};
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
  Exceptions : none
  Caller     : general

=cut


sub cdna2genomic {
  my ($self,$start,$end) = @_;

  if( !defined $end ) {
    throw("Must call with start/end");
  }

  my $mapper = $self->exon_coord_mapper();

  return $mapper->map_coordinates( 'cdna', $start, $end, 1, "cdna" );
}


=head2 genomic2cdna

  Arg [1]    : $start
               The start position in genomic coordinates
  Arg [2]    : $end
               The end position in genomic coordinates
  Arg [3]    : (optional) $strand
               The strand of the genomic coordinates (default value 1)
  Example    : @coords = $trans_mapper->genomic2cdna($start, $end, $strnd);
  Description: Converts genomic coordinates to cdna coordinates.  The
               return value is a list of coordinates and gaps.  Gaps
               represent intronic or upstream/downstream regions which do
               not comprise this transcripts cdna.  Coordinate objects
               represent genomic regions which map to exons (utrs included).
  Returntype : list of Bio::EnsEMBL::Mapper::Coordinate and
               Bio::EnsEMBL::Mapper::Gap objects
  Exceptions : none
  Caller     : general

=cut

sub genomic2cdna {
  my ($self, $start, $end, $strand) = @_;

  unless(defined $start && defined $end && defined $strand) {
    throw("start, end and strand arguments are required\n");
  }

  my $mapper = $self->exon_coord_mapper();

  return $mapper->map_coordinates($self->{'seqname'}, $start, $end, $strand,
                                  "genomic");

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
  Exceptions : none
  Caller     : general

=cut

sub pep2genomic {
  my ($self,$start,$end) = @_;

  if( !defined $end ) {
    throw("Must call with start/end");
  }

  # move start end into translate cDNA coordinates now.
  # much easier!
  $start = 3* $start-2 + ($self->{'cdna_coding_start'} - 1);
  $end   = 3* $end + ($self->{'cdna_coding_start'} - 1);

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
  Exceptions : throw if incorrect arguments
  Caller     : general

=cut

sub genomic2cds {
 my ($self, $start, $end, $strand) = @_;

 if(!defined($start) || !defined($end) || !defined($strand)) {
   throw("start, end and strand arguments are required");
 }

 if($start > $end) {
   throw("start arg must be less than or equal to end arg");
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
  Exceptions : none
  Caller     : general

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
