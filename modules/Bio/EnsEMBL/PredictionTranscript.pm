# EnsEMBL module for Transcript
# Copyright EMBL-EBI/Sanger center 2002
#
# Cared for by Arne Stabenau
#

#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

PredictionTranscript

=head1 SYNOPSIS


=head1 DESCRIPTION

Container for single transcript ab initio gene prediction ala GenScan or SNAP.
Is directly storable/retrievable in EnsEMBL using PredictionTranscript Adaptor.

Creation:

     my $tran = new Bio::EnsEMBL::PredictionTranscript();
     $tran->add_Exon( $exon );

     my $tran = new Bio::EnsEMBL::PredictionTranscript(@exons);

     The order of the exons has to be right, as PT cant judge how to sort them.
     ( no sort as in Bio::EnsEMBL::Transcript )

     PredictionTranscript is geared towards the partial retrieve case from db.
     Exons can be missing in the middle. For storage though its necessary to 
     have them all and in contig coord system. 

Manipulation:

     # Returns an array of Exon objects, might contain undef instead of exon.
     my @exons = @{$tran->get_all_Exons};  

     # Returns the peptide translation as string 
     my $pep   = $tran->translate;

     # phase padded Exons cdna sequence. Phases usually match.
     my $cdna = $trans->get_cdna() 


=head1 CONTACT

contact EnsEMBL dev <ensembl-dev@ebi.ac.uk> for information

=cut

package Bio::EnsEMBL::PredictionTranscript;
use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Translation;

@ISA = qw(Bio::EnsEMBL::Transcript);


=head2 coding_region_start

  Arg [1]    : none
  Example    : $coding_region_start = $pt->coding_region_start
  Description: Retrieves the start of the coding region of this transcript in
               slice coordinates.  For prediction transcripts this
               is always the start of the transcript (i.e. there is no UTR).
               By convention, the coding_region_start is always lower than
               the value returned by the coding_end method.
               The value returned by this function is NOT the biological
               coding start since on the reverse strand the biological coding
               start would be the higher genomic value.
  Returntype : int
  Exceptions : none
  Caller     : general

=cut

sub coding_region_start {
  my $self = shift;
  return $self->start();
}


=head2 coding_region_end

  Arg [1]    : none
  Example    : $coding_region_end = $transcript->coding_region_end
  Description: Retrieves the start of the coding region of this prediction
               transcript. For prediction transcripts this is always the same
               as the end since no UTRs are stored.
               By convention, the coding_region_end is always higher than the
               value returned by the coding_region_start method.
               The value returned by this function is NOT the biological
               coding start since on the reverse strand the biological coding
               end would be the lower genomic value.
  Returntype : int
  Exceptions : none
  Caller     : general

=cut

sub coding_region_end {
  my $self = shift;
  return $self->end();
}



=head2 get_all_translateable_Exons

  Arg [1]    : none
  Example    : $exons = $self->get_all_translateable_Exons
  Description: Retrieves the translateable portion of all exons in this
               transcript.  For prediction transcripts this means all exons
               since no UTRs are stored for them.
  Returntype : listref of Bio::EnsEMBL::PredictionExons
  Exceptions : none
  Caller     : general

=cut

sub get_all_translateable_Exons {
  my $self = shift;
  return $self->get_all_Exons();
}


sub get_all_DBEntries { return []; }

sub get_all_DBLinks { return []; }

sub add_DBEntry {}

sub external_db { return undef; }

sub external_status { return undef; }

sub external_name { return undef; }

sub is_known { return 0;}


sub translation {
  my $self = shift;

  #calculate translation on the fly
  my $strand = $self->strand();

  my $start_exon;
  my $end_exon;

  my @exons = $self->get_all_Exons();

  return undef if(!@$exons);

  if($strand == 1) {
    $start_exon = $exons[0];
    $end_exon = $exons[-1];
  } else {
    $start_exon = $exons[-1];
    $end_exon = $exons[0];
  }

  return
    Bio::EnsEMBL::Translation->new(-START_EXON => $start_exon,
                                   -END_EXON   => $end_exon,
                                   -SEQ_START  => 1,
                                   -SEQ_END    => $end_exon->length());
}



=head2 translate

  Args      : none
  Function  : Give a peptide translation of all exons currently in
              the PT. Gives empty string when none is in.
  Returntype: a Bio::Seq as in transcript->translate()
  Exceptions: if the exons come in two or more groups, with an undef exon
              in the middle, only the first group is translated.
  Caller    : general

=cut


sub translate {
  my ($self) = @_;

  my $dna = $self->get_cdna();
  $dna    =~ s/TAG$|TGA$|TAA$//i;
  # the above line will remove the final stop codon from the mrna
  # sequence produced if it is present, this is so any peptide produced
  # won't have a terminal stop codon
  # if you want to have a terminal stop codon either comment this line out
  # or call translatable seq directly and produce a translation from it

  my $bioseq = new Bio::Seq( -seq => $dna, -moltype => 'dna' );

  return $bioseq->translate();
}
	 


=head2 get_cdna

  Args      : none
  Function  : Give a concat cdna of all exons currently in
              the PT. Gives empty string when none is in. Pads between not  
              phase matching exons. Builds internal coord translation table.
  Returntype: txt
  Exceptions: if the exons come in two or more groups, with an undef exon
              in the middle, only the first groups cdna is returned.
  Caller    : general, $self->translate()

=cut

sub get_cdna {
  my $self = shift;

  my $exons = $self->get_all_Exons();
  
  my $cdna = undef;
  my $lastphase = 0;

  my ( $cdna_start, $cdna_end );
  my ( $pep_start, $pep_end );
  my ( $new_cdna, $pep_count );

  $cdna_start = 1;
  $pep_start = 1;

  $self->{'_exon_align'} = [];

  for my $exon ( @$exons ) {
    my $exon_align = {};
    if( ! defined $exon ) {
      if( ! defined $cdna ) {
	next;
      } else {
	last;
      }
    } 

    push( @{$self->{'_exon_align'}}, $exon_align );

    my $phase = 0;

    if (defined($exon->phase)) {
      $phase = $exon->phase;
    }

    if( $phase != $lastphase ) {

      if( $lastphase == 1 ) {
	$cdna .= 'NN';
	$cdna_start += 2;
	$pep_start++;
      } elsif( $lastphase == 2 ) {
	$cdna .= 'N';
	$cdna_start += 1;
	$pep_start++;
      }

      #startpadding for this exon
      $cdna .= 'N' x $phase;
      $cdna_start += $phase;
    }
    
    $new_cdna = $exon->seq->seq();
#    print $new_cdna."\n";
    $cdna .= $new_cdna;
    $cdna_end = $cdna_start + CORE::length( $new_cdna ) - 1;

    # how many peptides are added by this exon??

    $pep_count = int( ( CORE::length( $new_cdna ) + $phase + 2 ) / 3 );

    $pep_end = $pep_start + $pep_count - 1; 
    $lastphase = $exon->end_phase();
      
    $exon_align->{ 'cdna_start' } = $cdna_start;
    $exon_align->{ 'cdna_end' } =  $cdna_end;
    $exon_align->{ 'pep_start' } = $pep_start;
    $exon_align->{ 'pep_end' } = $pep_end;
    $exon_align->{ 'exon' } = $exon;

    if( $lastphase == 0 ) { 
      $pep_start = $pep_end + 1;
    } else {
      $pep_start = $pep_end;
    }
    $cdna_start = $cdna_end+1;

  }

  if( ! defined $cdna ) { $cdna = '' };
  return $cdna
}



=head1 pep2genomic

  Arg  1   : integer start - relative to peptide
  Arg  2   : integer end   - relative to peptide

  Function : Provides a list of Bio::EnsEMBL::SeqFeatures which
             is the genomic coordinates of this start/end on the peptide

  Returns  : list of Bio::EnsEMBL::SeqFeature

=cut

sub pep2genomic {
  my ($self,$start,$end) = @_;

  if( !defined $end ) {
    $self->throw("Must call with start/end");
  }

  # move start end into translate cDNA coordinates now.
  # much easier!
  $start = 3* $start-2;
  $end   = 3* $end;

  #
  # Adjust the phase
  #
  my $exons = $self->get_all_Exons;
  if($exons && (my $e = $exons->[0])) {
    $start -= $e->phase;
    $end   -= $e->phase;
  }

  return $self->cdna2genomic( $start, $end );
}



=head2 cdna2genomic

  Arg [1]    : $start
               The start position in genomic coordinates
  Arg [2]    : $end
               The end position in genomic coordinates
  Arg [3]    : (optional) $strand
               The strand of the genomic coordinates
  Example    : @coords = $transcript->cdna2genomic($start, $end);
  Description: Converts cdna coordinates to genomic coordinates.  The
               return value is a list of coordinates and gaps.
  Returntype : list of Bio::EnsEMBL::Mapper::Coordinate and
               Bio::EnsEMBL::Mapper::Gap objects
  Exceptions : none
  Caller     : general

=cut

sub cdna2genomic {
  my $self = shift;
  my $start = shift || 1;
  my $end = shift || $self->length;

  if( !defined $end ) {
    $self->throw("Must call with start/end");
  }

  my $mapper = $self->_get_cdna_coord_mapper();

  return $mapper->map_coordinates( $self, $start, $end, 1, "cdna" );
}



=head2 genomic2cdna

  Arg [1]    : $start
               The start position in genomic coordinates
  Arg [2]    : $end
               The end position in genomic coordinates
  Arg [3]    : (optional) $strand
               The strand of the genomic coordinates (default value 1)
  Arg [4]    : (optional) $contig
               The contig the coordinates are on.  This can be a slice
               or RawContig, but must be the same object in memory as
               the contig(s) of this transcripts exon(s), because of the
               use of object identity. If no contig argument is specified the
               contig of the first exon is used, which is fine for slice
               coordinates but may cause incorrect mappings in raw contig
               coords if this transcript spans multiple contigs.
  Example    : @coords = $transcript->genomic2cdna($start, $end, $strnd, $ctg);
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
  my ($self, $start, $end, $strand, $contig) = @_;

  unless(defined $start && defined $end) {
    $self->throw("start and end arguments are required\n");
  }
  $strand = 1 unless(defined $strand);

  #"ids" in mapper are contigs of exons, so use the same contig that should
  #be attached to all of the exons...
  unless(defined $contig) {
    my @exons = @{$self->get_all_translateable_Exons};
    return () unless(@exons);
    $contig = $exons[0]->contig;
  }
  my $mapper = $self->_get_cdna_coord_mapper;

  return $mapper->map_coordinates($contig, $start, $end, $strand, "genomic");
}


=head2 _get_cdna_coord_mapper

  Args       : none
  Example    : none
  Description: creates and caches a mapper from "cdna" coordinate system to 
               "genomic" coordinate system. Uses Exons to help with that. Only
               calculates in the translateable part. 
  Returntype : Bio::EnsEMBL::Mapper( "cdna", "genomic" );
  Exceptions : none
  Caller     : cdna2genomic, pep2genomic

=cut

sub _get_cdna_coord_mapper {
  my ( $self ) = @_;

  if( defined $self->{'_exon_coord_mapper'} ) {
    return $self->{'_exon_coord_mapper'};
  }

  #
  # the mapper is loaded with OBJECTS in place of the IDs !!!!
  #  the objects are the contigs in the exons
  #
  my $mapper;
  $mapper = Bio::EnsEMBL::Mapper->new( "cdna", "genomic" );
  my @exons = @{$self->get_all_translateable_Exons() };
  my $start = 1;
  for my $exon ( @exons ) {
    $exon->load_genomic_mapper( $mapper, $self, $start );
    $start += $exon->length;
  }
  $self->{'_exon_coord_mapper'} = $mapper;
  return $mapper;
}



=head2 type

  Arg [1]    : (optional) string $type 
  Example    : none
  Description: Getter setter for the type of this prediction transcript
  Returntype : string
  Exceptions : none
  Caller     : none

=cut

sub type {
  my ($self, $type) = @_;

  if(defined $type) {
    $self->{'_type'} = $type;
  }

  return $self->{'_type'};
}


=head2 transform

  Arg [1]    : (optional) $slice
               The slice to transform this transcript to.  If not provided
               the transcript is transformed to raw contig coords.
  Example    : $pt->transform($slice)
  Description: Transforms this prediction transcript to slice or chromosomal
               coordinates
  Returntype : none
  Exceptions : none
  Caller     : general

=cut

sub transform {
  my $self = shift;

  my @exons;

  foreach my $exon (@{$self->get_all_Exons()}) {
    push @exons, $exon->transform(@_);
  }

  #flush the exons and all related internal caches
  $self->flush_Exons();

  # attach the new list of exons to the transcript
  $self->{'exons'} = \@exons;
}




=head2 get_exon_count

  Description: DEPRECATED - use get_all_Exons instead

=cut

sub get_exon_count {
   my $self = shift;
   deprecate("Use scalar(@{get_all_Exons}) instead");
   return scalar( @{$self->get_all_Exons} );
}


=head2 set_exon_count

  Description: DEPRECATED - this method does nothing now

=cut

sub set_exon_count {
  deprecate('This method no longer does anything.');
}



1;
