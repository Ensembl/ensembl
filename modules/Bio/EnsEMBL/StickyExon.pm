#
# BioPerl module for Exon
#
# Cared for by Ewan Birney <birney@sanger.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::StickyExon - A Confirmed Exon which spans two contigs internally

=head1 SYNOPSIS

    $sticky = new Bio::EnsEMBL::Exon;

    # is a normal exon
    $sticky->start();
    $sticky->end();

    # has component_Exons
    foreach $sub ( $sticky->each_component_Exon ) {
       # $sub is an exon that ends on a contig
    }

=head1 DESCRIPTION

Sticky Exons represent Exons which internally span contigs. They are made during the
write back on virtual contigs, which writes the exons that span joins into the database.


=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::StickyExon;

use Bio::Seq;

use vars qw(@ISA $AUTOLOAD);
use strict;

# Object preamble - inherits from Bio::SeqFeature::Generic

use Bio::EnsEMBL::Exon;


@ISA = qw(Bio::EnsEMBL::Exon);
# new() is inherited from Bio::Root::Object

# _initialize is where the heavy stuff will happen when new is called

sub new {
  my($class,@args) = @_;
  
  my $self = Bio::EnsEMBL::Exon->new(@args);
  bless $self,$class;



  # Array to store exon tags
  $self->{_component_exons} = [];
  
  return $self;
}


=head2 id

 Title   : id
 Usage   : overrides id to get/set locally and sets all component exons
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub id{
   my ($self,$value) = @_;

   if( defined $value ) {
       $self->{'_sticky_id'} = $value;
       foreach my $c ( $self->each_component_Exon() ) {
	   $c->id($value);
       }
   }

   return $self->{'_sticky_id'};

}


=head2 each_component_Exons

 Title   : each_component_Exon
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub each_component_Exon{
   my ($self,@args) = @_;

   return @{$self->{'_component_exons'}};
}



=head1

  Arg  1   : integer start - relative to the exon
  Arg  2   : integer end   - relative to the exon

  Function : Provides a list of Bio::EnsEMBL::SeqFeatures which
             is the genomic coordinates of this start/end on the exon
             For simple exons this is one feature  for Stickies this
             is overridden and gives out a list of Bio::EnsEMBL::SeqFeatures

  Returns  : list of Bio::EnsEMBL::SeqFeature


=cut

sub contig_seqfeatures_from_relative_position {
  my ($self,$start,$end) = @_;

  if( !defined $end ) {
    $self->throw("Have not defined all the methods!");
  }

  # easy
  if( $start < 1 ) {
    $self->warn("Attempting to fetch start less than 1 ($start)");
    $start = 1;
  }

  if( $end > $self->length ) {
    $self->warn("Attempting to fetch end greater than end of exon ($end)");
    $end = $self->length;
  }

  my @out;
  my $sf;
  my @exons = $self->each_component_Exon();
  my $len = 0;
  while( scalar(@exons) > 0 ) {
    if( $exons[0]->length + $len > $start ) {
       last;
    } else {
       my $discard = shift @exons;
       $len += $discard;
    }
  }

  # handle the first component exon

  if( scalar(@exons) == 0 ) {
     return @out;
  }
  
  $sf = Bio::EnsEMBL::SeqFeature->new();
  $sf->seqname($exons[0]->contig->id);
  $sf->strand($exons[0]->strand);
  $sf->start($exons[0]->start + $start - $len);

  if( $end < $len + $exons[0]->length ) {
      $sf->end($exons[0]->start + $end - $len);
      return $sf;
  } else {
      $sf->end($exons[0]->end);
      push(@out,$sf);
  }


  while( scalar(@exons) ) {
     if( $exons[0]->length + $len > $end ) {
        last;
     }
     $sf = Bio::EnsEMBL::SeqFeature->new();
     $sf->seqname($exons[0]->contig->id);
     $sf->strand($exons[0]->strand);
     $sf->start($exons[0]->start);
     $sf->start($exons[0]->end);
     push(@out,$sf);
     $len += $exons[0]->length;
  }

  if( scalar(@exons) == 0 ) {
     return @out;
  }

  # handle the last exon

  $sf = Bio::EnsEMBL::SeqFeature->new();
  $sf->seqname($exons[0]->contig->id);
  $sf->strand($exons[0]->strand);
  $sf->start($exons[0]->start);
  $sf->start($exons[0]->start + $end - $len);



  return @out;
}

=head2 add_component_Exon

 Title   : add_component_Exon
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub add_component_Exon{
   my ($self,$exon) = @_;

   if( !ref $exon || ! $exon->isa('Bio::EnsEMBL::Exon')) {
       $self->throw("$exon is not an exon");
   }

   push(@{$self->{'_component_exons'}},$exon);
}

=head2 length

 Title   : length
 Usage   : length
 Function: calculate number of  nucleotides constituting the Exon
 Example :
 Returns : a number
 Args    : none

=cut


sub length {
    my $self = shift;

    my $len =0; 

    foreach my $subexon ( $self->each_component_Exon ) {
        $len += $subexon->length;
    }
    return $len;
}

=head2 _sort_by_sticky_rank

 Title   : _sort_by_sticky_rank
 Usage   : 
 Function: put the contained exons in the right order
 Example :
 Returns : 
 Args    : 

=cut


sub _sort_by_sticky_rank {
    my $self = shift;

    my @sorted;

    @sorted= sort {$a->sticky_rank <=> $b->sticky_rank } 
      @{$self->{'_component_exons'}};
    $self->{'_component_exons'} = \@sorted;
    return 1;
}




=head1 seq

  Arg [1]  : String $seq
             You can set the seq of a sticky Exon. Omit, if you want to retrieve.
             Attachseq is probably the better way (sigh)
  Function : retrieve sequence of sticky exon from db or return stored one.
             If sequence was retrieved once, its cached.

  Returns  : Has to return Bio::Seq


=cut


sub seq {
  my $self = shift;
  my $seq = shift;

  if( defined $seq ) {
    if( $seq ) {
      $self->{'_seq'} = $seq;
    } else {
      $self->{'_seq'} = undef;
    }
  } else {
    my $seqString = "";
    for my $cExon ( $self->each_component_Exon() ) {
      $seqString .= $cExon->seq()->seq();
    }
    $self->{'_seq'} = $seqString;
  }

  return Bio::Seq->new( -seq => $self->{'_seq'} );
}




=head2 transform

  Arg  1    : Bio::EnsEMBL::Slice $slice
              make this slice coords.
  Function  : make slice coords from raw contig coords or vice versa
  Returntype: Bio::EnsEMBL::Exon (Bio::EnsEMBL::StickyExon)
  Exceptions: none
  Caller    : Gene::transform()

=cut


sub transform {
  my $self = shift;
  my $slice = shift;

  # print "Calling transform on sticky exon\n";

  if( defined $self->contig() and 
      $self->contig()->isa( "Bio::EnsEMBL::RawContig" ) )  {

    my $mapper = $slice->adaptor->db->get_AssemblyMapperAdaptor->fetch_by_type
      ( $slice->assembly_type() );

    my $dna_seq = "";
    my $mapped_start = 0;
    my $mapped_end = -1;
    my $composite_exon_strand = 0;
    my $composite_exon_phase = -1;

    my @supporting_features;

    # sort the component exons
    $self->_sort_by_sticky_rank(); 

    # and now retrieve them
    my @component_exons = $self->each_component_Exon();

    foreach my $c_exon ( @component_exons ) {

      my @mapped = $mapper->map_coordinates_to_assembly
	(
	 $c_exon->contig()->dbID,
	 $c_exon->start(),
	 $c_exon->end(),
	 $c_exon->strand()
	);

  #    print STDERR "\nStickyExon.pm transform method:\n"; 
  #    print STDERR "[Mapped start  " . $mapped[0]->start . "\t";
  #    print STDERR "Mapped end    " . $mapped[0]->end . "\t";
  #    print STDERR "Mapped strand " . $mapped[0]->strand . "]\n";
  #    print STDERR "Sticky rank : " . $c_exon->sticky_rank() . "\n";

      # exons should always transform so in theory no error check
      # necessary
      if( ! @mapped ) {
	$self->throw( "Component Sticky Exon couldnt map" );
      }
    

      # should get a gap object returned if an exon lies outside of 
      # the current slice.  Simply return the exon as is - i.e. untransformed.
      # this untransformed exon will be distinguishable as it will still have
      # contig attached to it and not a slice.
      if( $mapped[0]->isa( "Bio::EnsEMBL::Mapper::Gap" )) {
	return $self;
      }

      # use empty Slice to map to chromosomal coords
      if( ! defined $slice->chr_name() ) {
	$slice->chr_name( $mapped[0]->id() );
	$slice->chr_start( 1 );
	$slice->strand( 1 );
      }

      # concatenate the raw sequence together
      $dna_seq .= $c_exon->seq();

      #add the supporting features from the exons
      push @supporting_features, $c_exon->get_all_supporting_features(); 

   #   print STDERR $c_exon->dbID . " " . $c_exon->seq . "\n";

      # now pull out the start and end points of the newly concatenated sequence
      # if we've got the first sticky exon, store the relevant info

      if ($c_exon->sticky_rank == 1 ) {
	# this assumes that the strand of the sticky exon is
	# set by the first one - is this correct?
	$composite_exon_strand = $mapped[0]->strand();

	if ( $composite_exon_strand == 1 ) {
	  $mapped_start = $mapped[0]->start();
	  $composite_exon_phase = $c_exon->phase();
	}
	else { # for the reverse strand case it is something different
	  $mapped_end = $mapped[0]->end();
	}
      }

      # now do the end point
      # keep storing as you iterate over the component exons
      # since the exons are previously sorted based on their sticky rank
      # then the last set of stored values will be the last component exon
      else {
	if ( $mapped[0]->strand == 1 ) {
	  $mapped_end = $mapped[0]->end;
	}
	else {
	  $mapped_start = $mapped[0]->start;
	  $composite_exon_phase = $c_exon->phase();
	}
      }
    }
    

    # now build the new composite exon
    my $newexon = Bio::EnsEMBL::Exon->new();
    $newexon->start( $mapped_start - $slice->chr_start() + 1 );
    $newexon->end( $mapped_end - $slice->chr_start() + 1);
    $newexon->strand( $composite_exon_strand * $slice->strand() );
    $newexon->dbID($component_exons[0]->dbID);
    $newexon->adaptor($component_exons[0]->adaptor);

    $newexon->contig( $slice );
    $newexon->add_supporting_features(@supporting_features);
    $newexon->phase( $composite_exon_phase );
#SMJS Hack

    $newexon->end_phase($self->end_phase);


    return $newexon;
  } else {
    $self->throw( "Unexpected StickyExon in Assembly coords ..." );
  }
}


1;
