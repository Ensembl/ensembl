#
# Ensembl module for Bio::EnsEMBL::Utils::EMBL::GeneWrapper
#
# Cared for by Ewan Birney <ensembl-dev@ebi.ac.uk>
#
# Copyright EMBL / EBI
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::Utils::EMBL::GeneWrapper

=head1 SYNOPSIS

   foreach $gene ( $slice->get_all_Genes ) {
     push @seq_feats, new Bio::EnsEMBL::Utils::EMBL::GeneWrapper($gene);
   }
       

=head1 DESCRIPTION

Allows genes to be EMBL dumped as though they were seq features (which in 
EnsEMBL they are not).  Essentially a workaround to force ensembl to fit 
the Bioperl embl dumping code. This wrapper is a replacement for some of the
old VirtualGene functionality.

=head1 AUTHOR - Graham McVicker

This modules is part of the Ensembl project http://www.ensembl.org

Email ensembl-dev@ebi.ac.uk


=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal 
methods are usually preceded with a _

=cut


use strict;

package Bio::EnsEMBL::Utils::EMBL::GeneWrapper;

use Bio::SeqIO::FTHelper;
use Bio::EnsEMBL::Root;
use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Root);


=head2 new

  Arg [1]    : Bio::EnsEMBL::Gene
  Example    : $gw = new Bio::EnsEMBL::Utils::EMBL::GeneWrapper($gene)
  Description: Creates a new GeneWrapper 'wrapped' around a gene object to
               allow it to be embl dumped. 
  Returntype : Bio::EnsEMBL::Utils::EMBL::GeneWrapper
  Exceptions : none
  Caller     : Embl dumping scripts

=cut

sub new {
  my ($class, $gene) = @_;

  my $self = $class->SUPER::new();

  unless($gene && ref $gene && $gene->isa('Bio::EnsEMBL::Gene')) {
    $self->throw("Gene argument must be a Bio::EnsEMBL::Gene");
  }

  $self->gene($gene);
  $self->{_strict_embl_dumping} = 0;

  return $self;
}


=head2 gene

  Arg [1]    : (optional) Bio::EnsEMBL::Gene $gene
  Example    : $gene = $gene_wrapper->gene();
  Description: getter/setter fot the gene contained inside this wrapper
  Returntype : Bio::EnsEMBL::Gene
  Exceptions : thrown if the $gene arg is not a Bio::EnsEMBL::Gene
  Caller     : internal

=cut

sub gene {
  my ($self, $gene) = @_;

  if($gene) {
    unless($gene->isa('Bio::EnsEMBL::Gene')) {
      $self->throw("[$gene] is not a Bio::EnsEMBL::Gene");
    }
    
    $self->{_gene} = $gene;
  }

  return $self->{_gene};
}


=head2 to_FTHelper

  Arg [1]    : none
  Example    : none
  Description: Wrapper function - 
               This is where all the embl dumping magic around genes takes 
               place. Most of the embl output for ensembl genes is defined 
               here.
  Returntype : list of Bio::SeqIO::FTHelper
  Exceptions : none
  Caller     : Bio::SeqIO::embl

=cut

sub to_FTHelper {
  my $self = shift;

  my @out;

  my @dblinks = @{$self->gene->get_all_DBLinks()};

  foreach my $trans (@{$self->gene->get_all_Transcripts()}) {
    my $join = "";
    
    foreach my $exon (@{$trans->get_all_translateable_Exons()}) {
      $join .= ',' if($join); #append a comma to the last coord set

      if($exon->strand() == 1) {
	$join .= $exon->start()."..".$exon->end();
      } else {
	$join .= "complement(".$exon->start()."..".$exon->end().")";
      }
    }

    my $ft = Bio::SeqIO::FTHelper->new();
    $ft->loc("join(".$join.")");
    $ft->key('CDS');

    $ft->add_field('translation', $trans->translate()->seq());
    $ft->add_field('cds', $trans->translation->stable_id());
    $ft->add_field('gene', $self->gene->stable_id());
    $ft->add_field('transcript', $trans->stable_id());
    foreach my $dbl (@dblinks) {
      $ft->add_field('db_xref', $dbl->database().":".$dbl->primary_id());
    }
    push(@out, $ft);
  }

  foreach my $exon (@{$self->gene->get_all_Exons()}) {
    my $ft = Bio::SeqIO::FTHelper->new();

    if( $exon->strand() == 1 ) {
      $ft->loc($exon->start."..".$exon->end());
    } else {
      $ft->loc("complement(".$exon->start()."..".$exon->end().")");
    }

    $ft->key("exon");
    
    if( $self->strict_EMBL_dumping()) {
      $ft->add_field('db_xref', 'ENSEMBL:HUMAN-Exon-'.$exon->stable_id());
    } else {
      $ft->add_field('exon_id', $exon->stable_id());
      $ft->add_field('start_phase', $exon->phase());
      $ft->add_field('end_phase', $exon->phase());
    }
    
    push (@out, $ft);
  }
  
  return @out;
}

  
=head2 strict_EMBL_dumping

  Arg [1]    : optional boolean $newval
  Example    : none
  Description: Getter/ setter for strict embl dumping attribute which 
               determines whether embl dumping will be 'strict'. Default is 0.
  Returntype : boolean
  Exceptions : none
  Caller     : to_FTHelper

=cut

sub strict_EMBL_dumping {
  my ($self, $newval) = @_;
  
  if(defined $newval) {
    $self->{_strict_embl_dumping} = $newval;
  }

  return $self->{_strict_embl_dumping};
}
    

=head2 source_tag

  Arg [1]    : none
  Example    : none
  Description: Wrapper function - returns a source tag for the gene
  Returntype : none
  Exceptions : none
  Caller     : Bio::SeqIO::embl

=cut

sub source_tag {
  my ($self, @args) = @_;

  return "ensembl";
}


=head2 primary_tag

  Arg [1]    : none
  Example    : none
  Description: Wrapper function returns a primary tag for thie gene
  Returntype : string
  Exceptions : none
  Caller     : Bio::SeqIO::embl

=cut

sub primary_tag {
  my ($self, @args) = @_;

  return "genefragment";
}


=head2 has_tag

  Arg [1]    : none
  Example    : none
  Description: Wrapper function - always returns 0
  Returntype : none
  Exceptions : none
  Caller     : Bio::SeqIO::embl

=cut

sub has_tag {
  my $self = shift;

  return 0;
}


=head2 all_tags

  Arg [1]    : none
  Example    : none
  Description: Wrapper function - does nothing
  Returntype : none
  Exceptions : none
  Caller     : Bio::SeqIO::embl

=cut

sub all_tags {
  my $self = shift;

  return;
}


=head2 each_tag_value

  Arg [1]    : none
  Example    : none
  Description: Wrapper function - returns empty list
  Returntype : empty list
  Exceptions : none
  Caller     : Bio::SeqIO::embl

=cut

sub each_tag_value {
  my $self = shift;

  return ();
}



1;
