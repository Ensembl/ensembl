
#
# BioPerl module for Bio::EnsEMBL::DB::WriteableVirtualContig
#
# Cared for by Ewan Birney <birney@ebi.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::DB::WriteableVirtualContig - Virtual Contig which has write_Gene methods

=head1 SYNOPSIS

    #As with VC, but add extra gene_obj

    $gene_obj= $obj->gene_Obj;

    $wvc= Bio::EnsEMBL::DB::WriteableVirtualContig->new( -focuscontig => $rawcontig,
					      -focusposition => 2,
					      -ori => 1,
					      -left => 5000,
					      -right => 5000,
					      -gene_obj => $gene_obj,
					      );

    # build a gene somehow
    # This call writes the gene to the database, mapping coordinates back to
    # rawcontig coordinates
    $wvc->write_Gene($gene);



=head1 DESCRIPTION

VirtualContig is the object that allows a region of the human genome to be viewed
as being a single piece of DNA when it is made from a set of individual contigs.

This allows genes build on virtual contigs to be written back to the database sensibly. 

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::DB::WriteableVirtualContig;
use vars qw(@ISA);
use strict;

# Object preamble - inheriets from Bio::Root::Object

use Bio::EnsEMBL::DB::VirtualContig;

@ISA = qw(Bio::EnsEMBL::DB::VirtualContig);


# _initialize is where the heavy stuff will happen when new is called

sub _initialize {
  my($self,@args) = @_;
  
  my ($gene_obj) = $self->_rearrange([qw( GENE_OBJ )],@args);
  
  if( !ref $gene_obj || !$gene_obj->isa('Bio::EnsEMBL::DBSQL::Gene_Obj') ) {
      $self->throw("must pass gene_obj argument to WriteableVirtualContig, got a [$gene_obj]");
  }

  my $make = $self->SUPER::_initialize(@args);
  $self->_gene_obj($gene_obj);

  
 
  return $make; # success - we hope!
}

=head2 write_Gene

 Title   : write_Gene
 Usage   : $wvc->write_Gene($gene)
 Function: Writes a gene built on this VirtualContig back to a set of 
           raw contig positions

           Internally this builds a copy of the genes,transcripts and translations
           in RC coordinate space, making heavy use of VirtualMap vcpos_to_rcpos
           function. Exons could be split into 
 Example :
 Returns : 
 Args    :


=cut

sub write_Gene{
   my ($self,$gene) = @_;

   if( !ref $gene || ! $gene->isa('Bio::EnsEMBL::Gene') ) {
       $self->throw("Got to write a gene, not a [$gene]");
   }

   # we need to map the exons back into RC coordinates.
   
   # in some cases the mapping will cause us to get a large number
   # of sticky exons.

   # the basically means we have to clone all the objects. :(

   my $clonedgene = Bio::EnsEMBL::Gene->new();
   my %translation;

   $clonedgene->id($gene->id);
   $clonedgene->version($gene->version);
   $clonedgene->created($gene->created);
   $clonedgene->modified($gene->modified);
   
   foreach my $trans ( $gene->each_Transcript ) {
       my $clonedtrans = Bio::EnsEMBL::Transcript->new();
       $clonedtrans->id($trans->id);
       $clonedtrans->version($trans->version);
       $clonedtrans->created($trans->created);
       $clonedtrans->modified($trans->modified);

       $clonedgene->add_Transcript($clonedtrans);

       foreach my $exon ( $trans->each_Exon ) {
	   my @clonedexons = $self->_reverse_map_Exon($exon);
	   foreach my $ce ( @clonedexons ) {
	       $clonedtrans->add_Exon($ce);
	   }
	   
	   # translations
	   if( exists $translation{$trans->translation->id} ) {
	       $clonedtrans->translation($translation{$trans->translation->id});
	   } else {
	       my $trl = $trans->translation(); 
	       my $clonedtrl = Bio::EnsEMBL::Translation->new();
	       $clonedtrl->id($trl->id);
	       $clonedtrl->created($trl->created);
	       $clonedtrl->modified($trl->modified);
	       $clonedtrl->start_exon_id($trl->start_exon_id);
	       $clonedtrl->end_exon_id($trl->end_exon_id);

	       my ($srawcontig,$start,$sstrand) = $self->_vmap->vcpos_to_rcpos($trl->start,1);
	       $clonedtrl->start($start);
	       my ($erawcontig,$end,$estrand) = $self->_vmap->vcpos_to_rcpos($trl->end,1);
	       $clonedtrl->end($end);
	       
	       $translation{$trl->id} = $clonedtrl;
	       $clonedtrans->translation($clonedtrl);
	   }
       }
   }

   $self->_gene_obj->write($clonedgene);
}

=head2 _reverse_map_Exon

 Title   : _reverse_map_Exon
 Usage   : (@exons) = $self->_reverse_map_Exon($exon)
 Function: Makes exons in RawContig coordinates from exon in VC coordinates.
           Multiple Exons might be returned when the Exons are made sticky
           due to exon crossing clone boundaries.
 Example :
 Returns : 
 Args    :


=cut

sub _reverse_map_Exon{
   my ($self,$exon) = @_;

   if( !ref $exon || !$exon->isa('Bio::EnsEMBL::Exon') ) {
       $self->throw("Must supply reverse map an exon not an [$exon]");
   }

   my ($scontig,$start,$sstrand) = $self->_vmap->vcpos_to_rcpos($exon->start,$exon->strand);
   my ($econtig,$end,$estrand)   = $self->_vmap->vcpos_to_rcpos($exon->end  ,$exon->strand);

   if( !ref $scontig || !ref $econtig || !$scontig->isa('Bio::EnsEMBL::DB::RawContigI') || !$econtig->isa('Bio::EnsEMBL::DB::RawContigI') ) {
       $self->throw("Exon on vc ".$exon->id." [".$exon->start.":".$exon->end."] is unmappable to rawcontig positions, probably being in a gap. Can't write");
   }

  
   if( $scontig->id eq $econtig->id ) {
       if( $sstrand != $estrand ) {
	   $self->throw("Bad internal error. Exon mapped to same contig but different strands!");
       }

       my $rmexon = Bio::EnsEMBL::Exon->new();
       $rmexon->id($exon->id);
       $rmexon->created($exon->created);
       $rmexon->modified($exon->modified);
       $rmexon->version($exon->version);
       $rmexon->phase($exon->phase);
       $rmexon->sticky_rank(1);
       $rmexon->start($start);
       $rmexon->end($end);
       $rmexon->strand($sstrand);
       $rmexon->contig_id($scontig->id);
       $rmexon->seqname($scontig->id);
       $rmexon->attach_seq($scontig->primary_seq);
       return ($rmexon);
   } else {
       # we are in the world of sticky-ness....


       my @mapcontigs = $self->_vmap->get_all_MapContigs();

       # walk to find scontig
       my $found = 0;
       foreach my $mc ( @mapcontigs ) {
	   print STDERR "Got ".$mc->contig->id.":".$scontig->id."\n";

	   if( $mc->contig->id eq $scontig->id ) {
	       unshift(@mapcontigs,$mc);
	       $found = 1;
	       last;
	   }
       }
       if( $found == 0 ) {
	   $self->throw("Internal error - unable to find map contig with this id");
       }


       my $vcstart = $exon->start;

       # ok. Move from start towards end, after we hit end.
       my @exported_exons;
       my $sticky = 1;

       foreach my $c ( @mapcontigs ) {	   
	   my $vcend;
	   print STDERR "***Looking at $c\n";

	   if( $c->contig->id eq $econtig->id ) {
	       # go to end position
	       $vcend = $exon->end();
	   } else {
	       $vcend = $c->end();
	   }

	   print STDERR "Going to call with $start:$end\n";
	   $self->_dump_map(\*STDERR);

	   my ($iscontig,$istart,$isstrand) = $self->_vmap->vcpos_to_rcpos($vcstart,$exon->strand);
	   my ($iecontig,$iend,$iestrand)   = $self->_vmap->vcpos_to_rcpos($vcend  ,$exon->strand);
  
	   if( $iscontig->id ne $iecontig->id || $isstrand != $iestrand) {
	       $self->throw("Bad internal error. Sticky Exon mapped to different contig/strand for a correct contig placement ".$iscontig->id.":".$iecontig->id);
	   }

	   my $rmexon = Bio::EnsEMBL::Exon->new();
	   $rmexon->id($exon->id);
	   $rmexon->created($exon->created);
	   $rmexon->modified($exon->modified);
	   $rmexon->version($exon->version);
	   $rmexon->phase($exon->phase);
	   $rmexon->sticky_rank($sticky++);
	   $rmexon->attach_seq($c->contig->primary_seq);	   
	   $rmexon->start($istart);
	   $rmexon->end($iend);
	   $rmexon->strand($isstrand);
	   $rmexon->contig_id($c->contig->id);
	   $rmexon->seqname($c->contig->id);
	   push(@exported_exons,$rmexon);

	   if( $c->contig->id eq $econtig->id ) {
	       last;
	   }
	   
       }

       return @exported_exons;
   }
       
   $self->throw("Internal error. Should not reach here!");
}



=head2 _gene_obj

 Title   : _gene_obj
 Usage   : $obj->_gene_obj($newval)
 Function: 
 Returns : value of _gene_obj
 Args    : newvalue (optional)


=cut

sub _gene_obj{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'_gene_obj'} = $value;
    }
    return $obj->{'_gene_obj'};

}
