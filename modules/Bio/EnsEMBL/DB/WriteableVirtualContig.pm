
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

   # sanity check 

   $self->_sanity_check($gene);

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
   foreach my $dbl ( $gene->each_DBLink() ) {
       $clonedgene->add_DBLink($dbl);
   }

   foreach my $trans ( $gene->each_Transcript ) {
       my $clonedtrans = Bio::EnsEMBL::Transcript->new();
       $clonedtrans->id($trans->id);
       $clonedtrans->version($trans->version);
       $clonedtrans->created($trans->created);
       $clonedtrans->modified($trans->modified);
       foreach my $dbl ( $trans->each_DBLink() ) {
	   $clonedtrans->add_DBLink($dbl);
       }

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
	       $clonedtrl->start_exon_id($trl->start_exon_id);
	       $clonedtrl->end_exon_id($trl->end_exon_id);
	       $clonedtrl->version($trl->version);

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

   my ($scontig,$start,$sstrand) = $self->_vmap->raw_contig_position($exon->start,$exon->strand);
   my ($econtig,$end,$estrand)   = $self->_vmap->raw_contig_position($exon->end  ,$exon->strand);

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
       foreach my $se ( $exon->each_Supporting_Feature ) {
	   my ($secontig,$sestart,$sestrand) = $self->_vmap->raw_contig_position($se->start,$se->strand);
	   my ($sncontig,$seend,$snstrand) = $self->_vmap->raw_contig_position($se->start,$se->strand);
	   if( !ref $secontig || !ref $sncontig || $secontig->id ne $sncontig->id ) {
	       $self->warn("supporting evidence spanning contigs. Cannot write");
	       next;
	   }
	   if( $sestart < $seend ) {
	       $se->start($sestart);
	       $se->end($seend);
	   } else {
	       $se->start($sestart);
	       $se->end($seend);
	   }
	   $se->strand($sestrand);
	   $se->seqname($secontig->id);
	   if( $se->can('attach_seq') ) {
	       $se->attach_seq($secontig->primary_seq);
	   }

	   $rmexon->add_Supporting_Feature($se);
       }

       # we could test on strand changes. This just assummes everything works
       # as it says on the tin ;)
       if( $start < $end ) {
	   $rmexon->start($start);
	   $rmexon->end($end);
       } else {
	   $rmexon->start($end);
	   $rmexon->end($start);
       }
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
       my $mc;
       while ( $mc = shift @mapcontigs ) { 
 
	   print STDERR "Got ".$mc->contig->id.":".$scontig->id."\n";

	   if( $mc->contig->id eq $scontig->id ) {
	       print STDERR "Unshifting ",$mc->contig->id,"\n";
	       unshift(@mapcontigs,$mc);
	       $found = 1;
	       last;
	   }
       }
       if( $found == 0 ) {
	   $self->throw("Internal error - unable to find map contig with this id");
       }


       my $vcstart = $exon->start;
       print STDERR "Looking from exon-wise",$exon->start,":",$exon->end,"\n";

       # ok. Move from start towards end, after we hit end.
       my @exported_exons;
       my $sticky = 1;

       foreach my $c ( @mapcontigs ) {	   
	   my $vcend;
	   print STDERR "***Looking at",$c->contig->id," - $vcstart...\n";

	   if( $c->contig->id eq $econtig->id ) {
	       # go to end position
	       print STDERR "Going for end...",$econtig->id,"\n";
	       $vcend = $exon->end();
	   } else {
	       print STDERR "Going for end of contig\n";
	       $vcend = $c->end();
	   }

	   print STDERR "....Going to call with $vcstart:$vcend\n";
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
	   $vcstart = $vcend+1;
       }

       return @exported_exons;
   }
       
   $self->throw("Internal error. Should not reach here!");
}

=head2 _sanity_check

 Title   : _sanity_check
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub _sanity_check{
   my ($self,$gene) = @_;
   my $error =0;
   my $message;
   if( !defined $gene->id ) {
       $error = 1;
       $message .= "Gene has no id;";
   }
   if( !defined $gene->version ) {
       $error = 1;
       $message .= "Gene has to have a version;";
   }
   foreach my $trans ( $gene->each_Transcript ) {
       if( !defined $trans->id ) {
	   $error = 1;
	   $message .= "Transcript has no id;";
       }
       if( !defined $trans->translation || !ref $trans->translation) {
	   $error = 1;
	   $message .= "Transcript has no translation;";
       } else {
	   if( !defined $trans->translation->id ) {
	       $error = 1;
	       $message .= "Translation has no id";
	   } 
	   if( !defined $trans->translation->start ) {
	       $error = 1;
	       $message .= "Translation has no start";
	   } 
	   if( !defined $trans->translation->start_exon_id ) {
	       $error = 1;
	       $message .= "Translation has no start exon id";
	   } 
	   if( !defined $trans->translation->end ) {
	       $error = 1;
	       $message .= "Translation has no end";
	   } 
	   if( !defined $trans->translation->end_exon_id ) {
	       $error = 1;
	       $message .= "Translation has no end exon id";
	   } 
       }
       foreach my $exon ( $trans->each_Exon ) {
	   if( !defined $exon->id ) {
	       $error = 1;
	       $message .= "Exon has no id";
	   } 
	   if( !defined $exon->created ) {
	       $error = 1;
	       $message .= "Exon has no id";
	   } 
	   if( !defined $exon->modified ) {
	       $error = 1;
	       $message .= "Exon has no id";
	   } 
	   if( !defined $exon->contig_id  ) {
	       $error = 1;
	       $message .= "Exon has no contig id";
	   } else {
	       if( $exon->contig_id ne $self->id ) {
		   $error = 1;
		   $message .= "Exon [".$exon->id."] does not seem to be on this VirtualContig";
	       }
	   }
	   if( !defined $exon->start || !defined $exon->end || !defined $exon->strand || !defined $exon->phase ) {
	       $error = 1;
	       $message .= "Exon has error in start/end/strand/phase";
	   } 
       }
   }

   if( $error == 1 ) {
       $self->throw("Cannot write gene due to: $message");
   }
	   

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
