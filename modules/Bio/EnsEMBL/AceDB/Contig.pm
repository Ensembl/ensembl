REVISION 1.3 WAS LOST.  JUNK LINE INSERTED IN DELTATEXT
REVISION 1.4 WAS LOST.  JUNK LINE INSERTED IN DELTATEXT
REVISION 1.5 WAS LOST.  JUNK LINE INSERTED IN DELTATEXT
REVISION 1.6 WAS LOST.  JUNK LINE INSERTED IN DELTATEXT

#
# BioPerl module for Contig
#
# Cared for by Ewan Birney <birney@sanger.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::AceDB::Contig - Handle onto a database stored contig

=head1 SYNOPSIS

Give standard usage here

=head1 DESCRIPTION

Describe the object here

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::AceDB::Contig;
use vars qw(@ISA);
use strict;

# Object preamble - inheriets from Bio::Root::Object

use Bio::Root::Object;
use Bio::SeqFeature::Generic;
use Bio::EnsEMBL::AceDB::Obj;
use Bio::EnsEMBL::DB::ContigI;
use Bio::Seq;
use Bio::EnsEMBL::Gene;

@ISA = qw(Bio::Root::Object Bio::EnsEMBL::DB::ContigI);
# new() is inherited from Bio::Root::Object

# _initialize is where the heavy stuff will happen when new is called

sub _initialize {
  my($self,@args) = @_;

  my $make = $self->SUPER::_initialize(@args);

  my ($dbobj,$id) = $self->_rearrange([qw(DBOBJ
					  ID
					  )],@args);

  $id || $self->throw("Cannot make contig db object without id");
  $dbobj || $self->throw("Cannot make contig db object without db object");
  $dbobj->isa('Bio::EnsEMBL::AceDB::Obj') || $self->throw("Cannot make contig db object with a $dbobj object");

  $self->id($id);
  $self->_dbobj($dbobj);

# set stuff in self from @args
  return $make; # success - we hope!
}

=head2 seq

 Title   : seq
 Usage   : $seq = $contig->seq();
 Function: Gets a Bio::Seq object out from the contig
 Example :
 Returns : Bio::Seq object
 Args    :


=cut

sub seq{
   my ($self) = @_;
   my $id = $self->id();

   my $seq = $self->_dbobj->fetch('Sequence',$id) || 
       $self->throw("Could not retrieve $id from acedb" . Ace->error());
  
   my $dna = $seq->asDNA || $self->throw("Could not retrieve DNA from $id");

   $dna =~ s/^>.*\n//g;
   $dna =~ s/\s//g;
   $dna =~ tr/[a-z]/[A-Z]/;
   my $out = Bio::Seq->new ( -seq => $dna , -id => $id, -type => 'DNA' ) ;
   return $out;
}

=head2 get_all_SeqFeatures

 Title   : get_all_SeqFeatures
 Usage   : foreach my $sf ( $contig->get_all_SeqFeatures ) 
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_all_SeqFeatures {
   my ($self,@args) = @_;
   my @array;

   my $id = $self->id();
   my $seq;

   ($seq) = $self->_dbobj->fetch( 'Sequence', $id) || $self->throw("Could not retrieve sequence object $id");

   # FIXME - should be a subroutine perhaps?

   # The maps convert the ACeDB objects into start/end
   # values

   foreach my $pep ($seq->at('Homol.Pep_homol[1]')) { 
       foreach my $hit ($pep->at('BLASTX[1]')) {
	   my $sf = &_from_ace_seqfeature(map "$_",$hit->row());
	   $sf->primary_tag('BLASTX');
	   $sf->source_tag('EnsEMBL');
	   push(@array,$sf);
       }
   }

   foreach my $est ($seq->at('Homol.EST_homol[1]')) { 
       foreach my $hit ($est->at('EST_eg[1]')) {
	   my $sf = &_from_ace_seqfeature(map "$_",$hit->row());
	   $sf->primary_tag('TBLASTX_ESTGENOME');
	   $sf->source_tag('EnsEMBL');
	   push(@array,$sf);
       }
   }

   return @array;
}

=head2 order

 Title   : order
 Usage   : $order = $contig->order()
 Function: Provides the order of the contig, starting at
         : zero.
 Example :
 Returns : 
 Args    :


=cut

sub order {
   my ($self) = @_;

   # All AceDB contigs are single contigs in a single clone.
   # easy!
   return 0; 

}

=head2 get_all_Genes

 Title   : get_all_Genes
 Usage   : @genes = $contig->get_all_Genes
 Function:
 Example : 
 Returns : 
 Args    :
 Note    : WARNING. At the moment this just
           gets genes from the Sequence object, not
           any potential genes in link objects above
           this sequence object. To be fixed!

=cut

sub get_all_Genes {
   my ($self) = @_;

   my @out;
   my $id = $self->id();
   my $seq;
   my $bioseq = $self->seq();

   ($seq) = $self->_dbobj->fetch( 'Sequence', $id) || $self->throw("Could not find sequence object $id");

   # here we don't look at Locus objects to determine whether
   # things are transcripts or not

   foreach my $sub ($seq->at('Structure.Subsequence')) {
       my $genename = "$sub";

       if( $sub->fetch->at("Method[1]") ne 'supported_CDS' ) {
	   next;
       }

       my ($start,$end) = $sub->row(1);
       my $strand;
       if( $start > $end ) {
	   $strand = -1;
       } else {
	   $strand = 1;
       }

       my $gene = new Bio::EnsEMBL::Gene;
       push(@out,$gene);
       my $trans = new Bio::EnsEMBL::Transcript;
       
       $gene->add_Transcript($trans);
       my $subseq = $sub->fetch();
       my $tag = $sub->fetch("Properties.Start_not_found[1]");
       my $total;
       if( $tag ) {
	   $total = $tag -1;
       } else {
	   $total = 0;
       }

       
       my $index=1;
       foreach my $hit ($subseq->at('Structure.From.Source_Exons[1]')) {

	   # we have to map acedb coordinates which are relative to the
	   # start/end in the subsequence to the exon coordinates, which
	   # are absolute.

	   my ($starte,$ende) = map("$_",$hit->row());
	   my $exon = new Bio::EnsEMBL::Exon;
	   print STDERR "Exon with $start and $starte and $end\n";
	   if( $strand == 1 ) {
	       $exon->start($start+$starte-1);
	       $exon->end($start+$ende-1);
	   } else {
	       $exon->start($start-$ende+1);
	       $exon->end($start-$starte+1);
	   }

	   # calculate phase - ensembl thinks of exons as being
	   # independent from other exons, so phase is important.
	   $exon->strand($strand);
	   $exon->phase($total % 3);
	   $total += ($end-$start+1);

	   # random bits and bobs

	   $exon->clone_id($id);
	   $exon->contig_id($id);
	   $exon->attach_seq($bioseq);

	   my $exonid;

	   if( $self->_dbobj->_exon_id_start() ) {
	       $exonid = $self->_dbobj->_exon_id_start();
	       my $nexte = $exonid++;
	       $self->_dbobj->_exon_id_start($exonid);
	   } else {
	       $exonid = "dummy_exon_id.$genename.$index";
	   }

	   $exon->id($exonid);
	   $exon->created("1999-07-12");
	   $exon->modified("1999-07-12");
	   $index++;
	   $trans->add_Exon($exon);
       }
   }

   return @out;
}
   

sub _from_ace_seqfeature {
    my ($score,$start,$end,$tstart,$tend) = @_;
    my $strand;
    my $out;

    if( $start !~ /^\d+$/ || $end !~ /^\d+$/ || $score !~ /^\d+\.?\d*$/ ) {
	&Bio::Root::Object::throw("start $start, end $end and score $score look dodgy");
    }

    if( $start > $end ) {
	my $temp =$end;
	$end = $start;
	$start = $temp;
	$strand = -1;
    } else {
	$strand = +1;
    }
    
    my $out = new Bio::SeqFeature::Generic;
    $out->start($start);
    $out->end($end);
    $out->strand($strand);
    $out->score($score);
    
    return $out;
}

=head2 id

 Title   : id
 Usage   : $obj->id($newval)
 Function: 
 Example : 
 Returns : value of id
 Args    : newvalue (optional)


=cut

sub id{
    my ($self,$value) = @_;
    if( defined $value) {
	$self->{'id'} = $value;
    }
    return $self->{'id'};
}

=head2 offset

 Title   : offset
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub offset{
   my ($self,@args) = @_;

   return 1;
}

=head2 orientation

 Title   : orientation
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub orientation{
   my ($self,@args) = @_;
   
   return 1;

}

=head2 _dbobj

 Title   : _dbobj
 Usage   : $obj->_dbobj($newval)
 Function: 
 Example : 
 Returns : value of _dbobj
 Args    : newvalue (optional)


=cut

sub _dbobj{
   my ($self,$value) = @_;
   if( defined $value) {
      $self->{'_dbobj'} = $value;
    }
    return $self->{'_dbobj'};

}

1;
