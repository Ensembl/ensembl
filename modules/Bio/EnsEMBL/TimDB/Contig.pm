#
# BioPerl module for Bio::EnsEMBL::TimDB::Contig
#
# Cared for by Ewan Birney <birney@sanger.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::TimDB::Contig - Perl wrapper over Tim's directories for Contigs

=head1 SYNOPSIS

    $contig = Bio::EnsEMBL::TimDB::Contig->new();
    
    # $sf is a Bio::SeqFeatureI type object. $seq is a Bio::Seq object

    $contig->add_SeqFeature($sf);
    $contig->seq($seq); 

=head1 DESCRIPTION

Tim's contigs

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...
package Bio::EnsEMBL::TimDB::Contig;
use vars qw($AUTOLOAD @ISA);
use strict;
use Bio::EnsEMBL::DB::ContigI;
use Bio::Seq;
use Bio::SeqIO::Fasta;
use Bio::EnsEMBL::Analysis::Genscan;
use FileHandle;

# Object preamble - inheriets from Bio::Root::Object
use Bio::Root::Object;

@ISA = qw(Bio::Root::Object Bio::EnsEMBL::DB::ContigI);
# new() is inherited from Bio::Root::Object

# _initialize is where the heavy stuff will happen when new is called

sub _initialize {
  my($self,@args) = @_;
  
  my $make = $self->SUPER::_initialize;
  my ($dbobj,$id,$disk_id,$cloneobj,$clone_dir)=$self->_rearrange([qw(DBOBJ
							   ID
							   DISK_ID
							   CLONEOBJ
							   CLONE_DIR
							   )],@args);
  $id || $self->throw("Cannot make contig object without id");
  $disk_id || $self->throw("Cannot make contig object without disk_id");
  $dbobj || $self->throw("Cannot make contig object without db object");
  $dbobj->isa('Bio::EnsEMBL::TimDB::Obj') || 
      $self->throw("Cannot make contig object with a $dbobj object");
  $cloneobj->isa('Bio::EnsEMBL::TimDB::Clone') || 
      $self->throw("Cannot make clone object with a $cloneobj object");
  # id of contig
  $self->id($id);
  $self->disk_id($disk_id);
  # db object
  $self->_dbobj($dbobj);
  # clone object
  $self->_cloneobj($cloneobj);
  $self->_clone_dir($clone_dir);



  # ok. Hell. We open the Genscan file using the Genscan object.
  # this is needed to remap the exons lower down

  my $gs = Bio::EnsEMBL::Analysis::Genscan->new($self->_clone_dir . "/" . $self->disk_id . ".gs",$self->seq());

  # we yank out each exon and build a hash on start position

  my %exhash;

  foreach my $t ( $gs->each_Transcript ) {
      foreach my $ex ( $t->each_Exon ) {
	  $exhash{$ex->start} = $ex;
      }
  }
  

  # build array of genes
  $self->{'_gene_array'} = [];
  {
      my $bioseq=$self->seq;
      # 1. loop over list of exons in this contig
      my %transcript_id;
      my %gene_id;
      foreach my $exon (@{$dbobj->{'_contig2exon'}->{$id}}){
	  my $exon_id=$exon->id;
	  $exon->attach_seq($bioseq);
	  
	  if( ! defined $exhash{$exon->start()} ) {
	      $self->warn("No exon in in genscan file. Ugh");
	      next;
	  } 
	  if( $exhash{$exon->start()}->end != $exon->end() ) {
	      $self->throw("Exons with same start but different end!");
	  }

	  $exon->phase($exhash{$exon->start()}->phase);

	  # 2. build list of transcripts containing these exons
	  foreach my $transcript (@{$dbobj->{'_exon2transcript'}->{$exon_id}}){
	      $transcript_id{$transcript->id()}=$transcript;
	  }
      }
      # 3. build list of genes containing these transcripts
      foreach my $transcript_id (keys %transcript_id){
	  foreach my $gene (@{$dbobj->{'_transcript2gene'}->{$transcript_id}}){
	      $gene_id{$gene->id()}=$gene;
	  }
      }
      foreach my $gene (values %gene_id){
	  push(@{$self->{'_gene_array'}},$gene);
      }
  }

  # FIXME
  # not implemented here or elsewhere
  $self->{'_sf_array'} = [];
 
  # set stuff in self from @args
  return $make; # success - we hope!
}


=head2 get_all_SeqFeatures

 Title   : get_all_SeqFeatures
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_all_SeqFeatures{
   my ($self) = @_;

   $self->throw("Tim has not reimplemented this function");

   return @{$self->{'_sf_array'}};
}

=head2 get_all_Genes

 Title   : get_all_Genes
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_all_Genes{
    my ($self) = @_;
    return @{$self->{'_gene_array'}};
}


=head2 add_SeqFeature

 Title   : add_SeqFeature
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :

=cut

sub add_SeqFeature{
   my ($self,$sf) = @_;
   $self->throw("SeqFeatures cannot be added in TimDB");
}


=head2 add_Gene

 Title   : add_Gene
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :

=cut

sub add_Gene{
    my ($self,$gene) = @_;
    $self->throw("Genes cannot be added to TimDB");
}


=head2 offset

 Title   : offset
 Usage   : $self->offset($newval)
 Function: 
 Returns : value of offset
 Args    : newvalue (optional)

=cut

sub offset{
    my $self = shift;

    # for now this only works if there is only one contig
    if(scalar($self->_cloneobj()->get_all_Contigs())!=1){
	$self->throw("Tim has not reimplemented this function for >1 contig");
    }
    1;
}


=head2 orientation

 Title   : orientation
 Usage   : $self->orientation($newval)
 Function: 
 Returns : value of orientation
 Args    : newvalue (optional)

=cut

sub orientation{
    my $self = shift;
    # for now this only works if there is only one contig
    if(scalar($self->_cloneobj()->get_all_Contigs())!=1){
	$self->throw("Tim has not reimplemented this function for >1 contig");
    }
    1;
}


=head2 seq

 Title   : seq
 Usage   : $self->seq($newval)
 Function: 
 Returns : value of seq
 Args    : newvalue (optional)

=cut

sub seq{
    my $self = shift;
    
    if( @_ ) {
	$self->throw("Cannot set a sequence in TimDB");
    }elsif(defined $self->{'seq'}){
	return $self->{'seq'};
    }

    my $id=$self->id;
    my $disk_id=$self->disk_id;

    # read from sequence file
    my $cloneobj=$self->_cloneobj();
    my $clonediskid=$cloneobj->disk_id;
    my $file=$cloneobj->{'_clone_dir'}."/$clonediskid.seq";

    my $seqin = Bio::SeqIO::Fasta->new( -file => $file);
    $self->{'seq'} = $seqin->next_seq();
    $self->{'seq'}->type('Dna');

    if( ! $self->{'seq'} ) {
	$self->throw("Could not read sequence in $file");
    }

    return $self->{'seq'};
}


=head2 id

 Title   : id
 Usage   : $self->id($newval)
 Function: 
 Returns : value of id
 Args    : newvalue (optional)

=cut

sub id{
    my $self = shift;
    if( @_ ) {
	my $value = shift;
	$self->{'id'} = $value;
    }
    return $self->{'id'};
}

sub disk_id{
    my $self = shift;
    if( @_ ) {
	my $value = shift;
	$self->{'disk_id'} = $value;
    }
    return $self->{'disk_id'};
}


=head2 _dbobj

 Title   : _dbobj
 Usage   : $obj->_dbobj($newval)
 Function: 
 Example : 
 Returns : value of _dbobj
 Args    : newvalue (optional)

=cut

sub _dbobj {
    my ($obj,$value) = @_;
    if( defined $value) {
	$obj->{'_dbobj'} = $value;
    }
    return $obj->{'_dbobj'};
}

=head2 _cloneobj

 Title   : _cloneobj
 Usage   : $obj->_cloneobj($newval)
 Function: 
 Example : 
 Returns : value of _cloneobj
 Args    : newvalue (optional)

=cut

sub _cloneobj {
    my ($obj,$value) = @_;
    if( defined $value) {
	$obj->{'_cloneobj'} = $value;
    }
    return $obj->{'_cloneobj'};
}

=head2 _clone_dir

 Title   : _clone_dir
 Usage   : $obj->_clone_dir($newval)
 Function: 
 Returns : value of _clone_dir
 Args    : newvalue (optional)


=cut

sub _clone_dir{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'_clone_dir'} = $value;
    }
    return $obj->{'_clone_dir'};

}

1;



