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
use FileHandle;

# Object preamble - inheriets from Bio::Root::Object
use Bio::Root::Object;

@ISA = qw(Bio::Root::Object Bio::EnsEMBL::DB::ContigI);
# new() is inherited from Bio::Root::Object

# _initialize is where the heavy stuff will happen when new is called

sub _initialize {
  my($self,@args) = @_;
  
  my $make = $self->SUPER::_initialize;
  my ($dbobj,$id,$cloneobj)=$self->_rearrange([qw(DBOBJ
						  ID
						  CLONEOBJ
						  )],@args);
  $id || $self->throw("Cannot make contig object without id");
  $dbobj || $self->throw("Cannot make contig object without db object");
  $dbobj->isa('Bio::EnsEMBL::TimDB::Obj') || 
      $self->throw("Cannot make contig object with a $dbobj object");
  $cloneobj->isa('Bio::EnsEMBL::TimDB::Clone') || 
      $self->throw("Cannot make clone object with a $cloneobj object");
  # id of contig
  $self->id($id);
  # db object
  $self->_dbobj($dbobj);
  # clone object
  $self->_cloneobj($cloneobj);


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

    # read from sequence file
    my $cloneobj=$self->_cloneobj();
    my $cloneid=$cloneobj->id;
    my $file=$cloneobj->{'_clone_dir'}."/$cloneid.seq";
    my $fh = new FileHandle;
    $fh->open($file) || 
	$self->throw("Could not open sequence file [", $file, "]");
    my $is = $fh->input_record_separator('>');
    my $flag;
    while(<$fh>){
	if(/^$id\s[^\n]+\n(.*)/s){
	    my $dna=$1;
	    $dna=~s/\n//g;
	    $dna=~s/\s//g;
	    $dna=~tr/[a-z]/[A-Z]/;
	    $self->{'seq'}=Bio::Seq->new ( -seq => $dna , -id => $id, -type => 'DNA' );
	    $flag=1;
	    last;
	}
    }
    unless($flag){
	$self->throw("Could not find contig $id in sequence file");
    }
    $fh->input_record_separator($is);
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

1;



