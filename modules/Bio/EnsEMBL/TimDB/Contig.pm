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

Bio::EnsEMBL::TimDB::Contig - Perl wrapper over Tim\'s directories for Contigs

=head1 SYNOPSIS

    $contig = Bio::EnsEMBL::TimDB::Contig->new();
    
    # $sf is a Bio::SeqFeatureI type object. $seq is a Bio::Seq object

    $contig->add_SeqFeature($sf);
    $contig->seq($seq); 

=head1 DESCRIPTION

Tim\'s contigs

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
use Bio::SeqIO;
use Bio::EnsEMBL::Analysis::Genscan;
use Bio::EnsEMBL::Analysis::FeatureParser;
use FileHandle;

# Object preamble - inheriets from Bio::Root::Object
use Bio::Root::Object;

@ISA = qw(Bio::Root::Object Bio::EnsEMBL::DB::ContigI);
# new() is inherited from Bio::Root::Object

# _initialize is where the heavy stuff will happen when new is called

sub _initialize {
  my($self,@args) = @_;
  
  my $make = $self->SUPER::_initialize;
  my ($dbobj,$id,$disk_id,$clone_dir,$order,$offset,$orientation,$length)=
      $self->_rearrange([qw(DBOBJ
			    ID
			    DISK_ID
			    CLONE_DIR
			    ORDER
			    OFFSET
			    ORIENTATION
			    LENGTH
			    )],@args);
  $id      || $self->throw("Cannot make contig object without id");
  $disk_id || $self->throw("Cannot make contig object without disk_id");
  $dbobj   || $self->throw("Cannot make contig object without db object");
  $dbobj->isa('Bio::EnsEMBL::TimDB::Obj') ||   $self->throw("Cannot make contig object with a $dbobj object");
  $order   || $self->throw("Cannot make contig object without order");
  $offset  || $self->throw("Cannot make contig object without offset");
  $orientation || $self->throw("Cannot make contig object without orientation");
  $length  || $self->throw("Cannot make contig object without length");
  
  # id of contig
  $self->id         ($id);
  $self->disk_id    ($disk_id);
  # db object
  $self->_dbobj     ($dbobj);
  # clone object
  $self->_clone_dir ($clone_dir);
  $self->order      ($order);
  $self->offset     ($offset);
  $self->orientation($orientation);
  $self->length     ($length);

  # ok. Hell. We open the Genscan file using the Genscan object.
  # this is needed to remap the exons lower down
  my $gs = Bio::EnsEMBL::Analysis::Genscan->new($self->_clone_dir . "/" . 
						$self->disk_id . ".gs",
						$self->seq());
  # save for later
  $self->_gs($gs);

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
	      $self->warn("No exon in in genscan file. Ugh [Exon $exon_id, Disk id ".
			  $self->disk_id);
	      next;
	  } 
	  if( $exhash{$exon->start()}->end != $exon->end() ) {
	      $self->throw("Exons with same start but different end!\n".
			   "Exon: $exon_id ".$exon->start."-".$exon->end.
			   " ".$exhash{$exon->start()}->id." ".
			   " ".$exhash{$exon->start()}->start.
			   "-".$exhash{$exon->start()}->end);
	  }

	  $exon->phase($exhash{$exon->start()}->phase);

	  # 2. build list of transcripts containing these exons
	  foreach my $transcript (@{$dbobj->{'_exon2transcript'}->{$exon_id}}){
	      $transcript_id{$transcript->id()}=$transcript;
	      # Now deal with adding translations!

	      #
	      # This puts in the Translation information
	      #
	      
	      my $fe = $transcript->first_exon();
	      my $le = $transcript->last_exon();
	      
	      if( !defined $fe ) {
		  $self->throw("Atempting to build a transcript with no Exons. problem!");
	      }
	      
	      my $trans = Bio::EnsEMBL::Translation->new();
	      $trans->start_exon_id($fe->id);
	      $trans->end_exon_id($le->id);
	      
	      if( $fe->strand == 1 ) {
		  $trans->start($fe->start + (3 -$fe->phase)%3 );
	      } else {
		  $trans->start($fe->end - (3 -$fe->phase)%3 );
	      }
	      
	      if( $le->strand == -1 ) {
		  $trans->end($le->start + $le->end_phase );
	      } else {
		  $trans->end($le->end - $le->end_phase );
	      }
	      my $tid = $transcript->id();
	      # horrible
	      $tid =~ s/ENST/ENSP/;
	      
	      $trans->id($tid);
	      $trans->version($transcript->version);
	      
	      $transcript->translation($trans);
	      
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

  # declared here as an array, but data is parsed in
  # method call to get features
  $self->{'_sf_array'} = [];
 
  # set stuff in self from @args
  return $make; # success - we hope!
}

=head2 get_all_SeqFeatures

 Title   : get_all_SeqFeatures
 Usage   : foreach my $sf ( $contig->get_all_SeqFeatures)
 Function: Gets all the sequence features on the whole contig.
 Example :
 Returns : 
 Args    :


=cut

sub get_all_SeqFeatures {
    my ($self) = @_;

    my @out;

    push(@out,$self->get_all_SimilarityFeatures);
    push(@out,$self->get_all_RepeatFeatures);

    return @out;
}

=head2 get_all_RepeatFeatures

 Title   : get_all_RepeatFeatures
 Usage   : foreach my $sf ( $contig->get_all_RepeatFeatures
 Function: Gets all the repeat features on the whole contig.
 Example :
 Returns : 
 Args    :


=cut

sub get_all_RepeatFeatures {
    my ($self) = @_;

    my $debug=0;

    my $sfobj=Bio::EnsEMBL::Analysis::FeatureParser->new($self->id,
							 $self->_clone_dir,
							 $self->disk_id,
							 $self->_gs,
							 $self->seq,
							 'repeat',
							 $debug);

    push(@{$self->{'_sf_array'}},$sfobj->each_Feature);
    
    # return array of objects
    return @{$self->{'_sf_array'}};
}

=head2 get_all_SimilarityFeatures

 Title   : get_all_SimilarityFeatures
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_all_SimilarityFeatures {
    my ($self) = @_;

    # get sf object
    my $debug=0;

    my $sfobj=Bio::EnsEMBL::Analysis::FeatureParser->new($self->id,
							 $self->_clone_dir,
							 $self->disk_id,
							 $self->_gs,
							 $self->seq,
							 'similarity',
							 $debug);

    push(@{$self->{'_sf_array'}},$sfobj->each_Feature);
    
    # return array of objects
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


=head2 length

 Title   : length
 Usage   : $obj->length($newval)
 Function: 
 Returns : value of length
 Args    : newvalue (optional)


=cut

sub length{
    my $self = shift;
    if( @_ ) {
	my $value = shift;
	$self->{'_length'} = $value;
    }
    return $self->{'_length'};
}


=head2 order

 Title   : order
 Usage   : $obj->order($newval)
 Function: 
 Returns : value of order
 Args    : newvalue (optional)


=cut

sub order{
    my $self = shift;
    if( @_ ) {
	my $value = shift;
	$self->{'_order'} = $value;
    }
    return $self->{'_order'};
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
    if( @_ ) {
	my $value = shift;
	$self->{'_offset'} = $value;
    }
    return $self->{'_offset'};
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
    if( @_ ) {
	my $value = shift;
	$self->{'_orientation'} = $value;
    }
    return $self->{'_orientation'};
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
    my $clonediskid=$disk_id;
    $clonediskid=~s/\.\d+$//;

    # read from sequence file
    my $file=$self->_clone_dir . "/$clonediskid.seq";

    local *IN;
    open(IN,$file) || die "cannot open $file";
    my $seqin = Bio::SeqIO->new( '-format' => 'Fasta', -fh => \*IN);
    my($seq,$seqid,$ffound);
    while($seq=$seqin->next_seq()){
	$seqid=$seq->id;
	#print STDERR "Read $seqid from $file\n";
	if($seqid eq $disk_id){
	    $ffound=1;
	    last;
	}
    }
    close(IN);
    if(!$ffound){
	$self->throw("Cannot find contig $id in $file");
    }
    
    $self->{'seq'}=$seq;
    $self->{'seq'}->type('Dna');
    $self->{'seq'}->id($self->id());
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

=head2 _gs

 Title   : _gs
 Usage   : $obj->_gs($newval)
 Function: 
 Returns : value of _gs
 Args    : newvalue (optional)


=cut

sub _gs{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'_gs'} = $value;
    }
    return $obj->{'_gs'};

}


=head2 created

 Title   : created
 Usage   :
 Function: created data - in TimDB this is based on the created data of the clone
 Example : $contig->created
 Returns : 
 Args    :


=cut

sub created{
    my ($self) = @_;
    my $id=$self->id;
    # get clone from contig
    $id=~s/\.\d+$//;
    return $self->_dbobj->get_Clone($id)->created;
}


=head2 seq_date

 Title   : seq_date
 Usage   : $contig->seq_date()
 Function: Gives the unix time value of the dna sequence file stored in TimDB
    Since DNA is stored by clone, there is just one time for all contigs
 Example : $contig->seq_date()
 Returns : unix time
 Args    : none


=cut

sub seq_date{
    my ($self) = @_;
    my $id=$self->id;
    # get clone from contig
    $id=~s/\.\d+$//;
    return $self->_dbobj->get_Clone($id)->seq_date;
}
1;
