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
use Bio::EnsEMBL::DB::RawContigI;
use Bio::Seq;
use Bio::SeqIO;
use Bio::EnsEMBL::Analysis::Genscan;
use Bio::EnsEMBL::Analysis::FeatureParser;
use Bio::EnsEMBL::Chromosome;
use FileHandle;

# Object preamble - inheriets from Bio::Root::Object
use Bio::Root::Object;

@ISA = qw(Bio::Root::Object Bio::EnsEMBL::DB::RawContigI);
# new() is inherited from Bio::Root::Object

# _initialize is where the heavy stuff will happen when new is called

sub _initialize {
    my($self,@args) = @_;
  
    my $make = $self->SUPER::_initialize;
    my ($dbobj,$id,$disk_id,$clone_dir,$orientation,$length,
	$chr,$species,$embl_offset,$embl_order,$international_id)=
	$self->_rearrange([qw(DBOBJ
			      ID
			      DISK_ID
			      CLONE_DIR
			      ORIENTATION
			      LENGTH
			      CHR
			      SPECIES
			      EMBL_OFFSET
			      EMBL_ORDER
			      INTERNATIONAL_ID
			      )],@args);
    
    $id          || $self->throw("Cannot make contig object without id");
    $dbobj       || $self->throw("Cannot make contig object without db object");
    $dbobj->isa('Bio::EnsEMBL::TimDB::Obj') || $self->throw("Cannot make contig object with a $dbobj object");
    
    $self->id         ($id);
    $self->disk_id    ($disk_id);
    $self->_dbobj     ($dbobj);
    $self->_clone_dir ($clone_dir);
    $self->orientation($orientation);
    $self->length     ($length);
    $self->chromosome ($chr,$species);
    $self->embl_offset($embl_offset);
    $self->embl_order ($embl_order);
    $self->international_id ($international_id);
    
    # declared here as an array, but data is parsed in
    # method call to get features
    
    $self->{'_sf_array'}        = [];   # Sequence features
    $self->{'_sf_repeat_array'} = [];   # Repeat features
    
    # set stuff in self from @args
    return $make; # success - we hope!
}

sub validate {

    my ($self) = @_;

    $self->id          || $self->throw("Cannot make contig object without id");
    $self->disk_id     || $self->throw("Cannot make contig object without disk_id");
    $self->_dbobj       || $self->throw("Cannot make contig object without db object");
#    $self->order       || $self->throw("Cannot make contig object without order");
#    $self->offset      || $self->throw("Cannot make contig object without offset");
#    $self->orientation || $self->throw("Cannot make contig object without orientation");
    $self->length      || $self->throw("Cannot make contig object without length");
    $self->chromosome  || $self->throw("Cannot make contig object without chromosome");
    $self->_clone_dir   || $self->throw("Cannot make contig object without clone_dir");
    $self->embl_offset || $self->throw("Cannot make contig object without embl_offset");
    $self->embl_order  || $self->throw("Cannot make contig object without embl_order");
#    $self->international_id  || $self->throw("Cannot make contig object without international_id");

    $self->_dbobj->isa('Bio::EnsEMBL::TimDB::Obj') ||   $self->throw("Cannot make contig object with a [" . $self->_dbobj ."] object");
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
    push(@out,$self->get_all_GenePredictions);
    push( @out, $self->get_all_StsFeatures );

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
    
    if (!$self->{_read_Repeats}) {
	$self->featureParser->read_Repeats;
	$self->{_read_Repeats} = 1;
    } 

    # return array of objects
    return $self->featureParser->each_Repeat;
}


sub featureParser {
    my ($self) = @_;

    my $debug = 0;

    if (!defined($self->{_featureParser})) {
	print STDERR "Looking at [",$self->primary_seq,"]\n";
	my $sfobj=Bio::EnsEMBL::Analysis::FeatureParser->new($self->id,
							     $self->_clone_dir,
							     $self->disk_id,
							     $self->_gs,
							     $self->primary_seq,
							     $debug);

	$self->{_featureParser} = $sfobj;
    }

    return $self->{_featureParser};
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

    if (!defined($self->{_read_Similarities})) {
	$self->featureParser->read_Similarities;
	$self->{_read_Similarities} = 1;
    }
    # return array of objects
    return $self->featureParser->each_Feature;
}

=head2 get_all_StsFeatures

 Title   : get_all_StsFeatures
 Usage   :
 Function:
 Example :
 Returns : List of FeaturePairs, describing sts-hits on the contig.
 Args    :


=cut

sub get_all_StsFeatures {
    my ($self) = @_;

    if (!defined($self->{_read_StsFeatures})) {
	$self->featureParser->read_StsFeatures;
	$self->{_read_StsFeatures} = 1;
    }
    # return array of objects
    return $self->featureParser->each_StsFeature;
}

=head2 get_all_GenePredictions

 Title   : get_all_GenePrediction
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_all_GenePredictions {
    my ($self) = @_;

    if (!defined($self->{_read_Genscan})) {
	$self->featureParser->read_Genscan;
	$self->{_read_Genscan} = 1;
    }
    # return array of objects
    my @ret= $self->featureParser->each_Genscan;
    # make sure they have the correct seqname
    my @out;
    foreach my $f ( @ret ) {
	$f->seqname($self->id());
	push(@out,$f);
    }
    
    return @out;
}


=head2 get_all_Genes

 Title   : get_all_Genes
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_all_Genes {
    my ($self) = @_;

    # map genes if not already mapped
    $self->_build_genes unless $self->{'_mapped'};
    
    return @{$self->{'_gene_array'}};
}


=head2 get_old_exons

 Title   : get_old_exons
 Usage   : 
 Function: 
 Example :
 Returns : 
 Args    :


=cut

sub get_old_exons {
    my ($self) = @_;
    
    # all old_exon stuff read at clone level on demand
    my $contig_id=$self->id;
    my $clone_id=$contig_id;
    $clone_id=~s/\.\d+$//;
    return $self->_dbobj->get_Clone($clone_id)->get_old_exons($contig_id);
}


=head2 add_SeqFeature

 Title   : add_SeqFeature
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :

=cut

sub add_SeqFeature {
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

sub length {
    my $self = shift;
    if( @_ ) {
	my $value = shift;
	$self->{'_length'} = $value;
    }
    return $self->{'_length'};
}


=head2 embl_order

 Title   : embl_order
 Usage   : $obj->embl_order($newval)
 Function: 
 Returns : value of order
 Args    : newvalue (optional)


=cut

sub embl_order{
    my ($self,$arg) = @_;

    if( defined($arg) ) {
	$self->{'_order'} = $arg;
    }

    return $self->{'_order'};
}

sub order {
    my ($self,$arg) = @_;

    $self->warn("Contig->order is deprecated in Bio::EnsEMBL::DB::ContigI. Use Contig->embl_order instead");

    return $self->embl_order($arg);
}

=head2 embl_offset

 Title   : embl_offset
 Usage   : $self->embl_offset($newval)
 Function: 
 Returns : value of offset
 Args    : newvalue (optional)

=cut

sub embl_offset{
    my $self = shift;
    if( @_ ) {
	my $value = shift;
	$self->{'_offset'} = $value;
    }
    return $self->{'_offset'};
}

sub offset {
    my $self = shift;

    $self->warn("Contig->offset is deprecated in Bio::EnsEMBL::DB::ContigI. Use Contig->embl_offset instead\n");

    return $self->embl_offset(@_);
}


=head2 international_id

 Title   : international_id
 Usage   : $obj->international_id($newval)
 Function: 
 Returns : value of order
 Args    : newvalue (optional)


=cut

sub international_id{
    my ($self,$arg) = @_;

    if( defined($arg) ) {
	$self->{'_international_id'} = $arg;
    }

    return $self->{'_international_id'};
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

sub seq {
    my ($self,$arg) = @_;

    $self->warn("Contig::seq is deprecated. Use primary_seq instead");

    return $self->primary_seq($arg);
}

=head2 primary_seq

 Title   : primary_seq
 Usage   : $self->primary_seq
 Function: 
 Returns : value of seq
 Args    : none

=cut

sub primary_seq {
    my ($self,$arg) = @_;
    
    if(defined($arg)) {
	$self->throw("Cannot set a sequence in TimDB [$arg]");
    } elsif(defined $self->{'seq'}){
	return $self->{'seq'};
    }

    my $id          = $self->id;
    my $disk_id     = $self->disk_id;
    my $clonediskid = $disk_id;
    $clonediskid    =~s/\.\d+$//;

    # read from sequence file
    my $file=$self->_clone_dir . "/$clonediskid.seq";

    local *IN;
    open(IN,$file) || die "cannot open $file";
    my $seqin = Bio::SeqIO->new( '-format' => 'Fasta', -fh => \*IN);
    my($seq,$seqid,$ffound);
    while($seq=$seqin->next_primary_seq()){
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
    $self->{'seq'}->moltype('dna');
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

=head2 chromosome

 Title   : chromosome
 Usage   : $chr = $self->chromosome([$chromosome[)
 Function: get/set chromosome for this contig. Defaults to the unknown
           human Chromosome.
 Returns : a Chromosome object
 Args    : 

=cut

sub chromosome {
    my ($self,$chr,$species) = @_;

    if (defined($chr) && defined($species)) {
#    if (!defined($self->{_chromosome})) {
	$self->{_chromosome} = Bio::EnsEMBL::Species->chromosome_by_name($species,$chr);
    } 

    return $self->{_chromosome};
}


=head2 internal_id

 Title   : internal_id
 Usage   : $self->internal_id($newval)
 Function: 
 Returns : value of internal_id
 Args    : newvalue (optional)

=cut

sub internal_id {
    my $self = shift;
    if ( @_ ) {
	my $value = shift;
	$self->{'internal_id'} = $value;
    }
    return $self->{'internal_id'};
}


=head2 disk_id

 Title   : disk_id
 Usage   : $self->disk_id($newval)
 Function: 
 Returns : value of disk_id
 Args    : newvalue (optional)

=cut

sub disk_id {
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

    return $obj->dbobj($value);

}

sub dbobj {
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
 Usage   : $self->_gs
 Function: 
 Returns : value of _gs
 Args    : newvalue (optional)


=cut

sub _gs{
    my ($self) = @_;

    if(!defined($self->{'_gs'})) {
	print STDERR "Passing in ",$self->primary_seq," to genscan\n";
	my $gs = Bio::EnsEMBL::Analysis::Genscan->new($self->_clone_dir . "/" . 
						      $self->disk_id . ".gs",
						      $self->primary_seq());
	
	
	$self->{'_gs'} = $gs;
    }

    return $self->{'_gs'};
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

=head2 _build_genes

 Title   : _build_genes
 Usage   : never called by user
 Function: deferred loading routine - builds etg map (if not yet called) then builds
    gene map for this contig
 Example : 
 Returns : none
 Args    : none


=cut

sub _build_genes{
    my $self=shift;

    # map if not already mapped
    $self->_dbobj->map_etg unless $self->_dbobj->{'_mapped'};
    
    # ok. Hell. We open the Genscan file using the Genscan object.
    # this is needed to remap the exons lower down
    # FIXME - sure we don't need this call
    $self->_gs;
  
    # we yank out each exon and build a hash on start position
    my %exhash = $self->_make_exon_hash;

    my $id=$self->id;
    my $dbobj=$self->_dbobj;

    # build array of genes
    $self->{'_gene_array'} = [];
    {
	my $bioseq=$self->primary_seq;
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
		$trans->end_exon_id  ($le->id);
		
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
    $self->{'_mapped'}=1;
}

sub _make_exon_hash {
    my ($self) = @_;

    my $gs = $self->_gs;

    my %exhash;

    foreach my $t ( $gs->each_Transcript ) {
	foreach my $ex ( $t->each_Exon ) {
	    $exhash{$ex->start} = $ex;
	}
    }
    return (%exhash);
}

sub checksum {
    my ($self,$arg) = @_;

    if (defined($arg)) {
	$self->{_checksum} = $arg;
    }

    return $self->{_checksum};
}

=head2 get_left_overlap

 Title   : get_left_overlap
 Usage   : $overlap_object = $contig->get_left_overlap();
 Function: Returns the overlap object of contig to the left.
           This could be undef, indicating no overlap
 Returns : A Bio::EnsEMBL::ContigOverlapHelper object
 Args    : None

=cut

sub get_left_overlap{
   my ($self,@args) = @_;

   return;
}


=head2 get_right_overlap

 Title   : get_right_overlap
 Usage   : $overlap_object = $contig->get_right_overlap();
 Function: Returns the overlap object of contig to the left.
           This could be undef, indicating no overlap
 Returns : A Bio::EnsEMBL::ContigOverlapHelper object
 Args    : None

=cut

sub get_right_overlap{
   my ($self,@args) = @_;

   return;
}

=head2 get_all_ExternalFeatures

 Title   : get_all_ExternalFeatures (Abstract)
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_all_ExternalFeatures{
   my ($self) = @_;
   
   return;
}



1;
