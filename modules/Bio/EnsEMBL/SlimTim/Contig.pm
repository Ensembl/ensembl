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

    $contig = Bio::EnsEMBL::SlimTim::Contig->new();
    
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
package Bio::EnsEMBL::SlimTim::Contig;
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
    my ($id,$diskid,$seq,$dir)=
	$self->_rearrange([qw(ID
			      DISK_ID
			      SEQ
			      DIR
			      )],@args);
    
    $id          || $self->throw("Cannot make contig object without id");
    $diskid      || $self->throw("Cannot make contig object without diskid");
    $seq         || $self->throw("Cannot make contig object without seq");
    $dir         || $self->throw("Cannot make contig object without dir");

    $self->id($id);
    $self->disk_id($diskid);
    $self->primary_seq($seq);
    $self->_dir($dir);
    
    $self->_parse_sequence_header();

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
	my $sfobj=Bio::EnsEMBL::Analysis::FeatureParser->new($self->id,
							     $self->_dir,
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
 Usage   : $obj->primary_seq($newval)
 Function: 
 Returns : value of primary_seq
 Args    : newvalue (optional)


=cut

sub primary_seq{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'primary_seq'} = $value;
    }
    return $obj->{'primary_seq'};

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
	my $gs = Bio::EnsEMBL::Analysis::Genscan->new($self->_dir . "/" . 
						      $self->disk_id . ".gs",
						      $self->primary_seq());
	
	
	$self->{'_gs'} = $gs;
    }

    return $self->{'_gs'};
}


=head2 created

 Title   : created
 Usage   : $obj->created($newval)
 Function: 
 Returns : value of created
 Args    : newvalue (optional)


=cut

sub created{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'created'} = $value;
    }
    return $obj->{'created'};

}

=head2 seq_date

 Title   : seq_date
 Usage   : $obj->seq_date($newval)
 Function: 
 Returns : value of seq_date
 Args    : newvalue (optional)


=cut

sub seq_date{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'seq_date'} = $value;
    }
    return $obj->{'seq_date'};

}

=head2 _parse_sequence_header

 Title   : _parse_sequence_header
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub _parse_sequence_header{
   my ($self) = @_;

   my $desc = $self->primary_seq->desc();
   my %hash;

   #Unfinished sequence: AC010481  Contig_ID: 00001 acc=AC010481 sv=AC010481.4 start=1 emblid=AC010481 div=HTGS_PHASE1 chr=5 freeze=1 int_id=AC010481.4~1 species=human Length: 6513 bp
   my @keyvalue = split(/ /,$desc);
   foreach my $keyvalue ( @keyvalue ) {
       if( $keyvalue =~ /(\S+)=(\S+)/ ) {
	   $hash{$1} = $2;
       }
   }

   $self->clone_accession($hash{'sv'});
   $self->embl_offset($hash{'start'});
   $self->clone_phase($hash{'div'});
   $self->international_id($hash{'int_id'});
   $hash{'int_id'} =~ /\S+\.\d+~(\d+)/;
   my $order = $1;
   $self->embl_order($order);
   my $id = $self->id();
   $id =~ /\S+\.(\S+)/;
   my $ext = $1;
   $self->id($hash{'acc'}.".".$ext);

}


=head2 clone_accession

 Title   : clone_accession
 Usage   : $obj->clone_accession($newval)
 Function: 
 Returns : value of clone_accession
 Args    : newvalue (optional)


=cut

sub clone_accession{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'clone_accession'} = $value;
    }
    return $obj->{'clone_accession'};

}


=head2 embl_offset

 Title   : embl_offset
 Usage   : $obj->embl_offset($newval)
 Function: 
 Returns : value of embl_offset
 Args    : newvalue (optional)


=cut

sub embl_offset{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'embl_offset'} = $value;
    }
    return $obj->{'embl_offset'};

}

=head2 embl_order

 Title   : embl_order
 Usage   : $obj->embl_order($newval)
 Function: 
 Returns : value of embl_order
 Args    : newvalue (optional)


=cut

sub embl_order{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'embl_order'} = $value;
    }
    return $obj->{'embl_order'};

}

=head2 _dir

 Title   : _dir
 Usage   : $obj->_dir($newval)
 Function: 
 Returns : value of _dir
 Args    : newvalue (optional)


=cut

sub _dir{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'_dir'} = $value;
    }
    return $obj->{'_dir'};

}

=head2 clone_phase

 Title   : clone_phase
 Usage   : $obj->clone_phase($newval)
 Function: 
 Returns : value of clone_phase
 Args    : newvalue (optional)


=cut

sub clone_phase{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'clone_phase'} = $value;
    }
    return $obj->{'clone_phase'};

}



1;





