
#
# BioPerl module for DB::Clone
#
# Cared for by EnsEMBL (www.ensembl.org)
#
# Copyright GRL and EBI
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::Clone - Object representing one clone

=head1 SYNOPSIS

    # $db is Bio::EnsEMBL::DB::Obj 

    @contig = $db->get_all_Contigs();

    $clone = $db->get_Clone();

    @genes    = $clone->get_all_Genes();

=head1 DESCRIPTION

Represents information on one Clone

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::Clone;
use vars qw(@ISA);
use strict;

# Object preamble - inheriets from Bio::EnsEMBL::Root

use Bio::EnsEMBL::Root;

@ISA = qw(Bio::EnsEMBL::Root);

sub new {
    my ($class,$adaptor,@args) = @_;

    my $self = {};
    bless $self,$class;


    my ($internal_id,$id,$embl_id,$version,$embl_version,$htg_phase,$created,$modified)=@args;

    $self->throw("Don't have a id [$id] for new clone") unless $id;
    $self->throw("Don't have a adaptor [$adaptor] for new clone $id") unless $adaptor;
    $self->throw("Don't have a internal id [$internal_id] for new clone $id") unless $internal_id;
    $self->throw("Don't have a embl id [$embl_id] for new clone $id") unless $embl_id;
    
    $self->throw("Don't have a version [$version] for new clone $id") unless defined $version;
    $self->throw("Don't have a embl verson [$embl_version] for new clone $id") unless defined $embl_version;
    

    if( $version == 0 ) {
	$self->warn("seq version $version and embl version $embl_version are 0 - this will not play nicely with external feature factories!");
    }

    # HACK htg_phase is NOT set correctly for many clones (th, 7/01)
   $self->warn("Don't have a htg phase [$htg_phase] for new clone $id") unless $htg_phase;  
    #$self->throw("Don't have a created [$created] for new clone") unless $created;
    #$self->throw("Don't have a modified [$modified] for new clone") unless $modified;

    $self->adaptor($adaptor);
    $self->dbID($internal_id);
    $self->id($id);
    $self->embl_id($embl_id);
    $self->version($version);
    $self->embl_version($embl_version);
    $self->htg_phase($htg_phase);
    $self->created($created);
    $self->modified($modified);

    return $self;
}





=head2 get_all_Genes

 Title   : get_all_Genes
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut



sub get_all_Genes
{
    my ($self,$supporting)=@_;

    return $self->adaptor->get_all_Genes($self->dbID,$supporting);
}




=head2 get_Contig

 Title   : get_Contig
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_Contig {
   my ($self,$contigid) = @_;

   my $contig = $self->adaptor->get_Contig($contigid);
   
   return $contig->fetch();
}

=head2 get_all_geneid

 Title   : get_all_geneid
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut



sub get_all_my_geneid
{
    my ($self)=shift;

    return $self->adaptor->get_all_my_geneid($self->dbID);
}






=head2 get_all_Contigs

 Title   : get_Contigs
 Usage   : foreach $contig ( $clone->get_all_Contigs ) 
 Function:
 Example :
 Returns : 
 Args    :


=cut




sub get_all_Contigs {
    my( $self ) = @_;

    my( $c_list );
    unless ($c_list = $self->{'_contig_list'}) {
        my $ra = $self->adaptor->db->get_RawContigAdaptor;
        $c_list = $ra->fetch_by_clone($self);
        $self->{'_contig_list'} = $c_list;
    }
    return @$c_list;
}



sub delete
{
    my ($self)=shift;
    $self->warn("delete is now deprecated, use delete_by_dbID instead");
    $self->delete_by_dbID;
}






=head2 delete_by_dbID

 Title   : delete_by_dbID
 Usage   : $clone->delete_by_dbID()
 Function: Deletes clone (itself), including contigs and features, but not its genes
 Example : 
 Returns : nothing
 Args    : none


=cut

sub delete_by_dbID {
    my ($self)=shift;
    return $self->adaptor->delete_by_dbID($self->dbID);
}






=head2 get_all_rawcontigs_by_position

 Title   : get_rawcontig_by_position
 Usage   : $obj->get_rawcontig_by_position($position)
 Function: 
 Example : 
 Returns : returns a raw contig object or undef on error
 Args    : a position (basepair) in clone


=cut

sub get_rawcontig_by_position {

    my ($self, $pos) = @_;

    if( !ref $self || ! $self->isa('Bio::EnsEMBL::DB::CloneI') ) {
        $self->throw("Must supply a clone to get_all_RawContigs: Bailing out...");
    }

    if ($pos < 1 ){
        $self->throw("get_rawcontig_by_position error: Position must be > 0");
    }
    
    my @contigs =  $self->get_all_Contigs();
    @contigs = sort { $a->embl_offset <=> $b->embl_offset } @contigs;
    
    foreach my $c (reverse @contigs ) {
        if ($pos > $c->embl_offset) {
            my $size = $c->embl_offset + $c->length;
            return $c;
        } else {
            my $size = $c->embl_offset + $c->length;
            next;
        }
    }
    
    return (undef);
}



=head2 is_golden

 Title   : is_golden
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub is_golden{
   my ($self,@args) = @_;
   
   foreach my $contig ($self->get_all_Contigs) {
       if ($contig->is_golden) {
	   return 1;
       }
   }
   return 0;
}





=head2 seq_date

 Title   : seq_date
 Usage   : $clone->seq_date()
 Function: loops over all $contig->seq_date, throws a warning if they are different and 
           returns the first unix time value of the dna created datetime field, which indicates
           the original time of the dna sequence data
 Example : $clone->seq_date()
 Returns : unix time
 Args    : none


=cut

sub seq_date {
   my ($self) = @_;

   my $id = $self->id();
   my ($seq_date,$old_seq_date);

   foreach my $contig ($self->get_all_Contigs) {
       $seq_date = $contig->seq_date;
       if ($old_seq_date) {
	   if ($seq_date != $old_seq_date) {
	       $self->warn ("The created date of the DNA sequence from contig 
                             $contig is different from that of the sequence 
                             from other contigs on the same clone!");
	   }
       }
       $old_seq_date = $seq_date;
   }
   
   return $seq_date;
}


=head2 htg_phase

 Title   : htg_phase
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub htg_phase {
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'htg_phase'} = $value;
    }
    return $obj->{'htg_phase'};
}

=head2 created

 Title   : created
 Usage   : $clone->created()
 Function: Gives the unix time value of the created datetime field, which indicates
           the first time this clone was put in ensembl
 Example : $clone->created()
 Returns : unix time
 Args    : none


=cut

sub created {
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'created'} = $value;
    }
    return $obj->{'created'};
}

=head2 modified

 Title   : modified
 Usage   : $clone->modified()
 Function: Gives the unix time value of the modified datetime field, which indicates
           the last time this clone was modified in ensembl
 Example : $clone->modified()
 Returns : unix time
 Args    : none


=cut


sub modified {
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'modified'} = $value;
    }
    return $obj->{'modified'};
}

=head2 version

 Title   : version
 Usage   : $clone->version()
 Function: Gives the value of version
           (Please note: replaces old sv method!!!)
 Example : $clone->version()
 Returns : version number
 Args    : none


=cut



sub version{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'version'} = $value;
    }
    return $obj->{'version'};

}

=head2 embl_version

 Title   : embl_version
 Usage   : $clone->embl_version()
 Function: Gives the value of the EMBL version, i.e. the data version
 Example : $clone->embl_version()
 Returns : version number
 Args    : none


=cut

sub embl_version {
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'embl_version'} = $value;
    }
    return $obj->{'embl_version'};
}



=head2 embl_id

 Title   : embl_id
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut




sub embl_id {
   my ($obj,$value) = @_;
   if( defined $value) {
      $obj->{'embl_id'} = $value;
    }
    return $obj->{'embl_id'};

}



=head2 id

 Title   : id
 Usage   : $obj->id($newval)
 Function: 
 Example : 
 Returns : value of id
 Args    : newvalue (optional)


=cut

sub id {
   my ($obj,$value) = @_;
   if( defined $value) {
      $obj->{'_clone_id'} = $value;
    }
    return $obj->{'_clone_id'};

}

=head2 dbID

 Title   : dbID
 Usage   : $obj->dbID($newval)
 Function: 
 Returns : value of dbID
 Args    : newvalue (optional)


=cut

sub dbID{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'dbID'} = $value;
    }
    return $obj->{'dbID'};

}


=head2 adaptor

 Title   : adaptor
 Usage   : $obj->adaptor($newval)
 Function: 
 Example : 
 Returns : value of adaptor
 Args    : newvalue (optional)


=cut

sub adaptor {
   my ($obj,$value) = @_;
   if( defined $value) {
      $obj->{'adaptor'} = $value;
    }
    return $obj->{'adaptor'};

}


1;






