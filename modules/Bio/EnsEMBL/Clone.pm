
#
# EnsEMBL module for Bio::EnsEMBL::Clone
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

@ISA = qw( Bio::EnsEMBL::Root );

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

  Args       : none
  Example    : none
  Description: gets all Genes that have coordinates on this Clone. They
               come in RawContig coords, but not all coords need to be on this 
               Clone 
  Returntype : list of Bio::EnsEMBL::Gene
  Exceptions : none
  Caller     : general

=cut


sub get_all_Genes
{
    my $self=shift;

    return $self->adaptor->get_all_Genes( $self->dbID );
}





=head2 get_Contig

  Args       : none
  Example    : none
  Description: deprecated, use ContigAdaptor to get Contig
  Returntype : none
  Exceptions : none
  Caller     : none

=cut


sub get_Contig {
   my ($self,$contigid) = @_;

   $self->throw("Clone::get_Contig is deprecated, " . 
                "use \$contig_adaptor->fetch_by_dbID(\$contig_id) instead\n");

   return undef;

#   my $contig = $self->adaptor->get_Contig($contigid);
   
#   return $contig->fetch();
}




=head2 get_all_Contigs

  Args       : none
  Example    : none
  Description: get RawContig objects from this Clone
  Returntype : list of Bio::EnsEMBL::RawContig
  Exceptions : none
  Caller     : general

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


=head2 delete

  Args       : none
  Example    : none
  Description: deprecated, use object adaptor for deletion
  Returntype : none
  Exceptions : none
  Caller     : none

=cut



sub delete
{
    my ($self)=shift;
    $self->warn("delete is now deprecated, use delete_by_dbID instead");
    $self->delete_by_dbID;
}




=head2 delete_by_dbID

  Args       : none
  Example    : none
  Description: Deletes the clone, contig, dna and features for this,
               Genes are not deleted
  Returntype : none
  Exceptions : none
  Caller     : general

=cut



sub delete_by_dbID {
    my ($self)=shift;
    return $self->adaptor->delete_by_dbID($self->dbID);
}




=head2 get_all_rawcontigs_by_position

  Arg   1    : int $base_pair
  Example    : none
  Description: returns the RawContig that contains that clone basepair
  Returntype : Bio:EnsEMBL::RawContig
  Exceptions : base_pair > 0, returns the last Contig if base_pair is outside
  Caller     : general

=cut



sub get_rawcontig_by_position {

    my ($self, $pos) = @_;


    if ($pos < 1 ){
        $self->throw("get_rawcontig_by_position error: Position must be > 0");
    }
    
    my @contigs =  $self->get_all_Contigs();
    @contigs = sort { $b->embl_offset <=> $a->embl_offset } @contigs;
    
    foreach my $c ( @contigs ) {
        if ($pos > $c->embl_offset) {
             return $c;
        } 
    }
    
    return (undef);
}



=head2 is_golden

  Args       : none
  Example    : none
  Description: deprecated, use assembly_mapper->in_assembly( $clone )
  Returntype : none
  Exceptions : none
  Caller     : none

=cut


sub is_golden{
   my ($self,@args) = @_;
   
   $self->warn("Clone::is_golden is deprecated. " .
	       "Use \$assembly_mapper->in_assembly(\$clone)");

   my $asma = $self->adaptor()->db()->get_AssemblyMapperAdaptor();
   my $am = $asma->fetch_by_type($self->db()->assembly_type());
   return $am->in_assembly($self);

#   foreach my $contig ($self->get_all_Contigs) {
#       if ($contig->is_golden) {
#	   return 1;
#       }
#   }
#   return 0;
}




=head2 seq_date

  Args       : none
  Example    : none
  Description: loops over all $contig->seq_date, throws a warning 
               if they are different and returns the first unix 
               time value of the dna created datetime field, which indicates
               the original time of the dna sequence data
  Returntype : a unix time ??
  Exceptions : none
  Caller     : general

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

  Arg [1]    : string $htg_phase
               0,1,2,3 representing how finished the clone is
  Example    : none
  Description: get/set for attribute htg_phase
               ( high throughput genome project phase ) 
  Returntype : string
  Exceptions : none
  Caller     : general

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

  Arg [1]    : string $created
  Example    : none
  Description: get/set for attribute created.
               Gives the unix time value of the created 
               datetime field, which indicates
               the first time this clone was put in ensembl
  Returntype : string
  Exceptions : none
  Caller     : general

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

  Arg [1]    : string $modified
  Example    : none
  Description: get/set for attribute modified
               Gives the unix time value of the modified 
               datetime field, which indicates
               the last time this clone was modified in ensembl
  Returntype : string
  Exceptions : none
  Caller     : general

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

  Arg [1]    : string $version
  Example    : none
  Description: get/set for attribute version
               this could contain an ensembl version for the clone.
               Usually we just use the EMBL one though. EnsEMBL version
               are currently not generated or maintained for clones.
  Returntype : string
  Exceptions : none
  Caller     : general

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

  Arg [1]    : string $embl_version
  Example    : none
  Description: get/set for attribute embl_version
  Returntype : string
  Exceptions : none
  Caller     : general

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

  Arg [1]    : string $embl_id
  Example    : none
  Description: get/set for attribute embl_id
  Returntype : string
  Exceptions : none
  Caller     : general

=cut




sub embl_id {
   my ($obj,$value) = @_;
   if( defined $value) {
      $obj->{'embl_id'} = $value;
    }
    return $obj->{'embl_id'};

}



=head2 id

  Args       : none
  Example    : none
  Description: should be deprecated, gives an optional ensembl name
               for the clone. Was used for non submitted clones and is
               probably useless now.
  Returntype : none
  Exceptions : none
  Caller     : none

=cut


sub id {
   my ($obj,$value) = @_;
   if( defined $value) {
      $obj->{'_clone_id'} = $value;
    }
    return $obj->{'_clone_id'};

}

=head2 dbID

  Arg [1]    : int $dbID
  Example    : none
  Description: get/set for the database internal id
  Returntype : int
  Exceptions : none
  Caller     : general, set from adaptor on store

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

  Arg [1]    : Bio::EnsEMBL::DBSQL::CloneAdaptor $adaptor
  Example    : none
  Description: get/set for this objects Adaptor
  Returntype : Bio::EnsEMBL::DBSQL::CloneAdaptor
  Exceptions : none
  Caller     : general, set from adaptor on store

=cut


sub adaptor {
   my ($obj,$value) = @_;
   if( defined $value) {
      $obj->{'adaptor'} = $value;
    }
    return $obj->{'adaptor'};

}


1;






