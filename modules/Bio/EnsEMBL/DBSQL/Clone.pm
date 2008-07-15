
#
# BioPerl module for DB::Clone
#
# Cared for by Ewan Birney <birney@sanger.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::DB::Clone - Object representing one clone

=head1 SYNOPSIS

    # $db is Bio::EnsEMBL::DB::Obj 

    @contig = $db->get_Contigs();

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


package Bio::EnsEMBL::DBSQL::Clone;
use vars qw(@ISA);
use strict;

# Object preamble - inheriets from Bio::Root::Object

use Bio::Root::Object;
use Bio::EnsEMBL::DBSQL::RawContig;
use Bio::EnsEMBL::DBSQL::Feature_Obj;
use Bio::EnsEMBL::DBSQL::Gene_Obj;
use Bio::EnsEMBL::DB::CloneI;

@ISA = qw(Bio::Root::Object Bio::EnsEMBL::DB::CloneI);

# new() is inherited from Bio::Root::Object

# _initialize is where the heavy stuff will happen when new is called

sub _initialize {
  my($self,@args) = @_;

  my $make = $self->SUPER::_initialize;

  my ($dbobj,$id) = $self->_rearrange([qw(DBOBJ
					  ID
					  )],@args);

  $id    || $self->throw("Cannot make clone db object without id");
  $dbobj || $self->throw("Cannot make clone db object without db object");
  $dbobj->isa('Bio::EnsEMBL::DBSQL::Obj') || $self->throw("Cannot make clone db object with a $dbobj object");

  $self->id($id);
  $self->_db_obj($dbobj);
  $self->fetch();

# set stuff in self from @args
  return $make; # success - we hope!
}

=head2 fetch

 Title   : fetch
 Usage   :
 Function:
 Example :
 Returns : nothing
 Args    :


=cut

sub fetch { 
    my ($self) = @_;
 
    my $id=$self->id();   
    my $sth = $self->_db_obj->prepare(
               "SELECT internal_id, id FROM clone WHERE id = \"$id\";");
    my $res = $sth ->execute();   
    my $rowhash = $sth->fetchrow_hashref;
    if( ! $rowhash ) {
	# make sure we deallocate sth - keeps DBI happy!
	$sth = 0;
	$self->throw("Clone $id does not seem to occur in the database!");
    }   
    $self->_internal_id($rowhash->{'internal_id'});
    return $self;
}


=head2 delete

 Title   : delete
 Usage   : $clone->delete()
 Function: Deletes clone (itself), including contigs and features, but not its genes
 Example : 
 Returns : nothing
 Args    : none


=cut

sub delete {
   my ($self) = @_;
   
 
   #(ref($clone_id)) && $self->throw ("Passing an object reference instead of a variable\n");

   my $internal_id = $self->_internal_id;

   my @contigs;
   my @dnas;

   # get a list of contigs to zap
  my $sth = $self->_db_obj->prepare(
      "SELECT internal_id, dna FROM contig WHERE clone = $internal_id");
   my $res = $sth->execute;

   while( my $rowhash = $sth->fetchrow_hashref) {
       push(@contigs,$rowhash->{'internal_id'});
       push(@dnas,$rowhash->{'dna'});
   }
   
   # Delete from DNA table, Contig table, Clone table
   
   foreach my $contig ( @contigs ) {
      my $sth = $self->_db_obj->prepare(
                      "DELETE FROM contig WHERE internal_id = $contig");
       my $res = $sth->execute;
   }


   foreach my $dna (@dnas) {
      $sth = $self->_db_obj->prepare("DELETE FROM dna WHERE id = $dna");
       $res = $sth->execute;

       # Mysql does not optimise or statements in where clauses
      $sth = $self->_db_obj->prepare(
                    "DELETE FROM contigoverlap WHERE dna_a_id = $dna;");
       $res = $sth ->execute;
      $sth = $self->_db_obj->prepare(
                    "DELETE FROM contigoverlap WHERE dna_b_id = $dna;");
       $res = $sth ->execute;

   }

  $sth = $self->_db_obj->prepare(
                  "DELETE FROM clone WHERE internal_id = $internal_id");
   $res = $sth->execute;
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
   my ($self, $supporting) = @_;
   my @out;
   my $clone_id = $self->_internal_id();   
   my %got;
    # prepare the SQL statement
   my $sth = $self->_db_obj->prepare("
        SELECT t.gene
        FROM transcript t,
             exon_transcript et,
             exon e,
             contig c
        WHERE e.contig = c.internal_id
          AND et.exon = e.id
          AND t.id = et.transcript
          AND c.clone = $clone_id
        ");

    my $res = $sth->execute();
   
    while (my $rowhash = $sth->fetchrow_hashref) { 
            
        if( ! exists $got{$rowhash->{'gene'}}) {  
            
           my $gene_obj = Bio::EnsEMBL::DBSQL::Gene_Obj->new($self->_db_obj);             
	   my $gene = $gene_obj->get($rowhash->{'gene'}, $supporting);
           if ($gene) {
	        push(@out, $gene);
           }
	   $got{$rowhash->{'gene'}} = 1;
        }       
    }
   
    if (@out) {
        return @out;
    }
    return;
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

   # should check this contig is in this clone?
   my $contig = $self->_db_obj->get_Contig($contigid);
   
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

sub get_all_my_geneid {
   my ($self) = @_;

   my $cloneid = $self->_internal_id;

  my $sth = $self->_db_obj->prepare(
                     "SELECT COUNT(*), cont.clone, ex.contig, tran.gene  "
                       . "FROM   contig          AS cont, "
                       . "       transcript      AS tran, "
                       . "       exon_transcript AS et, "
                       . "       exon            AS ex "
                       . "WHERE  ex.id            = et.exon "
                       . "AND    tran.id          = et.transcript "
                       . "AND    cont.clone       = $cloneid  "
                       . "AND    cont.internal_id = ex.contig "
                       . "GROUP BY tran.gene" );

   my @out;

   $sth->execute;
   while( my $rowhash = $sth->fetchrow_hashref) {
       push(@out,$rowhash->{'gene'});
   }
   return @out;
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
   my ($self) = @_;
   my $sth;
   my @res;
   my $internal_id = $self->_internal_id();

  my $sql =
    "SELECT id, internal_id FROM contig WHERE clone = $internal_id";
  # warn $sql;
   $sth= $self->_db_obj->prepare($sql);
   my $res  = $sth->execute();
   my $seen = 0;

   my $count   = 0;
   my $total   = 0;
   my $version = $self->embl_version();

   while( my $rowhash = $sth->fetchrow_hashref) {
       my $contig = $self->_db_obj->get_Contig( $rowhash->{'id'});
       $contig->internal_id($rowhash->{internal_id});
       $contig->seq_version($version);

       push(@res,$contig);
       $seen = 1;
   }

   if( $seen == 0  ) {
       $self->throw("Clone [$internal_id] has no contigs in the database. Should be impossible, but clearly isn't...");
   }

   return @res;   
}

=head2 get_all_RawContigs

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


=head2 get_all_ContigOverlaps 

 Title   : get_all_ContigOverlaps
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_all_ContigOverlaps {
    my ($self) = @_;
    
    
    my( %overlap );
    foreach my $contig ($self->get_all_Contigs) {
	foreach my $lap ($contig->get_all_Overlaps) {
            $overlap{$lap->hash_string} = $lap;
        }
    }
    return values %overlap;
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


=head2 htg_phase

 Title   : htg_phase
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub htg_phase {
   my $self = shift;
   if( defined $self->{'_htg_phase'} ) {
       return $self->{'_htg_phase'};
   }

   my $internal_id = $self->_internal_id();

  my $sth = $self->_db_obj->prepare(
        "SELECT htg_phase FROM clone WHERE internal_id = $internal_id");
   $sth->execute();
   my $rowhash = $sth->fetchrow_hashref();
   $self->{'_htg_phase'} = $rowhash->{'htg_phase'};
   return $self->{'_htg_phase'};
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
   my ($self) = @_;

   my $internal_id = $self->_internal_id();

  my $sth = $self->_db_obj->prepare( "SELECT UNIX_TIMESTAMP(created) "
                      . "FROM clone WHERE internal_id = $internal_id" );
   $sth->execute();
   my $rowhash = $sth->fetchrow_hashref();
   return $rowhash->{'UNIX_TIMESTAMP(created)'};
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
   my ($self) = @_;

   my $internal_id = $self->_internal_id();

  my $sth = $self->_db_obj->prepare( "SELECT UNIX_TIMESTAMP(modified) "
                      . "FROM clone WHERE internal_id = $internal_id" );
   $sth->execute();
   my $rowhash = $sth->fetchrow_hashref();
   return $rowhash->{'UNIX_TIMESTAMP(modified)'};
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

sub version {
   my $self = shift;

   my $internal_id = $self->_internal_id();

  my $sth = $self->_db_obj->prepare(
          "SELECT version FROM clone WHERE internal_id = $internal_id");
   $sth->execute();
   my $rowhash = $sth->fetchrow_hashref();
   return $rowhash->{'version'};
}

=head2 _stored

 Title   : _stored
 Usage   : $obj->_stored($newval)
 Function: Internal method should not really be needed
           stores the time of storage of the deleted object
 Returns : value of stored
 Args    : newvalue (optional)


=cut

sub _stored {
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'_stored'} = $value;
    }
    return $obj->{'_stored'};
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
   my $self = shift;
   my $internal_id = $self->_internal_id();

  my $sth = $self->_db_obj->prepare(
     "SELECT embl_version FROM clone WHERE internal_id = $internal_id");
   $sth->execute();
   my $rowhash = $sth->fetchrow_hashref();
   return $rowhash->{'embl_version'};
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
	       $self->warn ("The created date of the DNA sequence from contig $contig is different from that of the sequence from other contigs on the same clone!");
	   }
       }
       $old_seq_date = $seq_date;
   }
   
   return $seq_date;
}

=head2 sv

 Title   : sv
 Usage   : $clone->sv
 Function: old version method
 Example : $clone->sv
 Returns : version number
 Args    : none


=cut

sub sv{
   my ($self) = @_;
    
   $self->warn("Clone::sv - deprecated method. From now on you should use the Clone::version method, which is consistent with our nomenclature");

   return $self->embl_version();
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
   my ($self) = @_;

   my $internal_id = $self->_internal_id();

  my $sth = $self->_db_obj->prepare(
          "SELECT embl_id FROM clone WHERE internal_id = $internal_id");
   $sth->execute();
   my $rowhash = $sth->fetchrow_hashref();
   return $rowhash->{'embl_id'};
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

=head2 _internal_id

 Title   : _internal_id
 Usage   : $obj->_internal_id($newval)
 Function: 
 Returns : value of _internal_id
 Args    : newvalue (optional)


=cut

sub _internal_id{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'_internal_id'} = $value;
    }
    return $obj->{'_internal_id'};

}


=head2 _db_obj

 Title   : _db_obj
 Usage   : $obj->_db_obj($newval)
 Function: 
 Example : 
 Returns : value of _db_obj
 Args    : newvalue (optional)


=cut

sub _db_obj {
   my ($obj,$value) = @_;
   if( defined $value) {
      $obj->{'_db_obj'} = $value;
    }
    return $obj->{'_db_obj'};

}

1;
