
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
use Bio::EnsEMBL::DBSQL::Contig;
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

  $id || $self->throw("Cannot make contig db object without id");
  $dbobj || $self->throw("Cannot make contig db object without db object");
  $dbobj->isa('Bio::EnsEMBL::DBSQL::Obj') || $self->throw("Cannot make contig db object with a $dbobj object");

  $self->id($id);
  $self->_dbobj($dbobj);

# set stuff in self from @args
  return $make; # success - we hope!
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
   my ($self,@args) = @_;
   my @out;
   my $id = $self->id();
   my @genes;

   #
   # A quick trip to the database to pull out from the neighbourhood the genes
   # that might be on this clone.
   #

   my $sth = $self->_dbobj->prepare("select gene from geneclone_neighbourhood where clone = '$id'");
   $sth->execute();
   while( (my $hash = $sth->fetchrow_hashref()) ) {
       push(@genes,$hash->{'gene'});
   }


   #
   # A Gene/Clone Neighbourhood positive does not guarentee that a gene is
   # actually on this clone. The SQL statement needs to double check this
   #

   foreach my $geneid ( @genes ) {
       #
       # The aim here is to get all the information for constructing the genes one
       # juicy SQL statement, effectively removing multiple SQL statement gets from this
       # construction.
       #
       
       #
       # I know this SQL statement is silly.
       #
       
       $sth = $self->_dbobj->prepare("select p3.gene,p4.id,p3.id,p1.exon,p1.rank,p2.seq_start,p2.seq_end,p2.created,p2.modified,p2.strand,p2.phase,p5.seq_start,p5.start_exon,p5.seq_end,p5.end_exon,p5.id from contig as p4, transcript as p3, exon_transcript as p1, exon as p2,translation as p5 where p3.gene = '$geneid' and p4.clone = '$id' and p2.contig = p4.id and p1.exon = p2.id and p3.id = p1.transcript and p5.id = p3.translation order by p3.gene,p3.id,p1.rank");
   
       $sth->execute();
       my $current_gene_id = '';
       my $current_transcript_id = '';
       my ($gene,$trans);
       while( (my $arr = $sth->fetchrow_arrayref()) ) {
	   my ($geneid,$contigid,$transcriptid,$exonid,$rank,$start,$end,$exoncreated,$exonmodified,$strand,$phase,$trans_start,$trans_exon_start,$trans_end,$trans_exon_end,$translationid) = @{$arr};
	   if( ! defined $phase ) {
	       $self->throw("Bad internal error! Have not got all the elements in gene array retrieval");
	   }

	   if( $geneid ne $current_gene_id ) {
	       if( $transcriptid eq $current_transcript_id ) {
		   $self->throw("Bad internal error. Switching genes without switching transcripts");
	       } 
	       $gene = Bio::EnsEMBL::Gene->new();
	       $gene->id($geneid);
	       push(@out,$gene);
	   }
	   if( $transcriptid != $current_transcript_id ) {
	       $trans = Bio::EnsEMBL::Transcript->new();
	       $trans->id($transcriptid);
	       my $translation = Bio::EnsEMBL::Translation->new();
	       $translation->start($trans_start);
	       $translation->end($trans_end);
	       $translation->start_exon_id($trans_exon_start);
	       $translation->end_exon_id($trans_exon_end);
	       $translation->id($translationid);
	       $gene->add_Transcript($trans);
	   }
	   
	   my $exon = Bio::EnsEMBL::Exon->new();
	   $exon->clone_id($id);
	   $exon->contig_id($contigid);
	   $exon->id($exonid);
	   $exon->created($exoncreated);
	   $exon->modified($exonmodified);
	   $exon->start($start);
	   $exon->end($end);
	   $exon->strand($strand);
	   $exon->phase($phase);
	   
	   #
	   # Attach the sequence, cached if necessary...
	   #
	   
	   my $seq;
	   
	   if( $self->_dbobj->_contig_seq_cache($exon->contig_id) ) {
	       $seq = $self->_dbobj->_contig_seq_cache($exon->contig_id);
	   } else {
	       my $contig = $self->_dbobj->get_Contig($exon->contig_id());
	       $seq = $contig->seq();
	       $self->_dbobj->_contig_seq_cache($exon->contig_id,$seq);
	   }
	   
	   $exon->attach_seq($seq);
	   $trans->add_Exon($exon);
       }
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

sub get_all_Genes_slow{
   my ($self,@args) = @_;
   my @out;
   my $id = $self->id();
   my %got;

   my $sth = $self->_dbobj->prepare("select p3.gene from contig as p4, transcript as p3, exon_transcript as p1, exon as p2, geneclone_neighbourhood as p5 where p5.clone = '$id' and p5.gene = p3.gene and p4.clone = '$id' and p2.contig = p4.id and p1.exon = p2.id and p3.id = p1.transcript");
   
   my $res = $sth->execute();
   while( my $rowhash = $sth->fetchrow_hashref) {
       if( ! exists $got{$rowhash->{'gene'}} ) {
	   my $gene = $self->_dbobj->get_Gene($rowhash->{'gene'});
	   push(@out,$gene);
	   $got{$rowhash->{'gene'}} = 1;
       }    
   }
   return @out;
}

=head2 get_Contig

 Title   : get_Contig
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_Contig{
   my ($self,$contigid) = @_;

   # should check this contig is in this clone?
   return $self->_dbobj->get_Contig($contigid);
}

=head2 get_all_Contigs

 Title   : get_Contigs
 Usage   : foreach $contig ( $clone->get_all_Contigs ) 
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_all_Contigs{
   my ($self) = @_;
   my $sth;
   my @res;
   my $name = $self->id();

   my $sql = "select id from contig where clone = \"$name\" ";
#   print STDERR "Looking at $name. [$sql]\n";
   $sth= $self->_dbobj->prepare($sql);
   my $res = $sth->execute();
   my $seen = 0;

   my $count = 0;
   my $total = 0;

   while( my $rowhash = $sth->fetchrow_hashref) {
       my $contig = new Bio::EnsEMBL::DBSQL::Contig ( -dbobj => $self->_dbobj,
						   -id => $rowhash->{'id'} );
       $contig->order($count++);
       #$contig->offset($total);
       
       #$total += $contig->length();
       #$total += 400;

       push(@res,$contig);
       $seen = 1;
   }
   if( $seen == 0  ) {
       $self->throw("Clone $name has no contigs in the database. Should be impossible, but clearly isnt...");
   }

   return @res;   
}

=head2 htg_phase

 Title   : htg_phase
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub htg_phase{
   my $self = shift;

   my $id = $self->id();

   my $sth = $self->_dbobj->prepare("select htg_phase from clone where id = \"$id\" ");
   $sth->execute();
   my $rowhash = $sth->fetchrow_hashref();
   return $rowhash->{'htg_phase'};
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

sub created{
   my ($self) = @_;

   my $id = $self->id();

   my $sth = $self->_dbobj->prepare("select created from clone where id = \"$id\" ");
   $sth->execute();
   my $rowhash = $sth->fetchrow_hashref();
   my $datetime = $rowhash->{'created'};
   $sth = $self->_dbobj->prepare("select UNIX_TIMESTAMP('".$datetime."')");
   $sth->execute();
   $rowhash = $sth->fetchrow_arrayref();
   return $rowhash->[0];
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

sub modified{
   my ($self) = @_;

   my $id = $self->id();

   my $sth = $self->_dbobj->prepare("select modified from clone where id = \"$id\" ");
   $sth->execute();
   my $rowhash = $sth->fetchrow_hashref();
   my $datetime = $rowhash->{'modified'};
   $sth = $self->_dbobj->prepare("select UNIX_TIMESTAMP('".$datetime."')");
   $sth->execute();
   $rowhash = $sth->fetchrow_arrayref();
   return $rowhash->[0];
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
   my $self = shift;
   my $id = $self->id();

   my $sth = $self->_dbobj->prepare("select version from clone where id = \"$id\" ");
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

sub _stored{
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

sub embl_version{
   my $self = shift;
   my $id = $self->id();

   my $sth = $self->_dbobj->prepare("select embl_version from clone where id = \"$id\" ");
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

sub seq_date{
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

sub embl_id{
   my ($self) = @_;

   my $id = $self->id();

   my $sth = $self->_dbobj->prepare("select embl_id from clone where id = \"$id\" ");
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


=head2 _dbobj

 Title   : _dbobj
 Usage   : $obj->_dbobj($newval)
 Function: 
 Example : 
 Returns : value of _dbobj
 Args    : newvalue (optional)


=cut

sub _dbobj{
   my ($obj,$value) = @_;
   if( defined $value) {
      $obj->{'_dbobj'} = $value;
    }
    return $obj->{'_dbobj'};

}

1;
