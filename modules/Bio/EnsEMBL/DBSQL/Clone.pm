
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
    my $sth = $self->_db_obj->prepare("select id from clone where id = \"$id\";");    
    my $res = $sth ->execute();   
    my $rv  = $sth ->rows;    
    if( ! $rv ) {
	# make sure we deallocate sth - keeps DBI happy!
	$sth = 0;
        print STDERR "Clone $id does not seem to occur in the database!\n";     
	$self->throw("Clone $id does not seem to occur in the database!");
    }   
    return $self;
}

=head2 get_all_id

 Title   : get_all_id
 Usage   : @cloneid = $clone->get_all_id
 Function: returns all the valid (live) Clone ids in the database
 Example :
 Returns : 
 Args    :


=cut

sub get_all_id {
   my ($self) = @_;
   my @out;

   my $sth = $self->_db_obj->prepare("select id from clone");
   my $res = $sth->execute;

   while( my $rowhash = $sth->fetchrow_hashref) {
       push(@out,$rowhash->{'id'});
   }

   return @out;
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

   my $clone_id = $self->id;

   my @contigs;
   my @dnas;

   # get a list of contigs to zap
   my $sth = $self->_db_obj->prepare("select internal_id,dna from contig where clone = '$clone_id'");
   my $res = $sth->execute;

   while( my $rowhash = $sth->fetchrow_hashref) {
       push(@contigs,$rowhash->{'internal_id'});
       push(@dnas,$rowhash->{'dna'});
   }
   
   # Delete from DNA table, Contig table, Clone table
   
   foreach my $contig ( @contigs ) {
       my $sth = $self->_db_obj->prepare("delete from contig where internal_id = '$contig'");
       my $res = $sth->execute;

       my $feature_obj=Bio::EnsEMBL::DBSQL::Feature_Obj->new($self->_db_obj);
       $feature_obj->delete($contig);
   }


   foreach my $dna (@dnas) {
       $sth = $self->_db_obj->prepare("delete from dna where id = $dna");
       $res = $sth->execute;

       $sth = $self->_db_obj->prepare("delete from contigoverlap where dna_a_id = $dna or dna_b_id = $dna");
       $res = $sth ->execute;

   }

   $sth = $self->_db_obj->prepare("delete from clone where id = '$clone_id'");
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
   my ($self,$supporting) = @_;
   my @out;
   my $id = $self->id();
   my @genes;
   my @sup_exons;

   #
   # A quick trip to the database to pull out from the neighbourhood the genes
   # that might be on this clone.
   #

   my $sth = $self->_db_obj->prepare("select gene from geneclone_neighbourhood where clone = '$id'");

   $sth->execute();

   while( (my $hash = $sth->fetchrow_hashref()) ) {
       push(@genes,$hash->{'gene'});
   }
   
   #
   # A Gene/Clone Neighbourhood positive does not guarentee that a gene is
   # actually on this clone. The SQL statement needs to double check this
   #

   foreach my $geneid ( @genes ) {
       print(STDERR "Finding gene $geneid\n");
       #
       # The aim here is to get all the information for constructing the genes one
       # juicy SQL statement, effectively removing multiple SQL statement gets from this
       # construction.
       #
       
       #
       # I know this SQL statement is silly.
       #
       #  select p3.gene,
       #      p4.id,
       #      p3.id,
       #      p1.exon,p1.rank,
       #      p2.seq_start,p2.seq_end,
       #      UNIX_TIMESTAMP(p2.created),UNIX_TIMESTAMP(p2.modified),
       #      p2.strand,p2.phase,
       #      p5.seq_start,p5.start_exon,p5.seq_end,p5.end_exon,p5.id,
       #      p6.version,p3.version,p2.version,p5.version
       #  from   gene as p6,
       #      contig as p4,
       #      transcript as p3,
       #      exon_transcript as p1,
       #      exon as p2,
       #      translation as p5
       #  where  p6.id    = '$geneid'
       #  and    p3.gene  = '$geneid'
       #  and    p4.clone = '$id'
       #  and    p2.contig = p4.internal_id
       #  and    p1.exon  = p2.id
       #  and    p3.id    = p1.transcript
       #  and    p5.id    = p3.translation
       #  order by p3.gene,p3.id,p1.rank
        
       
       my $query = q{
         SELECT con.id
              , tscript.id
              , e_t.exon, e_t.rank
              , exon.seq_start, exon.seq_end
              , UNIX_TIMESTAMP(exon.created)
              , UNIX_TIMESTAMP(exon.modified)
              , exon.strand, exon.phase
              , transl.seq_start, transl.start_exon
              , transl.seq_end, transl.end_exon
              , transl.id
              , gene.version
              , tscript.version
              , exon.version
              , transl.version
            FROM gene
              , contig con
              , transcript tscript
              , exon_transcript e_t
              , exon
              , translation transl
            WHERE con.internal_id = exon.contig
              AND exon.id = e_t.exon
              AND e_t.transcript = tscript.id
              AND tscript.translation = transl.id
              AND tscript.gene = gene.id
              AND con.clone = ?
              AND gene.id = ?
            ORDER BY tscript.gene
              , tscript.id
              , e_t.rank
         };

       my $sth = $self->_db_obj->prepare($query);
       $sth->execute($id, $geneid);
       my $current_gene_id       = '';
       my $current_transcript_id = '';

       my ($gene,$trans);

       while( (my $arr = $sth->fetchrow_arrayref()) ) {
           my ($contigid,
               $transcriptid,
               $exonid, $rank,
               $start, $end,
               $exoncreated,
               $exonmodified,
               $strand, $phase,
               $trans_start, $trans_exon_start,
               $trans_end, $trans_exon_end,
               $translationid,
               $geneversion,
               $transcriptversion,
               $exonversion,
               $translationversion) = @{$arr};

           print STDERR "Got exon $exonid\n";

           if( ! defined $phase ) {
	       $self->throw("Bad internal error! Have not got all the elements in gene array retrieval");
           }

           if( $geneid ne $current_gene_id ) {

	       if( $transcriptid eq $current_transcript_id ) {
	           $self->throw("Bad internal error. Switching genes without switching transcripts");
	       } 

	       $gene = Bio::EnsEMBL::Gene->new();
	       $gene->id($geneid);

	    #   $sth = $self->_db_obj->prepare("select version from gene where id='".$gene->id."'");
	    #   $sth->execute();
	    #   my $rowhash = $sth->fetchrow_hashref();

	       $gene->version($geneversion);
	       $gene->add_cloneid_neighbourhood($id);

	       $current_gene_id = $geneid;
	       push(@out,$gene);
	       #print STDERR "Made new gene\n";
           }

           if( $transcriptid ne $current_transcript_id ) {
	       $trans = Bio::EnsEMBL::Transcript->new();
	       $trans->id($transcriptid);
	       $trans->version($transcriptversion);
	       $current_transcript_id = $transcriptid;

	       my $translation = Bio::EnsEMBL::Translation->new();
	       $translation->start        ($trans_start);
	       $translation->end          ($trans_end);
	       $translation->start_exon_id($trans_exon_start);
	       $translation->end_exon_id  ($trans_exon_end);
	       $translation->id           ($translationid);
	       $translation->version      ($translationversion);
	       $trans->translation        ($translation);
	       $gene ->add_Transcript     ($trans);
           }

           my $exon = Bio::EnsEMBL::Exon->new();
           #print(STDERR "Creating exon in clone  = $contigid\n");
           $exon->clone_id ($id);
           $exon->contig_id($contigid);
           $exon->id       ($exonid);
           $exon->created  ($exoncreated);
           $exon->modified ($exonmodified);
           $exon->start    ($start);
           $exon->end      ($end);
           $exon->strand   ($strand);
           $exon->phase    ($phase);
           $exon->version  ($exonversion);
           #
           # Attach the sequence, cached if necessary...
           #

           if ($supporting && $supporting eq 'evidence') {
	       push @sup_exons, $exon;
           }

           my $seq;

           if( $self->_db_obj->_contig_seq_cache($exon->contig_id) ) {
	       $seq = $self->_db_obj->_contig_seq_cache($exon->contig_id);
           } else {

	       my $contig      = $self->_db_obj->get_Contig(-id => $exon->contig_id );
	       $contig->fetch();

	       $seq = $contig->primary_seq();
	       $self->_db_obj->_contig_seq_cache($exon->contig_id,$seq);
           }

           $exon ->attach_seq(new Bio::Seq(-id  => $seq->id,
				           -seq => $seq->seq));
           $trans->add_Exon($exon);
       }
    }

    if ($supporting && $supporting eq 'evidence' && @sup_exons != 0) {

       my $gene_obj=Bio::EnsEMBL::DBSQL::Gene_Obj->new($self->_db_obj);
       $gene_obj->get_supporting_evidence(@sup_exons);
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

sub get_all_Genes_slow {
   my ($self,@args) = @_;
   my @out;
   my $id = $self->id();
   my %got;

   my $sth = $self->_db_obj->prepare("select p3.gene from contig as p4, transcript as p3, exon_transcript as p1, exon as p2, geneclone_neighbourhood as p5 where p5.clone = '$id' and p5.gene = p3.gene and p4.clone = '$id' and p2.contig = p4.internal_id and p1.exon = p2.id and p3.id = p1.transcript");
   
   my $res = $sth->execute();
   while( my $rowhash = $sth->fetchrow_hashref) {
       if( ! exists $got{$rowhash->{'gene'}} ) {
	   
	   my $gene_obj=Bio::EnsEMBL::DBSQL::Gene_Obj->new($self->_db_obj);
	   my $gene = $gene_obj->get($rowhash->{'gene'});

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

   my $cloneid = $self->id;

   my $sth = $self->_db_obj->prepare("select count(*),cont.clone ,ex.contig,tran.gene  " .
			    "from   contig          as cont, ".
			    "       transcript      as tran, " .
			    "       exon_transcript as et, " .
			    "       exon            as ex " .
			    "where  ex.id            = et.exon " .
			    "and    tran.id          = et.transcript " .
			    "and    cont.clone       = '$cloneid'  " .
			    "and    cont.internal_id = ex.contig " .
			    "group by tran.gene");

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
   my $name = $self->id();

   my $sql = "select id,internal_id from contig where clone = \"$name\" ";

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
       $self->throw("Clone $name has no contigs in the database. Should be impossible, but clearly isn't...");
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
    
    my @overlaps;

    foreach my $contig ($self->get_all_Contigs) {
	if (defined($contig->get_left_overlap)) {
	    
	    my $overlap    = $contig->get_left_overlap;
	    my $type;
	    
	    if ($overlap->sister_polarity == 1) {
		$type = 'left2right';
	    } elsif ($overlap->sister_polarity == -1) {
		$type = 'left2left';
	    } else {
		$self->throw("Invalid value [" .$overlap->sister_polarity . "] for polarity");
	    }
	    
	    my $tmpoverlap = new Bio::EnsEMBL::ContigOverlap(-contiga => $contig,
							     -contigb => $overlap->sister,
							     -positiona => $overlap->self_position,
							     -positionb => $overlap->sister_position,
							     -source    => $overlap->source,
							     -distance  => $overlap->distance,
							     -overlap_type => $type);
	    
	    push(@overlaps,$tmpoverlap);
	}

	if (defined($contig->get_right_overlap)) {
	    
	    my $overlap    = $contig->get_right_overlap;
	    my $type;
	    
	    if ($overlap->sister_polarity == 1) {
		$type = 'right2left';
	    } elsif ($overlap->sister_polarity == -1) {
		$type = 'right2right';
	    } else {
		$self->throw("Invalid value [" .$overlap->sister_polarity . "] for polarity");
	    }
	    
	    my $tmpoverlap = new Bio::EnsEMBL::ContigOverlap(-contiga => $contig,
							     -contigb => $overlap->sister,
							     -positiona => $overlap->self_position,
							     -positionb => $overlap->sister_position,
							     -source    => $overlap->source,
							     -distance  => $overlap->distance,
							     -overlap_type => $type);
	    
	    push(@overlaps,$tmpoverlap);
	    
	}
    }

    return (@overlaps);
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

   my $id = $self->id();

   my $sth = $self->_db_obj->prepare("select htg_phase from clone where id = \"$id\" ");
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

sub created {
   my ($self) = @_;

   my $id = $self->id();

   my $sth = $self->_db_obj->prepare("select UNIX_TIMESTAMP(created) from clone where id = \"$id\" ");
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

   my $id = $self->id();

   my $sth = $self->_db_obj->prepare("select modified from clone where id = \"$id\" ");
   $sth->execute();
   my $rowhash = $sth->fetchrow_hashref();
   my $datetime = $rowhash->{'modified'};
   $sth = $self->_db_obj->prepare("select UNIX_TIMESTAMP('".$datetime."')");
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

sub version {
   my $self = shift;
   my $id = $self->id();

   my $sth = $self->_db_obj->prepare("select version from clone where id = \"$id\" ");
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
   my $id = $self->id();

   my $sth = $self->_db_obj->prepare("select embl_version from clone where id = \"$id\" ");
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

   my $id = $self->id();

   my $sth = $self->_db_obj->prepare("select embl_id from clone where id = \"$id\" ");
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
