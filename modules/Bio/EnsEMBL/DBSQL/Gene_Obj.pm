#
# EnsEMBL module for Bio::EnsEMBL::DBSQL::Gene_Obj
#
# Cared for by Elia Stupka <elia@ebi.ac.uk>
#
# Copyright Elia Stupka
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::DBSQL::Gene_Obj - MySQL database adapter class for EnsEMBL genes, transcripts,
exons, etc.

=head1 SYNOPSIS

  $gene   = $gene_obj->get('HG45501');

  use Bio::EnsEMBL::Gene;
  use Bio::EnsEMBL::DBSQL::Gene_Obj;

  # Get a gene object from the database
  my $gene = $gene_obj->get('HG45501', $db_obj);

=head1 DESCRIPTION

This is one of the objects contained in Bio:EnsEMBL::DBSQL::Obj,
dealing with Gene methods, such as writing and getting genes,
transcripts, translations, and exons.

The Obj object represents a database that is implemented somehow (you
shouldn\'t care much as long as you can get the object).

=head1 CONTACT

Elia Stupka: elia@ebi.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...

package Bio::EnsEMBL::DBSQL::Gene_Obj;

use vars qw(@ISA);
use strict;

# Object preamble - inheriets from Bio::Root::Object

use Bio::Root::Object;
use Bio::EnsEMBL::DBSQL::Obj;
use Bio::EnsEMBL::DB::Gene_ObjI;
use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::Exon;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::DB::VirtualContig;
use DBI;
use Bio::EnsEMBL::StickyExon;

use Bio::EnsEMBL::DBSQL::DummyStatement;
use Bio::EnsEMBL::DB::Gene_ObjI;

@ISA = qw(Bio::EnsEMBL::DB::Gene_ObjI Bio::Root::Object);

# new() is inherited from Bio::Root::Object

# _initialize is where the heavy stuff will happen when new is called

sub _initialize {
  my($self,$db_obj) = @_;

  my $make = $self->SUPER::_initialize;
  
  $db_obj || $self->throw("Database Gene object must be passed a db obj!");
  $self->_db_obj($db_obj);
  $self->use_delayed_insert(1);
  return $make; # success - we hope!

}

=head2 delete

 Title   : delete
 Usage   : $Gene_Obj->delete_Gene($gene_id)
 Function: deletes a gene from the database, i.e. exons, transcripts, translations
 Example : $geneobj->delete_Gene('ENSG00000019482')
 Returns : nothing
 Args    : gene id

=cut

sub delete {
    my ($self, $gene_id) = @_;
    
    my $db = $self->_db_obj;
    
    # Get the transcript, translation and exon IDs for this gene
    my $sth = $db->prepare(q{
        SELECT t.id
          , t.translation
          , et.exon
        FROM gene g
          , transcript t
          , exon_transcript et
          , exon e
        WHERE g.id = t.gene
          AND t.id = et.transcript
          AND et.exon = e.id
          AND g.id = ?
        });
    $sth->execute($gene_id);
    
    my( %transcript, %translation, %exon );
    while (my $row = $sth->fetchrow_arrayref) {
         $transcript{$row->[0]} = 1;
        $translation{$row->[1]} = 1;
               $exon{$row->[2]} = 1;
    }
    
    # Deletes which use the gene ID
    my $gene_delete      = $db->prepare(q{DELETE FROM gene WHERE id = ?});
    my $gene_type_delete = $db->prepare(q{DELETE FROM genetype WHERE gene_id = ?});

    $gene_delete     ->execute($gene_id);
    $gene_type_delete->execute($gene_id);
    
    # Deletes which use the transcript ID
    my $transcript_delete      = $db->prepare(q{DELETE FROM transcript WHERE id = ?});
    my $exon_transcript_delete = $db->prepare(q{DELETE FROM exon_transcript WHERE transcript = ?});
    
    foreach my $trans_id (keys %transcript) {
        $transcript_delete     ->execute($trans_id);
        $exon_transcript_delete->execute($trans_id);
    }
    
    # Translation delete
    my $translation_delete = $db->prepare(q{DELETE FROM translation WHERE id = ?});
    
    foreach my $transl_id (keys %translation) {
        $translation_delete->execute($transl_id);
    }
    
    # Deletes which use the exon ID
    my $exon_delete       = $db->prepare(q{DELETE FROM exon WHERE id = ?});
    my $supporting_delete = $db->prepare(q{DELETE FROM supporting_feature WHERE exon = ?});
    
    foreach my $exon_id (keys %exon) {
        $exon_delete      ->execute($exon_id);
        $supporting_delete->execute($exon_id);
    }
}   



=head2 delete_Exon

 Title   : delete_Exon
 Usage   : $obj->delete_Exon($exon_id)
 Function: Deletes exon, including exon_transcript rows
 Example : $obj->delete_Exon(ENSE000034)
 Returns : nothing
 Args    : $exon_id

=cut

sub delete_Exon{
    my ($self,$exon_id) = @_;

    $exon_id || $self->throw ("Trying to delete an exon without an exon_id\n");
    
    #Delete exon_transcript rows
    my $sth = $self->_db_obj->prepare("delete from exon_transcript where exon = '$exon_id'");
    $sth->execute;

    #Delete exon rows
    $sth = $self->_db_obj->prepare("delete from exon where id = '$exon_id'");
    $sth->execute;

    $self->delete_Supporting_Evidence($exon_id);
}

=head2 delete_Supporting_Evidence

 Title   : delete_Supporting_Evidence
 Usage   : $obj->delete_Supporting_Evidence($exon_id)
 Function: Deletes exon\'s supporting evidence entries
 Example : $obj->delete_Supporting_Evidence(ENSE000034)
 Returns : nothing
 Args    : $exon_id


=cut

sub delete_Supporting_Evidence {
    my ($self,$exon_id) = @_;

    $exon_id || $self->throw ("Trying to delete supporting_evidence without an exon_id\n");

    my $sth = $self->_db_obj->prepare("delete from supporting_feature where exon = '" . $exon_id . "'");
    my $res = $sth->execute;
}


=head2 get_all_Gene_id

 Title   : get_all_Gene_id
 Usage   : $geneobj->get_all_Gene_id
 Function: Gets an array of ids for all genes in the current db
 Example : $geneobj->get_all_Gene_id
 Returns : array of ids
 Args    : none

=cut

sub get_all_Gene_id{
   my ($self) = @_;

   my @out;
   my $sth = $self->_db_obj->prepare("select id from gene");
   my $res = $sth->execute || $self->throw("Could not get any gene ids!");
   while (my ($id) = $sth->fetchrow) {
       push(@out, $id);
   }
   return @out;
}



=head2 get_all_Transcript_id

 Title   : get_all_Transcript_id
 Usage   : $geneobj->get_all_Transcript_id
 Function: Gets an array of ids for all genes in the current db
 Example : $geneobj->get_all_Transcript_id
 Returns : array of ids
 Args    : none

=cut

sub get_all_Transcript_id{
   my ($self) = @_;

   my @out;
   my $sth = $self->_db_obj->prepare("select id from transcript");
   my $res = $sth->execute || $self->throw("Could not get any transcript ids!");
   while( my $rowhash = $sth->fetchrow_hashref) {
       push(@out,$rowhash->{'id'});
   }
   return @out;
}




=head2 get_Gene_by_Transcript_id

 Title   : get_Gene_by_Transcript_id
 Usage   : $gene_obj->get_Gene_by_Transcript_id($transid, $supporting)
 Function: gets one gene out of the db with or without supporting evidence
 Returns : gene object (with transcripts, exons and supp.evidence if wanted)
 Args    : transcript id and supporting tag (if latter not specified, assumes without
           Note that it is much faster to get genes without supp.evidence!


=cut

sub get_Gene_by_Transcript_id {
    my $self = shift;
    my $transid = shift;
    my $supporting = shift;

    # this is a cheap SQL call
    my $sth = $self->_db_obj->prepare("select gene from transcript where id = '$transid'");
    $sth->execute;

    my ($geneid) = $sth->fetchrow_array();
    if( !defined $geneid ) {
        return undef;
    }
    return $self->get($geneid,$supporting);
}

=head2 get_geneids_by_hids

 Title   : get_geneids_by_hids
 Usage   : @geneids = $obj->get_geneids_by_hids(@hids)
 Function: gives back geneids with these hids as supporting evidence
 Example :
 Returns : 
 Args    :


=cut

sub get_geneids_by_hids{
   my ($self,@hids) = @_;

    my $inlist = join(',',map "'$_'", @hids);
       $inlist = "($inlist)";

   my $sth = $self->_db_obj->prepare("select transcript.gene from transcript as transcript, exon_transcript as exon_transcript, exon as exon, supporting_feature as supporting_feature where exon.id = supporting_feature.exon and exon_transcript.exon = exon.id and exon_transcript.transcript = transcript.id and supporting_feature.hid in $inlist");

   $sth->execute();
   my %gene;

   while( (my $arr = $sth->fetchrow_arrayref()) ) {
       my ($geneid) = @{$arr};
       $gene{$geneid} =1;
   }

   return keys %gene;
}

=head2 get

 Title   : get
 Usage   : $geneobj->get($geneid, $supporting)
 Function: gets one gene out of the db with or without supporting evidence
 Example : $obj->get('ENSG00000009151','evidence')
 Returns : gene object (with transcripts, exons and supp.evidence if wanted)
 Args    : gene id and supporting tag (if latter not specified, assumes without
	   Note that it is much faster to get genes without supp.evidence!

=cut

sub get {
    my ($self,$geneid, $supporting) = @_;
    
    my @out;
    
    if (!$supporting) {
        @out = $self->get_array_supporting('without', $geneid);
    }
    else {
        @out = $self->get_array_supporting($supporting, $geneid);
    }
    
    $self->throw("Error retrieving gene with ID: $geneid") unless $out[0]; 
    
    return $out[0];
}


=head2 get_array_supporting

    Title   : get_array_supporting
    Usage   : $obj->get_array_supporting($supporting,@geneid)
    Function: Gets an array of genes, with transcripts and exons. If $supporting
           equal to 'evidence' the supporting evidence for each exon is also read
    from the supporting evidence table
    Example : $obj->get_array_supporting ('evidence',@geneid)
    Returns : an array of gene objects
    Args    : 'evidence' and gene id array

    
=cut
    
sub get_array_supporting {
    my ($self,$supporting,@geneid) = @_;
    
    $supporting || $self->throw("You need to specify whether to retrieve supporting evidence or not!");

    if( @geneid == 0 ) {
	$self->throw("Attempting to create gene with no id");
    }
   
    my (@out, @sup_exons);
    
    my $inlist = join(',', map "'$_'", @geneid);
    
    # I know this SQL statement is silly.
    #    
     
    my $query = qq{
        SELECT tscript.gene
          , con.id
          , tscript.id
          , e_t.exon, e_t.rank
          , exon.seq_start, exon.seq_end
          , UNIX_TIMESTAMP(exon.created)
          , UNIX_TIMESTAMP(exon.modified)
          , exon.strand
          , exon.phase
	  , exon.sticky_rank
          , transl.seq_start, transl.start_exon
          , transl.seq_end, transl.end_exon
          , transl.id
          , gene.version
          , UNIX_TIMESTAMP(gene.created)
	  , UNIX_TIMESTAMP(gene.modified)
          , UNIX_TIMESTAMP(gene.stored)
	  , tscript.version
          , exon.version
          , transl.version
          , con.clone
	  , genetype.type
        FROM contig con
          , gene
          , transcript tscript
          , exon_transcript e_t
          , exon
          , translation transl
	  , genetype
        WHERE gene.id = tscript.gene
          AND tscript.id = e_t.transcript
          AND e_t.exon = exon.id
          AND exon.contig = con.internal_id
          AND tscript.translation = transl.id
          AND gene.id IN ($inlist)
          AND genetype.gene_id = gene.id
        ORDER BY tscript.gene
          , tscript.id
          , e_t.rank
          , exon.sticky_rank
        };

    #print STDERR "query [$query]\n"; 

    my $sth = $self->_db_obj->prepare($query);
    my $res = $sth ->execute();
   
    my $current_gene_id       = '';
    my $current_transcript_id = '';
    my $previous_exon = undef;
    my $sticky_exon = 0;
    
    my ($gene,$trans);
    my @transcript_exons;
    
    while( (my $arr = $sth->fetchrow_arrayref()) ) {
	#print STDERR "Getting into this row now....\n";

	my ($geneid,$contigid,$transcriptid,$exonid,$rank,$start,$end,
	    $exoncreated,$exonmodified,$strand,$phase,$exon_rank,$trans_start,
	    $trans_exon_start,$trans_end,$trans_exon_end,$translationid,
	    $geneversion,$genecreated,$genemodified,$genestored,$transcriptversion,$exonversion,$translationversion,$cloneid,$genetype) = @{$arr};

 	
	if( ! defined $phase ) {
	    $self->throw("Bad internal error! Have not got all the elements in gene array retrieval");
	}
	
	# I think this is a dirty hack 
	#if( exists $seen{"$exonid-$rank"} ) {
	#    next;
	#}

	
	
	# Create new gene if the id has changed
	if( $geneid ne $current_gene_id ) {
	    
	    if( $transcriptid eq $current_transcript_id ) {
		$self->throw("Bad internal error. Switching genes without switching transcripts");
	    } 
	    
	    $gene = Bio::EnsEMBL::Gene->new();
	    
	    $gene->id                       ($geneid);
	    $gene->type                     ($genetype);
	    $gene->version                  ($geneversion);
	    $gene->created                  ($genecreated);
	    $gene->modified                 ($genemodified);
	    $gene->_stored                  ($genestored);
	    $gene->add_cloneid_neighbourhood($cloneid);
	    
	    $current_gene_id = $geneid;
	    push(@out,$gene);
	    
	}
	
	# Create new transcript if the id has changed
	if( $transcriptid ne $current_transcript_id ) {

	    # put away old exons
             if( defined $trans ) {
		#print "Storing\n";
	        $self->_store_exons_in_transcript($trans,@transcript_exons);
            }
	    # put in new exons
	    
	    $trans = Bio::EnsEMBL::Transcript->new();
	    
	    $trans->id     ($transcriptid);
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

	#print(STDERR "Creating exon - contig id $contigid\n");

	$exon->clone_id ($cloneid);
	$exon->contig_id($contigid);
	$exon->id       ($exonid);
	$exon->created  ($exoncreated);
	$exon->modified ($exonmodified);
	$exon->start    ($start);
	$exon->end      ($end);
	$exon->phase   ($phase);
	$exon->strand    ($strand);
	$exon->version  ($exonversion);
	$exon->seqname  ($contigid);
	$exon->sticky_rank($exon_rank);
        
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
	    my $contig      = $self->_db_obj->get_Contig($exon->contig_id);
	    $contig->fetch(); 
	    $seq = $contig->primary_seq();
	    $self->_db_obj->_contig_seq_cache($exon->contig_id,$seq);
	}

	$exon ->attach_seq($seq);
	push(@transcript_exons,$exon);


    }
    if( $current_gene_id eq '' ) {
	return ();
    }
    
    $self->_store_exons_in_transcript($trans,@transcript_exons);
   
    if ($supporting && $supporting eq 'evidence') {
	$self->get_supporting_evidence(@sup_exons);
    }

    foreach my $g ( @out) {
	$self->_get_dblinks($g);
    }

    
    return @out;
}

=head2 _store_exons_in_transcript

 Title   : _store_exons_in_transcript
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub _store_exons_in_transcript{
   my ($self,$trans,@exons) = @_;

   if( !ref $trans || !$trans->isa('Bio::EnsEMBL::Transcript') ) {
       $self->throw(" $trans is not a transcript");
   }
   #print STDERR "Got ",scalar(@exons),"to store...\n";

   my $exon;
   while ( ($exon = shift @exons)) {
       #print STDERR "Handling exon",$exon->id,":",$exon->sticky_rank,"\n";

       if( $#exons >= 0 && $exons[0]->id eq $exon->id ) {
        
	   # sticky exons.
	   my @sticky_exons;
	   push(@sticky_exons,$exon);
	   while( my $newexon = shift @exons ) {
	       if( $newexon->id eq $exon->id ) {
                            
		   push(@sticky_exons,$newexon);
                   
	       } else {
               
		   unshift(@exons,$exon);
		   last;
	       }
	   }
           
	   my $sticky = $self->_make_sticky_exon(@sticky_exons);
	   #print STDERR "Added sticky exon... $sticky\n";
	   $trans->add_Exon($sticky);
           
       } else {
	   $trans->add_Exon($exon);
       }
   }

}

=head2 _make_sticky_exon

 Title   : _make_sticky_exon
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub _make_sticky_exon{
   my ($self,@exons) = @_;

   my $sticky = Bio::EnsEMBL::StickyExon->new();
   my $seq  = Bio::PrimarySeq->new();
   my $seqstr = "";

   @exons = sort { $a->sticky_rank <=> $b->sticky_rank } ( @exons );
   $seq->id("exon.sticky_contig.".$exons[0]->id);
   
   $sticky->id($exons[0]->id);
   $sticky->phase($exons[0]->phase);
   $sticky->contig_id($seq->id);
   $sticky->clone_id ($exons[0]->clone_id);
   $sticky->created  ($exons[0]->created);
   $sticky->modified ($exons[0]->modified);
	
   $sticky->version  ($exons[0]->version);
   $sticky->seqname  ($seq->id);

   foreach my $exon ( @exons ) {
       #print STDERR "Exon ",$exon->start," ",$exon->end," ",$exon->seqname,"\n";
       $seqstr .= $exon->seq->seq();
       #print STDERR "Sticking in ",$exon->id,":",$exon->sticky_rank," $seqstr\n";

       $sticky->add_component_Exon($exon);
   }

   $seq->seq($seqstr);
   $seq->display_id("sequence.join.".$exons[0]->id);
   $sticky->start    (1);
   $sticky->end      ($seq->length);
   $sticky->strand   (1);
   $sticky->seqname  ($seq->id);
   $sticky->attach_seq($seq);
   return $sticky;

}

=head2 _get_dblinks

 Title   : _get_dblinks
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub _get_dblinks{
   my ($self,$gene) = @_;

   if( !defined $gene || ! ref $gene || !$gene->isa('Bio::EnsEMBL::Gene') ) {
       $self->throw("no gene passed to get_dblinks");
   }
   my $geneid = $gene->id;

   my $query = "select external_db,external_id from genedblink where gene_id = '$geneid'";
   my $sth = $self->_db_obj->prepare($query);
   my $res = $sth ->execute();
   while( (my $hash = $sth->fetchrow_hashref()) ) {
       my $dblink = Bio::Annotation::DBLink->new();
       $dblink->database($hash->{'external_db'});
       $dblink->primary_id($hash->{'external_id'});
       $gene->add_DBLink($dblink);
   }

   foreach my $trans ( $gene->each_Transcript ) {
       my $transid = $trans->id;
       
       $query = "select external_db,external_id from transcriptdblink where transcript_id = '$transid'";
       $sth = $self->_db_obj->prepare($query);
       $res = $sth ->execute();
       while( (my $hash = $sth->fetchrow_hashref()) ) {
	   
	   my $dblink = Bio::Annotation::DBLink->new();
	   $dblink->database($hash->{'external_db'});
	   $dblink->primary_id($hash->{'external_id'});
	   $trans->add_DBLink($dblink);
       }
   }



}



=head2 get_Gene_by_DBLink

 Title   : get_Gene_by_DBLink
 Usage   : $gene_obj->get_Gene_by_DBLink($ext_id, $supporting)
 Function: gets one gene out of the db with or without supporting evidence
 Returns : gene object (with transcripts, exons and supp.evidence if wanted)
 Args    : transcript id and supporting tag (if latter not specified,
assumes without
           Note that it is much faster to get genes without supp.evidence!


=cut

sub get_Gene_by_DBLink {
    my $self = shift;
    my $external_id = shift;
    my $supporting = shift;
   
    my $sth = $self->_db_obj->prepare("select gene_id from genedblink where external_id = '$external_id'");
    $sth->execute;

    my ($geneid) = $sth->fetchrow_array();
    if( !defined $geneid ) {
        return undef;
    }
    return $self->get($geneid,$supporting);
}

=head2 get_Exon

 Title   : get_Exon
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_Exon{
   my ($self,$exonid) = @_;

   my $sth     = $self->_db_obj->prepare("select e.id as exonid,e.version,e.contig," .
 			        "       UNIX_TIMESTAMP(e.created),UNIX_TIMESTAMP(e.modified), " .
				"       e.seq_start,e.seq_end,e.strand,e.phase, " .
				"       c.id as contigid " .
				"from   exon as e," .
				"       contig as c " .
				"where  e.id = '$exonid'" . 
				"and    e.contig = c.internal_id");

   my $res     = $sth->execute;
   my $rowhash = $sth->fetchrow_hashref;

   if( ! defined $rowhash ) {
       $self->throw("No exon of this id $exonid");
   }
   my $exon = Bio::EnsEMBL::Exon->new();

   $exon->contig_id($rowhash->{'contigid'});
   $exon->seqname($rowhash->{'contigid'});
   $exon->version  ($rowhash->{'version'});
   
   my $contig_id = $exon->contig_id();

   # we have to make another trip to the database to get out the contig to clone mapping.
   my $sth2     = $self->_db_obj->prepare("select clone from contig where
id = '$contig_id'");
   my $res2     = $sth2->execute;
   my $rowhash2 = $sth2->fetchrow_hashref;

   $exon->clone_id($rowhash2->{'clone'});

   # rest of the attributes
   $exon->id      ($rowhash->{'exonid'});
   $exon->created ($rowhash->{'UNIX_TIMESTAMP(created)'});
   $exon->modified($rowhash->{'UNIX_TIMESTAMP(modified)'});
   $exon->start   ($rowhash->{'seq_start'});
   $exon->end     ($rowhash->{'seq_end'});
   $exon->strand  ($rowhash->{'strand'});
   $exon->phase   ($rowhash->{'phase'});
   
   # we need to attach this to a sequence. For the moment, do it the stupid
   # way perhaps?

   my $seq;

   if( $self->_db_obj->_contig_seq_cache($exon->contig_id) ) {
       $seq = $self->_db_obj->_contig_seq_cache($exon->contig_id);
   } else {
       my $contig  = $self->_db_obj->get_Contig($exon->contig_id());
       $contig->fetch(); 
       $seq = $contig->primary_seq();
       $self->_db_obj->_contig_seq_cache($exon->contig_id,$seq);
   }

   $exon->attach_seq($seq);

   return $exon;
}

=head2 get_supporting_evidence

 Title   : get_supporting_evidence
 Usage   : $obj->get_supporting_evidence
 Function: 
 Example :
 Returns : nothing
 Args    : array of exon objects, needed to know which exon to attach the evidence to


=cut

sub get_supporting_evidence {
    my ($self,@exons) = @_;

    my $instring = "'";
    my %exhash;

    if (@exons == 0) {
	$self->throw("No exon objects were passed on!");
    }

    foreach my $exon (@exons) {

	$exhash{$exon->id} = $exon;

	$instring = $instring . $exon->id . "','";
    }
    $instring = substr($instring,0,-2);
   
    my $statement = "select * from supporting_feature where exon in (" . $instring . ")";
    #print STDERR "going to execute... [$statement]\n";

    my $sth = $self->_db_obj->prepare($statement);
    $sth->execute || $self->throw("execute failed for supporting evidence get!");

    my %anahash;

    while (my $rowhash = $sth->fetchrow_hashref) {
	my $f1 = new Bio::EnsEMBL::SeqFeature;
	my $f2 = new Bio::EnsEMBL::SeqFeature;
	
	my $f = new Bio::EnsEMBL::FeaturePair(-feature1 => $f1,
					      -feature2 => $f2);

	my $exon = $rowhash->{exon};
#	$f1->seqname($rowhash->{contig});
	$f1->seqname("Supporting_feature");
	$f1->start  ($rowhash->{seq_start});
	$f1->end    ($rowhash->{seq_end});
	$f1->strand ($rowhash->{strand});
	$f1->source_tag($rowhash->{name});
	$f1->primary_tag('similarity');
	$f1->score  ($rowhash->{score});
	
	$f2->seqname($rowhash->{hid});
	$f2->start  ($rowhash->{hstart});
	$f2->end    ($rowhash->{hend});
	$f2->strand ($rowhash->{strand});
	$f2->source_tag($rowhash->{name});
	$f2->primary_tag('similarity');
	$f2->score  ($rowhash->{score});

	my $analysisid = $rowhash->{analysis};
	
	if ($anahash{$analysisid}) {
	    $f->analysis($anahash{$analysisid});

	} else {
	    my $feature_obj=Bio::EnsEMBL::DBSQL::Feature_Obj->new($self->_db_obj);
	    $f->analysis($feature_obj->get_Analysis($analysisid));

	    $anahash{$analysisid} = $f->analysis;
	}
	
	$f->validate;

	$exhash{$exon}->add_Supporting_Feature($f);
    }

}

=head2 get_Transcript
    
 Title   : get_Transcript
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut
    
sub get_Transcript{
    my ($self,$transid) = @_;

    my $seen = 0;
    my $trans = Bio::EnsEMBL::Transcript->new();

    my $sth = $self->_db_obj->prepare("select exon from exon_transcript where transcript = '$transid'");
    my $res = $sth->execute();

    while( my $rowhash = $sth->fetchrow_hashref) {
	my $exon = $self->get_Exon($rowhash->{'exon'});
	$trans->add_Exon($exon);
	$seen = 1;
    }

    if ($seen == 0 ) {
	$self->throw("transcript $transid is not present in db");
    }
    
    $trans->id($transid);

    return $trans;
}

=head2 get_Transcript_by_est
    
 Title   : get_Transcript_by_est
 Usage   : $db->get_Transcript_by_est($est_accession)
 Function: Gets a transcript object for a specific est id
 Example : 
 Returns : Bio::EnsEMBL::Transcript object
 Args    : est genbank id


=cut
    
sub get_Transcript_by_est{
    my ($self,$est_id) = @_;
    my @out;
    my $seen=0;
    $est_id || $self->throw("You need to provide the accession number of the est to get a transcript!\n");

    my $est = "gb|$est_id%";

    my $sth = $self->_db_obj->prepare("select distinct e_t.transcript from feature as f, exon as e,exon_transcript as e_t where f.hid like '".$est."' and e.seq_start<=f.seq_start and e.seq_end >= f.seq_end and e.contig = f.contig and e_t.exon = e.id;");
    my $res = $sth->execute();
    my $transcript;
    while( my $rowhash = $sth->fetchrow_hashref) {
	$transcript = $self->get_Transcript($rowhash->{'transcript'});
	push @out, $transcript;
	$seen = 1;
    }
    
    if ($seen == 0) {
	$self->throw("Could not get transcript for est $est!");
    }

    return @out;
}



=head2 get_Translation

 Title   : get_Translation
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_Translation{
   my ($self,$translation_id) = @_;

   my $sth     = $self->_db_obj->prepare("select version,seq_start,start_exon,seq_end,end_exon from translation where id = '$translation_id'");
   my $res     = $sth->execute();
   my $rowhash = $sth->fetchrow_hashref;

   if( !defined $rowhash ) {
       $self->throw("no translation of $translation_id");
   }

   my $out = Bio::EnsEMBL::Translation->new();

   $out->version      ($rowhash->{'version'});
   $out->start        ($rowhash->{'seq_start'});
   $out->end          ($rowhash->{'seq_end'});
   $out->start_exon_id($rowhash->{'start_exon'});
   $out->end_exon_id  ($rowhash->{'end_exon'});
   $out->id           ($translation_id);

   return $out;
}

=head2 get_Virtual_Contig
    
 Title   : get_Virtual_Contig
 Usage   : $gene_obj->get_Virtual_Contig($transcript,$max_length)
 Function: Gets a Bio::EnsEMBL::DB::Virtual Contig object which 
           spans the whole sequence on which the given 
           Bio::EnsEMBL::Transcript object lies, as long 
           as its length does not exceed max_length. If max_length
           is exceeded, undef is returned instead.
 Example : $gene_obj->get_Virtual_Contig($transcript,50000)
 Returns : VirtualContig Object (or undef)
 Args    : Bio::EnsEMBL::Transcript object and max_length int variable


=cut
    
sub get_Virtual_Contig{
    my ($self,$trans_id, $max_length) = @_;

    if ( !defined $trans_id) {
	$self->throw("Cannot create Virtual Contig without a transcript id!");
    }

    if (!defined $max_length) {
	$max_length=500000;
    }

    my $transcript=$self->get_Transcript($trans_id);

    #First of all, get out first exon, and create a 10000 bp Virtual Contig
    #starting from the contig on which the first exon is lying

    my $first_exon=$transcript->start_exon();
    my $first_contig=$self->_db_obj->get_Contig($first_exon->contig_id());
    my $first_ori=$first_exon->strand();
    my $vc=Bio::EnsEMBL::DB::VirtualContig->new( -focuscontig => $first_contig,
					      -focusposition => $first_exon->start,
					      -ori => $first_ori,
					      -left => 5000,
					      -right => 5000
					      );
    
    #Now get the last exon, and extend the virtual contig until the last exon is also
    #contained in the vc

    my $last_exon=$transcript->end_exon();
    my $last_contig=$self->_db_obj->get_Contig($last_exon->contig_id());
    my $last_ori=$last_exon->strand();
    my $not_finished=1;
    my $length;

    EXTEND:while (){
	foreach my $vcraw ($vc->rawcontig_ids()) {
	    if ($last_contig->id eq $vcraw) {		
		#extend full length of contig not 5000! 
		$vc=$vc->extend(-$last_contig->length,$last_contig->length);
		
		last EXTEND;
	    }
	}

	$vc=$vc->extend(-5000,5000);


	#If virtual contig longer than max. length, exit and return undef
	if ($vc->length > $max_length) {
	    print STDERR "Hit max. length!\n";
	    return undef;
	}
	# to cope with badly built genes
	if ($length==$vc->length){last EXTEND;}
	$length=$vc->length;
	
    }

   

    #Check that the Virtual Contig contains all exons of this transcript
    my $ok=undef;
    my $oldc;
    my $old_e_cont;
    my @notin;

    foreach my $exon ($transcript->each_Exon) {

	my $contig_id=$exon->contig_id();
	if ($contig_id ne $old_e_cont) {
	    foreach my $vcraw ($vc->rawcontig_ids) {
		if ($contig_id eq $vcraw) {
		    $ok=1;
		}
	    }

	    #If contigs found not contained in virtual contig, push them in an array
	    if (!$ok) {
		if ($oldc ne $contig_id){
		    push @notin,$contig_id;
		}
		$oldc=$contig_id;
	    }
	}
	$ok=undef;
	$old_e_cont=$contig_id;
    }
    
    #If missing contigs are found, extend again
    my $n=@notin+0;
    foreach my $notin (@notin) {
      EXTEND_MORE:while (){
	  foreach my $vcraw ($vc->rawcontig_ids()) {
	      if ($notin eq $vcraw) {
		  last EXTEND_MORE;
	      }
	  }
	  $vc=$vc->extend(-5000,5000);

	  #If virtual contig longer than max. length, exit and return undef
	  if ($vc->length > $max_length) {
	      print STDERR "Hit max. length!\n";
	      return undef;
	  }
	  # to cope with badly built genes
	  if ($length==$vc->length){last EXTEND_MORE;}
	  $length=$vc->length;


      }
    }
   
    return $vc;
}




=head2 get_Transcript_in_VC_coordinates
    
 Title   : get_Transcript_in_VC_coordinates
 Usage   : $gene_obj->get_Transcript_in_VC_coordinates($transcript_id)
 Function: Gets a Bio::EnsEMBL::Transcript object in vc coordinates
 Example : $gene_obj->get_Virtual_Contig($transcript_id)
 Returns : Bio::EnsEMBL::Transcript
 Args    : transcript id


=cut




sub get_Transcript_in_VC_coordinates
{
    
    my ($self,$transcript_id)=@_; 
 
    $transcript_id || $self->throw("need a transcript id");   
    my $vc=$self->get_Virtual_Contig($transcript_id);
    
    my $found;
  GENE: foreach my $gene ($vc->get_all_Genes){
      foreach my $transcript($gene->each_Transcript){
	
	  #print STDERR $transcript->id,"\n";

	  if ($transcript->id eq $transcript_id){$found=$transcript;last GENE;}
      }
  }   
    if (!defined $found ){$self->throw("Created VC around wrong transcript, sorry !")};
    return $found;    
    
}

=head2 write

 Title   : write
 Usage   : $Gene_obj->write_Gene($gene)
 Function: writes a particular gene into the database
 Example :
 Returns : nothing
 Args    : $gene object


=cut

sub write{
   my ($self,$gene) = @_;
   my $old_gene;
   my %done;
   
   if ( !defined $gene || ! $gene->isa('Bio::EnsEMBL::Gene') ) {
       $self->throw("$gene is not a EnsEMBL gene - not writing!");
   }

   # get out unique contig ids from gene to check against
   # database.

   my %contighash;


   foreach my $contig_id ( $gene->unique_contig_ids() ) {
       eval {
	   #print STDERR "Getting out contig for $contig_id\n";
	   my $contig      = $self->_db_obj->get_Contig($contig_id);
	   $contig->fetch();
	   
	   $contighash{$contig_id} = $contig;

	   # if there is no exception then it is there. Get rid of it
	   $contig = 0;
       };
       if( $@ ) {
	   $self->throw("In trying to write gene " . $gene->id(). " into the database, unable to find contig $contig_id. Aborting write\n\nFull Exception\n\n$@\n");
	   # done before locks, so we are ok.
       }
       
   }

   # gene is big daddy object

   foreach my $trans ( $gene->each_Transcript() ) {
       $self->write_Transcript($trans,$gene);
       my $c = 1;
       foreach my $exon ( $trans->each_Exon() ) {

	   my $sth = $self->_db_obj->prepare("insert into exon_transcript (exon,transcript,rank) values ('". $exon->id()."','".$trans->id()."',".$c.")");
	   $sth->execute();
	   $c++;


	   if( $done{$exon->id()} ) { 
	       next; 
	   }
	   $done{$exon->id()} = 1;

	   if( $exon->isa('Bio::EnsEMBL::StickyExon') ) {
	       $self->write_StickyExon($exon); 
	   } else {
	       $self->write_Exon($exon);
	   }

       }
   }

   !$gene->created() && $gene->created(0);
   !$gene->modified() && $gene->modified(0);
 
   my $sth2 = $self->_db_obj->prepare("insert into gene (id,version,created,modified,stored) values ('". 
			     $gene->id       . "','".
			     $gene->version  . "',FROM_UNIXTIME(".
			     $gene->created  . "),FROM_UNIXTIME(".
			     $gene->modified . "),now())");
   $sth2->execute();

   foreach my $dbl ( $gene->each_DBLink ) {
       my $sth3 = $self->_db_obj->prepare("insert into genedblink (gene_id,external_id,external_db) values ('". 
			     $gene->id        . "','".
			     $dbl->primary_id . "','".
			     $dbl->database   . "')");
       $sth3->execute();
   }
    
   my $id=$gene->id;
   my $type=$gene->type;

   my $sth4 = $self->_db_obj->prepare("insert into genetype (gene_id,type) values ('$id','$type')");

   $sth4->execute();





   
   return 1;
}

=head2 write_Exon

 Title   : write_Exon
 Usage   : $obj->write_Exon($exon)
 Function: writes a particular exon into the database
 Example :
 Returns : 
 Args    :


=cut

sub write_Exon {
    my ($self,$exon,$no_supporting) = @_;
    my $old_exon;
    
    if( ! $exon->isa('Bio::EnsEMBL::Exon') ) {
	$self->throw("$exon is not a EnsEMBL exon - not dumping!");
    }
    $exon->id() || $self->throw("Missing exon id");
    $exon->version() || $self->throw("Missing exon version number"); 
    $exon->contig_id() || $self->throw("Missing exon contig id");
    $exon->start || $self->throw("Missing exon start position"); 
    $exon->end || $self->throw("Missing exon end position");
    $exon->created || $self->throw("Missing exon created time");
    $exon->modified || $self->throw("Missing exon modified time");

    # got to convert contig_id to internal_id

    my $contig = $self->_db_obj->get_Contig($exon->contig_id);
    
    if( $exon->start > $exon->end ) {
	$self->throw("Start is greater than end for exon. Not writing it!");
    }

    my $exonst = q{
        insert into exon (id, version, contig, created, modified
          , seq_start, seq_end, strand, phase, stored, end_phase, sticky_rank) 
        values (?,?,?,FROM_UNIXTIME(?),FROM_UNIXTIME(?),?,?,?,?,NOW(),?,?)
        };
    
    my $sth = $self->_db_obj->prepare($exonst);
    $sth->execute(
        $exon->id(),
        $exon->version(),
        $contig->internal_id,
        $exon->created(),
        $exon->modified(),
        $exon->start,
        $exon->end,
        $exon->strand,
        $exon->phase,
        $exon->end_phase,
	$exon->sticky_rank
        );
    
    # Now the supporting evidence
    
    if( !defined $no_supporting || !$no_supporting ) {
	$self->write_supporting_evidence($exon);
    }

    return 1;
}


=head2 write_StickyExon

 Title   : write_StickyExon
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub write_StickyExon{
   my ($self,$exon) = @_;

   if( ! $exon->isa('Bio::EnsEMBL::StickyExon') ) {
       $self->throw("$exon is not a EnsEMBL exon - not dumping!");
   }

   foreach my $e ( $exon->each_component_Exon() ) {
       $self->write_Exon($e,1);
   }
   my ($f) = $exon->each_component_Exon();
   $self->write_supporting_evidence($f);
}



=head2 write_supporting_evidence

 Title   : write_supporting_evidence
 Usage   : $obj->write_supporting_evidence
 Function: Writes supporting evidence features to the database
 Example :
 Returns : nothing
 Args    : None


=cut

sub write_supporting_evidence {
    my ($self,$exon) = @_;

    $self->throw("Argument must be Bio::EnsEMBL::Exon. You entered [$exon]\n") unless $exon->isa("Bio::EnsEMBL::Exon");


    my $string;
    if( $self->use_delayed_insert == 1 ) {
	$string = 'DELAYED';
    } else {
	$string = '';
    }

    #$string = '';

    my $sth  = $self->_db_obj->prepare("insert $string into supporting_feature(id,exon,seq_start,seq_end,score,strand,analysis,name,hstart,hend,hid) values(?,?,?,?,?,?,?,?,?,?,?)");
    
    FEATURE: foreach my $f ($exon->each_Supporting_Feature) {

	eval {
	    $f->validate();
	};

	if ($@) {
	    print(STDERR "Supporting feature invalid. Skipping feature\n");
	    next FEATURE;
	}
	my $feature_obj=Bio::EnsEMBL::DBSQL::Feature_Obj->new($self->_db_obj);
  	my $analysisid = $feature_obj->write_Analysis($f->analysis);
	
	if ($f->isa("Bio::EnsEMBL::FeaturePairI")) {
	    $sth->execute('NULL',
			  $exon->id,
			  $f->start,
			  $f->end,
			  $f->score,
			  $f->strand,
			  $analysisid,
			  $f->source_tag,
			  $f->hstart,
			  $f->hend,
			  $f->hseqname
			  );
	} else {
	    $self->warn("Feature is not a Bio::EnsEMBL::FeaturePair");
	}
    }
}

=head2 write_Transcript

 Title   : write_Transcript
 Usage   : $obj->write_Transcript($trans,$gene)
 Function: writes a particular transcript *but not the exons* into
           the database
 Example :
 Returns : 
 Args    :


=cut

sub write_Transcript{
   my ($self,$trans,$gene) = @_;
   my $old_trans;

   if( ! $trans->isa('Bio::EnsEMBL::Transcript') ) {
       $self->throw("$trans is not a EnsEMBL transcript - not dumping!");
   }

   if( ! $gene->isa('Bio::EnsEMBL::Gene') ) {
       $self->throw("$gene is not a EnsEMBL gene - not dumping!");
   }

   $trans->id || $self->throw("You need a transcript id to write a transcript!");
   $gene->id || $self->throw("You need a gene id to write a transcript!");
   $trans->translation->id || $self->throw("You need a translation id to write a transcript!");
   $trans->version || $self->throw("You need a transcript version number to write a transcript!");

   # ok - now load this line in
   my $tst = $self->_db_obj->prepare("
        insert into transcript (id, gene, translation, version) 
        values (?, ?, ?, ?)
        ");
                
   $tst->execute(
        $trans->id,
        $gene->id, 
        $trans->translation->id,
        $trans->version   
        );

   #print STDERR "Going to look at gene links\n";

   foreach my $dbl ( $trans->each_DBLink ) {
       #print STDERR "Going to insert for",$trans->id," ",$dbl->primary_id," ",$dbl->database,"\n";
       my $sth3 = $self->_db_obj->prepare("insert into transcriptdblink (transcript_id,external_id,external_db) values ('". 
					  $trans->id        . "','".
					  $dbl->primary_id . "','".
					  $dbl->database   . "')");
       $sth3->execute();
       
   }

   $self->write_Translation($trans->translation());
   return 1;
}

=head2 write_Translation

 Title   : write_Translation
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub write_Translation{
    my ($self,$translation) = @_;
    my $old_transl;
    
    if( !$translation->isa('Bio::EnsEMBL::Translation') ) {
	$self->throw("Is not a translation. Cannot write!");
    }
    
    if ( !defined $translation->version  ) {
	$self->throw("No version number on translation");
    }
    
    my $tst = $self->_db_obj->prepare("insert into translation (id,version,seq_start,start_exon,seq_end,end_exon) values ('" 
			     . $translation->id . "',"
			     . $translation->version . ","
			     . $translation->start . ",'"  
			     . $translation->start_exon_id. "',"
			     . $translation->end . ",'"
			     . $translation->end_exon_id . "')");
    $tst->execute();
    return 1;
}




sub get_NewId {
    my ($self,$table,$stub) = @_;

    $table || $self->throw("Need to provide a table name to get a new id!\n");
    
    my $query = "select max(id) from $table";

    my $sth   = $self->_db_obj->prepare($query);
    $sth->execute;
    my ($id)   = $sth->fetchrow;

    #print(STDERR "max id is '$id'\n");
    
    if (!defined $id || $id eq "") {
	$id = $stub . "00000000000";
    }
    
    #print(STDERR "max id is '$id'\n");

    if ($id =~ /$stub(\d+)$/) {
	my $newid  = $1;
	$newid++;
	
	if (length($newid) > 11) {
	    if ($newid =~ /^0/) {
		$newid =~ s/^0//;
	    } else {
		$self->throw("Can't truncate number string to generate new id [$newid]");
	    }
	}
	$newid = $stub . $newid;
	#print STDERR ("New id is $newid\n");
	return $newid;
    } else {
	$self->throw("[$id] does not look like an object id (e.g. ENST00000019784)");
    }
    

}



=head2 get_new_GeneID

 Title   : get_new_GeneID
 Usage   : my $id = $geneobj->get_new_GeneID
 Function: 
 Example : 
 Returns : Gets the next unused gene id from the database
 Args    : none


=cut

sub get_new_GeneID {
    my ($self,$stub) = @_;

    $stub = "ENSG" unless defined($stub);

    return $self->get_NewId("gene",$stub);
    
}

=head2 get_new_TranscriptID

 Title   : get_new_TranscriptID
 Usage   : my $id = $geneobj->get_new_TranscriptID
 Function: 
 Example : 
 Returns : Gets the next unused transcript id from the database
 Args    : none


=cut

sub get_new_TranscriptID {
    my ($self,$stub) = @_;

    $stub = "ENST" unless defined($stub);

    return $self->get_NewId("transcript",$stub);

}

=head2 get_new_ExonID

 Title   : get_new_ExonID
 Usage   : my $id = $geneobj->get_new_ExonID
 Function: 
 Example : 
 Returns : Gets the next unused exon id from the database
 Args    : none


=cut

sub get_new_ExonID {
    my ($self,$stub) = @_;

    $stub = "ENSE" unless defined($stub);

    return $self->get_NewId("exon",$stub);
}


sub get_new_TranslationID {
    my ($self,$stub) = @_;

    $stub  = "ENSP" unless defined($stub);

    return $self->get_NewId("translation",$stub);

}

=head2 _db_obj

 Title   : _db_obj
 Usage   : $obj->_db_obj($newval)
 Function: 
 Example : 
 Returns : value of _db_obj
 Args    : newvalue (optional)


=cut

sub _db_obj{
   my ($self,$value) = @_;
   if( defined $value) {
      $self->{'_db_obj'} = $value;
    }
    return $self->{'_db_obj'};

}

=head2 use_delayed_insert

 Title   : use_delayed_insert
 Usage   : $obj->use_delayed_insert($newval)
 Function: 
 Returns : value of use_delayed_insert
 Args    : newvalue (optional)


=cut

sub use_delayed_insert{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'use_delayed_insert'} = $value;
    }
    return $obj->{'use_delayed_insert'};

}
