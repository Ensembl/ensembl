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

This is one of the objects contained in Bio:EnsEMBL::DBSQL::Obj, dealing with
Gene methods, such as writing and getting genes, transcripts, translations, and exons.

The Obj object represents a database that is implemented somehow (you shouldn\'t
care much as long as you can get the object).

=head1 CONTACT

Elia Stupka: elia@ebi.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are 
usually preceded with a _

=cut


# Let the code begin...

package Bio::EnsEMBL::DBSQL::Gene_Obj;

use vars qw(@ISA);
use strict;

# Object preamble - inheriets from Bio::Root::Object

use Bio::Root::Object;
use Bio::EnsEMBL::DBSQL::Obj;

use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::Exon;
use Bio::EnsEMBL::Transcript;
use DBI;

use Bio::EnsEMBL::DBSQL::DummyStatement;

@ISA = qw(Bio::Root::Object);
# new() is inherited from Bio::Root::Object

# _initialize is where the heavy stuff will happen when new is called

sub _initialize {
  my($self,$db_obj) = @_;

  my $make = $self->SUPER::_initialize;
  
  $db_obj || $self->throw("Database Gene object must be passed a db obj!");
  $self->_db_obj($db_obj);

  #
  # This needs to be rethought. We are caching sequences
  # here to allow multiple exons to be retrieved fine
  # And now more cache's. I think cache's might be a fact of life...
  # 

  $self->{'_contig_seq_cache'} = {};

  return $make; # success - we hope!

}

=head2 each_cloneid

 Title   : each_cloneid
 Usage   : @each_cloneid = $gene_obj->each_cloneid($geneid);
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub each_cloneid{
   my ($self,$geneid) = @_;

   my $sth = $self->_db_obj->prepare("select clone from geneclone_neighbourhood where gene = '$geneid'");

   my @out;

   $sth->execute;
   while( my $rowhash = $sth->fetchrow_hashref) {
       push(@out,$rowhash->{'clone'});
   }
   return @out;
}

=head2 delete

 Title   : delete
 Usage   : $Gene_Obj->delete_Gene($gene_id)
 Function: deletes a gene from the database, i.e. exons, transcripts, translations
 Example : $geneobj->delete_Gene('ENSG00000019482')
 Returns : nothing
 Args    : gene id

=cut

sub delete{
   my ($self,$geneid) = @_;
   my @trans;
   my %exon;
   my $translation;
   # get out exons, transcripts for gene. 

   my $sth = $self->_db_obj->prepare("select id,translation from transcript where gene = '$geneid'");
   $sth->execute;
   while( my $rowhash = $sth->fetchrow_hashref) {
       push(@trans,$rowhash->{'id'});
       $translation = $rowhash->{'translation'};
   }

   foreach my $trans ( @trans ) {
       my $sth = $self->_db_obj->prepare("select exon from exon_transcript where transcript = '$trans'");
       $sth->execute;
       while( my $rowhash = $sth->fetchrow_hashref) {
	   $exon{$rowhash->{'exon'}} =1;
       }
   }

   my $sth2 = $self->_db_obj->prepare("delete from translation where id = '$translation'");
   $sth2->execute;

   # delete exons, transcripts, gene rows

   foreach my $exon ( keys %exon ) {
       my $sth = $self->_db_obj->prepare("delete from exon where id = '$exon'");
       $sth->execute;

       $sth = $self->_db_obj->prepare("delete from supporting_feature where exon = '$exon'");
       $sth->execute;
   }

   foreach my $trans ( @trans ) {
       my $sth= $self->_db_obj->prepare("delete from transcript where id = '$trans'");
       $sth->execute;
       $sth= $self->_db_obj->prepare("delete from exon_transcript where transcript = '$trans'");
       $sth->execute;
   }

   $sth = $self->_db_obj->prepare("delete from gene where id = '$geneid'");
   $sth->execute;
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
    my $sth = $self->_db_obj->prepare("delete from exon_transcript where transcript = '".$exon_id."'");
    my $res = $sth ->execute;

    #Delete exon rows
    $sth = $self->_db_obj->prepare("delete from exon where id = '".$exon_id."'");
    $res = $sth->execute;

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

   if ($supporting && $supporting eq 'evidence') {
       @out = $self->get_array_supporting('evidence',$geneid);
   }   else {
       @out = $self->get_array_supporting('without',$geneid);
   }
   return $out[0];
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
   my $res = $sth->execute;

   while( my $rowhash = $sth->fetchrow_hashref) {
       push(@out,$rowhash->{'id'});
   }

   return @out;
}

=head2 get_array_supporting

    Title   : get_Gene_array_supporting
    Usage   : $obj->get_Gene_array_supporting($supporting,@geneid)
    Function: Gets an array of genes, with transcripts and exons. If $supporting
           equal to 'evidence' the supporting evidence for each exon is also read
    from the supporting evidence table
    Example : $obj->get_Gene_array_supporting ('evidence',@geneid)
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
    
    my $inlist = join(',',map "'$_'", @geneid);
    $inlist = "($inlist)";
    
    # I know this SQL statement is silly.
    #
    
    
    my $query =              "select p3.gene," .
	                             "p4.id," .
			             "p3.id," .
			             "p1.exon," .
			             "p1.rank," .
			             "p2.seq_start,p2.seq_end," .
			             "UNIX_TIMESTAMP(p2.created),UNIX_TIMESTAMP(p2.modified)," .
			             "p2.strand,p2.phase," .
            			     "p5.seq_start,p5.start_exon,p5.seq_end,p5.end_exon,p5.id," .
			             "p6.version,p3.version,p2.version,p5.version,p4.clone " .
			     "from   gene as p6," .
			             "contig as p4, " .
			             "transcript as p3, " .
			             "exon_transcript as p1, " .
			             "exon as p2," .
			             "translation as p5," .
			             "geneclone_neighbourhood as p7  " .
			     "where  p6.id in $inlist " .
			     "and    p3.gene = p6.id " .
			     "and    p4.clone = p7.clone " .
			     "and    p7.gene = p6.id  " .
			     "and    p2.contig = p4.internal_id " .
			     "and    p1.exon = p2.id " .
			     "and    p3.id = p1.transcript " .
			     "and    p5.id = p3.translation " .
			     "order by p3.gene,p3.id,p1.rank";
    
    my $sth = $self->_db_obj->prepare($query);
    my $res = $sth ->execute();
    
    my $current_gene_id       = '';
    my $current_transcript_id = '';
    
    my ($gene,$trans);
    
    while( (my $arr = $sth->fetchrow_arrayref()) ) {
	
	my ($geneid,$contigid,$transcriptid,$exonid,$rank,$start,$end,
	    $exoncreated,$exonmodified,$strand,$phase,$trans_start,
	    $trans_exon_start,$trans_end,$trans_exon_end,$translationid,
	    $geneversion,$transcriptversion,$exonversion,$translationversion,$cloneid) = @{$arr};
	
	if( ! defined $phase ) {
	    $self->throw("Bad internal error! Have not got all the elements in gene array retrieval");
	}
	
	# Create new gene if the id has changed
	if( $geneid ne $current_gene_id ) {
	    
	    if( $transcriptid eq $current_transcript_id ) {
		$self->throw("Bad internal error. Switching genes without switching transcripts");
	    } 
	    
	    $gene = Bio::EnsEMBL::Gene->new();
	    
	    $gene->id                       ($geneid);
	    $gene->version                  ($geneversion);
	    $gene->add_cloneid_neighbourhood($cloneid);
	    
	    $current_gene_id = $geneid;
	    push(@out,$gene);
	    
	}
	
	# Create new transcript if the id has changed
	if( $transcriptid ne $current_transcript_id ) {

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
	$exon->strand   ($strand);
	$exon->phase    ($phase);
	$exon->version  ($exonversion);
	$exon->seqname  ($contigid);
	
	#
	# Attach the sequence, cached if necessary...
	#
	if ($supporting && $supporting eq 'evidence') {
	    push @sup_exons, $exon;
	}
	
	my $seq;
	
	if( $self->_contig_seq_cache($exon->contig_id) ) {
	    $seq = $self->_contig_seq_cache($exon->contig_id);
	} else {
	    my $contig = $self->_db_obj->get_Contig($exon->contig_id());
	    $seq = $contig->primary_seq();
	    $self->_contig_seq_cache($exon->contig_id,$seq);
	}
	
	$exon ->attach_seq($seq);
	$trans->add_Exon($exon);

    }
    
    if ($supporting && $supporting eq 'evidence') {
	$self->get_supporting_evidence(@sup_exons);
    }
    
    return @out;
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
      $exon->version  ($rowhash->{'version'});

   my $contig_id = $exon->contig_id();

   # we have to make another trip to the database to get out the contig to clone mapping.
   my $sth2     = $self->_db_obj->prepare("select clone from contig where internal_id = '$contig_id'");
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

   if( $self->_contig_seq_cache($exon->contig_id) ) {
       $seq = $self->_contig_seq_cache($exon->contig_id);
   } else {
       my $contig = $self->_db_obj->get_Contig($exon->contig_id());
       $seq = $contig->primary_seq();
       $self->_contig_seq_cache($exon->contig_id,$seq);
   }

   $exon->attach_seq($seq);

   return $exon;
}

=head2 get_supporting_evidence

 Title   : get_supporting_evidence
 Usage   : $obj->get_supporting_evidence
 Function: Writes supporting evidence features to the database
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
   
    my $sth = $self->_db_obj->prepare("select * from supporting_feature where exon in (" . $instring . ")");
    $sth->execute;

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
	    $f->analysis($self->_db_obj->get_Analysis($analysisid));
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
	   my $contig = $self->_db_obj->get_Contig($contig_id);

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
       $self->_db_obj->write_Transcript($trans,$gene);
       my $c = 1;
       foreach my $exon ( $trans->each_Exon() ) {
	   my $sth = $self->_db_obj->prepare("insert into exon_transcript (exon,transcript,rank) values ('". $exon->id()."','".$trans->id()."',".$c.")");
	   $sth->execute();
	   if( $done{$exon->id()} ) { next; }
	   $done{$exon->id()} = 1;

	   my $internal_contig_id = $contighash{$exon->contig_id}->internal_id;

	   if (!defined($internal_contig_id)) {
	       $self->throw("Internal id not found for contig [" . $exon->contig_id . "]");
	   }
	   my $tmpid = $exon->contig_id;
	   $exon->contig_id($internal_contig_id);
	   $self->_db_obj->write_Exon($exon);
	   $exon->contig_id($tmpid);
	   
	   $c++;
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

   foreach my $cloneid ($gene->each_cloneid_neighbourhood) {

       my $sth = $self->_db_obj->prepare("select gene,clone from geneclone_neighbourhood where gene='".$gene->id."' && clone='$cloneid'");
       $sth->execute();
       my $rowhash =  $sth->fetchrow_arrayref();
       my  $rv = $sth->rows;
       if( ! $rv ) {
	   $sth = $self->_db_obj->prepare("insert into geneclone_neighbourhood (gene,clone) values ('" . 
				 $gene->id . "','". 
				 $cloneid ."')");
	   $sth->execute();
       }
   }
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
    my ($self,$exon) = @_;
    my $old_exon;
    
    if( ! $exon->isa('Bio::EnsEMBL::Exon') ) {
	$self->throw("$exon is not a EnsEMBL exon - not dumping!");
    }
    
    my $exonst = "insert into exon (id,version,contig,created,modified,seq_start,seq_end,strand,phase,stored,end_phase) values ('" .
	
	$exon->id()        . "'," .
	$exon->version()   . ",'".
	$exon->contig_id() . "', FROM_UNIXTIME(" .
	$exon->created()   . "), FROM_UNIXTIME(" .
	$exon->modified()  . ")," .
	$exon->start       . ",".
	$exon->end         . ",".
	$exon->strand      . ",".
	$exon->phase       . ",now(),".
        $exon->end_phase . ")";
    
    my $sth = $self->_db_obj->prepare($exonst);
    $sth->execute();
    
    # Now the supporting evidence
    
    $self->write_supporting_evidence($exon);
    return 1;
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

    my $sth  = $self->_db_obj->prepare("insert into supporting_feature(id,exon,seq_start,seq_end,score,strand,analysis,name,hstart,hend,hid) values(?,?,?,?,?,?,?,?,?,?,?)");
    
    FEATURE: foreach my $f ($exon->each_Supporting_Feature) {

	eval {
	    $f->validate();
	};

	if ($@) {
	    print(STDERR "Supporting feature invalid. Skipping feature\n");
	    next FEATURE;
	}

	my $analysisid = $self->_db_obj->write_Analysis($f->analysis);
	
	if ($f->isa("Bio::EnsEMBL::FeaturePair")) {
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
	    #$self->warn("Feature is not a Bio::EnsEMBL::FeaturePair");
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

   # ok - now load this line in
   my $tst = $self->_db_obj->prepare("insert into transcript (id,gene,translation,version) values ('" . $trans->id . "','" . $gene->id . "','" . $trans->translation->id() . "',".$trans->version.")");
   $tst->execute();
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

=head2 _contig_seq_cache

 Title   : _contig_seq_cache
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub _contig_seq_cache{
   my ($self,$id,$seq) = @_;

   if( $seq ) {
       
       #
       # Every 100 hits, flush the cache
       #
       if( $self->{'_contig_seq_cnt'} > 100 ) {
	   $self->_flush_seq_cache;
	   $self->{'_contig_seq_cnt'} = 0;
       }

       $self->{'_contig_seq_cnt'}++;
       $self->{'_contig_seq_cache'}->{$id} = $seq;
   }

   return $self->{'_contig_seq_cache'}->{$id};
}
