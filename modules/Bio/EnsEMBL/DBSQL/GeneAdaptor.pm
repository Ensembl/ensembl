# EnsEMBL Gene reading writing adaptor for mySQL
#
# Copyright EMBL-EBI 2001
#
# Author: Arne Stabenau
# based on 
# Elia Stupkas Gene_Obj
# 
# Date : 20.02.2001
#

=head1 NAME

Bio::EnsEMBL::DBSQL::GeneAdaptor - MySQL Database queries to generate and store gens.

=head1 SYNOPSIS

=head1 CONTACT

  Arne Stabenau: stabenau@ebi.ac.uk
  Elia Stupka  : elia@ebi.ac.uk
  Ewan Birney  : birney@ebi.ac.uk

=head1 APPENDIX

=cut

;

package Bio::EnsEMBL::DBSQL::GeneAdaptor;

use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::Exon;
use Bio::EnsEMBL::Transcript;

=head2 remove_by_dbID

 Title   : remove_by_dbID
 Usage   : $geneAdaptor->remove_by_dbID($gene_id)
 Function: deletes a gene from the database, i.e. exons, transcripts, translations
 Example : $geneAdaptor->remove_by_dbID('ENSG00000019482')
 Returns : nothing
 Args    : gene id

=cut

sub remove_by_dbID {
   my ($self,$geneid) = @_;
   my @trans;
   my %exon;
   my @translation;
   # get out exons, transcripts for gene. 

   my $sth = $self->prepare("select id,translation from transcript where gene = '$geneid'");
   $sth->execute;
   while( my $rowhash = $sth->fetchrow_hashref) {
       push(@trans,$rowhash->{'id'});
       push(@translation,$rowhash->{'translation'});
   }

   foreach my $trans ( @trans ) {
       my $sth = $self->prepare("select exon from exon_transcript where transcript = '$trans'");
       $sth->execute;
       while( my $rowhash = $sth->fetchrow_hashref) {
	   $exon{$rowhash->{'exon'}} =1;
       }
   }

   foreach my $translation (@translation) {
       my $sth2 = $self->prepare("delete from translation where id = '$translation'");
       $sth2->execute;
   }
   # delete exons, transcripts, gene rows

   foreach my $exon ( keys %exon ) {
       my $sth = $self->prepare("delete from exon where id = '$exon'");
       $sth->execute;

       $sth = $self->prepare("delete from supporting_feature where exon = '$exon'");
       $sth->execute;
   }

   foreach my $trans ( @trans ) {
       my $sth= $self->prepare("delete from transcript where id = '$trans'");
       $sth->execute;
       $sth= $self->prepare("delete from exon_transcript where transcript = '$trans'");
       $sth->execute;
   }

   $sth = $self->prepare("delete from gene where id = '$geneid'");
   $sth->execute;

   $sth = $self->prepare("delete from genetype where gene_id = '$geneid'");
   $sth->execute;
}   


=head2 list_geneIds

 Title   : list_geneIds
 Usage   : $geneAdaptor->list_geneIds
 Function: Gets an array of ids for all genes in the current db
 Example : 
 Returns : array of ids
 Args    : none

=cut

sub list_geneIds {
   my ($self) = @_;

   my @out;
   my $sth = $self->prepare("select id from gene");
   my $res = $sth->execute;

   while( my $rowhash = $sth->fetchrow_hashref) {
       push(@out,$rowhash->{'id'});
   }

   return @out;
}



=head2 list_geneIds_by_hids

 Title   : list_geneIds_by_hids
 Usage   : @geneids = $obj->list_geneIds_by_hids(@hids)
 Function: gives back geneids with these hids as supporting evidence
 Example :
 Returns : 
 Args    :


=cut

sub list_geneIds_by_hids{
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


=head2 fetch_by_dbID

 Title   : fetch_by_dbID
 Usage   : $geneobj->fetch_by_dbID( $geneid, $supporting )
 Function: gets one gene out of the db with or without supporting evidence
 Example : $obj->get('ENSG00000009151','evidence')
 Returns : gene object (with transcripts, exons and supp.evidence if wanted)
 Args    : gene id and supporting tag if latter not specified, assumes without
	   Note that it is much faster to get genes without supp.evidence!

=cut


sub fetch_by_dbID {
    my ($self,$geneId ) = @_;
    
     my $query = qq{
        SELECT tscript.gene
          , tscript.id
          , e_t.exon, e_t.rank
          , transl.seq_start, transl.start_exon
          , transl.seq_end, transl.end_exon
          , transl.id
          , gene.version
          , tscript.version
          , transl.version
	  , genetype.type
          , gene.analysisId
        FROM gene
          , transcript tscript
          , exon_transcript e_t
          , translation transl
          , genetype genetype
        WHERE gene.id = tscript.gene
          AND tscript.id = e_t.transcript
          AND tscript.translation = transl.id
          AND genetype.gene_id = gene.id
          AND gene.id = $geneId

        ORDER BY tscript.gene
          , tscript.id
          , e_t.rank
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




=head2 get_Gene_by_Transcript_id

 Title   : get_Gene_by_Transcript_id
 Usage   : $gene_obj->get_Gene_by_Transcript_id($transid, $supporting)
 Function: gets one gene out of the db with or without supporting evidence
 Returns : gene object (with transcripts, exons and supp.evidence if wanted)
 Args    : transcript id and supporting tag (if latter not specified,
assumes without
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
    my @genes=$self->get_Gene_array_by_DBLink($external_id,$supporting);

    my $biggest;
    my $max=0;
    my $size=scalar(@genes);
    if ($size > 0) {
	foreach my $gene (@genes) {
	    my $size = (scalar($gene->each_unique_Exon));
	    if ($size > $max) {
		$biggest = $gene;
		$max=$size;
	    }
	}
	return $biggest;
    }
    return;
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

sub get_Gene_array_by_DBLink {
    my $self = shift;
    my $external_id = shift;
    my $supporting = shift;

    my @genes;
    my $sth = $self->_db_obj->prepare("select gene_id from genedblink where external_id = '$external_id'");
    $sth->execute;
    my $seen=0;
    while (my ($geneid)=$sth->fetchrow_array()) {
	push (@genes,$self->get($geneid,$supporting));
	$seen=1;
    }
    if( !$seen ) {
	return;
    }
    else {
	return @genes;
    }
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
	    my $feature_obj=Bio::EnsEMBL::DBSQL::Feature_Obj->new($self->_db_obj);
	    $f->analysis($feature_obj->get_Analysis($analysisid));

	    $anahash{$analysisid} = $f->analysis;
	}
	
	$f->validate;

	$exhash{$exon}->add_Supporting_Feature($f);
    }
}

=head2 get_supporting_evidence_direct

 Title   : get_supporting_evidence_direct
 Usage   : $obj->get_supporting_evidence_driect
 Function: Gets supporting evidence features from the feature table
 Example :
 Returns : nothing
 Args    : array of exon objects, needed to know which exon to attach the evidence to


=cut

sub get_supporting_evidence_direct {
    my ($self,@exons) = @_;

    my %exhash;
    my %analhash;
    if (@exons == 0) {
	$self->throw("No exon objects were passed on!");
    }
    my $list = "";

    foreach my $exon (@exons) {
	$list .= "'".$exon->id."',";
	$exhash{$exon->id} = $exon;
    }
    $list =~ s/\,$//;

    my $sth2=$self->_db_obj->prepare("select f.seq_start,f.seq_end,f.score,f.strand,f.analysis,f.name,f.hstart,f.hend,f.hid,f.evalue,f.perc_id,e.id,c.id from feature f,exon e,contig c where c.internal_id = f.contig and f.contig = e.contig and e.id in ($list) and !(f.seq_end < e.seq_start or f.seq_start > e.seq_end) and f.strand = e.strand and f.analysis < 5");
    $sth2->execute;

    while (my $arrayref = $sth2->fetchrow_arrayref) {
	my ($start,$end,$f_score,$strand,$analysisid,$name,$hstart,$hend,$hid,$evalue,$perc_id,$exonid,$contig) = @{$arrayref};
	my $analysis;
	if (!$analhash{$analysisid}) {
	    my $feature_obj=Bio::EnsEMBL::DBSQL::Feature_Obj->new($self->_db_obj);
	    $analysis = $feature_obj->get_Analysis($analysisid);
	    $analhash{$analysisid} = $analysis;	   
	} 
	else {
	    $analysis = $analhash{$analysisid};
	}
	
	
	if( !defined $name ) {
	    $name = 'no_source';
	}
	
	my $out = Bio::EnsEMBL::FeatureFactory->new_feature_pair();   
	$out->set_all_fields($start,$end,$strand,$f_score,$name,'similarity',$contig,$hstart,$hend,1,$f_score,$name,'similarity',$hid);
	$out->analysis($analysis);

	#$out->validate();
	$exhash{$exonid}->add_Supporting_Feature($out);
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
    $sth->execute();

    while( my $rowhash = $sth->fetchrow_hashref) {
	my $exon = $self->get_Exon($rowhash->{'exon'});
	$trans->add_Exon($exon);
	$seen = 1;
    }
    $sth = $self->_db_obj->prepare("select version,translation from transcript where id = '$transid'");
    $sth->execute();

    while( my $rowhash = $sth->fetchrow_hashref) {
	my $translation = $self->get_Translation($rowhash->{'translation'});
	$trans->translation($translation);
	$trans->version($rowhash->{'version'});
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
    
    
}




=head2 store

 Title   : store
 Usage   : $geneAdaptor->store($gene)
 Function: writes a particular gene into the database
 Example :
 Returns : nothing
 Args    : $gene object


=cut

sub store {
   my ($self,$gene) = @_;
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
    
    my $exonst = q{
        insert into exon (id, version, contig, created, modified
          , seq_start, seq_end, strand, phase, stored, end_phase, rank) 
        values (?,?,?,FROM_UNIXTIME(?),FROM_UNIXTIME(?),?,?,?,?,NOW(),?,?)
        };
    
    my $contig = $self->_db_obj->get_Contig($exon->contig_id);

    #print STDERR $exon->id . " " . 
	#$exon->version . " " .
	#$contig->internal_id . " " . 
	#$exon->created . " " . 
	#$exon->modified . " " . 
	#$exon->start . " " . 
	#$exon->end . " " . 
	#$exon->strand . " " . 
	#$exon->phase . " " . 
	#$exon->end_phase . " " . 
	#$exon->sticky_rank . "\n";

    my $sth = $self->_db_obj->prepare($exonst);
    $sth->execute(
        $exon->id(),
        $exon->version(),
        $contig->internal_id(),
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
    #print STDERR "Writing supporting evidence for exon ".$exon->id."\n";
    my $sth  = $self->_db_obj->prepare("insert into supporting_feature(id,exon,seq_start,seq_end,score,strand,analysis,name,hstart,hend,hid) values(?,?,?,?,?,?,?,?,?,?,?)");
    
  FEATURE: foreach my $f ($exon->each_Supporting_Feature) {
	#print STDERR "Writing supporting feature ".$f->source_tag."\n";
	eval {
	    $f->validate();
	};

	if ($@) {
	    print(STDERR "Supporting feature invalid. Skipping feature\n");
	    next FEATURE;
	}
	my $feature_obj=Bio::EnsEMBL::DBSQL::Feature_Obj->new($self->_db_obj);
  	my $analysisid = $feature_obj->write_Analysis($f->analysis);
	
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



sub get_New_external_id {
    my ($self,$table,$stub,$number) = @_;

    $table .= "_external";
    if( !defined $number ) {
	$number = 1;
    }

    my @out;


    my $lsth   = $self->_db_obj->prepare("lock table $table write");
    $lsth->execute;

    # wrap critical region in an eval so we can catch errors and release table

    eval {

	my $query = "select max(external_id) as id from $table where external_id like '$stub%'";
	
	my $sth   = $self->_db_obj->prepare($query);
	my $res   = $sth->execute;
	my $row   = $sth->fetchrow_hashref;
	my $id    = $row->{id};
	
	if (!defined($id) || $id eq "") {
	    $id = $stub . "00000000000";
	}
	
	if ($id =~ /\D+(\d+)$/) {
	    
	    my $newid  = $1;
	    my $i;
	    
	    foreach $i ( 1..$number ) {

		$newid++;
		
		
		if (length($newid) > 11) {
		    if ($newid =~ /^0/) {
			$newid =~ s/^0//;
		    } else {
			$self->throw("Can't truncate number string to generate new id [$newid]");
		    }
		}
		my $c = $stub . $newid;
		my $query = "insert into $table (internal_id,external_id) values (NULL,'$c')";
		my $sth   = $self->_db_obj->prepare($query);
		my $res   = $sth->execute;
		
		push(@out,$c);
	    }
	    
	    
	} else {
	    $self->throw("[$id] does not look like an object id (e.g. ENST00000019784)");
	}
    };

    my $error = undef;

    if( $@ ) {
	$error = $@;
    }


    my $usth   = $self->_db_obj->prepare("unlock tables");
    $usth->execute;


    if( defined $error ) {
	$self->throw("Problem in making IDs. Unlocked tables. \n\n Error $@");
    }

    return @out;
    
}

sub get_NewId {
    my ($self,$table,$stub) = @_;

    my $query = "select max(id) as id from $table where id like '$stub%'";
    my $sth   = $self->_db_obj->prepare($query);
    my $res   = $sth->execute;
    my $row   = $sth->fetchrow_hashref;
    my $id    = $row->{id};

    if (!defined($id) || $id eq "") {
	$id = $stub . "00000000000";
    }

    if ($id =~ /\D+(\d+)$/) {

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

