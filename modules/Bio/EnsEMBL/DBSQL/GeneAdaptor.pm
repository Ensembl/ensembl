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
use Bio::EnsEMBL::NewGene;
use Bio::EnsEMBL::NewExon;
use Bio::EnsEMBL::NewTranscript;

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
  my ( $self, $geneId ) = @_;
  
  my $exonAdaptor = $self->db->get_ExonAdaptor();
  my @exons = $exonAdaptor->fetch_by_geneId( $geneId );

  foreach my $exon ( @exons ) {
    $exonIds{$exon->dbID} = $exon;
  }

  my $transcriptAdaptor = $self->db->get_TranscriptAdaptor();
  
  # fetching all exons by gene
  # fetching all transcripts
  # adding the exons
  # adding the transcripts
  my %exonIds;
  my %transcriptExons;

  my $query = qq{
    SELECT tscript.gene
      , tscript.id
      , e_t.exon, e_t.rank
      , gene.version
      , gene.analysisId
    FROM gene
      , transcript tscript
      , exon_transcript e_t
    WHERE gene.id = tscript.gene
      AND tscript.id = e_t.transcript
      AND gene.id = $geneId

    ORDER BY tscript.gene
      , tscript.id
      , e_t.rank
    }

  my $sth = $self->prepare( $query );
  $sth->execute();

  my $first = 1;
  while( my @arr = $sth->fetchrow_array() ) {
    # building a gene
    if( $first ) {
      $gene = Bio::EnsEMBL::NewGene->new();
      $gene->id( $geneId );
      $gene->version( $arr[4] );
      $first = 0;
    } 
    push( @{$transcriptExons{$arr[1]}}, $arr[2] );
  }

  if( $first ) {
    return undef;
  }
  
  foreach my $transcriptId ( keys %transcripts ) {
    # should be fetch_by_geneId ..
    my $transcript = $transcriptAdaptor->fetch_by_dbID( $transcriptId );
    foreach my $exonId ( @{$transcriptExons{$transcriptId}} ) {
      $transcript->add_Exon( $exonIds{$exonId } );
    }
    $gene->add_Transcript( $transcript );
  }
  
  # now add links
  my $sth3 = $self->prepare("
        SELECT gene_id,external_id,external_db 
        FROM genedblink 
        WHERE gene_id = '$geneId'
    ");
  $sth3->execute();
  while( my @arr = $sth3->fetchrow_array ) {
       my $dblink = Bio::Annotation::DBLink->new();
       $dblink->database($arr[2]);
       $dblink->primary_id($arr[1]);
       $gene->add_DBLink( $dblink );
     }

  # and add type
  my $sth4 = $self->prepare("
     SELECT gene_id,type 
     FROM genetype
     WHERE gene_id = '$geneId'
  " );
  $sth4->execute();
  
  if( $arr = $sth4->fetchrow_array ) {
    $gene->type( $arr[1] );
  }
  # we dont fetch the evidence by default
  # $exonAdaptor->fetch_evidence_by_gene( $gene );

}





=head2 fetch_by_Transcript_id

 Title   : fetch_by_Transcript_id
 Usage   : $gene_obj->get_Gene_by_Transcript_id($transid, $supporting)
 Function: gets one gene out of the db with or without supporting evidence
 Returns : gene object (with transcripts, exons and supp.evidence if wanted)
 Args    : transcript id and supporting tag (if latter not specified,
assumes without
           Note that it is much faster to get genes without supp.evidence!


=cut

sub fetch_by_Transcript_id {
    my $self = shift;
    my $transid = shift;

    # this is a cheap SQL call
    my $sth = $self->prepare("select gene from transcript where id = '$transid'");
    $sth->execute;

    my ($geneid) = $sth->fetchrow_array();
    if( !defined $geneid ) {
        return undef;
    }
    my $gene = $self->fetch_by_dbID( $geneid );
    return $gene;
}



=head2 fetch_by_Peptide_id 

 Title   : fetch_by_Peptide_id
 Usage   : $geneAdaptor->fetch_by_Peptide_id($peptideid)
 Function: gets one gene out of the db with or without supporting evidence
 Returns : gene object (with transcripts, exons and supp.evidence if wanted)
 Args    : peptide id and supporting tag (if latter not specified,
assumes without
           Note that it is much faster to get genes without supp.evidence!


=cut

sub fetch_by_Peptide_id {
    my $self = shift;
    my $peptideid = shift;

    # this is a cheap SQL call
    my $sth = $self->prepare("select gene from transcript where translation = '$peptideid'");
    $sth->execute;

    my ($geneid) = $sth->fetchrow_array();
    if( !defined $geneid ) {
        return undef;
    }
    return $self->fetch_by_dbID($geneid);
}


=head2 fetch_by_maximum_DBLink

 Title   : fetch_by_maximum_DBLink
 Usage   : $geneAdaptor->fetch_by_maximum_DBLink($ext_id, $supporting)
 Function: gets one gene out of the db with or without supporting evidence
 Returns : gene object (with transcripts, exons and supp.evidence if wanted)
 Args    : transcript id and supporting tag (if latter not specified,
assumes without
           Note that it is much faster to get genes without supp.evidence!


=cut

sub fetch_by_maximum_DBLink {
    my $self = shift;
    my $external_id = shift;

    my @genes=$self->fetch_by_DBLink( $external_id );

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

=head2 fetch_by_DBLink

 Title   : fetch_by_DBLink
 Usage   : $geneAdptor->fetch_by_DBLink($ext_id)
 Function: gets one gene out of the db with or without supporting evidence
 Returns : gene object (with transcripts, exons and supp.evidence if wanted)
 Args    : transcript id and supporting tag (if latter not specified,
assumes without
           Note that it is much faster to get genes without supp.evidence!


=cut

sub fetch_by_DBLink {
    my $self = shift;
    my $external_id = shift;

    my @genes;
    my $sth = $self->prepare("select gene_id from genedblink where external_id = '$external_id'");
    $sth->execute;
    my $seen=0;
    while (my ($geneid)=$sth->fetchrow_array()) {
	push (@genes,$self->fetch_by_dbID($geneid));
	$seen=1;
    }
    if( !$seen ) {
	return;
    }
    else {
	return @genes;
    }
}


sub get_evidence_by_gene {
  my ($self,$gene) = @_;
  
  $self->get_supporting_evidence_direct( $gene->each_unique_Exon() );
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

    # gets all the features overlapping with exon coordinates

    my $query = qq{
         SELECT f.seq_start,f.seq_end,
            f.score,f.strand,
            f.analysis,f.name,
            f.hstart,f.hend,f.hid,
            f.evalue,f.perc_id,e.id,c.id 
         FROM feature f, exon e, contig c 
         WHERE c.internal_id = f.contig 
            AND f.contig = e.contig 
           AND e.id in ($list) 
           AND !(f.seq_end < e.seq_start OR f.seq_start > e.seq_end) 
	       AND f.strand = e.strand
              AND f.analysis != 3}; #hack for genscan

 # PL: query not checked thoroughly
    my $sth2=$self->prepare($query);
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
    # All features overlapping are added to supporting feature 
    # now the supporting features

    $query = qq{
         SELECT sf.seq_start,sf.seq_end,
            sf.score,sf.strand,
            sf.analysis,sf.name,
            sf.hstart,sf.hend,sf.hid,
            sf.evalue,sf.perc_id,sf.exon 
         FROM supporting_feature sf ,
              exon e
         WHERE e.id  = sf.exon 
           AND e.id in ($list) 
           AND !(sf.seq_end < e.seq_start OR sf.seq_start > e.seq_end) 
	       AND sf.strand = e.strand
              AND sf.analysis != 3}; #hack for genscan

 # PL: query not checked thoroughly
    $sth2=$self->prepare($query);
    $sth2->execute;

    while (my $arrayref = $sth2->fetchrow_arrayref) {
	my ($start,$end,$f_score,$strand,$analysisid,$name,$hstart,$hend,$hid,$evalue,$perc_id,$exonid) = @{$arrayref};
	my $analysis;
	if (!$analhash{$analysisid}) {
	    my $feature_obj=Bio::EnsEMBL::DBSQL::Feature_Obj->new($self->db);
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
	$out->set_all_fields($start,$end,$strand,$f_score,$name,'similarity',1,$hstart,$hend,1,$f_score,$name,'similarity',$hid);
	$out->analysis($analysis);

	#$out->validate(); 
	$exhash{$exonid}->add_Supporting_Feature($out);
    }
  # supporting features are added as well, wired...
}                                       # get_supporting_evidence_direct



=head2 store

 Title   : store
 Usage   : $geneAdaptor->store($gene)
 Function: writes a particular gene into the database. Assumes that everything 
           has dbIDs ....
 Example :
 Returns : nothing
 Args    : $gene object


=cut

sub store {
   my ($self,$gene) = @_;
   my %done;
   my $transcriptAdaptor = $self->db->get_TranscriptAdaptor();
   my $exonAdaptor = $self->db->get_ExonAdaptor();
   my $stub = 'ENSG';

   if ( !defined $gene || ! $gene->isa('Bio::EnsEMBL::Gene') ) {
       $self->throw("$gene is not a EnsEMBL gene - not writing!");
   }

   !$gene->created() && $gene->created(0);
   !$gene->modified() && $gene->modified(0);
 
   my $lsth   = $self->prepare("lock table gene write");
   $lsth->execute;
   
   eval {
     
     my $query = "select max(id) from gene";
	
     my $sth   = $self->prepare($query);
     my $res   = $sth->execute;
     my $row   = $sth->fetchrow_hashref;
     my $id    = $row->{id};
     
     if (!defined($id) || $id eq "") {
       $id = $stub . "00000000001";
     } else {
       $id++;
     }


     my $sth2 = $self->prepare("insert into gene (id,version,created,modified,stored) values ('". 
			       $id       . "','".
			       $gene->version  . "',FROM_UNIXTIME(".
			       $gene->created  . "),FROM_UNIXTIME(".
			       $gene->modified . "),now())");
     $sth2->execute();
     my $usth   = $self->prepare("unlock tables");
     $usth->execute;
     $gene->id( $id );
     $gene->adaptor( $self );
     $gene->dbID( $id );
   };
   my $error;
   
   if( $@ ) {
     $error = $@;
     my $usth   = $self->prepare("unlock tables");
     $usth->execute;
     $self->throw( "Gene store failed" );
   }

   foreach my $dbl ( $gene->each_DBLink ) {
       my $sth3 = $self->prepare("insert into genedblink (gene_id,external_id,external_db) values ('". 
			     $gene->id        . "','".
			     $dbl->primary_id . "','".
			     $dbl->database   . "')");
       $sth3->execute();
   }

   my $id=$gene->id;
   my $type=$gene->type;

   my $sth4 = $self->prepare("insert into genetype (gene_id,type) values ('$id','$type')");

   $sth4->execute();

   # write exons transcripts and exon_transcript table
   foreach my $trans ( $gene->each_Transcript() ) {
       $transcriptAdaptor->store($trans,$gene);
       my $c = 1;
       foreach my $exon ( $trans->each_Exon() ) {
	   
	   if( $done{$exon->id()} ) { 
	       next; 
	   }
	   $done{$exon->id()} = 1;
	   
	   $exonAdaptor->store($exon); 
	   my $sth = $self->prepare("insert into exon_transcript (exon,transcript,rank) values ('". $exon->id()."','".$trans->id()."',".$c.")");
	   $sth->execute();
	   $c++;
	   
       }
   }
   
   return $gene->id;
}

1;



