# EnsEMBL Exon reading writing adaptor for mySQL
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

Bio::EnsEMBL::DBSQL::ExonAdaptor - MySQL Database queries to generate and store gens.

=head1 SYNOPSIS

=head1 CONTACT

  Arne Stabenau: stabenau@ebi.ac.uk
  Elia Stupka  : elia@ebi.ac.uk
  Ewan Birney  : 

=head1 APPENDIX

=cut



package Bio::EnsEMBL::DBSQL::ExonAdaptor;

use vars qw( @ISA );
use strict;


use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
#use Bio::EnsEMBL::NewGene;
#use Bio::EnsEMBL::NewExon;
#use Bio::EnsEMBL::NewTranscript;

@ISA = qw( Bio::EnsEMBL::DBSQL::BaseAdaptor );

=head2 fetch_by_dbID

 Title   : fetch_by_dbID
 Usage   : $exonAdaptor->fetch_by_dbID($exon_id)
 Function: 
 Example : $obj->remove_by_dbID(ENSE000034)
 Returns : nothing
 Args    : $exon_id

=cut


sub fetch_by_dbID {
  my ( $self, $dbID ) = @_;
  my ( $hashref, $arrRef );

  my $query = qq {
    SELECT id, version, contig, 
           UNIX_TIMESTAMP(created), 
           UNIX_TIMESTAMP(modified)
	   , seq_start, seq_end, strand, 
           phase, 
           UNIX_TIMESTAMP(stored), 
           end_phase, sticky_rank
    FROM exon 
    WHERE id = $dbID
    ORDER BY sticky_rank
  };

  my $sth = $self->prepare( $query );
  $sth->execute();
  
  my( $exon );
  if( $arrRef = $sth->fetchrow_arrayref() ) {
    # exons are made from start end strand context_name context_type
    $exon = Bio::EnsEMBL::NewExon->new( $arrRef->[5..7], $arrRef->[2], "RawContig" );
    $exon->id       ( $arrRef->[0]  ),
    $exon->version  ( $arrRef->[1]  ),
    $exon->created  ( $arrRef->[3]  ),
    $exon->modified ( $arrRef->[4]  ),
    $exon->phase    ( $arrRef->[8]  ),
    $exon->end_phase( $arrRef->[10] ),
  } else {
    # exon id not known
    return undef;
  }

  while( my $arrRef = $sth->fetchrow_arrayref() ) {
    $exon->add_Range( $arrRef->[5..7], $arrRef->[2], "RawContig" );
  }
  return $exon;
}


# returns list of exons or maybe empty list (gene not known)
sub fetch_by_geneId {
  my ( $self, $geneId ) = @_;
  my $hashRef;
  my %exons;
  my @out;
  my ( $currentId, $currentTranscript );

  my $query = qq {
    SELECT t.id tid
      , e.id eid
      , e.version
      , e.contig
      , UNIX_TIMESTAMP(e.created)
      , UNIX_TIMESTAMP(e.modified)
      , e.seq_start
      , e.seq_end
      , e.strand
      , e.phase
      , UNIX_TIMESTAMP(e.stored)
      , e.end_phase
      , e.sticky_rank
    FROM exon e
      , exon_transcript et
      , transcript t
    WHERE t.gene = '$geneId'
      AND et.transcript = t.id
      AND e.id = et.exon
    ORDER BY t.id
      , e.id
      , e.sticky_rank
  };

  my $sth = $self->prepare( $query );
  $sth->execute();

  $hashRef = $sth->fetchrow_hashref();

    my( $exon );
  while( 1 ) {
    my( $loading );
    if( $hashRef ) {
      if( ! exists $exons{ $hashRef->{eid} } ) {
        # exons are made from start end strand context_name context_type
        $exon = Bio::EnsEMBL::NewExon->new
          ( $hashRef->{start}, $hashRef->{end},
            $hashRef->{strand}, $hashRef->{contig}, "RawContig" );
        $exon->id( $hashRef->{eid}),
        $currentId = $hashRef->{eid};
        $currentTranscript = $hashRef->{tid};
        $exon->version( $hashRef->{version} ),
        $exon->created( $hashRef->{created} ), 
        $exon->modified( $hashRef->{modified} ),
        $exon->phase( $hashRef->{phase} ),
        $exon->end_phase( $hashRef->{end_phase} ),
 
        $exons{$currentId} = $exon;
        $loading = 1;
      } else {
        # exon known
        $loading = 0;
      }
    } else {
      return values %exons;
    }

    while( $hashRef = $sth->fetchrow_hashref() ) {
      if( $hashRef->{eid} == $currentId &&
          $hashRef->{tid} == $currentTranscript ) {
        if( $loading ) {
          $exon->add_Range( $hashRef->{start}, $hashRef->{end},
            $hashRef->{strand}, $hashRef->{contig}, "RawContig" );
        }
      } else {
        last;
      }
    }
  }
  return @out;
}

=head2 fetch_by_Transcript

 Title   : fetch_by_Transcript
 Usage   : $exonAdaptor->fetch_by_Transcript( $transcript )
 Function: 
 Example : $obj->fetch_by_Transcript( $transcript )
 Returns : all exons in that transcript
 Args    : $exon_id

=cut


sub fetch_by_Transcript {
  my ( $self, $transcript ) = @_;
  my $transcriptId = $transcript->dbID;
  return $self->fetch_by_transcriptId( $transcriptId );
}

sub fetch_by_transcriptId {
  my ( $self, $transcriptId ) = @_;

}


sub fetch_evidence_by_Exons {
  my ( $self, @exons ) = @_;

  my $instring = "'";
  my %exhash;
  my $anaAdaptor = $self->db->get_AnalysisAdaptor();

  foreach my $exon (@exons) {
     $exhash{$exon->id} = $exon;
     $instring = $instring . $exon->id . "','";
  }
    
  $instring = substr($instring,0,-2);
   
  my $query = "select * from supporting_feature where exon in (" . $instring . ")";

  my $sth = $self->prepare($query);
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
        $f1->percent_id($rowhash->{perc_id});
        $f1->p_value($rowhash->{p_value});
        $f1->phase  ($rowhash->{phase});
	$f1->score  ($rowhash->{score});
	$f2->seqname($rowhash->{hid});
	$f2->start  ($rowhash->{hstart});
	$f2->end    ($rowhash->{hend});
	$f2->strand ($rowhash->{strand});
	$f2->source_tag($rowhash->{name});
	$f2->primary_tag('similarity');
	$f2->score  ($rowhash->{score});
        $f2->percent_id($rowhash->{perc_id});
        $f2->p_value($rowhash->{p_value});
        $f2->phase  ($rowhash->{phase});

        my $analysis = $anaAdaptor->fetch_by_dbID( $rowhash->{analysis} );

        $f->analysis($analysis);

	$f->validate;
 
	$exhash{$exon}->add_Supporting_Feature($f);
    }

}


=head2 fetch_evidence_by_Gene

 Title   : fetch_evidence_by_Gene
 Usage   : $exonAdaptor->fetch_evidence_by_Gene($geneObject)
 Function: Fetch evidence for all exons in the Gene.
 Example : $exonAdaptor->fetch_evidence_by_Gene( $geneObject );
 Returns : nothing
 Args    : 

=cut

sub fetch_evidence_by_Gene {
  my ( $self, $gene )  = @_;
  my $instring;

  my @exons = $gene->each_unique_exon();
  foreach my $exon ( @exons ) {
    if( $exon->get_context_type() ne "RawContig" ) {
      $self->warn( "Exons not in right coord system" );       
      return;
    }
    if( $exon->is_split() ) {
      for my $range ( $exon->each_range() ) {
	$instring = $instring . $range->[3] . ","; # range->[3] context_name
      }
    } else {
        $instring = $exon->get_context_name();
    }
  }
    
  $instring = substr($instring,0,-1);

   
  my $statement = "select * from feature f where contig in (" . $instring . ")";
  
  #my $statement = "select * from supporting_feature where exon in (" . $instring . ")";
  # print STDERR "going to execute... [$statement]\n";
  
  my $sth = $self->prepare($statement);
  $sth->execute || $self->throw("execute failed for supporting evidence get!");
  
  my @features;
    
  my $anaAdaptor = $self->db->get_AnalysisAdaptor;
  
  while (my $rowhash = $sth->fetchrow_hashref) {
      my $f1 = new Bio::EnsEMBL::SeqFeature;
      my $f2 = new Bio::EnsEMBL::SeqFeature;
      
      my $f = new Bio::EnsEMBL::FeaturePair(-feature1 => $f1,
					    -feature2 => $f2);
      
      #	    my $exon = $rowhash->{exon};
      
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
      
      my $analysis = $anaAdaptor->fetch_by_dbID( $rowhash->{analysis} );
      
      $f->analysis($analysis);
	
      $f->validate;
      push(@features,$f);
    }

  # dodgy, start end in which coord system?

  foreach my $exon (@exons) {
    foreach my $f (@features) {
      if ($f->start == $exon->start && $f->end == $exon->end) {
	$exon->add_Supporting_Feature($f);
      }
    }
  }
  
}



=head2 store

 Title   : store
 Usage   : $exonAdaptor->store($exonObject)
 Function: Stores the exon. Be careful to call it just once per gene.
 Example : $dbID = $exonAdaptor->store( $exon );
 Returns : nothing
 Args    : $exon_id

=cut

sub store {
  my ( $self, $exon ) = @_;
  my $stub = 'ENSE';

  if( ! $exon->isa('Bio::EnsEMBL::NewExon') ) {
    $self->throw("$exon is not a EnsEMBL exon - not dumping!");
  }

  if( $exon->get_context_type() ne "RawContig" ) {
    $self->throw( "Exon not stored in the right coord system" );
  }

  my $lsth   = $self->prepare("lock table exon write");
  $lsth->execute;
   
  eval {
     
     my $query = "select max(id) from exon";
	
     my $sth   = $self->prepare($query);
     my $res   = $sth->execute;
     my $row   = $sth->fetchrow_hashref;
     my $id    = $row->{id};
     
     if (!defined($id) || $id eq "") {
       $id = $stub . "00000000001";
     } else {
       $id++;
     }


     my $exonst = q{
       INSERT into exon (id, version, contig, created, modified
			 , seq_start, seq_end, strand, phase, 
			 stored, end_phase, rank) 
	 VALUES (?,?,?,FROM_UNIXTIME(?),FROM_UNIXTIME(?),?,?,?,?,NOW(),?,?)
       };

     my @loc = $exon->each_range_context();

     my $sticky_rank = 1;

     foreach my $range ( @loc ) {
       my $contigId = $range->[4]->dbID; # dbID from context / RawContig
       
       my $sth = $self->prepare($exonst);
       $sth->execute(
		     $id,
		     $exon->version(),
		     $contigId,
		     $exon->created(),
		     $exon->modified(),
		     $range->[0], # start
		     $range->[1], # end
		     $range->[2], # strand
		     $exon->phase,
		     $exon->end_phase,
		     $sticky_rank
		    );
       $sticky_rank++;
     }
     my $usth   = $self->prepare("unlock tables");
     $usth->execute;
     $exon->id( $id );
     $exon->adaptor( $self );
     $exon->dbID( $id );
   };
  my $error;
  
  if( $@ ) {
    $error = $@;
    $self->throw( "Exon store failed" );
  }

  # Now the supporting evidence
  # should be stored from featureAdaptor

  my $sth  = $self->prepare("
     INSERT INTO supporting_feature( 
        id,exon,seq_start,seq_end,score,
        strand,analysis,name,hstart,hend,hid) 
     VALUES(?,?,?,?,?,?,?,?,?,?,?) 
   ");
    
#  my $feature_obj=Bio::EnsEMBL::DBSQL::Feature_Obj->new($self->dbadaptor);
  my $anaAdaptor = $self->dbadaptor->get_AnalysisAdaptor();

  FEATURE: foreach my $f ($exon->each_Supporting_Feature) {
	#print STDERR "Writing supporting feature ".$f->source_tag."\n";
	eval {
	    $f->validate();
	};

	if ($@) {
	    print(STDERR "Supporting feature invalid. Skipping feature\n");
	    next FEATURE;
	}

  	# my $analysisid = $feature_obj->write_Analysis($f->analysis);
	my $analysisid;
	if( !$f->analysis->adaptor == $anaAdaptor ) {
	  $analysisid = $f->analysis->dbID();
	} else {
	  $analysisid = $anaAdaptor->store( $f->analysis );
	  $f->analysis->dbID( $analysisid );
	}
	
	if ($f->isa("Bio::EnsEMBL::FeaturePair")) {
	    $sth->execute('NULL',
			  $exon->dbID(),
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


=head2 update

 Title   : update
 Usage   : $exonAdaptor->update($exonObject)
 Function: Updates the exon. When you change the exon (eg by adding 
           supporting features) then you do an update.
 Example : $dbID = $exonAdaptor->store( $exon );
 Returns : nothing
 Args    : $exon_id

=cut

sub update {
  my ( $self, $exon ) = @_;
  if( ! $exon->isa('Bio::EnsEMBL::Exon') ) {
    $self->throw("$exon is not a EnsEMBL exon - not dumping!");
  }

  if( $exon->get_context_type() ne "RawContig" ) {
    $self->throw( "Exon not stored in the right coord system" );
  }

  if( $exon->adaptor != $self ) {
    $self->throw( "You have to store before you can update" );
  }

  my $sth = $self->prepare( "
    DELETE FROM exon
    WHERE id = ?
  " );
 
  $sth->execute( $exon->dbID );
  $sth = $self->prepare( "
    DELETE FROM supporting_evidence
    WHERE exon = ?
  " );
  $sth->execute( $exon->dbID );

  my $exonst = q{
        INSERT into exon (id, version, contig, created, modified
          , seq_start, seq_end, strand, phase, stored, end_phase, rank) 
        VALUES (?,?,?,FROM_UNIXTIME(?),FROM_UNIXTIME(?),?,?,?,?,NOW(),?,?)
        };

  my @loc = $exon->each_range_context();

  my $sticky_rank = 1;

  foreach my $range ( @loc ) {
    my $contigId = $range->[4]->dbID; # dbID from context / RawContig

    my $sth = $self->prepare($exonst);
    $sth->execute(
        $exon->id(),
        $exon->version(),
        $contigId,
        $exon->created(),
        $exon->modified(),
        $range->[0], # start
        $range->[1], # end
        $range->[2], # strand
        $exon->phase,
        $exon->end_phase,
	$sticky_rank
        );
    $sticky_rank++;
  }

  # Now the supporting evidence
  # should be stored from featureAdaptor

  $sth  = $self->prepare("
     INSERT INTO supporting_feature( 
        id,exon,seq_start,seq_end,score,
        strand,analysis,name,hstart,hend,hid) 
     VALUES(?,?,?,?,?,?,?,?,?,?,?) 
   ");
    
#  my $feature_obj=Bio::EnsEMBL::DBSQL::Feature_Obj->new($self->dbadaptor);
  my $anaAdaptor = $self->dbadaptor->get_AnalysisAdaptor();

  FEATURE: foreach my $f ($exon->each_Supporting_Feature) {
	#print STDERR "Writing supporting feature ".$f->source_tag."\n";
	eval {
	    $f->validate();
	};

	if ($@) {
	    print(STDERR "Supporting feature invalid. Skipping feature\n");
	    next FEATURE;
	}

  	# my $analysisid = $feature_obj->write_Analysis($f->analysis);
	my $analysisid;
	if( !$f->analysis->adaptor == $anaAdaptor ) {
	  $analysisid = $f->analysis->dbID();
	} else {
	  $analysisid = $anaAdaptor->store( $f->analysis );
	}
	
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

=head2 remove_by_dbID

 Title   : remove_by_dbID
 Usage   : $exonAdaptor->remove_by_dbID($exon_id)
 Function: Deletes exon, including exon_transcript rows
 Example : $obj->remove_by_dbID(ENSE000034)
 Returns : nothing
 Args    : $exon_id

=cut

sub remove_by_dbID {
    my ($self,$exon_id) = @_;

    $exon_id || $self->throw ("Trying to delete an exon without an exon_id\n");
    
    # Delete exon_transcript rows
    my $sth = $self->prepare("delete from exon_transcript where exon = '".$exon_id."'");
    my $res = $sth ->execute;

    #Delete exon rows
    $sth = $self->prepare("delete from exon where id = '".$exon_id."'");
    $res = $sth->execute;

    # delete the supporting evidence
    $sth = $self->prepare("delete from supporting_feature where exon = '" . $exon_id . "'");
    $res = $sth->execute;
}



sub create_tables {
  my $self = shift;

  my $sth = $self->prepare( "drop table if exists exon" );
  $sth->execute();

  $sth = $self->prepare( qq{
    CREATE TABLE exon (
      id            varchar(40) NOT NULL,
      contig        int(10) unsigned NOT NULL,
      version       int(10) DEFAULT '1' NOT NULL,
      created       datetime NOT NULL,
      modified      datetime NOT NULL,
      stored        datetime NOT NULL,
      seq_start     int(10) NOT NULL,
      seq_end       int(10) NOT NULL,
      strand        int(2) NOT NULL,
      phase         int(11) NOT NULL,
      end_phase     int(11) NOT NULL,
      sticky_rank   int(10) DEFAULT '1' NOT NULL,
  
      PRIMARY KEY (id,sticky_rank),
      KEY id_contig (id,contig),
      KEY contig (contig)
    )
  } );
  $sth->execute();
}
