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

Bio::EnsEMBL::DBOLD::ExonAdaptor - MySQL Database queries to generate and store gens.

=head1 SYNOPSIS

=head1 CONTACT

  Arne Stabenau: stabenau@ebi.ac.uk
  Elia Stupka  : elia@ebi.ac.uk
  Ewan Birney  : 

=head1 APPENDIX

=cut



package Bio::EnsEMBL::DBOLD::ExonAdaptor;

use vars qw( @ISA );
use strict;


use Bio::EnsEMBL::DBOLD::BaseAdaptor;
use Bio::EnsEMBL::DBOLD::DBAdaptor;
use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::Exon;
use Bio::EnsEMBL::Transcript;

@ISA = qw( Bio::EnsEMBL::DBOLD::BaseAdaptor );

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
  $self->throw( "Not implemented yet" );
}

=head2 fetch_by_Transcript

 Title   : fetch_by_Transcript
 Usage   : $exonAdaptor->fetch_by_Transcript( $transcript )
 Function: 
 Example : $obj->remove_by_dbID(ENSE000034)
 Returns : nothing
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


  foreach my $exon (@exons) {
     $exhash{$exon->id} = $exon;
     $instring = $instring . $exon->id . "','";
  }
    
  $instring = substr($instring,0,-2);
   
  my $sth = $self->prepare("select * from supporting_feature where exon in (" . $instring . ")");
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
	    my $feature_obj=Bio::EnsEMBL::DBOLD::Feature_Obj->new($self->_db_obj);
	    $f->analysis($feature_obj->get_Analysis($analysisid));

	    $anahash{$analysisid} = $f->analysis;
	}
	
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
   
  foreach my $exon ($gene->each_unique_exon()) {
    if( $exon->get_context_type() ne "RawContig" ) {
      $self->warn( "Exons not in right coord system" );       
      return;
    }
    if( $exon->is_split() ) {
      for my $range ( $exon->each_range() ) {
	$instring = $instring . $range->[3] . "','"; # range->[3] context_name
      }
    } else {
        $instring = $instring . $exon->get_context_name() . "','";
    }
  }
    
  $instring = substr($instring,0,-2);

   
  my $statement = "select * from feature f,contig c where c.id in (" . $instring . ")";

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
	
	$f->analysis($analysis)
	    
	$f->validate;
	push(@features,$f);
    }
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
  if( ! $exon->isa('Bio::EnsEMBL::Exon') ) {
    $self->throw("$exon is not a EnsEMBL exon - not dumping!");
  }

  if( $exon->get_context_type() ne "RawContig" ) {
    $self->throw( "Exon not stored in the right coord system" );
  }

  my $exonst = q{
        INSERT into exon (id, version, contig, created, modified
          , seq_start, seq_end, strand, phase, stored, end_phase, rank) 
        VALUES (?,?,?,FROM_UNIXTIME(?),FROM_UNIXTIME(?),?,?,?,?,NOW(),?,?)
        };

  my @loc = $exon->each_range_context();

  my $sticky_rank = 1;

  foreach my $range ( @loc ) {
    my $contigId = $range->[4]->dbID;

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

  $sth  = $self->prepare("insert into supporting_feature(id,exon,seq_start,seq_end,score,strand,analysis,name,hstart,hend,hid) values(?,?,?,?,?,?,?,?,?,?,?)");
    
#  my $feature_obj=Bio::EnsEMBL::DBOLD::Feature_Obj->new($self->dbadaptor);
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
  	my $analysisid = $anaAdaptor->store( $f->analysis );
	
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
