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


package Bio::EnsEMBL::DBSQL::GeneAdaptor;

use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use vars '@ISA';
#use Bio::EnsEMBL::NewGene;
#use Bio::EnsEMBL::NewExon;
#use Bio::EnsEMBL::NewTranscript;

@ISA = ('Bio::EnsEMBL::DBSQL::BaseAdaptor');


=head2 list_geneIds

 Title   : list_geneIds
 Usage   : $geneAdaptor->list_geneIds
 Function: Gets an array of internal ids for all genes in the current db
 Example : 
 Returns : array of ids
 Args    : none

=cut

sub list_geneIds {
   my ($self) = @_;

   my @out;
   my $sth = $self->prepare("SELECT gene_id FROM gene");
   $sth->execute;

   while (my ($id) = $sth->fetchrow) {
       push(@out, $id);
   }

   return @out;
}

=head2 list_stable_geneIds

 Title   : list_stable_geneIds
 Usage   : $geneAdaptor->list_stable_geneIds
 Function: Gets an array of stable ids for all genes in the current db
 Example : 
 Returns : array of ids
 Args    : none

=cut

sub list_stable_geneIds {
   my ($self) = @_;

   my @out;
   my $sth = $self->prepare("SELECT stable_id FROM gene_stable_id");
   $sth->execute;

   while (my ($id) = $sth->fetchrow) {
       push(@out, $id);
   }

   return @out;
}


=head2 fetch_by_dbID

 Title   : fetch_by_dbID
 Usage   : $geneobj->fetch_by_dbID( $geneid)
 Function: gets one gene out of the db
 Example : $obj->get($dbID)
 Returns : gene object (with transcripts, exons and supp.evidence if wanted)
 Args    : gene id and supporting tag

=cut


sub fetch_by_dbID {
  my ( $self, $geneId ) = @_;
  
  my $exonAdaptor = $self->db->get_ExonAdaptor();
  my @exons = $exonAdaptor->fetch_by_geneId( $geneId );

  my %exonIds;

  #
  # We assumme that if we get no exons this is a throwing situation
  #
  if( scalar( @exons ) == 0 ) {
    $self->throw("No exons for gene $geneId, assumming no gene");
  }


  foreach my $exon ( @exons ) {
    $exonIds{$exon->dbID} = $exon;
  }

  my $transcriptAdaptor = $self->db->get_TranscriptAdaptor();
  
  # fetching all exons by gene
  # fetching all transcripts
  # adding the exons
  # adding the transcripts
  my %transcriptExons;
  my %transcripts;

  my $query = qq{
    SELECT tscript.gene_id
      , tscript.transcript_id
      , e_t.exon_id, e_t.rank
      , gene.analysisId
      , gene.type
      , tscript.translation_id
    FROM gene
      , transcript tscript
      , exon_transcript e_t
    WHERE gene.gene_id = tscript.gene_id
      AND tscript.transcript_id = e_t.transcript_id
      AND gene.gene_id = $geneId

    ORDER BY tscript.gene_id
      , tscript.transcript_id
      , e_t.rank
    };

  my $sth = $self->prepare( $query );
  $sth->execute();

  my $first = 1;
  while( my @arr = $sth->fetchrow_array() ) {
    # building a gene
    if( $first ) {
      $gene = Bio::EnsEMBL::Gene->new();
      $gene->adaptor($self);
      $gene->dbID( $geneId );
      my $ana = $self->db->get_AnalysisAdaptor->fetch_by_dbID($arr[4]);
      $gene->analysis($ana);
      $gene->type($arr[5]);
      $first = 0;
    }

    # store an array of exon ids for each transcript
    if( !exists $transcriptExons{$arr[1]} ) {
      $transcriptExons{$arr[1]} = [];
    }

    push( @{$transcriptExons{$arr[1]}}, $arr[2] );
    $transcripts{$arr[1]} = $arr[6];
  }

  if( $first ) {
    return undef;
  }
  
  foreach my $transcriptId ( keys %transcripts ) {
    # should be fetch_by_geneId ..
    my $transcript = Bio::EnsEMBL::Transcript->new();
    $transcript->dbID( $transcriptId );
    $transcript->adaptor( $self->db->get_TranscriptAdaptor() );
    $transcript->_translation_id($transcripts{$transcriptId} );

    foreach my $exonId ( @{$transcriptExons{$transcriptId}} ) {
      $transcript->add_Exon( $exonIds{$exonId} );
    }
    $gene->add_Transcript( $transcript );
  }
  
  return $gene;
}

=head2 fetch_by_stable_id

 Title   : fetch_by_stable_id
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub fetch_by_stable_id{
   my ($self,$id) = @_;

   my $sth = $self->prepare("select gene_id from gene_stable_id where stable_id = '$id'");
   $sth->execute;

   my ($dbID) = $sth->fetchrow_array();

   if( !defined $dbID ) {
       $self->throw("No stable id with $id, cannot fetch");
   }

   return $self->fetch_by_dbID($dbID);
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
    my $sth = $self->prepare("select gene_id from transcript where transcript_id = '$transid'");
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
    my $sth = $self->prepare("select gene_id from transcript where translation_id = '$peptideid'");
    $sth->execute;

    my ($geneid) = $sth->fetchrow_array();
    if( !defined $geneid ) {
        return undef;
    }
    return $self->fetch_by_dbID($geneid);
}


=head2 fetch_by_maximum_DBLink

 Title   : fetch_by_maximum_DBLink
 Usage   : $geneAdaptor->fetch_by_maximum_DBLink($ext_id)
 Function: gets one gene out of the db with 
 Returns : gene object (with transcripts, exons)
 Args    : 


=cut

sub fetch_by_maximum_DBLink {
    my $self = shift;
    my $external_id = shift;

    my @genes=$self->fetch_by_DBEntry( $external_id );

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


=head2 get_description

 Title   : get_description
 Usage   : $geneAdptor->get_description($dbID)
 Function: gets gene description line 
 Returns : a string
 Args    : 


=cut

sub get_description {
  my ($self,$dbID) = @_;

  if( !defined $dbID ) {
      $self->throw("must call with dbID");
  }


  my $sth = $self->prepare("select description from gene_description where gene_id = $dbID");

  $sth->execute;
  my @array = $sth->fetchrow_array();
  return $array[0];
}

=head2 get_stable_entry_info

 Title   : get_stable_entry_info
 Usage   : $geneAdptor->get_stable_entry_info($gene)
 Function: gets stable info for gene and places it into the hash
 Returns : 
 Args    : 


=cut

sub get_stable_entry_info {
  my ($self,$gene) = @_;

  if( !defined $gene || !ref $gene || !$gene->isa('Bio::EnsEMBL::Gene') ) {
     $self->throw("Needs a gene object, not a $gene");
  }

  my $sth = $self->prepare("select stable_id,UNIX_TIMESTAMP(created),UNIX_TIMESTAMP(modified),version from gene_stable_id where gene_id = ".$gene->dbID);
  $sth->execute();

  my @array = $sth->fetchrow_array();
  $gene->{'_stable_id'} = $array[0];
  $gene->{'_created'}   = $array[1];
  $gene->{'_modified'}  = $array[2];
  $gene->{'_version'}   = $array[3];
  

  return 1;
}

=head2 fetch_by_DBEntry

 Title   : fetch_by_DBLink
 Usage   : $geneAdptor->fetch_by_DBLink($ext_id)
 Function: gets one gene out of the db with or without supporting evidence
 Returns : gene object (with transcripts, exons and supp.evidence if wanted)
 Args    : transcript id and supporting tag (if latter not specified,
assumes without
           Note that it is much faster to get genes without supp.evidence!


=cut

sub fetch_by_DBEntry {
  my $self = shift;
  $self->throw( "Not implemented yet" );
}




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

   if( !defined $gene || !ref $gene || !$gene->isa('Bio::EnsEMBL::Gene') ) {
       $self->throw("Must store a gene object, not a $gene");
   }

   if( !defined $gene->analysis ) {
       $self->throw("Genes must have an analysis object!");
   }

   my $analysisId = $self->db->get_AnalysisAdaptor()->store( $gene->analysis );


   if ( !defined $gene || ! $gene->isa('Bio::EnsEMBL::Gene') ) {
       $self->throw("$gene is not a EnsEMBL gene - not writing!");
   }

 
   my $sth2 = $self->prepare("insert into gene set analysisId=$analysisId, type='".
			    $gene->type()."'" );
   $sth2->execute();
   
   $gene->adaptor( $self );
   $gene->dbID( $sth2->{'mysql_insertid'} );

   my $dbEntryAdaptor = $self->db->get_DBEntryAdaptor();

   foreach my $dbl ( $gene->each_DBLink ) {
     $dbEntryAdaptor->store( $dbl, $gene->dbID, "Gene" );
   }

   # write exons transcripts and exon_transcript table
   foreach my $trans ( $gene->each_Transcript() ) {
     $transcriptAdaptor->store($trans,$gene);
   }
   
   return $gene->dbID;
}


sub remove {
  my $self = shift;
  my $gene = shift;

  if( ! defined $gene->dbID() ) {
    return;
  }

  my $sth= $self->prepare( "delete from gene where gene_id = ? " );
  $sth->execute( $gene->dbID );
  $sth= $self->prepare( "delete from gene_stable_id where gene_id = ? " );
  $sth->execute( $gene->dbID );
  my $transcriptAdaptor = $self->db->get_TranscriptAdaptor();
  # my $dbEntryAdaptor = $self->db->get_DBEntryAdaptor();
  foreach my $trans ( $gene->each_Transcript() ) {
    $transcriptAdaptor->remove($trans,$gene);
  }

  $gene->{'_dbID'} = undef;
}





sub create_tables {
# read sql from geneAdaptor.sql
}



1;
__END__

