#EnsEMBL Exon reading writing adaptor for mySQL
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

Bio::EnsEMBL::DBSQL::ExonAdaptor - MySQL Database queries to generate and store exons (including supporting evidence)

=head1 SYNOPSIS

$exon_adaptor = $database_adaptor->get_ExonAdaptor();
$exon = $exon_adaptor->fetch_by_dbID

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
use Bio::EnsEMBL::Exon;
use Bio::EnsEMBL::StickyExon;

@ISA = qw( Bio::EnsEMBL::DBSQL::BaseAdaptor );


=head2 fetch_by_dbID

  Arg [1]    : int $exon_id
               the unique internal id of the exon to retrieve
  Example    : $exon = $exon_adaptor->fetch_by_dbID($exon_id);
  Description: Retrieves an exon from the database via its internal id
  Returntype : Bio::EnsEMBL::Exon in contig coordinates
  Exceptions : none
  Caller     : general

=cut

sub fetch_by_dbID {
  my $self = shift;
  my $dbID = shift;

  my $query = qq {
    SELECT  e.exon_id
      , e.contig_id
      , e.contig_start
      , e.contig_end
      , e.contig_strand
      , e.phase
      , e.end_phase
      , e.sticky_rank
    FROM exon e
    WHERE e.exon_id = $dbID
    ORDER BY e.sticky_rank DESC  };

  
  my $sth = $self->prepare($query);

  $sth->execute();


  my $hashRef;
  my $exon;

  if( $hashRef = $sth->fetchrow_hashref() ) {
    $exon = $self->_exon_from_sth( $sth, $hashRef );
  }

  delete $self->{rchash};
  return $exon;
}


=head2 fetch_by_stable_id

  Arg [1]    : string $stable_id
               the stable id of the exon to retrieve
  Example    : $exon = $exon_adaptor->fetch_by_stable_id('ENSE0000988221');
  Description: Retrieves an Exon from the database via its stable id
  Returntype : Bio::EnsEMBL::Exon in contig coordinates
  Exceptions : none
  Caller     : general

=cut

sub fetch_by_stable_id {
  my $self = shift;
  my $stable_id = shift;

  my $sth = $self->prepare( "SELECT exon_id from exon_stable_id
                             WHERE stable_id = ?" );
  $sth->execute( $stable_id );
  if( my $arr = $sth->fetchrow_arrayref ) {
    my $exon = $self->fetch_by_dbID( $arr->[0] );
    return $exon;
  } else {
    $self->warn( "No Exon with this stable_id in the database!" );
    return undef;
  }
}



=head2 fetch_all_by_gene_id

  Arg [1]    : int $id
               The identifier of the gene whose exons will be retrieved 
  Example    : @exons = $exon_adaptor->fetch_all_by_gene_id(1234); 
  Description: Retrieves all exons from the gene specified by $geneId
  Returntype : listref of Bio::EnsEMBL::Exon in contig coordinates
  Exceptions : thrown if $geneId is not defined  
  Caller     : general

=cut

sub fetch_all_by_gene_id {
  my ( $self, $gene_id ) = @_;
  my %exons;
  my $hashRef;
  my ( $currentId, $currentTranscript );

  if( !$gene_id ) {
      $self->throw("Gene dbID not defined");
  }
  $self->{rchash} = {};
  my $query = qq {
    SELECT 
      STRAIGHT_JOIN 
	e.exon_id
      , e.contig_id
      , e.contig_start
      , e.contig_end
      , e.contig_strand
      , e.phase
      , e.end_phase
      , e.sticky_rank
    FROM transcript t
      , exon_transcript et
      , exon e
    WHERE t.gene_id = ?
      AND et.transcript_id = t.transcript_id
      AND e.exon_id = et.exon_id
    ORDER BY t.transcript_id,e.exon_id
      , e.sticky_rank DESC
  };

  my $sth = $self->prepare( $query );
  $sth->execute($gene_id);

  while( $hashRef = $sth->fetchrow_hashref() ) {
    if( ! exists $exons{ $hashRef->{exon_id} } ) {

      my $exon = $self->_exon_from_sth( $sth, $hashRef );

      $exons{$exon->dbID} = $exon;
    }
  }
  delete $self->{rchash};
  
  my @out = ();

  push @out, values %exons;

  return \@out;
}


=head2 _exon_from_sth

  Arg [1]    : DBI statement handle $sth
               the prepared SQL statement used to generate the hashref
               this arg is required to retrieve component exons from
               the database if the exon happens to be a sticky exon 
  Arg [2]    : DBI row hashref $hashRef
               a hash reference representing a single exon or a portion
               of a sticky exon
  Example    : my $exon = $self->_exon_from_sth($sth, $hashref); 
  Description: PROTECTED
               Creates an exon (normal or sticky) from its statement handle
               and hashreference.  
  Returntype : Bio::EnsEMBL::Exon
  Exceptions : none
  Caller     : Bio::EnsEMBL::DBSQL::ExonAdaptor

=cut


# build an exon (possibly sticky exon) from given statement handle and
# retrieved row in hashref. return it. uses $self->{rchash} to store contig
# information. This should be cleared after leaving ExonAdaptor

sub _exon_from_sth {

  my ( $self, $sth, $hashRef ) = @_;
  my $sticky_length = 0;
  my $exon;
  if( $hashRef->{'sticky_rank'} >1 ) {	
    
    # sticky exon
    $exon = Bio::EnsEMBL::StickyExon->new();
    $exon->dbID($hashRef->{'exon_id'});
    # make first component exon
    my $component = $self->_new_Exon_from_hashRef($hashRef);
    
    $exon->add_component_Exon($component);
    $sticky_length += $component->length;

    $exon->phase($component->phase);
    $exon->end_phase($component->end_phase);
    $exon->adaptor($self);

    # continue while loop until we hit sticky_rank 1
    while( $hashRef = $sth->fetchrow_hashref() ) {
      my $component = $self->_new_Exon_from_hashRef($hashRef);

      $exon->add_component_Exon($component);
      $sticky_length += $component->length;

      if( $component->sticky_rank == 1 ) {
	$exon->contig( $component->contig );
	last;
      }
    }

    $exon->_sort_by_sticky_rank();

    # set start = 1 and end = length of sticky exon
    # build a minature sequence representing the sticky region and
    # attach

    $exon->start(1);
    $exon->end($sticky_length);
    $exon->strand( 1 );

  } 
  else {
    $exon = $self->_new_Exon_from_hashRef($hashRef);
  }

  return $exon;
}



=head2 _new_Exon_from_hashref

  Arg [1]    : DBI row hashref $hashRef
  Example    : my $exon = $hashref
  Description: PROTECTED
               creates a normal (non-sticky) exon from a row hash reference 
  Returntype : Bio::EnsEMBL::Exon
  Exceptions : none
  Caller     : _exon_from_sth

=cut

sub _new_Exon_from_hashRef {
   my $self = shift;
   my $hashRef = shift;

   my $exon = Bio::EnsEMBL::Exon->new();
   $exon->start( $hashRef->{'contig_start'} );
   $exon->end( $hashRef->{'contig_end'} );
   $exon->strand( $hashRef->{'contig_strand'} );
   $exon->phase( $hashRef->{phase} );
   $exon->end_phase( $hashRef->{end_phase} );
   $exon->dbID($hashRef->{'exon_id'});
   $exon->sticky_rank($hashRef->{'sticky_rank'});
   $exon->adaptor($self);

   my $rc = 
     $self->db->get_RawContigAdaptor->fetch_by_dbID($hashRef->{'contig_id'});

   $exon->contig( $rc );
   
  return $exon;
}





=head2 store

  Arg [1]    : Bio::EnsEMBL::Exon $exon
               the exon to store in this database
  Example    : $exon_adaptor->store($exon);
  Description: Stores an exon in the database
  Returntype : none
  Exceptions : thrown if exon (or component exons) do not have a contig_id
               or if $exon->start, $exon->end, $exon->strand, or $exon->phase 
               are not defined or if $exon is not a Bio::EnsEMBL::Exon 
  Caller     : general

=cut

sub store {
  my ( $self, $exon ) = @_;

  if( ! $exon->isa('Bio::EnsEMBL::Exon') ) {
    $self->throw("$exon is not a EnsEMBL exon - not storing!");
  }

  if( $exon->dbID && $exon->adaptor && $exon->adaptor == $self ) {
    return $exon->dbID();
  }

  if( ! $exon->start || ! $exon->end ||
      ! $exon->strand || ! defined $exon->phase ) {
    $self->throw("Exon does not have all attributes to store");
  }

  # trap contig_id separately as it is likely to be a common mistake

  my $exon_sql = q{
    INSERT into exon ( exon_id, contig_id, contig_start, 
		       contig_end, contig_strand, phase, 
		       end_phase, sticky_rank )
    VALUES ( ?, ?, ?, ?, ?, ?, ?,? ) 
  };
  my $exonst = $self->prepare($exon_sql);

  my $exonId = undef;

  if( $exon->isa( 'Bio::EnsEMBL::StickyExon' )) {
    # sticky storing. Sticky exons contain normal exons ...

    my $componentExons = $exon->get_all_component_Exons;
    for my $componentExon ( @$componentExons ) {
      my $contig = $componentExon->contig();

      unless($contig && ref $contig && $contig->dbID()) {
	$self->throw("Component Exon does not have an attached contig " .
		     "with a valid set database id. " .
		     "Needs to have one set");
      }

      $exonst->execute( $exonId, $contig->dbID,
			$componentExon->start(),
			$componentExon->end(),
			$componentExon->strand(),
			$componentExon->phase(),
			$componentExon->end_phase(),
			$componentExon->sticky_rank() );
      if( ! $exonId ) {
	$exonId = $exonst->{'mysql_insertid'};
      }
    }
  } else {
    # normal storing

    my $contig = $exon->contig();

    unless( $contig && ref $contig && $contig->dbID() ) {
      $self->throw("Exon does not have an attached contig with a valid " . 
		   "database id.  Needs to have one set");
    }

    $exonst->execute( undef,$contig->dbID,
		      $exon->start(),
		      $exon->end(),
		      $exon->strand(),
		      $exon->phase(),
		      $exon->end_phase(),
		      $exon->sticky_rank() );
    $exonId = $exonst->{'mysql_insertid'};
  }

  if ($exon->stable_id) {
    if (!$exon->created || 
        !$exon->modified ||
        !$exon->version) {
      $self->throw("Trying to store incomplete stable id information for exon");
    }

    my $statement = "INSERT INTO exon_stable_id(exon_id," .
	"version, stable_id, created, modified)".
                    " VALUES(" . $exonId . "," .
                               $exon->version . "," .
                               "'" . $exon->stable_id . "'," .
                               "FROM_UNIXTIME(".$exon->created."),".
                               "FROM_UNIXTIME(". $exon->modified . "))";
    #print $statement . "\n";

     my $sth = $self->prepare($statement);
     $sth->execute();
   }


  # Now the supporting evidence
  # should be stored from featureAdaptor
  my $sql = "insert into supporting_feature (exon_id, feature_id, feature_type)
             values(?, ?, ?)";  
  
  my $sf_sth = $self->db->prepare($sql);

  my $anaAdaptor = $self->db->get_AnalysisAdaptor();
  my $dna_adaptor = $self->db->get_DnaAlignFeatureAdaptor();
  my $pep_adaptor = $self->db->get_ProteinAlignFeatureAdaptor();
  my $type;

  my @exons = ();
  if($exon->isa('Bio::EnsEMBL::StickyExon')) {
    @exons = @{$exon->get_all_component_Exons};
  } else {
    @exons = ($exon);
  }

  foreach my $e (@exons) {
    foreach my $sf (@{$e->get_all_supporting_features}) {
      unless($sf->isa("Bio::EnsEMBL::BaseAlignFeature")){
	$self->throw("$sf must be an align feature otherwise" .
		     "it can't be stored");
      }

      #sanity check
      eval { $sf->validate(); };
      if ($@) {
        $self->warn("Supporting feature invalid. Skipping feature\n");
	next;
      }

      $sf->contig($e->contig);

      if($sf->isa("Bio::EnsEMBL::DnaDnaAlignFeature")){
	$dna_adaptor->store($sf);
	$type = 'dna_align_feature';
      }elsif($sf->isa("Bio::EnsEMBL::DnaPepAlignFeature")){
	$pep_adaptor->store($sf);
	$type = 'protein_align_feature';
      } else {
	$self->warn("Supporting feature of unknown type. Skipping : [$sf]\n");
	next;
      }

      $sf_sth->execute($exonId, $sf->dbID, $type);
    }
  }

  #
  # Finally, update the dbID and adaptor of the exon (and any component exons)
  # to point to the new database
  #
  foreach my $e (@exons) {
    $e->dbID($exonId);
    $e->adaptor($self);
  }

  $exon->adaptor($self);
  $exon->dbID($exonId);
}


=head2 get_stable_entry_info

  Arg [1]    : Bio::EnsEMBL::Exon $exon
  Example    : $exon_adaptor->get_stable_entry_info($exon);
  Description: gets stable info for an exon. this is not usually done at
               creation time for speed purposes, and can be lazy-loaded later
               if it is needed..
  Returntype : none
  Exceptions : none
  Caller     : Bio::EnsEMBL::Exon

=cut

sub get_stable_entry_info {
  my ($self,$exon) = @_;

  if( !$exon || !ref $exon || !$exon->isa('Bio::EnsEMBL::Exon') ) {
     $self->throw("Needs a exon object, not a $exon");
  }
  if(!$exon->dbID){
    #$self->throw("can't fetch stable info with no dbID");
    return;
  }
  my $sth = $self->prepare("SELECT stable_id, UNIX_TIMESTAMP(created),
                                   UNIX_TIMESTAMP(modified), version 
                            FROM   exon_stable_id 
                            WHERE  exon_id = " . $exon->dbID);

  $sth->execute();

  # my @array = $sth->fetchrow_array();
  if( my $aref = $sth->fetchrow_arrayref() ) {
    $exon->{'_stable_id'} = $aref->[0];
    $exon->{'_created'}   = $aref->[1];
    $exon->{'_modified'}  = $aref->[2];
    $exon->{'_version'}   = $aref->[3];
  }

  return 1;
}


=head2 remove

  Arg [1]    : Bio::EnsEMBL::Exon $exon
               the exon to remove from the database 
  Example    : $exon_adaptor->remove($exon);
  Description: Removes an exon from the database
  Returntype : none
  Exceptions : none
  Caller     : general

=cut

sub remove {
  my $self = shift;
  my $exon = shift;
  
  if ( ! $exon->dbID() ) {
    return;
  }
  #print "have ".$self->db."\n";
  my $sth = $self->prepare( "delete from exon where exon_id = ?" );
  $sth->execute( $exon->dbID );
  #print "have deleted ".$exon->dbID."\n";
  $sth = $self->prepare( "delete from exon_stable_id where exon_id = ?" );
  $sth->execute( $exon->dbID );
  
  my $sql = "select feature_type, feature_id from supporting_feature where exon_id = ".$exon->dbID." ";

  #    print STDERR "sql = ".$sql."\n";
  $sth = $self->prepare($sql);
  
  $sth->execute;
  
  my $prot_adp = $self->db->get_ProteinAlignFeatureAdaptor;
  my $dna_adp = $self->db->get_DnaAlignFeatureAdaptor;
  
  while(my ($type, $feature_id) = $sth->fetchrow){
    
    if($type eq 'protein_align_feature'){
      my $f = $prot_adp->fetch_by_dbID($feature_id);
      $prot_adp->remove($f);
      #print "have removed ".$f->dbID."\n";
    }
    elsif($type eq 'dna_align_feature'){
      my $f = $dna_adp->fetch_by_dbID($feature_id);
      #print "have removed ".$f->dbID."\n";
      $dna_adp->remove($f);
    }
  }

  $sth = $self->prepare( "delete from supporting_feature where exon_id = ?" );
  $sth->execute( $exon->dbID );

  # uhh, didnt know another way of resetting to undef ...
  $exon->{dbID} = undef;
}

=head2 list_dbIDs

  Arg [1]    : none
  Example    : @exon_ids = @{$exon_adaptor->list_dbIDs()};
  Description: Gets an array of internal ids for all exons in the current db
  Returntype : list of ints
  Exceptions : none
  Caller     : ?

=cut

sub list_dbIDs {
   my ($self) = @_;

   return $self->_list_dbIDs("exon");
}

=head2 list_stable_ids

  Arg [1]    : none
  Example    : @stable_exon_ids = @{$exon_adaptor->list_stable_dbIDs()};
  Description: Gets an array of stable ids for all exons in the current db
  Returntype : list of ints
  Exceptions : none
  Caller     : ?

=cut

sub list_stable_ids {
   my ($self) = @_;

   return $self->_list_dbIDs("exon_stable_id", "stable_id");
}

1;
