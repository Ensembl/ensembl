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


use Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor;
use Bio::EnsEMBL::Exon;
use Bio::EnsEMBL::Utils::Exception qw( warning throw deprecate );
 
@ISA = qw( Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor );



=head2 _tablename

  Arg [1]    : none
  Example    : none
  Description: PROTECTED implementation of superclass abstract method
               returns the names, aliases of the tables to use for queries
  Returntype : list of listrefs of strings
  Exceptions : none
  Caller     : internal

=cut

sub _tables {
  my $self = shift;
  
  return ([ 'exon_transcript', 'et' ],[ 'exon', 'e' ], 
	  [ 'exon_stable_id', 'esi' ] );

}



=head2 _default_where_clause

  Arg [1]    : none
  Example    : none
  Description: PROTECTED implementation of superclass abstract method
               Makes join between et and e table
  Returntype : string
  Exceptions : none
  Caller     : generic_fetch

=cut

sub _default_where_clause {
  return "et.exon_id = e.exon_id";
}


=head2 _columns

  Arg [1]    : none
  Example    : none
  Description: PROTECTED implementation of superclass abstract method
               returns a list of columns to use for queries
  Returntype : list of strings
  Exceptions : none
  Caller     : internal

=cut

sub _columns {
  my $self = shift;

  return qw( e.exon_id e.seq_region_id e.seq_region_start e.seq_region_end 
	     e.seq_region_strand e.phase e.end_phase
	     esi.stable_id esi.version );
}


sub _left_join {
  return ( [ 'exon_stable_id', "esi.exon_id = e.exon_id" ] );
}



=head2 _final_clause

  Arg [1]    : none
  Example    : none
  Description: PROTECTED implementation of superclass abstract method
               returns a default end for the SQL-query (ORDER BY)
  Returntype : string
  Exceptions : none
  Caller     : internal

=cut

sub _final_clause {
  return "ORDER BY et.transcript_id, et.rank";
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
  my ( $self, $stable_id, $cs_name, $cs_version ) = @_;

  my $constraint = "esi.stable_id = \"$stable_id\"";

  # should be only one :-)
  my $exons = $self->SUPER::generic_fetch( $constraint );

  if( ! @$exons ) { return undef }

  my @new_exons = map { $_->transform( $cs_name, $cs_version ) } @$exons;

  return $new_exons[0];
}


=head2 fetch_all_by_Transcript

  Arg [1]    : Bio::EnsEMBL::Transcript $transcript
  Example    : none
  Description: Retrieves all Exons for the Transcript in 5-3 order
  Returntype : listref Bio::EnsEMBL::Exon on Transcript slice 
  Exceptions : none
  Caller     : Transcript->get_all_Exons()

=cut

sub fetch_all_by_Transcript {
  my ( $self, $transcript ) = @_;

  my $constraint = "et.transcript_id = ".$transcript->dbID();

  my $exons = $self->SUPER::generic_fetch( $constraint );

  if( ! @$exons ) { return [] }

  my @new_exons = map { $_->transfer( $transcript->slice() ) } @$exons;

  return \@new_exons;
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
    $self->throw("$exon is not a EnsEMBL exon - not dumping!");
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
    INSERT into exon ( seq_region_id, seq_region_start, 
		       seq_region_end, seq_region_strand, phase, 
		       end_phase )
    VALUES ( ?, ?, ?, ?, ?, ? ) 
  };
  my $exonst = $self->prepare($exon_sql);

  my $exonId = undef;
  # normal storing

  my $slice = $exon->slice();

  unless( $slice && ref $slice ) {
    $self->throw("Exon does not have an attached slice with a valid " . 
		 "database id.  Needs to have one set");
  }

  my $seq_region_id = $self->db()->get_SliceAdaptor()->
    get_seq_region_id( $exon->slice() );

  if( ! $seq_region_id ) {
    throw( "Attached slice is not valid in database" );
  }

  $exonst->execute( $seq_region_id,
		    $exon->start(),
		    $exon->end(),
		    $exon->strand(),
		    $exon->phase(),
		    $exon->end_phase());
  $exonId = $exonst->{'mysql_insertid'};

  
  if ($exon->stable_id && $exon->version()) {

    my $statement = "
      INSERT INTO exon_stable_id
      SET version = ?,
          stable_id = ?,
          exon_id = ? ";
    
    my $sth = $self->prepare($statement);
    $sth->execute( $exon->version, $exon->stable_id, $exonId );
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

  foreach my $sf (@{$exon->get_all_supporting_features}) {
    unless($sf->isa("Bio::EnsEMBL::BaseAlignFeature")){
      $self->throw("$sf must be an align feature otherwise" .
		   "it can't be stored");
    }

    #sanity check
    
    $sf->slice($exon->slice());

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

  #
  # Finally, update the dbID and adaptor of the exon (and any component exons)
  # to point to the new database
  #

  $exon->adaptor($self);
  $exon->dbID($exonId);
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

=head2 _obj_from_hashref

  Arg [1]    : Hashreference $hashref
  Example    : none 
  Description: PROTECTED implementation of abstract superclass method.
               responsible for the creation of Genes 
  Returntype : listref of Bio::EnsEMBL::Genes in target coordinate system
  Exceptions : none
  Caller     : internal

=cut

sub _objs_from_sth {
  my ($self, $sth, $mapper, $dest_slice) = @_;

  #
  # This code is ugly because an attempt has been made to remove as many
  # function calls as possible for speed purposes.  Thus many caches and
  # a fair bit of gymnastics is used.
  #

  my $sa = $self->db()->get_SliceAdaptor();

  my @exons;
  my %rc_hash;
  my %analysis_hash;
  my %slice_hash;
  my %sr_name_hash;
  my %sr_cs_hash;


  my ( $exon_id, $seq_region_id, $seq_region_start,
       $seq_region_end, $seq_region_strand, $phase,
       $end_phase, $stable_id, $version );


  $sth->bind_columns(  \$exon_id, \$seq_region_id, 
        \$seq_region_start,
        \$seq_region_end, \$seq_region_strand, \$phase,
        \$end_phase, \$stable_id, \$version );

  my $asm_cs;
  my $cmp_cs;
  my $asm_cs_vers;
  my $asm_cs_name;
  my $cmp_cs_vers;
  my $cmp_cs_name;
  if($mapper) {
    $asm_cs = $mapper->assembled_CoordSystem();
    $cmp_cs = $mapper->component_CoordSystem();
    $asm_cs_name = $asm_cs->name();
    $asm_cs_vers = $asm_cs->version();
    $cmp_cs_name = $cmp_cs->name();
    $asm_cs_vers = $cmp_cs->version();
  }

  my $dest_slice_start;
  my $dest_slice_end;
  my $dest_slice_strand;
  my $dest_slice_length;
  if($dest_slice) {
    $dest_slice_start  = $dest_slice->start();
    $dest_slice_end    = $dest_slice->end();
    $dest_slice_strand = $dest_slice->strand();
    $dest_slice_length = $dest_slice->length();
  }

  FEATURE: while($sth->fetch()) {
    my $slice = $slice_hash{"ID:".$seq_region_id};

    if(!$slice) {
      $slice = $sa->fetch_by_seq_region_id($seq_region_id);
      $slice_hash{"ID:".$seq_region_id} = $slice;
      $sr_name_hash{$seq_region_id} = $slice->seq_region_name();
      $sr_cs_hash{$seq_region_id} = $slice->coord_system();
    }

    #
    # remap the feature coordinates to another coord system 
    # if a mapper was provided
    #
    if($mapper) {
      my $sr_name = $sr_name_hash{$seq_region_id};
      my $sr_cs   = $sr_cs_hash{$seq_region_id};

      ($sr_name,$seq_region_start,$seq_region_end,$seq_region_strand) =
        $mapper->fastmap($sr_name, $seq_region_start, $seq_region_end,
			 $seq_region_strand, $sr_cs);

      #skip features that map to gaps or coord system boundaries
      next FEATURE if(!defined($sr_name));

      #get a slice in the coord system we just mapped to
      if($asm_cs == $sr_cs || ($asm_cs != $sr_cs && $asm_cs->equals($sr_cs))) {
        $slice = $slice_hash{"NAME:$sr_name:$cmp_cs_name:$cmp_cs_vers"} ||=
          $sa->fetch_by_region($cmp_cs_name, $sr_name,undef, undef, undef,
                               $cmp_cs_vers);
      } else {
        $slice = $slice_hash{"NAME:$sr_name:$asm_cs_name:$asm_cs_vers"} ||=
          $sa->fetch_by_region($asm_cs_name, $sr_name, undef, undef, undef,
                               $asm_cs_vers);
      }
    }

    #
    # If a destination slice was provided convert the coords
    # If the dest_slice starts at 1 and is foward strand, nothing needs doing
    #
    if($dest_slice && ($dest_slice_start != 1 || $dest_slice_strand != 1)) {
      if($dest_slice_strand == 1) {
        $seq_region_start = $seq_region_start - $dest_slice_start + 1;
        $seq_region_end   = $seq_region_end   - $dest_slice_start + 1;
      } else {
        my $tmp_seq_region_start = $seq_region_start;
        $seq_region_start = $dest_slice_end - $seq_region_end + 1;
        $seq_region_end   = $dest_slice_end - $tmp_seq_region_start + 1;
        $seq_region_strand *= -1;
      }

      $slice = $dest_slice;

      #throw away features off the end of the requested slice
      if($seq_region_end < 1 || $seq_region_start > $dest_slice_length) {
        next FEATURE;
      }
    }


    #finally, create the new repeat feature
    push @exons, Bio::EnsEMBL::Exon->new
      ( '-start'         =>  $seq_region_start,
	'-end'           =>  $seq_region_end,
	'-strand'        =>  $seq_region_strand,
	'-adaptor'       =>  $self,
	'-slice'         =>  $slice,
	'-dbID'          =>  $exon_id,
        '-stable_id'     =>  $stable_id,
        '-version'       =>  $version,
        '-phase'         =>  $phase,
        '-end_phase'     =>  $end_phase )

  }

  return \@exons;
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

  deprecated( "This method call shouldnt be necessary" );

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

  deprecated( "Hopefully this method is not needed any more. Exons should be fetched by Transcript" );

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


1;
