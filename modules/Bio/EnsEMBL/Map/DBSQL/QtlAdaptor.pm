=head1 LICENSE

  Copyright (c) 1999-2012 The European Bioinformatics Institute and
  Genome Research Limited.  All rights reserved.

  This software is distributed under a modified Apache license.
  For license details, please see

    http://www.ensembl.org/info/about/code_licence.html

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <dev@ensembl.org>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk@ensembl.org>.

=cut

=head1 NAME

Bio::EnsEMBL::Map::DBSQL::QtlAdaptor

=head1 SYNOPSIS

=head1 DESCRIPTION

This module is responsible of retrieving QTLs from the database.

=head1 METHODS

=cut

package Bio::EnsEMBL::Map::DBSQL::QtlAdaptor;

use strict;

use Bio::EnsEMBL::Map::Qtl;
use Bio::EnsEMBL::DBSQL::BaseAdaptor;


use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::DBSQL::BaseAdaptor);


=head2 fetch_by_dbID

  Arg  1     : int $dbID
  Example    : none
  Description: get by database internal identifier
  Returntype : Bio::EnsEMBL::Map::Qtl
  Exceptions : none
  Caller     : general
  Status     : stable

=cut

sub fetch_by_dbID {
  my $self = shift;
  my $dbID = shift;

  return unless $dbID;
  
  my $res =  $self->_generic_fetch( [ "q.qtl_id = $dbID" ] );
  return $res->[0];
}


=head2 fetch_all

  Example    : none
  Description: get all the qtl's
  Returntype : listref Bio::EnsEMBL::Map::Qtl
  Exceptions : none
  Caller     : general
  Status     : stable

=cut

sub fetch_all {
  my $self = shift;
  $self->_generic_fetch( [] );
}




=head2 fetch_all_by_trait

  Arg [1]    : string $trait
               The phenotype we are looking for
  Example    : none
  Description: get by phenotype/trait string
  Returntype : listref Bio::EnsEMBL::Map::Qtl
  Exceptions : none
  Caller     : general
  Status     : stable

=cut


sub fetch_all_by_trait {
  my $self = shift;
  my $trait = shift;

  return [] unless $trait;

  return $self->_generic_fetch( [ "q.trait = '$trait'" ] );
}




=head2 fetch_all_by_source_database

  Arg  1     : string $database_name
               Name of the database that provides the Qtl information
  Arg [2]    : string $database_primary_id
               The primary id of the qtl in that database
  Example    : none
  Description: retrieve Qtl by given information 
  Returntype : listref Bio::EnsEMBL::Map::Qtl 
  Exceptions : none
  Caller     : general
  Status     : stable

=cut

sub fetch_all_by_source_database {

  my $self = shift;
  my $database_name = shift;
  my $database_primary_id = shift;
  
  return [] unless $database_name;

  my @conditions;

  if( $database_name ) {
    push( @conditions, "q.source_database=\"$database_name\"" );
  }

  if( $database_primary_id ) {
    push( @conditions, "q.source_primary_id=\"$database_primary_id\"" );
  }

  return $self->_generic_fetch( \@conditions );
}


sub _generic_fetch {
  my $self = shift;
  my $conditions = shift;

  my $where = '';

  if( @$conditions ) {
    $where = "WHERE ".join( " and ", @$conditions );
  }

  my $query = "SELECT ".
    join( ", ", $self->_columns() ). 
      " FROM qtl q LEFT JOIN qtl_synonym qs ON q.qtl_id = qs.qtl_id ".
	$where;

  my $sth = $self->prepare( $query );
  $sth->execute();
  
  return $self->_obj_from_sth( $sth );
}


sub _columns {
  return ( 'q.qtl_id','qs.source_database','qs.source_primary_id',
	   'q.trait','q.lod_score','q.flank_marker_id_1',
	   'q.flank_marker_id_2','q.peak_marker_id' );
}


sub _obj_from_sth {
  my $self = shift;
  my $sth = shift;
  
  my ( $qtl_id, $source_database,
       $source_primary_id, $trait, $lod_score, $flank_marker_id_1,
       $flank_marker_id_2, $peak_marker_id );

  #warning: ordering depends on _columns function implementation
  $sth->bind_columns( \$qtl_id, 
		      \$source_database, \$source_primary_id, \$trait, 
		      \$lod_score, \$flank_marker_id_1,
		      \$flank_marker_id_2, \$peak_marker_id );

  my @out = ();
  my %already_seen;

  while( $sth->fetch()) {

    #multiple columns with same qtl are multiple synonyms
    if(my $qtl = $already_seen{$qtl_id}) {
      $qtl->add_synonym($source_database, $source_primary_id);
      next;
    }

    my $mad = $self->db()->get_MarkerAdaptor();

    my $flank_marker_1 = $flank_marker_id_1 ? $mad->fetch_by_dbID( $flank_marker_id_1 ) : undef ;
    my $flank_marker_2 = $flank_marker_id_2 ? $mad->fetch_by_dbID( $flank_marker_id_2 ) : undef;
    my $peak_marker = $peak_marker_id ? $mad->fetch_by_dbID( $peak_marker_id ) : undef;
    
    my $qtl = Bio::EnsEMBL::Map::Qtl->new
      (
       $qtl_id,
       $self->db->get_QtlAdaptor(),
       $flank_marker_1,
       $peak_marker,
       $flank_marker_2,
       $trait, 
       $lod_score,
       {$source_database => $source_primary_id}
      );
    
    push @out, $qtl;
    $already_seen{$qtl_id} = $qtl;
  }

  return \@out;
}



=head2 list_traits

  Args       : none
  Example    : none
  Description: list of all the different traits
  Returntype : listref string
  Exceptions : none
  Caller     : general
  Status     : stable

=cut



sub list_traits {
  my $self = shift;
  
  my $sth = $self->prepare( "
   SELECT DISTINCT trait
             FROM  qtl q
  " );
  
  my $res = []; 

  $sth->execute();
  push ( @$res ,
	 map { $_->[0] } @{$sth->fetchall_arrayref()}
       );

  return $res;
}


 
			  


1;
