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

Bio::EnsEMBL::Map::DBSQL::QtlFeatureAdaptor

=head1 SYNOPSIS


=head1 DESCRIPTION

This module is responsible of retrieving QtlFeatures (and their
associated Qtls) from the database.

The bulk of this objects' methods are inherited from
Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor

=head1 METHODS

=cut

package Bio::EnsEMBL::Map::DBSQL::QtlFeatureAdaptor;

use strict;

use Bio::EnsEMBL::Map::Qtl;
use Bio::EnsEMBL::Map::QtlFeature;
use Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor;

use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor);



=head2 fetch_all_by_Qtl

  Arg [1]    : Bio::EnsEMBL::Map::Qtl
  Example    : none
  Description: Retrieves a list of QtlFeatures for a given Qtl
  Returntype : listref of Bio::EnsEMBL::QtlFeatures
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub fetch_all_by_Qtl {
  my $self = shift;
  my $qtl = shift;

  my $constraint = 'q.qtl_id = ' . $qtl->dbID;

  return $self->generic_fetch($constraint, @_);
}




sub _columns {
  my $self = shift;

  return ( 'qf.seq_region_id', 'qf.seq_region_start', 'qf.seq_region_end',
           'q.qtl_id',
           'qf.analysis_id',
           'qs.source_database', 'qs.source_primary_id',
           'q.trait', 'q.lod_score', 'q.flank_marker_id_1',
           'q.flank_marker_id_2', 'q.peak_marker_id' );
}

sub _tables {
  my $self = shift;

  return (['qtl_feature', 'qf'], #primary table
	        ['qtl', 'q'],
          ['qtl_synonym', 'qs']);
}

sub _left_join {
  return ( [ 'qtl_synonym', 'q.qtl_id = qs.qtl_id' ] );
}

sub _default_where_clause {
  my $self = shift;

  return ('qf.qtl_id = q.qtl_id');
}


sub _objs_from_sth {
  my $self = shift;
  my $sth  = shift;

  my ( $seq_region_id, $seq_region_start, $seq_region_end, $qtl_id,
       $analysis_id, $source_database,
       $source_primary_id, $trait, $lod_score, $flank_marker_id_1,
       $flank_marker_id_2, $peak_marker_id );

  #warning: ordering depends on _columns function implementation
  $sth->bind_columns( \$seq_region_id, \$seq_region_start, \$seq_region_end, 
                      \$qtl_id, \$analysis_id,
                      \$source_database, \$source_primary_id, \$trait,
                      \$lod_score, \$flank_marker_id_1,
                      \$flank_marker_id_2, \$peak_marker_id );

  my @out = ();
  my %already_seen;

  my $mad = $self->db()->get_MarkerAdaptor();
  my $aad = $self->db()->get_AnalysisAdaptor();
  my $sad = $self->db()->get_SliceAdaptor();

  while( $sth->fetch()) {

    my $flank_marker_1 = $flank_marker_id_1 ?
                         $mad->fetch_by_dbID( $flank_marker_id_1 ) :
                         undef;
    my $flank_marker_2 = $flank_marker_id_2 ?
                         $mad->fetch_by_dbID( $flank_marker_id_2 ) :
                         undef;
    my $peak_marker = $peak_marker_id ?
                      $mad->fetch_by_dbID( $peak_marker_id ) :
                      undef;

    my $analysis = $aad->fetch_by_dbID( $analysis_id );

    my $slice = $sad->fetch_by_seq_region_id($seq_region_id);

    #rows with the same qtl contain additional synonyms of the qtl
    if(my $qtl = $already_seen{$qtl_id}) {
      $qtl->add_synonym($source_database, $source_primary_id);
      next;
    }

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

    $already_seen{$qtl_id} = $qtl;

    #now create a new marker_feature using the marker
    push @out, Bio::EnsEMBL::Map::QtlFeature->new
      ($self,
       $slice,
       $seq_region_start,
       $seq_region_end,
       $qtl,
       $analysis);
  }

  return \@out;
}


1;
