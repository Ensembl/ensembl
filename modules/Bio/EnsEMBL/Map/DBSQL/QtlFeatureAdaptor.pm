# EnsEMBL module for QtlFeatureAdaptor
# Copyright EMBL-EBI/Sanger center 2003
#
#
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::Map::DBSQL::QtlFeatureAdaptor

=head1 SYNOPSIS


=head1 DESCRIPTION

This module is responsible of retrieving QtlFeatures (and their associated Qtls)
from the database.

The bulk of this objects methods are inherited from 
Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor


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

=cut

sub fetch_all_by_Qtl {
  my $self = shift;
  my $qtl = shift;
  
  my $constraint = 'q.qtl_id = ' . $qtl->dbID;
  
  return $self->generic_fetch($constraint, @_);
}




sub _columns {
  my $self = shift;

  return ( 'c.name', 'qf.start', 'qf.end', 'q.qtl_id', 'qf.analysis_id',
	   'qs.source_database', 'qs.source_primary_id',
	   'q.trait', 'q.lod_score', 'q.flank_marker_id_1',
	   'q.flank_marker_id_2', 'q.peak_marker_id' );
}

sub _tables {
  my $self = shift;

  return (['qtl_feature', 'qf'], #primary table
	  ['qtl', 'q'],
	  ['chromosome' ,'c' ],
          ['qtl_synonym', 'qs']);
}

sub _left_join {
  return ( [ 'qtl_synonym', 'q.qtl_id = qs.qtl_id' ] );
}
          
sub _default_where_clause {
  my $self = shift;

  return ('qf.qtl_id = q.qtl_id AND qf.chromosome_id = c.chromosome_id');
}


sub _objs_from_sth {
  my $self = shift;
  my $sth  = shift;
  my $mapper = shift;
  my $slice = shift; # optional Slice

  
    
  my ( $chromosome_name, $chr_start, $chr_end, $qtl_id, $analysis_id, 
       $source_database,
       $source_primary_id, $trait, $lod_score, $flank_marker_id_1,
       $flank_marker_id_2, $peak_marker_id );

  #warning: ordering depends on _columns function implementation
  $sth->bind_columns( \$chromosome_name, \$chr_start, \$chr_end, \$qtl_id, 
		      \$analysis_id,
		      \$source_database, \$source_primary_id, \$trait, 
		      \$lod_score, \$flank_marker_id_1,
		      \$flank_marker_id_2, \$peak_marker_id );

  my @out = ();
  my %already_seen;
  while( $sth->fetch()) {

    my $mad = $self->db()->get_MarkerAdaptor();

    my $flank_marker_1 = $flank_marker_id_1 ? $mad->fetch_by_dbID( $flank_marker_id_1 ) : undef;
    my $flank_marker_2 = $flank_marker_id_2 ? $mad->fetch_by_dbID( $flank_marker_id_2 ) : undef;
    my $peak_marker = $peak_marker_id ? $mad->fetch_by_dbID( $peak_marker_id ) : undef;
    
    my $analysis = $self->db()->get_AnalysisAdaptor()->fetch_by_dbID( $analysis_id );

    if( ! $slice ) {
      $slice = $self->db->get_SliceAdaptor()->
	fetch_by_chr_name( $chromosome_name );
    }

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
    if( $slice->strand == 1 ) {
      push @out, Bio::EnsEMBL::Map::QtlFeature->new
	(
	 $self,
	 $slice,
	 $chr_start - $slice->chr_start() + 1,
	 $chr_end - $slice->chr_start() + 1,
	 $qtl,
	 $analysis
	);
    } else {
      push @out, Bio::EnsEMBL::Map::QtlFeature->new
	(
	 $self,
	 $slice,
	 $slice->chr_end() - $chr_end + 1,
	 $slice->chr_end() - $chr_start + 1,
	 $qtl,
	 $analysis
	);
    }
  }

  return \@out;
}

=head2 fetch_all_by_Slice_constraint

  Arg [1]    : Bio::EnsEMBL::Slice $slice
               the slice from which to obtain features
  Arg [2]    : (optional) string $constraint
               An SQL query constraint (i.e. part of the WHERE clause)
  Arg [3]    : (optional) string $logic_name
               the logic name of the type of features to obtain
  Example    : $fs = $a->fetch_all_by_Slice_constraint($slc, 'perc_ident > 5');
  Description: Returns a listref of features created from the database which 
               are on the Slice defined by $slice and fulfill the SQL 
               constraint defined by $constraint. If logic name is defined, 
               only features with an analysis of type $logic_name will be 
               returned. 
  Returntype : listref of Bio::EnsEMBL::SeqFeatures in Slice coordinates
  Exceptions : thrown if $slice is not defined
  Caller     : Bio::EnsEMBL::Slice

=cut

sub fetch_all_by_Slice_constraint {
  my($self, $slice, $constraint, $logic_name) = @_;

  unless(defined $slice && ref $slice && $slice->isa("Bio::EnsEMBL::Slice")) {
    $self->throw("Slice arg must be a Bio::EnsEMBL::Slice not a [$slice]\n");
  }

  $logic_name = '' unless $logic_name;
  $constraint = '' unless $constraint;

  #check the cache and return if we have already done this query
  my $key = join($slice->name, $constraint, $logic_name);
  return $self->{'_slice_feature_cache'}{$key} 
    if $self->{'_slice_feature_cache'}{$key};
    
  my $slice_start  = $slice->chr_start();
  my $slice_end    = $slice->chr_end();

  #get the synonym of the primary_table
  my @tabs = $self->_tables;
  my $syn = $tabs[0]->[1];
  
  my $chromosome_id = $slice->get_Chromosome->dbID();

  $constraint .= ' AND ' if($constraint);

  $constraint .= " ${syn}.chromosome_id = $chromosome_id " .
    "AND ${syn}.end >= $slice_start AND ${syn}.start <= $slice_end";
  
  #for speed the remapping to slice may be done at the time of object creation
  my $features = 
    $self->generic_fetch($constraint, $logic_name, undef, $slice); 
  
  if(@$features && (!$features->[0]->can('contig') || 
		    $features->[0]->contig == $slice)) {
    #features have been converted to slice coords already, cache and return
    return $self->{'_slice_feature_cache'}{$key} = $features;
  }
}




1;
