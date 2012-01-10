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

Bio::EnsEMBL::DBSQL::SlicingEventFeatureAdaptor - Database adaptor for the retrieval and
storage of SplicingEventFeature objects

=head1 SYNOPSIS

  use Bio::EnsEMBL::Registry;

  Bio::EnsEMBL::Registry->load_registry_from_db(
    -host => 'ensembldb.ensembl.org',
    -user => 'anonymous',
  );

  $se_adaptor =
    Bio::EnsEMBL::Registry->get_adaptor( "human", "core", "SplicingEventFeature" );


=head1 DESCRIPTION

This is a database aware adaptor for the retrieval and storage of SlicingEventFeatures
objects.

=head1 METHODS

=cut
package Bio::EnsEMBL::DBSQL::SplicingEventFeatureAdaptor;

use strict;

use Bio::EnsEMBL::Utils::Exception qw( deprecate throw warning );
use Bio::EnsEMBL::DBSQL::SliceAdaptor;
use Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::SplicingEventFeature;

use vars '@ISA';
@ISA = qw(Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor);



sub fetch_all_by_SplicingEvent{
  my ($self, $splicing_event) = @_;
  

  my ($splicing_event_feature_id, $splicing_event_id, $exon_id, $transcript_id, $feature_order, $transcript_association, $type, $start, $end);

  $splicing_event_id = $splicing_event->dbID;
  
  my $sql = "select splicing_event_feature_id, exon_id, transcript_id, feature_order, transcript_association, type, start, end from splicing_event_feature where splicing_event_id = $splicing_event_id";

  my $sth = $self->prepare($sql);

  $sth->execute();
  $sth->bind_columns(\$splicing_event_feature_id, \$exon_id, \$transcript_id, \$feature_order, \$transcript_association, \$type, \$start, \$end);

  my @features;

  while($sth->fetch()){
    push( @features,
          $self->_create_feature_fast( 'Bio::EnsEMBL::SplicingEventFeature', {
                                    'start'     => $start,
                                    'end'       => $end,
                                    'adaptor'   => $self,
                                    'dbID'      => $splicing_event_feature_id,
                                    'exon_id'   => $exon_id,
                                    'transcript_id' => $transcript_id,
                                    'slice'     => $splicing_event->slice,
                                    'type'      => $type,
				    'feature_order'            => $feature_order,
                                    'transcript_association'   => $transcript_association 
                                  } ) );
    
  }
  $sth->finish;
  return \@features;

}




# _tables
#  Arg [1]    : none
#  Description: PROTECTED implementation of superclass abstract method.
#               Returns the names, aliases of the tables to use for queries.
#  Returntype : list of listrefs of strings
#  Exceptions : none
#  Caller     : internal
#  Status     : At Risk

sub _tables {
  my $self = shift;

  return ([ 'splicing_event_feature', 'sef' ]);
}

# _columns
#  Arg [1]    : none
#  Example    : none
#  Description: PROTECTED implementation of superclass abstract method.
#               Returns a list of columns to use for queries.
#  Returntype : list of strings
#  Exceptions : none
#  Caller     : internal
#  Status     : At Risk

sub _columns {
  my $self = shift;

#  my $created_date = $self->db->dbc->from_date_to_seconds("gsi.created_date");
#  my $modified_date = $self->db->dbc->from_date_to_seconds("gsi.modified_date");

  return ( 'sef.splicing_event_id','sef.exon_id', 'sef.feature_order', 'sef.transcript_association', 'sef.type', 'sef.start', 'sef.end' );

}

sub list_dbIDs {
  my ($self,$ordered) = @_;

  return $self->_list_dbIDs("splicing_event_feature", undef, $ordered);
}


1;


