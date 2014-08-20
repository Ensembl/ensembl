=head1 LICENSE

Copyright [1999-2014] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=cut


=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=cut

=head1 NAME

Bio::EnsEMBL::DBSQL::SlicingTranscriptPairAdaptor - Database adaptor for the retrieval and
storage of SplicingTranscriptPair objects

=head1 SYNOPSIS

  use Bio::EnsEMBL::Registry;

  Bio::EnsEMBL::Registry->load_registry_from_db(
    -host => 'ensembldb.ensembl.org',
    -user => 'anonymous',
  );

  $se_adaptor =
    Bio::EnsEMBL::Registry->get_adaptor( "human", "core", "SplicingTranscriptPair" );


=head1 DESCRIPTION

This is a database aware adaptor for the retrieval and storage of SplicingTranscriptPairs
objects.

=head1 METHODS

=cut
package Bio::EnsEMBL::DBSQL::SplicingTranscriptPairAdaptor;

use strict;

use Bio::EnsEMBL::Utils::Exception qw( deprecate throw warning );
use Bio::EnsEMBL::DBSQL::SliceAdaptor;
use Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::SplicingTranscriptPair;

use vars '@ISA';
@ISA = qw(Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor);



sub fetch_all_by_SplicingEvent{
  my ($self, $splicing_event) = @_;
  

  my ($splicing_transcript_pair_id, $transcript_id_1, $transcript_id_2);

  my $splicing_event_id = $splicing_event->dbID;
  
  my $sql = "select splicing_transcript_pair_id, transcript_id_1, transcript_id_2 from splicing_transcript_pair where splicing_event_id = $splicing_event_id";

  my $sth = $self->prepare($sql);

  $sth->execute();
  $sth->bind_columns(\$splicing_transcript_pair_id, \$transcript_id_1, \$transcript_id_2);

  my @pairs;

  while($sth->fetch()){
    push( @pairs,
          $self->_create_feature_fast( 'Bio::EnsEMBL::SplicingTranscriptPair', {
                                    'adaptor'   => $self,
                                    'dbID'      => $splicing_transcript_pair_id,
                                    'transcript_id_1'   => $transcript_id_1,
                                    'transcript_id_2'   => $transcript_id_2,
				    'start'             => $splicing_event->start,
                                    'end'               => $splicing_event->end
                                  } ) );
    
  }
  $sth->finish;
  return \@pairs;

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

  return ([ 'splicing_transcript_pair', 'stp' ]);
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

  return ( 'stp.splicing_transcript_pair_id','stp.transcript_id_1', 'stp.transcript_id_2');

}

sub list_dbIDs {
  my ($self,$ordered) = @_;

  return $self->_list_dbIDs("splicing_transcript_pair", undef, $ordered);
}


1;


