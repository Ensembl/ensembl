# EnsEMBL Gene reading writing adaptor for mySQL
#
# Copyright EMBL-EBI 2001
#
# Author: Graham McVicker
# 
# Date : 22.07.2002
#

=head1 NAME

Bio::EnsEMBL::DBSQL::ProxyGeneAdaptor

=head1 SYNOPSIS

Designed as an abstraction over the database specific GeneAdaptors.  The proxy
gene adaptor normally behaves just as a normal core GeneAdaptor, however, for
certain requests it may decide to instead forward the request to another
database (such as the lite database if it is available).

=head1 CONTACT

  Arne Stabenau: stabenau@ebi.ac.uk
  Ewan Birney  : birney@ebi.ac.uk
  Graham McVicker : mcvicker@ebi.ac.uk

=head1 APPENDIX

=cut

use strict;

package Bio::EnsEMBL::DBSQL::ProxyGeneAdaptor;

use Bio::EnsEMBL::DBSQL::ProxyAdaptor;

use vars '@ISA';

@ISA = qw(Bio::EnsEMBL::DBSQL::ProxyAdaptor);

#new inherited from ProxyAdaptor


=head2 fetch_by_Slice

  Arg [1]    : generic list @args
  Example    : @genes = $gene_adaptor->fetch_by_Slice($slice);
  Description: Overrides the fetch_by_Slice method and makes a decision
               about what database to use for the request.  If the Lite
               database is available then the lite database will be
               used for this request - it is faster and will also
               fill the transcript and exon objects contained by the gene.
               If the lite database is not available then the primary database
               will be used.
  Returntype : Bio::EnsEMBL::Gene
  Exceptions : none
  Caller     : general

=cut

sub fetch_all_by_Slice {
  my ($self, @args) = @_;

  my $lite_db = $self->db()->get_db_adaptor('lite');

  if(defined $lite_db) {
    #use the Lite database if it is available

    return $lite_db->get_GeneAdaptor()->fetch_all_by_Slice(@args);
  }

  #otherwise use the core database
  return $self->{'_primary_adaptor'}->fetch_all_by_Slice(@args);
}



=head2 fetch_by_transcript_stable_id

  Arg [1]    : list of arbitrary args @args
  Example    : none
  Description: The proxy overrides this method and automatically calls
               the lite GeneAdaptor if it is available for greater speed.
                If it is not available than the core adaptor is used
  Returntype : Bio::EnsEMBL::Gene
  Exceptions : none
  Caller     : general

=cut

sub fetch_by_transcript_stable_id {
  my ($self, @args) = @_;

  my $lite_db = $self->db()->get_db_adaptor('lite');
  
  if(defined $lite_db) {
    #use the Lite database if it is available
    return $lite_db->get_GeneAdaptor()->fetch_by_transcript_stable_id(@args);
  }
  
  #otherwise use the core database
  return $self->{'_primary_adaptor'}->fetch_by_transcript_stable_id(@args);
}


=head2 fetch_by_stable_id

  Arg [1]    : list of arbitrary args @args
  Example    : none
  Description: The proxy overrides this method and automatically calls
               the lite GeneAdaptor if it is available for greater speed.
                If it is not available than the core adaptor is used
  Returntype : Bio::EnsEMBL::Gene
  Exceptions : none
  Caller     : general

=cut

sub fetch_by_stable_id{
  my ($self, @args) = @_;

  my $lite_db = $self->db()->get_db_adaptor('lite');
  
  if(defined $lite_db) {
    #use the Lite database if it is available
    return $lite_db->get_GeneAdaptor()->fetch_by_stable_id(@args);
  }
  
  #otherwise use the core database
  return $self->{'_primary_adaptor'}->fetch_by_stable_id(@args);
}

1;
