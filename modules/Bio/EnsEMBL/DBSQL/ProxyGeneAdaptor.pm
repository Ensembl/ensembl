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

sub fetch_by_Slice {
  my ($self, @args) = @_;

  my $lite_db = $self->db()->get_db_adaptor('lite');
 
  if(defined $lite_db) {
    #use the Lite database if it is available
    return $lite_db->get_GeneAdaptor()->fetch_by_Slice(@args);
  }

  #otherwise use the core database
  return $self->{'_primary_adaptor'}->fetch_by_Slice(@args);
}


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


sub fetch_by_stable_id{
  my ($self, @args) = @_;
  print STDERR ( "ProxyGeneAdaptor is called.\n" );

  my $lite_db = $self->db()->get_db_adaptor('lite');
  
  if(defined $lite_db) {
    #use the Lite database if it is available
    print STDERR "Lite Database used for fetch_by_stable_id.\n" ; 
    return $lite_db->get_GeneAdaptor()->fetch_by_stable_id(@args);
  }
  
  #otherwise use the core database
  return $self->{'_core_adaptor'}->fetch_by_stable_id(@args);
}

