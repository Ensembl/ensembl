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

use Bio::EnsEMBL::DBSQL::GeneAdaptor;
use Bio::EnsEMBL::DBSQL::DBAdaptor;

use vars '@ISA';

@ISA = ('Bio::EnsEMBL::DBSQL::GeneAdaptor');


#use inherited superclass constructor


sub fetch_by_Slice {
  my ($self, @args) = @_;

  my $lite_db = $self->db()->lite_DBAdaptor();
 
  if(defined $lite_db) {
    #use the Lite database if it is available
    return $lite_db->get_GeneAdaptor()->fetch_by_Slice(@args);
  }

  #otherwise use the core database
  return $self->SUPER::fetch_by_Slice(@args);
}



