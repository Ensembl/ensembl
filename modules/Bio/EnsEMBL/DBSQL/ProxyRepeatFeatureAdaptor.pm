# EnsEMBL Proxy reading writing adaptor for mySQL
#
# Copyright EMBL-EBI 2001
#
# Author: Graham McVicker
# 
# Date : 12.08.2002
#

=head1 NAME

Bio::EnsEMBL::DBSQL::ProxyRepeatFeatureAdaptor

=head1 SYNOPSIS

Designed as an abstraction over the database specific RepeatFeatureAdaptors.  
The proxy repeat feature adaptor normally behaves just as a normal core 
RepeatFeatureAdaptor, however, for certain requests it may decide to 
instead forward the request to another database (such as the lite database 
if it is available).

=head1 CONTACT

  Arne Stabenau: stabenau@ebi.ac.uk
  Ewan Birney  : birney@ebi.ac.uk
  Graham McVicker : mcvicker@ebi.ac.uk

=head1 APPENDIX

=cut

use strict;


package Bio::EnsEMBL::DBSQL::ProxyRepeatFeatureAdaptor;

use Bio::EnsEMBL::DBSQL::ProxyAdaptor;
use vars '@ISA';

@ISA = qw(Bio::EnsEMBL::DBSQL::ProxyAdaptor);



=head2 fetch_all_by_Slice

  Arg [1]    : arbitrary list of args @args 
  Example    : none
  Description: Forwards requests for fetch_all_by_Slice to the lite database if
               is available.  This is done for improved performance.
  Returntype : listref of Bio::EnsEMBL::Gene
  Exceptions : none
  Caller     : general

=cut

sub fetch_all_by_Slice {
  my ($self, @args) = @_;

  my $lite_db = $self->db()->get_db_adaptor('lite');

  if(defined $lite_db) {
    #use the lite database if it is available
    return $lite_db->get_RepeatFeatureAdaptor()->fetch_all_by_Slice(@args);
  } 

  #otherwise use the core database
  return $self->{'_primary_adaptor'}->fetch_all_by_Slice(@args);
}


1;

