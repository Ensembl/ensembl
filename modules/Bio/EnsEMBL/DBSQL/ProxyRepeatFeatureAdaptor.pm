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
The proxy gene adaptor normally behaves just as a normal core 
RepeatFeatureAdaptor, however, for certain requests it may decide to 
instead forward the request to another
database (such as the lite database if it is available).

=head1 CONTACT

  Arne Stabenau: stabenau@ebi.ac.uk
  Ewan Birney  : birney@ebi.ac.uk
  Graham McVicker : mcvicker@ebi.ac.uk

=head1 APPENDIX

=cut

use strict;

use Bio::EnsEMBL::DBSQL::RepeatFeatureAdaptorI;
use Bio::EnsEMBL::DBSQL::BaseAdaptor;

package Bio::EnsEMBL::DBSQL::ProxyRepeatFeatureAdaptor;

use vars '@ISA';

@ISA = qw(Bio::EnsEMBL::DBSQL::RepeatFeatureAdaptorI 
	  Bio::EnsEMBL::DBSQL::BaseAdaptor);

#implement the interface ProxyRepeatFeatureAdaptorI
use implements qw(Bio::EnsEMBL::DBSQL::RepeatFeatureAdaptorI);

sub new {
  my($class, $db, $core_adaptor) = @_;

  #call superclass constructor
  my $self = $class->SUPER::new($db);
  $self->{'_core_adaptor'} = $core_adaptor;
  return $self;
}

sub fetch_by_Slice {
  my ($self, @args) = @_;

  my $lite_db = $self->db()->lite_DBAdaptor();

  if(defined $lite_db) {
    #use the lite database if it is available
    return $lite_db->get_RepeatFeatureAdaptor()->fetch_by_Slice(@args);
  } 

  #otherwise use the core database
  return $self->{'_core_adaptor'}->fetch_by_Slice(@args);
}
  
sub fetch_by_RawContig {
  my ($self, @args) = @_;

  return $self->{'_core_adaptor'}->fetch_by_RawContig(@args);
}

sub fetch_by_contig_id {
  my ($self, @args) = @_;

  return $self->{'_core_adaptor'}->fetch_by_contig_id(@args);
}

sub fetch_by_dbID {
  my ($self, @args) = @_;

  return $self->{'_core_adaptor'}->fetch_by_dbID(@args);
}

sub fetch_by_assembly_location {
  my ($self, @args) = @_;

  return $self->{'_core_adaptor'}->fetch_by_assembly_location(@args);
}

sub fetch_by_assembly_location_constraint {
  my ($self, @args) = @_;

 return $self->{'_core_adaptor'}->fetch_by_assembly_location_constraint(@args);
}

sub store {
  my ($self, @args) = @_;

  return $self->{'_core_adaptor'}->store(@args);
}

