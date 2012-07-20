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

Bio::EnsEMBL::DBSQL::Support::FullIdCache - ID based caching using all available values  

=head1 SYNOPSIS

  my $cache = Bio::EnsEMBL::DBSQL::Support::FullIdCache->new($adaptor);
  my $obj = $cache->get(21);

=head1 DESCRIPTION

An implementation of caching which uses a raw hash to hold all available
values from an adaptor. Useful for working with a controlled vocabulary
table where cardinality is low.

=head1 METHODS

=cut

package Bio::EnsEMBL::DBSQL::Support::FullIdCache;

use strict;
use warnings;
use base qw/Bio::EnsEMBL::DBSQL::Support::BaseCache/;
use Bio::EnsEMBL::Utils::Exception qw/throw/;

=head2 build_cache

  Description: Builds a cache keyed by dbID and populated from a call
               using C<generic_fetch()>
  Returntype : Hash
  Exceptions : None
  Caller     : BaseAdaptors
  Status     : Beta

=cut

sub build_cache {
  my ($self) = @_;
  my $adaptor = $self->adaptor();
  my %cache; 
  my $objs = $adaptor->generic_fetch();
  foreach my $obj (@{$objs}) {
    $cache{$obj->dbID()} = $obj;
  }
  return \%cache;
}

=head2 clear_cache

  Description: Delegates to C<delete_cache()> in order to clear all values
               and on the next cache request will force a C<build_cache()> 
               call
  Returntype : None
  Exceptions : None
  Caller     : BaseAdaptors
  Status     : Beta

=cut

sub clear_cache {
  my ($self) = @_;
  $self->delete_cache();
  return;
}

=head2 put

  Description: Unsupported operation since this cache is read only apart from
               during the build process
  Returntype : None
  Exceptions : Thrown if ever called
  Caller     : BaseAdaptors
  Status     : Beta

=cut

sub put {
  my ($self) = @_;
  throw 'Unsupported operation';
  return;
}

=head2 remove

  Description: Unsupported operation since this cache is read only apart from
               during the build process
  Returntype : None
  Exceptions : Thrown if ever called
  Caller     : BaseAdaptors
  Status     : Beta

=cut

sub remove {
  my ($self) = @_;
  throw 'Unsupported operation';
  return;
}

1;