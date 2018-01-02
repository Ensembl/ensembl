=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

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

Bio::EnsEMBL::DBSQL::Support::FullIdCache - ID based caching using all available values  

=head1 SYNOPSIS

  my $cache = Bio::EnsEMBL::DBSQL::Support::FullIdCache->new($adaptor);
  my $obj = $cache->get(21);

=head1 DESCRIPTION

An implementation of caching which uses a raw hash to hold all available
values from an adaptor. Useful for working with a controlled vocabulary
table where cardinality is low.

Provides extra functionality to compute additional lookup keys.

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
  my $objects = $adaptor->generic_fetch();
  foreach my $object (@{$objects}) {
    my $key = $object->dbID();
    $cache{$key} = $object;
    #Add to additional lookup
    $self->add_to_additional_lookups($key, $object);
  }
  return \%cache;
}

########### Additional lookup code

sub put {
  my ($self, $key, $object) = @_;
  my $old = $self->SUPER::put($key, $object);
  #Add to additional lookup
  $self->remove_from_additional_lookup($key, $old) if $old;
  $self->add_to_additional_lookups($key, $object);
  return $old if $old;
  return;
}

sub remove {
  my ($self, $key) = @_;
  my $old = $self->SUPER::remove($key);
  if($old) {
    #Remove it from the additional lookup
    $self->remove_from_additional_lookup($key, $old);
    return $old;
  }
  return;
}

sub clear_cache {
  my ($self) = @_;
  $self->delete_cache();
  #Remove the additional lookup hash contents
  delete $self->{_additional_lookup};
  return;
}

=head2 get_by_additional_lookup

  Arg [1]     : String key of the lookup to search for the value in
  Arg [2]     : String value to search for. We expect exact lookups in the hash
  Description : Returns the object linked to the value in the specified lookup.
  Example     : my $analysis = $cache->get_by_additional_lookup('logic_name', 'xrefchecksum');
  Returntype  : Object a single object
  Exceptions  : Throws an exception if there are more than one ID linked to the
                value lookup. Also thrown if additional lookups are not supported 
  Caller      : BaseAdaptors
  Status      : Beta

=cut

sub get_by_additional_lookup {
  my ($self, $key, $value) = @_;
  $self->cache(); # trigger cache building
  my $additional_lookup = $self->_additional_lookup();
  if(exists $additional_lookup->{$key}) {
    if(exists $additional_lookup->{$key}->{$value}) {
      my $ids = $additional_lookup->{$key}->{$value};
      my $size = scalar(@{$ids});
      if($size > 1) {
        throw "The lookup $key and search value $value has more than one value attached. Use get_all_by_additional_lookup() instead to fetch";
      }
      elsif($size == 1) {
        return $self->get($ids->[0]);
      }
    }
  }
  return;
}

=head2 get_all_by_additional_lookup

  Arg [1]     : String key of the lookup to search for the value in
  Arg [2]     : String value to search for. We expect exact lookups in the hash
  Description : Returns an array of all the objects linked to the value
                in the specified lookup.
  Example     : my $array = $cache->get_all_by_additional_lookup('logic_name', 'xrefchecksum');
  Returntype  : ArrayRef of objects keyed agains the second argument
  Exceptions  : Throws an exception if there are more than one ID linked to the
                value lookup. Also thrown if additional lookups are not supported
  Caller      : BaseAdaptors
  Status      : Beta

=cut

sub get_all_by_additional_lookup {
  my ($self, $key, $value) = @_;
  $self->cache(); # trigger cache building
  my $additional_lookup = $self->_additional_lookup();
  if(exists $additional_lookup->{$key}) {
    if(exists $additional_lookup->{$key}->{$value}) {
      my $ids = $additional_lookup->{$key}->{$value};
      return $self->get_by_list($ids);
    }
  }
  return [];
}

=head2 remove_from_additional_lookup

  Arg [1]     : String The lookup key to remove from the additional lookup hash
  Arg [2]     : Object The object to remove from the additional lookup hash
  Description : Re-computes the additional keys for this object 
  Example     : $cache->remove_Object_from_additional_lookup($lookup_key, $object);
  Returntype  : None
  Exceptions  : Thrown if we do not support additional lookups
  Caller      : BaseAdaptors
  Status      : Beta

=cut

sub remove_from_additional_lookup {
  my ($self, $lookup_key, $object) = @_;

  # Compute the keys
  my $keys = $self->compute_keys($object);
  return if scalar(keys %{$keys}) == 0;

  my $additional_lookup = $self->_additional_lookup();

  foreach my $key (keys %{$keys}) {
    my $value = $keys->{$key};

    #Only remove if we had originally stored this as an
    #additional lookup
    if(exists $additional_lookup->{$key}) {
      if(exists $additional_lookup->{$key}->{$value}) {

        #Get the object ID & lookup the array of DBIDs
        my $lookup_keys = $additional_lookup->{$key}->{$value};
        my $length = scalar(@{$lookup_keys});
        for(my $i = 0; $i < $length; $i++) {
          if($lookup_keys->[$i] == $lookup_key) {
            #remove the 1 lookup key from the array and then terminate the
            #loop as we found our value
            splice(@{$lookup_keys}, $i, 1); 
            last;
          }
        }

        #If the size has hit 0 then delete the array
        if(scalar(@{$lookup_keys}) == 0) {
          delete $additional_lookup->{$key}->{$value};
        }
      }
    }
  }

  return;
}

=head2 compute_keys

  Arg [1]     : Object The object to compute keys from
  Description : Override to provide support for additional key lookup. The
                keys of the hash should represent the lookup name and the
                value is the computed key.
  Example     : Example of returning hash not of its usage. Proposed Analysis encoding
                { logic_name => 'xrefalignment', display_label => 'Xref Alignment'}
  Returntype  : HashRef key is the lookup name and value is the computed key
  Exceptions  : none
  Caller      : BaseAdaptors
  Status      : Beta

=cut

sub compute_keys {
  my ($self, $object) = @_;
  return {};
}

=head2 add_to_additional_lookups

  Arg [1]     : String The key used in the primary lookup hash. Normally
                a DB identifier
  Arg [2]     : Object The object to add to the additional lookups
  Description : Internally calls the C<compute_keys()> method and adds
                the object to the C<_additional_lookup()> hash.
  Returntype  : None
  Exceptions  : Thrown if additional lookups are not supported
  Caller      : BaseAdaptors
  Status      : Beta

=cut

sub add_to_additional_lookups {
  my ($self, $lookup_key, $object) = @_;
  my $keys = $self->compute_keys($object);
  return if scalar(keys %{$keys}) == 0;
  my $additional_lookup = $self->_additional_lookup();
  foreach my $key (keys %{$keys}) {
    my $value = $keys->{$key};
    push(@{$additional_lookup->{$key}->{$value}}, $lookup_key);
  }
  return;
}

=head2 _additional_lookup

  Description : Returns the additional lookup hash
  Example     : Example of additional hash structure (key is 
                lookup name, second key is value to search for
                and value is an array of dbIDs)
                {
                  logic_name => {
                    xrefalignment => [1]
                  },
                  display_label => {
                    'Xref Alignment' => [1]
                  }
                }
  Returntype  : HashRef
  Exceptions  : none
  Caller      : BaseAdaptors
  Status      : Beta

=cut

sub _additional_lookup {
  my ($self) = @_;
  $self->{_additional_lookup} ||= {};
  return $self->{_additional_lookup};
}

1;
