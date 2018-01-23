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

Bio::EnsEMBL::IdMapping::TinyFeature - lightweight feature object

=head1 SYNOPSIS

This object isn't instantiated. See objects which inherit from it
(TinyGene, TinyTranscript, etc.) for examples.

=head1 DESCRIPTION

This is the base class for the lightweight feature objects used by the
stable Id maping application. For performance reasons, these objects
are instantiated using a new_fast() method. The internal implementation
is an arrayref (rather than the more common hashref), which optimises
memory usage.

There are no adaptors to fetch TinyFeatures from the database. You
rather use the normal feature adaptors and then create the TinyFeatures
from the heavy objects you get. The memory saving will therefore mainly
take effect when serialising and reloading these objects.

Also note that TinyFeatures don't have a slice attached to them - all
location information (where required) is stored on the feature object
directly.

=head1 METHODS

  new_fast
  id
  stable_id
  version
  created_date
  modified_date
  to_string

=cut

package Bio::EnsEMBL::IdMapping::TinyFeature;

# internal data structure (array indices):
#
#  0  dbID
#  1  stable_id
#  2  version
#  3  created_date
#  4  modified_date
#
# other instance variables differ by subclass implementation, so look there.


use strict;
use warnings;
no warnings 'uninitialized';

use Bio::EnsEMBL::Utils::Exception qw(throw warning);


=head2 new_fast

  Arg[1]      : Arrayref $array_ref - the arrayref to bless into the new object 
  Description : Constructor.
  Return type : Bio::EnsEMBL::IdMapping::TinyFeature implementing class
  Exceptions  : none
  Caller      : Bio::EnsEMBL::IdMapping::Cache
  Status      : At Risk
              : under development

=cut

sub new_fast {
  my $class = shift;
  my $array_ref = shift;
  return bless $array_ref, $class;
}


=head2 id

  Arg[1]      : (optional) Int - the feature's internal Id ("dbID")
  Description : Getter/setter for the feature's internal Id.
  Return type : Int
  Exceptions  : none
  Caller      : Bio::EnsEMBL::IdMapping::Cache
  Status      : At Risk
              : under development

=cut

sub id {
  my $self = shift;
  $self->[0] = shift if (@_);
  return $self->[0];
}


=head2 stable_id

  Arg[1]      : (optional) String - the feature's stable Id
  Description : Getter/setter for the feature's stable Id.
  Return type : String
  Exceptions  : none
  Caller      : Bio::EnsEMBL::IdMapping::Cache
  Status      : At Risk
              : under development

=cut

sub stable_id {
  my $self = shift;
  $self->[1] = shift if (@_);
  return $self->[1];
}


=head2 version

  Arg[1]      : (optional) Int - the feature's stable Id version
  Description : Getter/setter for the feature's stable Id version.
  Return type : Int
  Exceptions  : none
  Caller      : Bio::EnsEMBL::IdMapping::Cache
  Status      : At Risk
              : under development

=cut

sub version {
  my $self = shift;
  $self->[2] = shift if (@_);
  return $self->[2];
}


=head2 created_date

  Arg[1]      : (optional) String - the feature's stable Id creation date
  Description : Getter/setter for the feature's stable Id creation date.
  Return type : String
  Exceptions  : none
  Caller      : Bio::EnsEMBL::IdMapping::Cache
  Status      : At Risk
              : under development

=cut

sub created_date {
  my $self = shift;
  $self->[3] = shift if (@_);
  return $self->[3];
}


=head2 modified_date

  Arg[1]      : (optional) String - the feature's stable Id modification date
  Description : Getter/setter for the feature's stable Id modification date.
  Return type : String
  Exceptions  : none
  Caller      : Bio::EnsEMBL::IdMapping::Cache
  Status      : At Risk
              : under development

=cut

sub modified_date {
  my $self = shift;
  $self->[4] = shift if (@_);
  return $self->[4];
}


=head2 to_string

  Example     : print LOG "Created ", $f->to_string, "\n";
  Description : Prints a string representation of the feature for debug
                purposes.
  Return type : String
  Exceptions  : none
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub to_string {
  my $self = shift;
  return $self->id.':'.$self->stable_id.'.'.$self->version;
}


1;

