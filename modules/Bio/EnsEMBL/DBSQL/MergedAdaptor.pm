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

Bio::EnsEMBL::DBSQL::MergedAdaptor

=head1 SYNOPSIS

  # load all available adaptors of the given type for the species
  my $merged_adaptor = Bio::EnsEMBL::DBSQL::MergedAdaptor->new(
    -species => "human",
    -type    => "gene"
  );

  # only load adaptors from the given groups of the given type for the species
  my $merged_adaptor = Bio::EnsEMBL::DBSQL::MergedAdaptor->new(
    -species => "human",
    -type    => "gene",
    -groups  => ['core','otherfeatures']
  );

=head1 DESCRIPTION

The MergedAdaptor object is merely a list of adaptors. AUTOLOAD is used
to call a subroutine on each adaptor and merge the results. This object structure
allows you to treat a set of adaptors as a logical single entity. The end result
is that disparate database source data sets are accessible through a single
adaptor call.

This code will convert single object return calls into ArrayRef returning calls
and so is only safe to use with the C<fetch_all_XXX> or C<get_all_XXX> methods.

=head1 METHODS

=cut


package Bio::EnsEMBL::DBSQL::MergedAdaptor;

use strict;
use warnings;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Utils::Scalar qw(wrap_array assert_ref);
require Bio::EnsEMBL::Registry;
my $registry = "Bio::EnsEMBL::Registry";


=head2 new

  Arg [SPECIES]: String species name to get adaptors for
  Arg [TYPE]   : String type to get adaptors for
  Arg [GROUPS] : (optional) ArrayRef of groups to load
  Example      : my $adaptor = Bio::EnsEMBL::DBSQL::MergedAdaptor->new(
                    -species=> 'human', -type =>'Population', -groups => ['Sanger','Ensembl']);
                 my $alL_adaptor = Bio::EnsEMBL::DBSQL::MergedAdaptor->new(
                    -species=> 'human', -type =>'Population');
  Description: Creates a new MergedAdaptor.
  Returntype : Bio::EnsEMBL::DBSQL::MergedAdaptor
  Exceptions : throws if species or type not specified
  Caller     : general
  Status     : At Risk
             : Under development

=cut

sub new {
  my ($class, @args) = @_;
  $class = ref($class) || $class;
  my $self = bless({}, $class);
  my ($species, $type, $groups) = rearrange([qw(SPECIES TYPE GROUPS)], @args);

  if(!defined($species)|| !defined($type)){
    throw "-SPECIES and -TYPE must be specified";
  }
  $self->_populate_adaptors($species, $type, $groups);
  return $self;
}

=head2 _populate_adaptors

  Arg [1]      : String species name to get adaptors for
  Arg [2]      : String type to get adaptors for
  Arg [3]      : (optional) ArrayRef of groups to load
  Description  : Auto-populates the current MergedAdaptor with the
                 adaptors linked to this species, type and optional set of groups
  Caller       : general
  Status       : At Risk

=cut

sub _populate_adaptors {
  my ($self, $species, $type, $groups) = @_;
  if($groups) {
    assert_ref($groups, 'ARRAY', '-GROUPS');
    my @adaptors = map { $registry->get_adaptor($species, $_, $type) } @{$groups};
    $self->add_list(@adaptors);
  }
  else {
    my $adaptors = $registry->get_all_adaptors(-SPECIES => $species, -TYPE => $type);
    $self->add_list(@{$adaptors});
  }
  return;
}

=head2 add_list

  Arg [n]     : Adaptors to add into this instance
  Description : Adds the given adaptors to the internal adaptor list

=cut

sub add_list {
  my ($self, @adaptors) = @_;
  push(@{$self->{adaptors}}, @adaptors);
  return;
}

=head2 add_list

  Arg [1]     : Adaptor to add into this instance
  Description : Adds the given adaptor to the internal adaptor list. For
                multiple adaptor addition use C<add_list()>.

=cut

sub add_adaptor {
  my ($self, $adaptor) = @_;
  $self->add_list($adaptor);
  return;
}

=head2 can

  Arg [1]     : String method name to be called
  Description : Implementation of UNIVERSAL::can(). We loop through the
                available adaptors and return true if any will respond
                to the given method name
  Returntype  : Boolean indicating if any delegating object will respond to this method

=cut

sub can {
  my ($self, $method) = @_;
  foreach my $adaptor (@{$self->{adaptors}}) {
    return 1 if $adaptor->can($method);
  }
  return 0;
}


=head2 isa

  Arg [1]     : String method name to be called
  Description : Implementation of UNIVERSAL::isa(). We loop through the
                available adaptors and return true if any inherited from
                the given class
  Returntype  : Boolean indicating if any delegating object inherits from the given class

=cut

sub isa {
  my ($self, $isa) = @_;
  foreach my $adaptor (@{$self->{adaptors}}) {
    return 1 if $adaptor->isa($isa);
  }
  return 0;
}

use vars '$AUTOLOAD';

=head2 AUTOLOAD

  Description : Internal override of AUTLOAD. The code will detect the requested
                method, loop through all available adaptors and will 
  Returntype  : Boolean indicating if any delegating object inherits from the given class

=cut

sub AUTOLOAD {
  my ($self, @args) = @_;
  my @return;

  #Detect required method
  $AUTOLOAD =~ /^.*::(\w+)+$/ ;
  my $sub = $1;

  foreach my $adaptor (@{$self->{adaptors}}) {
    my @local_return;
    if(my $method_ref = $adaptor->can($sub)) {
      my $ref = $method_ref->($adaptor, @args);
      push(@return, @{wrap_array($ref)});
    }
    else{
      warning("In Merged Adaptor $adaptor cannot call sub $sub");
    }
  }
  return \@return;
}

sub DESTROY {

}

1;
