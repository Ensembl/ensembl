=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2020] EMBL-European Bioinformatics Institute

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

Bio::EnsEMBL::Biotype

=head1 SYNOPSIS

    my $biotype = new Bio::EnsEMBL::Biotype(
      -name          => 'new_biotype,
      -object_type   => 'gene',
      -biotype_group => 'a_biotype_group',
      -so_acc        => 'SO::1234567',
      -description   => 'New biotype'
    );

    my $name = $biotype->name();
    my $biotype_group = $biotype->biotype_group();
    my $so_acc = $biotype->so_acc();

=head1 DESCRIPTION

  This is the Biotype object class.
  Gene and Transcript objects used to have a biotype() method that returned the biotype name
  (the biotype field in the gene and transcript tables).
  From e93 a new biotype table was added. However because of legacy code using direct sql
  queries on the biotype column of gene and transcript tables, that column that contains the
  biotype name was not replaced by biotype_id containing a foreign key to the new biotype table.
  Gene and Transcripts can still link to a Biotype through the key (name, object_type).

=head1 METHODS

=cut


package Bio::EnsEMBL::Biotype;

use strict;
use warnings;

use Bio::EnsEMBL::Storable;
use Bio::EnsEMBL::Utils::Exception qw(throw deprecate warning);
use Bio::EnsEMBL::Utils::Scalar qw(check_ref assert_ref);
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Scalar::Util qw(weaken isweak);

use parent qw(Bio::EnsEMBL::Storable);

=head2 new

  Arg [-BIOTYPE_ID]  :
      int - dbID of the biotype
  Arg [-NAME]    :
      string - the name of the biotype (for ensembl)
  Arg [-OBJECT_TYPE] :
      string - the object type this biotype applies to (gene or transcript)
  Arg [-BIOTYPE_GROUP]  :
      string - the name of the biotype group (for ensembl)
  Arg [-SO_ACC] :
      string - the Sequence Ontology accession of this biotype
  Arg [-SO_TERM] :
      string - the Sequence Ontology term for the SO accession of this biotype
  Arg [-DESCRIPTION] :
      string - the biotype description
  Arg [-DB_TYPE] :
      string - the database type for this biotype
  Arg [-ATTRIB_TYPE_ID] :
      int - attrib_type_id

  Example    : $biotype = Bio::EnsEMBL::Biotype->new(...);
  Description: Creates a new biotype object
  Returntype : Bio::EnsEMBL::Biotype
  Exceptions : none

=cut

sub new {
  my ( $caller, @args ) = @_;

  my $class = ref($caller) || $caller;

  my $self = $class->SUPER::new();

  my($dbID, $name, $object_type, $biotype_group, $so_acc, $so_term, $description, $db_type, $attrib_type_id) =
    rearrange([qw(BIOTYPE_ID NAME OBJECT_TYPE BIOTYPE_GROUP SO_ACC SO_TERM DESCRIPTION DB_TYPE ATTRIB_TYPE_ID)], @args);

  $self->{'dbID'} = $dbID;
  $self->{'name'} = $name;
  $self->{'object_type'} = $object_type;
  $self->{'biotype_group'} = $biotype_group;
  $self->{'so_acc'} = $so_acc;
  $self->{'so_term'} = $so_term;
  $self->{'description'} = $description;
  $self->{'db_type'} = $db_type;
  $self->{'attrib_type_id'} = $attrib_type_id;

  return $self;
}

=head2 new_fast

  Arg [1]    : hashref to be blessed
  Description: Construct a new Bio::EnsEMBL::Biotype using the hashref.
  Exceptions : none
  Returntype : Bio::EnsEMBL::Biotype

=cut


sub new_fast {
  my ( $class, $hashref ) = @_;

  my $self = bless $hashref, $class;

  if ( !isweak($self->{adaptor}) ) {
    weaken($self->{adaptor})
  }

  return $self;
}

=head2 name

  Arg [1]    : (optional) string $name
               The name of this biotype according to ensembl.
  Example    : $name = $biotype->name()
  Description: Getter/Setter for the name of this biotype.
  Returntype : string
  Exceptions : none

=cut

sub name {
  my ( $self, $name ) = @_;

  if ( defined($name) ) {
    $self->{'name'} = $name;
  }

  return $self->{'name'};
}

=head2 biotype_group

  Arg [1]    : (optional) string $biotype_group
  Example    : $biotype_group = $biotype->biotype_group();
  Description: Getter/Setter for the biotype_group of this biotype.
               Biotype groups are used internally at ensembl pipelines
               and consist on few defined categories.
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub biotype_group {
  my ( $self, $biotype_group ) = @_;

  if ( defined($biotype_group) ) {
    $self->{'biotype_group'} = $biotype_group;
  }

  return $self->{'biotype_group'};
}

=head2 so_acc

  Arg [1]    : (optional) string $so_acc
  Example    : $feat->so_acc();
  Description: Getter/Setter for the Sequence Ontology accession of this biotype.
               It must be a SO like accession.
  Returntype : string
  Exceptions : thrown if an invalid so_acc argument is passed

=cut

sub so_acc {
  my ( $self, $so_acc ) = @_;

  if ( defined($so_acc) ) {
    # throw an error if setting something that does not look like an SO acc
    unless ( $so_acc =~ m/\ASO:\d+\z/x ) {
      throw("so_acc must be a Sequence Ontology accession. '$so_acc' does not look like one.")
    }

    $self->{'so_acc'} = $so_acc;
  }

  return $self->{'so_acc'};
}

=head2 so_term

  Arg [1]    : (optional) string $so_term
  Example    : $feat->so_term();
  Description: Getter/Setter for the Sequence Ontology term of this biotype.
  Returntype : string
  Exceptions : none

=cut

sub so_term {
  my ( $self, $so_term ) = @_;

  if ( defined($so_term) ) {
    $self->{'so_term'} = $so_term;
  }

  return $self->{'so_term'};
}

=head2 object_type

  Arg [1]    : (optional) string $object_type
  Example    : $object_type = $biotype->object_type();
  Description: Getter/Setter for the object_type of this biotype.
               Biotypes can be assigned to either genes or transcripts,
               object_type refers to which of them.
  Returntype : string
  Exceptions : thrown if an invalid object_type argument is passed (not gene or transcript)

=cut

sub object_type {
  my ( $self, $object_type ) = @_;

  if ( defined($object_type) ) {
    $object_type = lc $object_type;
    # throw an error if setting something that does not look like an SO acc
    unless ( $object_type eq 'gene' || $object_type eq 'transcript' ) {
      throw("object_type must be gene or transcript. Got '$object_type'.")
    }

    $self->{'object_type'} = $object_type;
  }

  return $self->{'object_type'};
}

1;
