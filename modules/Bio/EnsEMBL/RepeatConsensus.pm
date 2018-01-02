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

package Bio::EnsEMBL::RepeatConsensus;

use strict;

use Bio::EnsEMBL::Storable;
use Bio::EnsEMBL::Utils::Argument qw(rearrange);

use vars qw(@ISA);
@ISA = qw(Bio::EnsEMBL::Storable);


=head2 new

  Arg [NAME] : string (optional)
               The name of this repeat consensus
  Arg [LENGTH]: int (optional)
               The length of the repeat consensus sequence
  Arg [REPEAT_CLASS]: string (optional)
               The type of repeat consensus
  Arg [REPEAT_CONSENSUS]: string (optional)
               The sequence of this repeat consensus
  Arg [REPEAT_TYPE]: string 
               Its like class only more general
  Arg [...]: Named arguments to superclass constructor
             (see Bio::EnsEMBL::Storable)
  Example    : $rc = Bio::EnsEMBL::RepeatConsensus->new
                       (-REPEAT_CONSENSUS => 'AATG'
                        -NAME => '(AATG)n',
                        -REPEAT_CLASS => 'Simple_repeat',
                        -LENGTH => '4',
                        -DBID => 1023,
                        -ADAPTOR => $rc_adaptor);
  Description: Creates a new Bio::EnsEMBL::RepeatConsensus object
  Returntype : Bio::EnsEMBL::RepeatConsensus
  Exceptions : none
  Caller     : RepeatFeatureAdaptors
  Status     : Stable

=cut

sub new {
  my $caller = shift;

  my $class = ref($caller) || $caller;

  my $self = $class->SUPER::new(@_);

  my ($name, $length, $repeat_class, $repeat_consensus, $repeat_type ) =
    rearrange(['NAME', 'LENGTH', 'REPEAT_CLASS', 'REPEAT_CONSENSUS', 'REPEAT_TYPE'], @_);

  $self->{'name'} = $name;
  $self->{'length'} = $length;
  $self->{'repeat_class'} = $repeat_class;
  $self->{'repeat_consensus'} = $repeat_consensus;
  $self->{'repeat_type'} = $repeat_type;

  return $self;
}


=head2 name

  Arg [1]    : string $name (optional)
  Example    : $name = $repeat_consensus->name()
  Description: Getter/Setter for the name of this repeat_consensus
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub name {
  my $self = shift;
  $self->{'name'} = shift if(@_);
  return $self->{'name'};
}


=head2 length

  Arg [1]    : int $length (optional)
  Example    : $length = $repeat_consensus->length()
  Description: Getter/Setter for the length of this repeat_consensus
  Returntype : int
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub length {
  my $self = shift;
  $self->{'length'} = shift if(@_);
  return $self->{'length'};
}


=head2 repeat_class

  Arg [1]    : string $class (optional)
               The class of 
  Example    : $class = $repeat_consensus->repeat_class()
  Description: Getter/Setter for the class of this repeat_consensus
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub repeat_class {
  my $self = shift;
  $self->{'repeat_class'} = shift if(@_);
  return $self->{'repeat_class'};
}

=head2 repeat_type

  Arg [1]    : string $type (optional)
               The type of the consensus 
  Example    : $type = $repeat_consensus->repeat_type()
  Description: Getter/Setter for the type of this repeat_consensus
  Returntype : string
  Exceptions : none
  Caller     : general

=cut

sub repeat_type {
  my $self = shift;
  $self->{'repeat_type'} = shift if(@_);
  return $self->{'repeat_type'};
}


=head2 desc

  Arg [1]    : none
  Example    : $desc = $repeat_consensus->desc()
  Description: Getter for the description of this repeat consensus as extracted
               from the repeat_class.  This method is probably useless.
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Medium risk

=cut

sub desc {
  my $self = shift;
  my $class = $self->repeat_class or return;
  return "class=$class";
}



=head2 repeat_consensus

  Arg [1]    : string $consensus_seq (optional)
               The sequence of this repeat consensus
  Example    : $consensus = $repeat_consensus->repeat_consensus();
  Description: Getter/Setter for the sequence of this repeat_consensus.
  Returntype : string
  Exceptions : none
  Caller     : general

=cut

sub repeat_consensus {
  my $self = shift;
  $self->{'repeat_consensus'} = shift if(@_);
  return $self->{'repeat_consensus'};
}



=head2 seq

  Arg [1]    : none
  Example    : none
  Description: Returns the repeat consensus.  This method is useless - Use
               repeat_consensus() instead.
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub seq {
  my( $self ) = @_;
  return $self->repeat_consensus;
}

1;

__END__

=head1 NAME - Bio::EnsEMBL::RepeatConsensus

=head1 DESCRIPTION

This object represents an entry in the
repeat_consensus table.

It can contain the consensus sequence for a
repeat such as a particular Alu, or "cag" for a
simple triplet repeat.

