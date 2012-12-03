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

package Bio::EnsEMBL::RepeatConsensus;

use strict;

use Bio::EnsEMBL::Storable;
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Scalar::Util qw(weaken isweak);

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


=head2 new_fast

  Arg [1] : hashref to bless as a new RepeatConsensus 

  Description: Creates a new Bio::EnsEMBL::RepeatConsensus object
  Returntype : Bio::EnsEMBL::RepeatConsensus
  Exceptions : none
  Caller     : internal
  Status     : Stable

=cut

sub new_fast {
  my $class = shift;
  my $hashref = shift;
  my $self = bless $hashref, $class;
  weaken($self->{adaptor})  if ( ! isweak($self->{adaptor}) );
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

