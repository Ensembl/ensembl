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

Bio::EnsEMBL::UTR - A UTR feature with a location and a type (five prime/3 prime)

=head1 SYNOPSIS

  use Bio::EnsEMBL::UTR;

  $feature = Bio::EnsEMBL::UTR->new(
    -start         => 100,
    -end           => 220,
    -strand        => -1,
    -slice         => $slice,
    -type          => 'five_prime_UTR',
    -transcript    => $transcript
  );

=head1 DESCRIPTION

This is a UTR feature within the Ensembl system.
A UTR repsents the non-coding (untranslated) regions 
of a transcript. It can be 5' or 3'

=head1 METHODS

=cut

use strict;

package Bio::EnsEMBL::UTR;

use vars qw(@ISA);

use Bio::EnsEMBL::Feature;
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Scalar::Util qw(weaken isweak);

@ISA = qw(Bio::EnsEMBL::Feature);


=head2 new

  Arg [...]  : Named arguments passed to superclass
  Example    : $feature = Bio::EnsEMBL::UTR->new
                        (-start   => 1,
                         -end     => 100,
                         -strand  => 1,
                         -slice   => $slice,
                         -dbID    => 10,
                         -transcript  => $transcript,
                         -type    => 'five_prime_UTR');
  Description: Constructs a new Bio::EnsEMBL::UTR.
  Returntype : Bio::EnsEMBL::UTR
  Exceptions : Thrown on invalid -SLICE, -STRAND arguments
  Caller     : general, subclass constructors
  Status     : Stable

=cut

sub new {
  my $caller = shift;

  #allow this to be called as class or object method
  my $class = ref($caller) || $caller;
  my $self = $class->SUPER::new(@_);

  my ($transcript, $type, $seq_region_start, $seq_region_end) = rearrange(['TRANSCRIPT','TYPE', 'SEQ_REGION_START', 'SEQ_REGION_END'],@_);

  $self->{'transcript'} = $transcript;
  $self->{'type'} = $type;
  $self->{'seq_region_start'} = $seq_region_start;
  $self->{'seq_region_end'} = $seq_region_end;

  return $self;
}


=head2 transcript

  Arg [1]    : (optional) Bio::EnsEMBL::Transcript
  Example    : $transcript = $utr->transcript();
  Description: Getter/Setter for the transcript associated with this
               UTR.
  Returntype : Bio::EnsEMBL::Transcript
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub transcript {
  my $self = shift;
  $self->{'transcript'} = shift if(@_);
  return $self->{'transcript'};
}

=head2 translation

    Description: Fetch the translation associated with
                 this transcript, if it exists. Return undef
                 if there is no translation, ie. a pseudogene
    Returntype : Bio::EnsEMBL::Translation or undef
    Caller     : general
    Status     : Stable

=cut

sub translation {
    my $self = shift;
    return $self->transcript()->translation();
}

=head2 seq_region_start

  Arg [1]    : (optional) string $seq_region_start
  Example    : $seq_region_start = $cds->seq_region_start();
  Description: Getter/Setter for the seq_region_start for this UTR.
               Overwrite default method from Feature as UTR does not have
               a table
  Returntype : String
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub seq_region_start {
  my $self = shift;
  $self->{'seq_region_start'} = shift if(@_);
  return $self->{'seq_region_start'};
}

=head2 seq_region_end

  Arg [1]    : (optional) string $seq_region_end
  Example    : $seq_region_end = $cds->seq_region_end();
  Description: Getter/Setter for the seq_region_end for this UTR.
               Overwrite default method from Feature as UTR does not have
               a table
  Returntype : String
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub seq_region_end {
  my $self = shift;
  $self->{'seq_region_end'} = shift if(@_);
  return $self->{'seq_region_end'};
}

=head2 get_Gene

  Description: Get the gene associated with the ExonTranscript,
               if a transcript has been set
  Returntype : Bio::EnsEMBL::Gene or undef
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_Gene {
    my $self = shift;

    if($self->{'transcript'}) {
	return $self->{'transcript'}->get_Gene();
    }

    return;
}

=head2 type

  Arg [1]    : (optional) string $type
  Example    : print $utr->type();
  Description: This method returns a string that describes
               the type of UTR feature (5 prime or 3 prime)
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub type {
  my $self = shift;
  $self->{'type'} = shift if( @_ );
  return ( $self->{'type'} || 'UTR' );
}


=head2 summary_as_hash

  Example    : my $hash = $utr->summary_as_hash();
  Description: Generates a HashRef compatible with GFFSerializer. Adds
               Parent, source and type to the basic Feature hash
  Returntype : Hash
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub summary_as_hash {
  my ($self) = @_;
  my $hash = $self->SUPER::summary_as_hash();
  $hash->{'type'} = $self->type() if $self->type();
  my $parent = $self->transcript();
  $hash->{'Parent'} = $parent->display_id if defined $parent;
  $hash->{'source'} = $self->transcript->source() if $self->transcript();
  delete $hash->{id};
  return $hash;
}


1;
