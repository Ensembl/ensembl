=head1 LICENSE

Copyright [2018] EMBL-European Bioinformatics Institute

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

Bio::EnsEMBL::RNAProduct - A class representing the mature RNA product
of a transcript

=head1 DESCRIPTION

TODO

=head1 SYNOPSIS

  my $rnaproduct = Bio::EnsEMBL::RNAProduct->new(
    -SEQ_START => 36,
    -SEQ_END   => 58
  );

  # Stable-ID setter
  $rnaproduct->stable_id('ENSM00090210');

  # Get start and end position in the precursor transcript
  my $start = $rnaproduct->start();
  my $end = $rnaproduct->end();

=cut


package Bio::EnsEMBL::RNAProduct;

use vars qw($AUTOLOAD);
use strict;

use Bio::EnsEMBL::Utils::Exception qw(throw warning );
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Bio::EnsEMBL::Utils::Scalar qw( assert_ref wrap_array );
use Scalar::Util qw(weaken);

use Bio::EnsEMBL::Storable;

use parent qw(Bio::EnsEMBL::Storable);


=head2 new

  Arg [-SEQ_START]    : The offset in the Translation indicating the start
                        position of the product sequence.
  Arg [-SEQ_END]      : The offset in the Translation indicating the end
                        position of the product sequence.
  Arg [-STABLE_ID]    : The stable identifier for this RNAPRoduct
  Arg [-VERSION]      : The version of the stable identifier
  Arg [-DBID]         : The internal identifier of this RNAProduct
  Arg [-ADAPTOR]      : The TranslationAdaptor for this RNAProduct
  Arg [-SEQ]          : Manually sets the nucleotide sequence of this
                        rnaproduct. May be useful if this rnaproduct is not
                        stored in a database.
  Arg [-CREATED_DATE] : the date the rnaproduct was created
  Arg [-MODIFIED_DATE]: the date the rnaproduct was modified
  Example    : my $rp = Bio::EnsEMBL::RNAProduct->new(
                 -SEQ_START => 36,
                 -SEQ_END   => 58
               );
  Description: Constructor.  Creates a new RNAProduct object
  Returntype : Bio::EnsEMBL::RNAProduct
  Exceptions : none
  Caller     : general
  Status     : In Development

=cut

sub new {
  my $caller = shift;

  my $class = ref($caller) || $caller;

  my ($seq_start, $seq_end, $stable_id, $version, $dbID, $adaptor, $seq,
      $created_date, $modified_date ) =
	rearrange(["SEQ_START", "SEQ_END", "STABLE_ID", "VERSION", "DBID",
		   "ADAPTOR", "SEQ", "CREATED_DATE", "MODIFIED_DATE"], @_);

  # Default version
  $version //= 1;

  my $self = bless {
    'start'      => $seq_start,
    'end'        => $seq_end,
    'stable_id'  => $stable_id,
    'version'    => $version,
    'dbID'       => $dbID,
    'seq'        => $seq,
    'created_date' => $created_date,
    'modified_date' => $modified_date
  }, $class;

  $self->adaptor($adaptor);

  return $self;
}


=head2 created_date

  Arg [1]    : (optional) string $created_date - created date to set
  Example    : $rnaproduct->created_date('2007-01-10 20:52:00');
  Description: Getter/setter for attribute created_date
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub created_date {
  my $self = shift;
  $self->{'created_date'} = shift if (@_);
  return $self->{'created_date'};
}


=head2 end

  Arg [1]    : (optional) int $end - end position to set
  Example    : $rnaproduct->end(39);
  Description: Getter/setter for the value of end, which is a position within
               the precursor Transcript.
  Returntype : int
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub end {
  my $self = shift;
  $self->{'end'} = shift if (@_);
  return $self->{'end'};
}


=head2 modified_date

  Arg [1]    : (optional) string $modified_date - modification date to set
  Example    : $rnaproduct->modified_date('2007-01-10 20:52:00');
  Description: Getter/setter for attribute modified_date
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub modified_date {
  my $self = shift;
  $self->{'modified_date'} = shift if (@_);
  return $self->{'modified_date'};
}


=head2 stable_id

  Arg [1]    : (optional) string $stable_id - stable ID to set
  Example    : $rnaproduct->stable_id('ENSM00090210');
  Description: Getter/setter for attribute stable_id
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub stable_id {
  my $self = shift;
  $self->{'stable_id'} = shift if (@_);
  return $self->{'stable_id'};
}


=head2 start

  Arg [1]    : (optional) int $start - start position to set
  Example    : $rnaproduct->start(17);
  Description: Getter/setter for the value of start, which is a position within
               the precursor Transcript.
  Returntype : int
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub start {
  my $self = shift;
  $self->{'start'} = shift if (@_);
  return $self->{'start'};
}


=head2 version

  Arg [1]    : (optional) string $version - version to set
  Example    : $rnaproduct->version(2);
  Description: Getter/setter for attribute version
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub version {
  my $self = shift;
  $self->{'version'} = shift if (@_);
  return $self->{'version'};
}

1;
