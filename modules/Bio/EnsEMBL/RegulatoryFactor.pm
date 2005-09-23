use strict;

package Bio::EnsEMBL::RegulatoryFactor;

use Bio::EnsEMBL::Storable;
use Bio::EnsEMBL::Utils::Argument qw(rearrange);

use vars qw(@ISA);
@ISA = qw(Bio::EnsEMBL::Storable);


=head2 new

  Arg [NAME] : string (optional)
               The name of this regulatory factor
  Arg [TYPE]: string (optional)
               The type of repeat consensus
  Arg [...]: Named arguments to superclass constructor
             (see Bio::EnsEMBL::Storable)
  Example    : $rc = Bio::EnsEMBL::RegulatoryFactor->new
                       (-NAME    => 'SomeFactorName'
                        -TYPE    => 'promoter',
                        -DBID    => 1023,
                        -ADAPTOR => $rc_adaptor);
  Description: Creates a new Bio::EnsEMBL::RegulatoryFactor object
  Returntype : Bio::EnsEMBL::RegulatoryFactor
  Exceptions : none
  Caller     : RegulatoryFactor
  Status     : At Risk
             : under development

=cut

sub new {
  my $caller = shift;

  my $class = ref($caller) || $caller;

  my $self = $class->SUPER::new(@_);

  my ($name, $type ) =
    rearrange(['NAME', 'TYPE'], @_);

  $self->{'name'} = $name;
  $self->{'type'} = $type;

  return $self;
}


sub new_fast {
  my $class = shift;
  my $hashref = shift;
  return bless $hashref, $class;
}


=head2 name

  Arg [1]    : none
  Example    : $trans = $regulatory_factor->get_coding_transcript()
  Description: Returns the transcript that codes for this regulatory factor,
               if it is defined.
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : At Risk
             : under development

=cut

sub name {
  my $self = shift;
  $self->{'name'} = shift if(@_);
  return $self->{'name'};
}


=head2 coding_transcript

  Arg [1]    : none
  Example    : $transcript = $regulatory_factor->coding_transcript()
  Description: Get the transcript which codes for this regulatory factor.
  Returntype : Bio::EnsEMBL::Transcript, or undef
  Exceptions : none
  Caller     : general
  Status     : At Risk
             : under development

=cut

sub coding_transcript {

  my $self = shift;

  my ($transcript_id);

  my $sth = $self->adaptor()->db()->dbc()->prepare("SELECT transcript_id
			                            FROM regulatory_factor_coding
			                            WHERE regulatory_factor_id=?");

  $sth->execute($self->dbID());
  $sth->bind_columns(\$transcript_id);

  my $transcript_adaptor = $self->adaptor()->db()->get_TranscriptAdaptor();
  my $transcript;
  if ($sth->fetch) {
    $transcript = $transcript_adaptor->fetch_by_dbID($transcript_id);
  }

  return $transcript;
}


=head2 coding_gene

  Arg [1]    : none
  Example    : $gene = $regulatory_factor->coding_gene()
  Description: Get the gene which codes for this regulatory factor.
  Returntype : Bio::EnsEMBL::Gene, or undef
  Exceptions : none
  Caller     : general
  Status     : At Risk
             : under development

=cut

sub coding_gene {

  my $self = shift;

  my ($gene_id);

  my $sth = $self->adaptor()->db()->dbc()->prepare("SELECT gene_id
			                            FROM regulatory_factor_coding
			                            WHERE regulatory_factor_id=?");

  $sth->execute($self->dbID());
  $sth->bind_columns(\$gene_id);

  my $gene_adaptor = $self->adaptor()->db()->get_GeneAdaptor();
  my $gene;
  if ($sth->fetch) {
    $gene = $gene_adaptor->fetch_by_dbID($gene_id);
  }

  return $gene;
}



=head2 type

  Arg [1]    : string $type
  Example    : $type = $regulatory_factor->type()
  Description: Getter/Setter for the type of this regulatory factor
  Returntype : int
  Exceptions : none
  Caller     : general
  Status     : At Risk
             : under development

=cut

sub type {
  my $self = shift;
  $self->{'type'} = shift if(@_);
  return $self->{'type'};
}

1;

__END__

=head1 NAME - Bio::EnsEMBL::RegulatoryFactor

=head1 DESCRIPTION

This object represents an entry in the
regulatory_factor table.

=head1 AUTHOR

Glenn Proctor< email> glenn@ebi.ac.uk

