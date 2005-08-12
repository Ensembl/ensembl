#
# Ensembl module for Bio::EnsEMBL::RegulatoryFeature
#
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::RegulatoryFeature - A feature representing a regulatory region

=head1 SYNOPSIS

use Bio::EnsEMBL::RegulatoryFeature;

$feature = Bio::EnsEMBL::RegulatoryFeature->new(-start     => 100,
                                                -end       => 220,
                                                -strand    => -1,
                                                -slice     => $slice,
                                                -analysis  => $analysis,
                                                -factor    => $factor,
                                                -dbID      => 1230,
                                                -adaptor   => $adaptor);

=head1 DESCRIPTION

This is a feature which represents a regulatory region.

=head1 AUTHOR - Glenn Proctor

This module is part of the Ensembl project http://www.ensembl.org

=head1 CONTACT

Post comments/questions to the ensembl development list: ensembl-dev@ebi.ac.uk

=head1 METHODS

=cut

use strict;

package Bio::EnsEMBL::RegulatoryFeature;

use vars qw(@ISA);

use Bio::EnsEMBL::Feature;
use Bio::EnsEMBL::Utils::Argument qw(rearrange);

@ISA = qw(Bio::EnsEMBL::Feature);


=head2 new

  Arg [...]  : Named arguments passed to superclass
  Example    : $feature = Bio::EnsEMBL::RegulatoryFeature->new
                                               (-name      => 'hsa-miR-108',
						-start     => 100,
                                                -end       => 220,
                                                -strand    => -1,
                                                -slice     => $slice,
                                                -analysis  => $analysis,
                                                -factor     => $factor,
                                                -dbID      => 1230,
                                                -adaptor   => $adaptor);
  Description: Constructs a new Bio::EnsEMBL::RegulatoryFeature.
  Exceptions : Thrown on invalid -SLICE, -ANALYSIS, -STRAND arguments
  Caller     : general, subclass constructors
  Status     : At Risk
             : under development

=cut

sub new {
  my $caller = shift;

  #allow this to be called as class or object method
  my $class = ref($caller) || $caller;
  my $self = $class->SUPER::new(@_);

  my ($name, $factor) = rearrange(['NAME', 'FACTOR'],@_);

  $self->{'name'} = $name;
  $self->{'factor'} = $factor;

  return $self;
}


sub new_fast {
  my $class = shift;
  my $hashref = shift;
  return bless $hashref, $class;
}


=head2 factor

  Arg [1]    : (optional) Bio::EnsEMBL::RegulatoryFactor
  Example    : 
  Description: Getter/Setter for the regulatory factor that this regulatory feature represents
  Returntype : Bio::EnsEMBL::RegulatoryFactor
  Exceptions : none
  Caller     : general
  Status     : At Risk
             : under development

=cut

sub factor {
  my $self = shift;

  if(@_) {
    my $rm = shift;
    if(defined($rm)) {
      if(!ref($rm) || !$rm->isa('Bio::EnsEMBL::RegulatoryFactor')) {
        throw('RegulatoryFactor arg must be a Bio::EnsEMBL::RegulatoryFactor');
      }
    }
    $self->{'factor'} = $rm;

  }

  return $self->{'factor'};
}


=head2 name

  Arg [1]    : none
  Example    : print $rf->name();
  Description: Getter/setter for this features name type.
  Returntype : string
  Exceptions : none
  Caller     : web drawing code
  Status     : At Risk
             : under development

=cut

sub name {
  my $self = shift;
  $self->{'name'} = shift if(@_);
  return $self->{'name'};
}

=head2 regulated_transcripts

  Arg [1]    : none
  Example    : my $transcipts = $rf->regulated_transcripts();
  Description: Return a list of all transcripts that are regulated by 
               this RegulatoryFeature.
  Returntype : listREF of Bio::EnsEMBL::RegulatoryFeatures
  Exceptions : none
  Caller     : ?
  Status     : At Risk
             : under development

=cut

sub regulated_transcripts {

  my $self = shift;

  my ($transcript_id);

  my $sth = $self->adaptor()->db()->dbc()->prepare("SELECT ensembl_object_id
			                            FROM regulatory_feature_object
			                            WHERE ensembl_object_type='Transcript'
			                            AND regulatory_feature_id=?");

  $sth->execute($self->dbID());
  $sth->bind_columns(\$transcript_id);

  my $transcript_adaptor = $self->adaptor()->db()->get_TranscriptAdaptor();
  my @transcripts;
  while ($sth->fetch) {
    push @transcripts, $transcript_adaptor->fetch_by_dbID($transcript_id);
  }

  return \@transcripts;

}



=head2 transcripts_regulated_by_same_factor

  Arg [1]    : none
  Example    : my $transcipts = $rf->transcripts_regulated_by_same_factor();
  Description: Return a list of all transcripts that are regulated by any
               RegulatoryFeature with the same factor as this one. If a transcript
               is regulated by more than one feature, it will only appear
               once in the returned list.
  Returntype : listREF of Bio::EnsEMBL::RegulatoryFeatures
  Exceptions : none
  Caller     : ?
  Status     : At Risk
             : under development

=cut

sub transcripts_regulated_by_same_factor {

  my $self = shift;


  # fetch all regulatory features by factor
  my $features = $self->adaptor()->fetch_all_by_factor($self->factor());

  # build cumulative list of transcripts, removing duplicates
  my %transcripts;
  foreach my $feature (@$features) {
    foreach my $feature_transcript (@{$feature->regulated_transcripts()}) {
      bless $feature_transcript, "Bio::EnsEMBL::Transcript";
      if (!exists($transcripts{$feature_transcript->dbID()})) {
	$transcripts{$feature_transcript->dbID()} = $feature_transcript;
      }
    }

  }

  my @tr = values %transcripts;

  return \@tr;

}

1;
