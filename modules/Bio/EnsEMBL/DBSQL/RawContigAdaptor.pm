#
# EnsEMBL module for Bio::EnsEMBL::DBSQL::RawContigAdaptor
#
#
# Copyright Imre Vastrik
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME Bio::EnsEMBL::DBSQL::RawContigAdaptor
- DEPRECATED use Bio::EnsEMBL::DBSQL::SliceAdaptor instead

=head1 CONTACT

Post questions to the EnsEMBL developer list <ensembl-dev@ebi.ac.uk>

=cut

use strict;

package Bio::EnsEMBL::DBSQL::RawContigAdaptor;

use vars qw(@ISA);

use Bio::EnsEMBL::DBSQL::SliceAdaptor;
use Bio::EnsEMBL::RawContig;
use Bio::EnsEMBL::Utils::Exception qw(throw deprecate warning);

@ISA = qw(Bio::EnsEMBL::DBSQL::SliceAdaptor);


=head2 fetch_by_dbID

  Description: DEPRECATED use SliceAdaptor::new instead

=cut

sub new {
  my $class = shift;

  deprecate("Bio::EnsEMBL::DBSQL::RawContigAdaptor is a deprecated class.\n" .
            "Use the Bio::EnsEMBL::DBSQL::SliceAdaptor class instead");
  return $class->SUPER::new(@_);
}


=head2 fetch_by_dbID

  Description: DEPRECATED use SliceAdaptor::fetch_by_seq_region_id instead

=cut

sub fetch_by_dbID {
  my $self = shift;
  my $dbID = shift;

  deprecate('Use SliceAdaptor::fetch_by_seq_region_id instead');

  my $contig = $self->fetch_by_seq_region_id($dbID);

  #change blessing to contig so old contig functions still work
  return bless($contig, 'Bio::EnsEMBL::RawContig');
}


=head2 fetch_all

  Description: DEPRECATED use SliceAdaptor::fetch_all instead

=cut

sub fetch_all {
  my $self  = shift;

  deprecate('Use SliceAdaptor::fetch_all instead');

  #assume that by 'contig' what is actually meant is sequence level
  my $csa = $self->db->get_CoordSystemAdaptor();
  my $cs = $csa->fetch_sequence_level();

  my $result = $self->SUPER::fetch_all($cs->name,$cs->version);

  #change blessings to RawContig so old contig functions still work
  foreach my $contig (@$result) {
    bless $contig, 'Bio::EnsEMBL::RawContig';
  }

  return $result;
}

=head2 fetch_by_name

  Description: DEPRECATED use SliceAdaptor::fetch_by_region instead

=cut

sub fetch_by_name {
  my $self = shift;
  my $name = shift;

  deprecate('Use SliceAdaptor::fetch_by_region instead');

  #assume that by 'contig' what is actually meant is seq level
  my $csa = $self->db()->get_CoordSystemAdaptor();
  my $cs = $csa->fetch_sequence_level();

  my $contig = $self->fetch_by_region($cs->name,
                                      $name,
                                      undef,
                                      undef,
                                      undef,
                                      $cs->version);

  #change blessing to RawContig so old contig functions still work
  return bless $contig, 'Bio::EnsEMBL::RawContig';
}


=head2 fetch_filled_by_dbIDs

  Description: DEPRECATED use SliceAdaptor::fetch_by_seq_region_id instead

=cut

sub fetch_filled_by_dbIDs {
  my $self = shift;
  my @contig_ids = @_;

  deprecate('Use SliceAdaptor::fetch_by_seq_region_id instead');

  my %result = ();

  foreach my $contig_id (@contig_ids) {
    my $contig = $self->fetch_by_seq_region_id($contig_id);
    #change blessing to Bio::EnsEMBL::RawContig so old functions still work
    bless $contig, 'Bio::EnsEMBL::RawContig';
    $result{$contig_id} = $contig;
  }

  return \%result;
}


=head2 fetch_all_by_Clone

  Description: DEPRECATED use Slice::project instead

=cut

sub fetch_all_by_Clone {
  my $self = shift;
  my $clone = shift;

  deprecate('Use Slice::project instead');

  #project this clone to contig coord system
  #assume that by 'contig' what is actually meant is seq level
  my $csa = $self->db()->get_CoordSystemAdaptor();
  my $cs = $csa->fetch_sequence_level();
  my $projection = $clone->project($cs->name, $cs->version);

  my @out;
  foreach my $seg (@$projection) {
    my $contig = $seg->[2];
    #rebless into contig so old method calls still work
    bless $contig, 'Bio::EnsEMBL::RawContig';
    push @out, $contig;
  }

  return \@out;
}



=head2 fetch_attributes

  Description: DEPRECATED This method is no longer needed and does nothing

=cut

sub fetch_attributes {
  deprecate('This method is no longer needed and does nothing');
  return;
}



=head2 store

  Description: DEPRECATED XXX Replacement may need to be written

=cut

sub store {
  ### XXX Does a replacement need to be implemented?
  deprecate('No replacement has been implemented');
  return 0;
}



=head2 remove

  Description: DEPRECATED XXX Replacement may need to be written

=cut

sub remove {
  ### XXX Does a replacement need to be implemented?
  deprecate('No replacement has been implemened');
  return 0;
}



sub fetch_all_names{
  my ($self) = @_;

  deprecate('Use SliceAdaptor::fetch_all instead');
  my @names = map {$_->seq_region_name} @{$self->fetch_all()};
  return \@names;
}

1;
