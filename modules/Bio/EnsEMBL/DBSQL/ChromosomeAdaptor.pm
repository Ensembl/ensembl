#
# Ensembl module for Bio::EnsEMBL::DBSQL::ChromosomeAdaptor
#
# Cared for by Ewan Birney <birney@ebi.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::DBSQL::ChromosomeAdaptor - DEPRECATED
use Bio::EnsEMBL::SliceAdator instead

=head1 DESCRIPTION

This class is deprecated.  SliceAdaptor should be used instead.

=head1 AUTHOR - Ewan Birney

This modules is part of the Ensembl project http://www.ensembl.org

Email birney@ebi.ac.uk

=cut

package Bio::EnsEMBL::DBSQL::ChromosomeAdaptor;
use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Utils::Exception qw(deprecate);
use Bio::EnsEMBL::DBSQL::SliceAdaptor;
use Bio::EnsEMBL::Chromosome;

@ISA = qw(Bio::EnsEMBL::DBSQL::SliceAdaptor);


=head2 fetch_by_dbID

  Description: DEPRECATED Use SliceAdaptor::fetch_by_seq_region_id instead

=cut

sub fetch_by_dbID {
  my ($self,$id) = @_;

  deprecate('Use SliceAdaptor::fetch_by_seq_region_id instead.');

  my $chr = $self->SUPER::fetch_by_seq_region_id($id);

  #change blessing to chromosome and reset adaptor
  bless $chr, 'Bio::EnsEMBL::Chromosome';
  $chr->adaptor($self);
  return $chr;
}


=head2 fetch_by_chr_name

  Description: DEPRECATED Use SliceAdaptor::fetch_by_region instead

=cut

sub fetch_by_chr_name{
  my ($self,$chr_name) = @_;

  deprecate('Use SliceAdaptor::fetch_by_region instead.');

  my $top_level = $self->db->get_CoordSystemAdaptor->fetch_top_level();

  my $chr = $self->SUPER::fetch_by_region($top_level->name(),
                                          $chr_name,
                                          undef,
                                          undef,
                                          undef,
                                          $top_level->version);
  #change blessing and adaptor
  bless $chr, 'Bio::EnsEMBL::Chromosome';
  $chr->adaptor($self);
  return $chr;
}


=head2 fetch_all

  Description: DEPRECATED Use SliceAdaptor::fetch_all instead

=cut

sub fetch_all {
  my($self) = @_;

  deprecate('Use SliceAdaptor::fetch_all instead.');

  my $top_level = $self->db->get_CoordSystemAdaptor->fetch_top_level();
  my $chrs = $self->SUPER::fetch_all($top_level->name, $top_level->version);

  foreach my $chr (@$chrs) {
    bless $chr, 'Bio::EnsEMBL::Chromosome';
    $chr->adaptor($self);
  }

  return $chrs;
}


=head2 store

  Description: DEPRECATED no replacement exists

=cut

sub store{
  ### XXX Should there be a replacement for this?
  deprecate('No replacement has been implemented.');
}

1;
