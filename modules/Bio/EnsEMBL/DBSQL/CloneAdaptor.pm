#
# BioPerl module for DB::Clone
#
# Cared for by EnsEMBL (www.ensembl.org)
#
# Copyright GRL and EBI
#
# You may distribute this module under the same terms as perl itself

=head1 NAME

Bio::EnsEMBL::DBSQL::CloneAdaptor - DEPRECATED use Bio::EnsEMBL::SliceAdaptor
instead

=head1 DESCRIPTION

DEPRECATED - Use Bio::EnsEMBL::DBSQL::SliceAdaptor to create slices on clones
instead

=head1 AUTHOR - Ewan Birney

This modules is part of the Ensembl project http://www.ensembl.org

Email birney@ebi.ac.uk

=cut


package Bio::EnsEMBL::DBSQL::CloneAdaptor;
use vars qw(@ISA);
use strict;
use Bio::EnsEMBL::DBSQL::SliceAdaptor;
use Bio::EnsEMBL::Clone;
use Bio::EnsEMBL::Utils::Exception qw( deprecate throw warning );

@ISA = qw(Bio::EnsEMBL::DBSQL::SliceAdaptor);

=head2 new

  Description: DEPRECATED.  Use Bio::EnsEMBL::DBSQL::SliceAdaptor instead

=cut

sub new {
  my $class = shift;

  deprecate("Bio::EnsEMBL::DBSQL::CloneAdaptor is a deprecated class\n" .
            "Use Bio::EnsEMBL::DBSQL::SliceAdaptor instead");

  return $class->SUPER::new(@_);
}


=head2 fetch_by_accession

  Description: DEPRECATED. Use SliceAdaptor::fetch_by_region instead

=cut

sub fetch_by_accession {
  my ($self,$acc) = @_;

  deprecate('Use SliceAdaptor::fetch_by_region instead');

  #this unfortunately needs a version on the end to work
  if(! ($acc =~ /\./)) {
    my $sth = $self->prepare("SELECT sr.name " .
                             "FROM   seq_region sr, coord_system cs " .
                             "WHERE  cs.name = 'clone' " .
                             "AND    cs.coord_system_id = sr.coord_system_id ".
                             "AND    sr.name LIKE '$acc.%'");
    $sth->execute();
    if(!$sth->rows()) {
      $sth->finish();
      throw("Clone $acc not found in database");
    }

    ($acc) = $sth->fetchrow_array();

    $sth->finish();
  }

  my $clone = $self->fetch_by_region('clone', $acc);

  #rebless slice into clone so old methods still work
  return bless($clone, 'Bio::EnsEMBL::Clone');
}


=head2 fetch_by_accession_version

  Description: DEPRECATED. Use SliceAdaptor::fetch_by_region instead

=cut

sub fetch_by_accession_version { 
  my ($self, $acc, $ver) = @_;

  deprecate('Use SliceAdaptor::fetch_by_region instead');

  my $clone = $self->fetch_by_region('clone', "$acc.$ver");

  #rebless slice into clone so old methods still work
  return bless($clone, 'Bio::EnsEMBL::Clone');
}


=head2 fetch_by_name

  Description: DEPRECATED. Use SliceAdaptor::fetch_by_region instead

=cut

sub fetch_by_name {
  my ($self, $name) = @_;
  deprecate('Use SliceAdaptor::fetch_by_region instead');
  return $self->fetch_by_accession($name);
}


=head2 fetch_by_dbID

  Description: DEPRECATED. Use SliceAdaptor::fetch_by_seq_region_id instead

=cut

sub fetch_by_dbID {
  my ($self,$id) = @_;

  my $clone = $self->fetch_by_seq_region_id($id);

  return bless $clone, 'Bio::EnsEMBL::Clone';
}


=head2 fetch_all

  Description: DEPRECATED. Use SliceAdaptor::fetch_all instead

=cut

sub fetch_all {
  my $self = shift;

  my $clones = $self->SUPER::fetch_all('clone');

  map {bless $_, 'Bio::EnsEMBL::Clone'} @$clones;

  return $clones;
}



=head2 list_embl_version_by_accesssion

  Description: DEPRECATED There is currently no replacement for this method
               Perhaps there should be?

=cut

sub list_embl_version_by_accession {
  deprecate('There is currently no replacement for this method');

  ### XXX Perhaps there should be a replacement?
  return ();
}


=head2 remove

  Description: DEPRECATED There is currently no replacement for this method
               Perhaps there should be?


=cut

sub remove {
  deprecate('There is currently no replacement for this method');

  ### XXX Perhaps there should be a replacement?
  return 0;
}



=head2 store

  Description: DEPRECATED There is currently no replacement for this method
               Perhaps there should be?


=cut

sub store{
  deprecate('There is currently no replacement for this method');

  ### XXX Perhaps there should be a replacement?
  return 0;
}


1;
