#
# Ensembl module for Bio::EnsEMBL::DBSQL::AssemblyExceptionFeatureAdaptor
#
#
# Copyright Ensembl
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::DBSQL::AssemblyExceptionFeatureAdaptor

=head1 SYNOPSIS

my $assembly_exception_feature_adaptor = $database_adaptor->get_AssemblyExceptionFeatureAdaptor();
@assembly_exception_features = $assembly_exception_feature_adaptor->fetch_by_Slice($slice);

=head1 DESCRIPTION

Assembly Exception Feature Adaptor - database access for assembly exception features

=head1 AUTHOR - Glenn Proctor

=head1 METHODS


=cut

package Bio::EnsEMBL::DBSQL::AssemblyExceptionFeatureAdaptor;
use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor;
use Bio::EnsEMBL::AssemblyExceptionFeature;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);

@ISA = qw(Bio::EnsEMBL::DBSQL::BaseAdaptor);



=head2 fetch_all

  Arg [1]    : none
  Example    : my @axfs = @{$axfa->fetch_all()};
  Description: Retrieves all assembly exception features which are in the
               database and builds internal caches of the features.
  Returntype : reference to list of Bio::EnsEMBL::AssemblyExceptionFeatures
  Exceptions : none
  Caller     : fetch_by_dbID, fetch_by_Slice

=cut

sub fetch_all {
  my $self = shift;

  if(defined($self->{'aexc_cache'})) {
    return $self->{'aexc_cache'};
  }

  my $sth = $self->prepare
    ("SELECT assembly_exception_id, seq_region_id, seq_region_start,
             seq_region_end, exc_type, exc_seq_region_id, exc_seq_region_start,
             exc_seq_region_end, ori
      FROM assembly_exception");

  $sth->execute();

  my ($ax_id, $sr_id, $sr_start, $sr_end,
      $x_type, $x_sr_id, $x_sr_start, $x_sr_end, $ori);

  $sth->bind_columns(\$ax_id, \$sr_id, \$sr_start, \$sr_end,
                     \$x_type, \$x_sr_id, \$x_sr_start, \$x_sr_end, \$ori);

  my @features;
  my $sa = $self->db()->get_SliceAdaptor();

  $self->{'aexc_dbID_cache'} = {};

  while($sth->fetch()) {
    my $slice   = $sa->fetch_by_seq_region_id($sr_id);
    my $x_slice = $sa->fetch_by_seq_region_id($x_sr_id);

    # each row creates TWO features, each of which has alternate_slice
    # pointing to the "other" one

    my $a = Bio::EnsEMBL::AssemblyExceptionFeature->new
          ('-dbID'            => $ax_id,
           '-start'           => $sr_start,
           '-end'             => $sr_end,
           '-strand'          => 1,
           '-adaptor'         => $self,
           '-slice'           => $slice,
           '-alternate_slice' => $x_slice->sub_Slice($x_sr_start, $x_sr_end),
           '-type'            => $x_type);

    push @features, $a;
    $self->{'aexc_dbID_cache'}->{$ax_id} = $a;

    push @features, Bio::EnsEMBL::AssemblyExceptionFeature->new
          ('-dbID'            => $ax_id,
           '-start'           => $x_sr_start,
           '-end'             => $x_sr_end,
           '-strand'          => 1,
           '-adaptor'         => $self,
           '-slice'           => $x_slice,
           '-alternate_slice' => $slice->sub_Slice($sr_start, $sr_end),
           '-type'            => "$x_type REF" );
  }

  $sth->finish();

  $self->{'aexc_cache'} = \@features;

  return \@features;
}

#
# cleans up internal caches during garbage collection
#
sub deleteObj {
  my $self = shift;

  $self->SUPER::deleteObj(@_);

  delete $self->{'aexc_cache'};
  delete $self->{'aexc_slice_cache'};
  delete $self->{'aexc_dbID_cache'};

  return;
}



=head2 fetch_by_dbID

  Arg [1]    : int $dbID
  Example    : my $axf = $axfa->fetch_by_dbID(3);
  Description: Retrieves a single assembly exception feature via its internal
               identifier.  Note that this only retrieves one of the two
               assembly exception features which are represented by a single
               row in the assembly_exception table.
  Returntype : Bio::EnsEMBL::AssemblyExceptionFeature
  Exceptions : none
  Caller     : general

=cut

sub fetch_by_dbID {
  my $self = shift;
  my $dbID = shift;

  if(!exists($self->{'aexc_dbID_cache'})) {
    # force loading of cache
    $self->fetch_all();
  }

  return $self->{'aexc_dbID_cache'}->{$dbID};
}


=head2 fetch_all_by_Slice

  Arg [1]    : Bio::EnsEMBL::Slice $slice
  Example    : my @axfs = @{$axfa->fetch_all_by_Slice($slice)};
  Description: Retrieves all assembly exception features which overlap the
               provided slice.  The returned features will be in coordinate
               system of the slice.
  Returntype : reference to list of Bio::EnsEMBL::AssemblyException features
  Exceptions : none
  Caller     : Feature::get_all_alt_locations, general

=cut

sub fetch_all_by_Slice {
  my $self = shift;
  my $slice = shift;

  my $key= uc($slice->name());

  if(exists($self->{'aexc_slice_cache'}->{$key})) {
    return $self->{'aexc_slice_cache'}->{$key};
  }

  my $all_features = $self->fetch_all();

  my $mcc = $self->db()->get_MetaCoordContainer();
  my $css = $mcc->fetch_all_CoordSystems_by_feature_type('assembly_exception');

  my @features;

  my $remap = \&Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor::_remap;
  my $ma = $self->db()->get_AssemblyMapperAdaptor();

  foreach my $cs (@$css) {
    my $mapper;
    if($cs->equals($slice->coord_system)) {
      $mapper = undef;
    } else {
      $mapper = $ma->fetch_by_CoordSystems($cs,$slice->coord_system());
    }

    push @features, @{&$remap($all_features, $mapper, $slice)};
  }

  $self->{'aexc_slice_cache'}->{$key} = \@features;

  return \@features;
}





1;
