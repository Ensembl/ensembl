=head1 LICENSE

  Copyright (c) 1999-2009 The European Bioinformatics Institute and
  Genome Research Limited.  All rights reserved.

  This software is distributed under a modified Apache license.
  For license details, please see

    http://www.ensembl.org/info/about/code_licence.html

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <ensembl-dev@ebi.ac.uk>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk@ensembl.org>.

=cut

=head1 NAME

Bio::EnsEMBL::DBSQL::AssemblyExceptionFeatureAdaptor

=head1 SYNOPSIS

  my $assembly_exception_feature_adaptor =
    $database_adaptor->get_AssemblyExceptionFeatureAdaptor();

  @assembly_exception_features =
    $assembly_exception_feature_adaptor->fetch_all_by_Slice($slice);

=head1 DESCRIPTION

Assembly Exception Feature Adaptor - database access for assembly
exception features.

=head1 METHODS

=cut

package Bio::EnsEMBL::DBSQL::AssemblyExceptionFeatureAdaptor;

use strict;
use warnings;
no warnings qw(uninitialized);

use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor;
use Bio::EnsEMBL::AssemblyExceptionFeature;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Cache;

our @ISA = qw(Bio::EnsEMBL::DBSQL::BaseAdaptor);

# set the number of slices you'd like to cache
our $ASSEMBLY_EXCEPTION_FEATURE_CACHE_SIZE = 100;

=head2 new

  Arg [1]    : list of args @args
               Superclass constructor arguments
  Example    : none
  Description: Constructor which just initializes internal cache structures
  Returntype : Bio::EnsEMBL::DBSQL::AssemblyExceptionFeatureAdaptor
  Exceptions : none
  Caller     : implementing subclass constructors
  Status     : Stable

=cut

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;

  my $self = $class->SUPER::new(@_);

  # initialize an LRU cache for slices
  my %cache;
  tie(%cache, 'Bio::EnsEMBL::Utils::Cache',
    $ASSEMBLY_EXCEPTION_FEATURE_CACHE_SIZE);

  $self->{'_aexc_slice_cache'} = \%cache;

  return $self;
}

=head2 fetch_all

  Arg [1]    : none
  Example    : my @axfs = @{$axfa->fetch_all()};
  Description: Retrieves all assembly exception features which are in the
               database and builds internal caches of the features.
  Returntype : reference to list of Bio::EnsEMBL::AssemblyExceptionFeatures
  Exceptions : none
  Caller     : fetch_by_dbID, fetch_by_Slice
  Status     : Stable

=cut

sub fetch_all {
  my $self = shift;

  # this is the "global" cache for all assembly exception features in the db
  if(defined($self->{'_aexc_cache'})) {
    return $self->{'_aexc_cache'};
  }

  my $statement = qq(
  SELECT    ae.assembly_exception_id,
            ae.seq_region_id,
            ae.seq_region_start,
            ae.seq_region_end,
            ae.exc_type,
            ae.exc_seq_region_id,
            ae.exc_seq_region_start,
            ae.exc_seq_region_end,
            ae.ori
  FROM      assembly_exception ae,
            coord_system cs,
            seq_region sr
  WHERE     cs.species_id = ?
    AND     sr.coord_system_id = cs.coord_system_id
    AND     sr.seq_region_id = ae.seq_region_id);

  my $sth = $self->prepare($statement);

  $sth->bind_param( 1, $self->species_id(), SQL_INTEGER );

  $sth->execute();

  my ($ax_id, $sr_id, $sr_start, $sr_end,
      $x_type, $x_sr_id, $x_sr_start, $x_sr_end, $ori);

  $sth->bind_columns(\$ax_id, \$sr_id, \$sr_start, \$sr_end,
                     \$x_type, \$x_sr_id, \$x_sr_start, \$x_sr_end, \$ori);

  my @features;
  my $sa = $self->db()->get_SliceAdaptor();

  $self->{'_aexc_dbID_cache'} = {};

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
    $self->{'_aexc_dbID_cache'}->{$ax_id} = $a;

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

  $self->{'_aexc_cache'} = \@features;
  
  return \@features;
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
  Status     : Stable

=cut

sub fetch_by_dbID {
  my $self = shift;
  my $dbID = shift;

  if(!exists($self->{'_aexc_dbID_cache'})) {
    # force loading of cache
    $self->fetch_all();
  }

  return $self->{'_aexc_dbID_cache'}->{$dbID};
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
  Status     : Stable

=cut

sub fetch_all_by_Slice {
  my $self = shift;
  my $slice = shift;

  my $key= uc($slice->name());

  # return features from the slice cache if present
  if(exists($self->{'_aexc_slice_cache'}->{$key})) {
    return $self->{'_aexc_slice_cache'}->{$key};
  }

  my $all_features = $self->fetch_all();

  my $mcc = $self->db()->get_MetaCoordContainer();
  my $css = $mcc->fetch_all_CoordSystems_by_feature_type('assembly_exception');

  my @features;

  my $ma = $self->db()->get_AssemblyMapperAdaptor();

  foreach my $cs (@$css) {
    my $mapper;
    if($cs->equals($slice->coord_system)) {
      $mapper = undef;
    } else {
      $mapper = $ma->fetch_by_CoordSystems($cs,$slice->coord_system());
    }

    push @features, @{ $self->_remap($all_features, $mapper, $slice) };
  }

  $self->{'_aexc_slice_cache'}->{$key} = \@features;

  return \@features;
}


#
# Given a list of features checks if they are in the correct coord system
# by looking at the first features slice.  If they are not then they are
# converted and placed on the slice.
#
# Note that this is a re-implementation of a method with the same name in
# BaseFeatureAdaptor, and in contrast to the latter which maps features in
# place, this method returns a remapped copy of each feature. The reason for
# this is to get around conflicts with caching.
#
sub _remap {
  my ($self, $features, $mapper, $slice) = @_;

  # check if any remapping is actually needed
  if(@$features && (!$features->[0]->isa('Bio::EnsEMBL::Feature') ||
                    $features->[0]->slice == $slice)) {
    return $features;
  }

  # remapping has not been done, we have to do our own conversion from
  # to slice coords

  my @out;

  my $slice_start = $slice->start();
  my $slice_end   = $slice->end();
  my $slice_strand = $slice->strand();
  my $slice_cs    = $slice->coord_system();

  my ($seq_region, $start, $end, $strand);

  my $slice_seq_region = $slice->seq_region_name();

  foreach my $f (@$features) {
    # since feats were obtained in contig coords, attached seq is a contig
    my $fslice = $f->slice();
    if(!$fslice) {
      throw("Feature does not have attached slice.\n");
    }
    my $fseq_region = $fslice->seq_region_name();
    my $fcs = $fslice->coord_system();

    if(!$slice_cs->equals($fcs)) {
      # slice of feature in different coord system, mapping required
      ($seq_region, $start, $end, $strand) =
        $mapper->fastmap($fseq_region,$f->start(),$f->end(),$f->strand(),$fcs);

      # undefined start means gap
      next if(!defined $start);
    } else {
      $start      = $f->start();
      $end        = $f->end();
      $strand     = $f->strand();
      $seq_region = $f->slice->seq_region_name();
    }

    # maps to region outside desired area
    next if ($start > $slice_end) || ($end < $slice_start) || 
      ($slice_seq_region ne $seq_region);

    # create new copies of successfully mapped feaatures with shifted start,
    # end and strand
    my ($new_start, $new_end);
    if($slice_strand == -1) {
      $new_start = $slice_end - $end + 1;
      $new_end = $slice_end - $start + 1;
    } else {
      $new_start = $start - $slice_start + 1;
      $new_end = $end - $slice_start + 1;
    }
    
    push @out, Bio::EnsEMBL::AssemblyExceptionFeature->new(
            '-dbID'            => $f->dbID,
            '-start'           => $new_start,
            '-end'             => $new_end,
            '-strand'          => $strand * $slice_strand,
            '-adaptor'         => $self,
            '-slice'           => $slice,
            '-alternate_slice' => $f->alternate_slice,
            '-type'            => $f->type,
    );
  }

  return \@out;
}




1;
