
#
# Ensembl module for Bio::EnsEMBL::DBSQL::AssemblyMapperAdaptor
#
# Written by Ewan Birney <birney@ebi.ac.uk>
#
# Copyright EnsEMBL
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=Head1 NAME

Bio::EnsEMBL::DBSQL::AssemblyMapperAdaptor

=head1 SYNOPSIS

  my $db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(...);

  my $asma = $dba->get_AssemblyMapperAdaptor();
  my $csa  = $dba->get_CoordSystemAdaptor();

  my $chr33_cs   = $csa->fetch_by_name('chromosome', 'NCBI33');
  my $chr34_cs   = $csa->fetch_by_name('chromosome', 'NCBI34');
  my $ctg_cs     = $csa->fetch_by_name('contig');
  my $clone_cs   = $csa->fetch_by_name('clone');

  my $chr_ctg_mapper =
    $asma->fetch_by_CoordSystems($chr33_cs, $ctg_cs);

  my $ncbi33_ncbi34_mapper =
    $asm_adptr->fetch_by_CoordSystems($chr33,$chr34);

  my $ctg_clone_mapper =
    $asm_adptr->fetch_by_CoordSystems($ctg_cs,$clone_cs);


=head1 DESCRIPTION

Adaptor for handling Assembly mappers.  This is a
I<Singleton> class.  ie: There is only one per
database (C<DBAdaptor>).

This is used to retrieve mappers between any two coordinate systems whose
makeup is described by the assembly table.  Currently only pairwise mapping
is supported (i.e there must be an explicit relationship between the coordinate
systems in the assembly table), but in the future 'chained' mapping between
coordinate systems with indirect relationships may be possible.

=head1 CONTACT

Post general queries to B<ensembl-dev@ebi.ac.uk>

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


package Bio::EnsEMBL::DBSQL::AssemblyMapperAdaptor;
use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::AssemblyMapper;

use Bio::EnsEMBL::Utils::Exception qw(deprecate throw);

@ISA = qw(Bio::EnsEMBL::DBSQL::BaseAdaptor);

=head2 new

  Arg [1]    : Bio::EnsEMBL::DBAdaptor $dbadaptor the adaptor for
               the database this assembly mapper is using.
  Example    : my $asma = new Bio::EnsEMBL::AssemblyMapperAdaptor($dbadaptor);
  Description: Creates a new AssemblyMapperAdaptor object
  Returntype : Bio::EnsEMBL::DBSQL::AssemblyMapperAdaptor
  Exceptions : none
  Caller     : Bio::EnsEMBL::DBSQL::DBAdaptor

=cut

sub new {
  my($class, $dbadaptor) = @_;

  my $self = $class->SUPER::new($dbadaptor);

  $self->{'_asm_mapper_cache'} = {};
  $self->{'_sr_id_cache'} = {};

  return $self;
}


=head2 fetch_by_CoordSystems

  Arg [1]    : Bio::EnsEMBL::CoordSystem $cs1
               One of the coordinate systems to retrieve the mapper
               between
  Arg [2]    : Bio::EnsEMBL::CoordSystem $cs2
               The other coordinate system to map between
  Description: Retrieves an Assembly mapper for two coordinate
               systems whose relationship is described in the
               assembly table.

               The ordering of the coodinate systems is arbitrary.
               The following two statements are equivalent:
               $mapper = $asma->fetch_by_coord_systems($cs1,$cs2);
               $mapper = $asma->fetch_by_coord_systems($cs2,$cs1);
  Returntype : Bio::EnsEMBL::AssemblyMapper
  Exceptions : none
  Caller     : general

=cut

sub fetch_by_CoordSystems {
  my $self = shift;
  my $cs1  = shift;
  my $cs2  = shift;

  if(!$cs1 || !ref($cs1) || !$cs1->isa('Bio::EnsEMBL::CoordSystem')) {
    throw("cs1 argument must be a Bio::EnsEMBL::CoordSystem");
  }
  if(!$cs2 || !ref($cs2) || !$cs2->isa('Bio::EnsEMBL::CoordSystem')) {
    throw("cs2 argument must be a Bio::EnsEMBL::CoordSystem");
  }

  if($cs1->equals($cs2)) {
    throw("Cannot create mapper between same coord systems: " .
          $cs1->name . " " . $cs1->version);
  }

  my $csa = $self->db->get_CoordSystemAdaptor();

  #retrieve the shortest possible mapping path between these systems
  my @mapping_path = @{$csa->get_mapping_path($cs1,$cs2)};

  if(!@mapping_path) {
    throw("There is no mapping defined between these coord systems:\n" .
          $cs1->name() . " " . $cs1->version() . " and " . $cs1->name() . " " .
          $cs2->version());
  }

  my $key = join(':', map({$_->dbID()} @mapping_path));

  my $asm_mapper = $self->{'_asm_mapper_cache'}->{$key};

  return $asm_mapper if($asm_mapper);

  if(@mapping_path > 2) {
    throw("Only explicit (1 step) coordinate system mapping is currently\n" .
          "supported.  Mapping between \n" .
          $cs1->name() . " " . $cs1->version() . " and " .
          $cs2->name() . " " . $cs2->version() .
          "\nrequires ". (scalar(@mapping_path)-1) . " steps.");
  }

  $asm_mapper = Bio::EnsEMBL::AssemblyMapper->new(@mapping_path);

  $self->{'_asm_mapper_cache'}->{$key} = $asm_mapper;

  return $asm_mapper;
}



=head2 register_assembled

  Arg [1]    : Bio::EnsEMBL::AssemblyMapper $asm_mapper
	       A valid AssemblyMapper object
  Arg [2]    : string $asm_seq_region
               The name of the seq_region to be registered
  Arg [3]    : int $asm_start
               The start of the region to be registered
  Arg [4]    : int $asm_end
               The end of the region to be registered
  Description: Declares an assembled region to the AssemblyMapper.
               This extracts the relevant data from the assembly
               table and stores it in Mapper internal to the $asm_mapper.
               It therefore must be called before any mapping is
               attempted on that region. Otherwise only gaps will
               be returned.  Note that the AssemblyMapper automatically
               calls this method when the need arises.
  Returntype : none
  Exceptions : throw if the seq_region to be registered does not exist
               or if it associated with multiple assembled pieces (bad data
               in assembly table)
  Caller     : Bio::EnsEMBL::AssemblyMapper

=cut

my $CHUNKFACTOR = 20;  # 2^20 = approx. 10^6

sub register_assembled {
  my $self = shift;
  my $asm_mapper = shift;
  my $asm_seq_region = shift;
  my $asm_start      = shift;
  my $asm_end        = shift;


  my $asm_cs_id = $asm_mapper->assembled_CoordSystem->dbID();
  my $cmp_cs_id = $asm_mapper->component_CoordSystem->dbID();

  #split up the region to be registered into fixed chunks
  #this allows us to keep track of regions that have already been
  #registered and also works under the assumption that if a small region
  #is requested it is likely that other requests will be made in the
  #vicinity (the minimum size registered the chunksize (2^chunkfactor)

  my @chunk_regions;
  {
    #determine span of chunks
    #bitwise shift right is fast and easy integer division
    my $start_chunk = $asm_start >> $CHUNKFACTOR;
    my $end_chunk   = $asm_end   >> $CHUNKFACTOR;

    #find regions of continuous unregistered chunks
    my $i;
    my ($begin_chunk_region,$end_chunk_region);
    for ($i = $start_chunk; $i <= $end_chunk; $i++) {
      if($asm_mapper->have_registered_assembled($asm_seq_region, $i)) {
        if(defined($begin_chunk_region)) {
          #this is the end of an unregistered region.
          my $region = [$begin_chunk_region   << $CHUNKFACTOR,
                        $end_chunk_region     << $CHUNKFACTOR];
          push @chunk_regions, $region;
          $begin_chunk_region = $end_chunk_region = undef;
        }
      } else {
        $begin_chunk_region = $i if(!defined($begin_chunk_region));
        $end_chunk_region   = $i+1;
        $asm_mapper->register_assembled($asm_seq_region,$i);
      }
    }

    #the last part may have been an unregistered region too
    if(defined($begin_chunk_region)) {
      my $region = [$begin_chunk_region << $CHUNKFACTOR,
                    $end_chunk_region   << $CHUNKFACTOR];
      push @chunk_regions, $region;
    }
  }

  return if(!@chunk_regions);

  my $asm_seq_region_id = $self->{'_sr_id_cache'}->{"$asm_seq_region:$asm_cs_id"};

  if(!$asm_seq_region_id) {
    # Get the seq_region_id via the name.  This would be quicker if we just
    # used internal ids instead but stored but then we lose the ability
    # the transform accross databases with different internal ids

    my $sth = $self->prepare("SELECT seq_region_id" .
                             "FROM   seq_region" .
                             "WHERE  name = ? AND coord_system_id = ?");

    $sth->execute($asm_seq_region, $asm_cs_id);

    if(!$sth->rows() == 1) {
      throw("Ambiguous or non-existant seq_region [$asm_seq_region]" .
            "in coord system $asm_cs_id");
    }

    ($asm_seq_region_id) = $sth->fetchrow_array();
    $self->{'_sr_id_cache'}->{"$asm_seq_region:$asm_cs_id"} =
      $asm_seq_region_id;

    $sth->finish();
  }

  # Retrieve the description of how the assembled region is made from
  # component regions for each of the continuous blocks of unregistered,
  # chunked regions

  my $q = qq{
      SELECT
         asm.cmp_start,
         asm.cmp_end,
         asm.cmp_seq_region_id,
         sr.name,
         asm.ori,
         asm.asm_start,
         asm.asm_end
      FROM
         assembly asm, seq_region sr
      WHERE
         asm.asm_seq_region_id = ? AND
         ? <= asm.asm_end AND
         ? >= asm.asm_start AND
         asm.cmp_seq_region_id = sr.seq_region_id AND
	 sr.coord_system_id = ?
   };

  my $sth = $self->prepare($q);

  foreach my $region (@chunk_regions) {
    my($region_start, $region_end) = @$region;
    $sth->execute($asm_seq_region_id, $region_start, $region_end, 
		  $asm_mapper->component_CoordSystem->dbID());

    my($cmp_start, $cmp_end, $cmp_seq_region_id, $cmp_seq_region, $ori);
    $sth->bind_columns(\$cmp_start, \$cmp_end, \$cmp_seq_region_id,
                       \$cmp_seq_region, \$ori);

    #
    # Load the unregistered regions of the mapper
    #
    while($sth->fetch()) {
      next if($asm_mapper->have_registered_component($cmp_seq_region));
      $asm_mapper->register_component($cmp_seq_region);
      $asm_mapper->mapper->add_map_coordinates(
                 $cmp_seq_region, $cmp_start, $cmp_end,
                 $ori,
                 $asm_seq_region, $region_start, $region_end);
      $self->{'_sr_id_cache'}->{"$cmp_seq_region:$cmp_cs_id"} =
        $cmp_seq_region_id;
    }
  }

  $sth->finish();
}


=head2 register_component

  Arg [1]    : Bio::EnsEMBL::AssemblyMapper $asm_mapper
	       A valid AssemblyMapper object
  Arg [2]    : string $cmp_seq_region
               The name of the seq_region to be registered
  Description: Declares a component region to the AssemblyMapper.
               This extracts the relevant data from the assembly
               table and stores it in Mapper internal to the $asm_mapper.
               It therefore must be called before any mapping is
               attempted on that region. Otherwise only gaps will
               be returned.  Note that the AssemblyMapper automatically
               calls this method when the need arises.
  Returntype : none
  Exceptions : throw if the seq_region to be registered does not exist
               or if it associated with multiple assembled pieces (bad data
               in assembly table)
  Caller     : Bio::EnsEMBL::AssemblyMapper

=cut

sub register_component {
  my $self = shift;
  my $asm_mapper = shift;
  my $cmp_seq_region = shift;

  my $cmp_cs_id = $asm_mapper->component_CoordSystem()->dbID();

  #do nothing if this region is already registered
  return if($asm_mapper->have_registered_component($cmp_seq_region));


  my $cmp_seq_region_id = 
    $self->{'_sr_id_cache'}->{"$cmp_seq_region:$cmp_cs_id"};


  if(!$cmp_seq_region_id) {
    my $sth = $self->prepare("SELECT seq_region_id" .
                             "FROM   seq_region" .
                             "WHERE  name = ? AND coord_system_id = ?");

    $sth->execute($cmp_seq_region, $cmp_cs_id);

    if(!$sth->rows() == 1) {
      throw("Ambiguous or non-existant seq_region [$cmp_seq_region] " .
            "in coord system $cmp_cs_id");
    }

    ($cmp_seq_region_id) = $sth->fetchrow_array();
    $self->{'_sr_id_cache'}->{"$cmp_seq_region:$cmp_cs_id"} =
      $cmp_seq_region_id;

    $sth->finish();
  }

  # Determine what part of the assembled region this component region makes up

  my $q = qq{
      SELECT
         asm.asm_start,
         asm.asm_end,
         asm.seq_region_id,
         sr.name
      FROM
         assembly asm, seq_region sr
      WHERE
         asm.cmp_seq_region_id = ? AND
         asm.asm_seq_region_id = sr.seq_region_id
   };

  my $sth = $self->prepare($q);
  $sth->execute($cmp_seq_region_id);

  if($sth->rows() == 0) {
    #this component is not used in the assembled part i.e. gap
    $asm_mapper->register_component($cmp_seq_region_id);
    return;
  }

  #something is wonky if there are multiple assembled regions for a component
  if($sth->rows() != 1) {
    throw("Multiple assembled regions for single " .
          "component region cmp_seq_region_id=[$cmp_seq_region_id]");
  }

  my ($asm_start, $asm_end, $asm_seq_region_id, $asm_seq_region) =
    $sth->fetchrow_array();

  $self->{'_sr_id_cache'}->{$asm_seq_region} = $asm_seq_region_id;

  $sth->finish();

  # Register the corresponding assembled region. This allows a us to
  # register things in assembled chunks which allows us to:
  # (1) Keep track of what assembled regions are registered
  # (2) Use locality of reference (if they want something in same general
  #     region it will already be registered).

  $self->register_assembled($asm_mapper,$asm_seq_region,$asm_start,$asm_end);
}



=head2 seq_regions_to_ids

  Arg [1]    : Bio::EnsEMBL::CoordSystem $coord_system
  Arg [2]    : listref of strings $seq_regions
  Example    : my @ids = @{$asma->seq_regions_to_ids($coord_sys, \@seq_regs)};
  Description: Converts a list of seq_region names to internal identifiers
               using the internal cache that has accumulated while registering
               regions for AssemblMappers.  Note that this only works for
               regions which have been registered and is only intended to be
               called from the AssemblMapper objects created by this adaptor.
  Returntype : listref of ints
  Exceptions : none
  Caller     : general

=cut

sub seq_regions_to_ids {
  my $self = shift;
  my $coord_system = shift;
  my $seq_regions = shift;

  my $cs_id = $coord_system->dbID();

  my @out;

  foreach my $sr (@$seq_regions) {
    my $id = $self->{'_sr_id_cache'}->{"$sr:$cs_id"};
    throw("Seq_region [$sr] with coord_system [$cs_id] not found") if(!$id);
    push @out, $id;
  }

  return \@out;
}



=head2 register_region

  Description: DEPRECATED use register_assembled instead

=cut

sub register_region{
  my ($self, $assmapper, $type, $chr_name, $start, $end) = @_;

  deprecate('Use register_assembled instead');

  $self->register_assembled($assmapper, $chr_name, $start, $end);
}


=head2 register_contig

  Description: DEPRECATED use register_component instead

=cut

sub register_contig {
   my ($self, $assmapper, $type, $contig_id ) = @_;

   deprecate('Use register_component instead');

   #not sure if the use is passing in a seq_region_name or a
   #seq_region_id...
   register_component($assmapper, $contig_id);
}


=head2 fetch_by_type

  Description: DEPRECATED use fetch_by_CoordSystems instead

=cut

sub fetch_by_type{
  my ($self,$type) = @_;

  deprecate('Use fetch_by_coord_systems instead');

  #assume that what the user wanted was a mapper between the sequence coord
  #level and the top coord level

  my $csa = $self->db()->get_CoordSystemAdaptor();

  my $cs1 = $csa->fetch_top_level($type);
  my $cs2 = $csa->fetch_sequence_level();

  return $self->fetch_by_coord_systems($cs1,$cs2);
}



1;
