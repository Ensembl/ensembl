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

Bio::EnsEMBL::DBSQL::AssemblyMapperAdaptor

=head1 SYNOPSIS

  use Bio::EnsEMBL::Registry;

  Bio::EnsEMBL::Registry->load_registry_from_db(
    -host => 'ensembldb.ensembl.org',
    -user => 'anonymous'
  );

  $asma = Bio::EnsEMBL::Registry->get_adaptor( "human", "core",
    "assemblymapper" );

  $csa = Bio::EnsEMBL::Registry->get_adaptor( "human", "core",
    "coordsystem" );

  my $chr33_cs = $csa->fetch_by_name( 'chromosome', 'NCBI33' );
  my $chr34_cs = $csa->fetch_by_name( 'chromosome', 'NCBI34' );
  my $ctg_cs   = $csa->fetch_by_name('contig');
  my $clone_cs = $csa->fetch_by_name('clone');

  my $chr_ctg_mapper =
    $asma->fetch_by_CoordSystems( $chr33_cs, $ctg_cs );

  my $ncbi33_ncbi34_mapper =
    $asm_adptr->fetch_by_CoordSystems( $chr33, $chr34 );

  my $ctg_clone_mapper =
    $asm_adptr->fetch_by_CoordSystems( $ctg_cs, $clone_cs );


=head1 DESCRIPTION

Adaptor for handling Assembly mappers.  This is a I<Singleton> class.
ie: There is only one per database (C<DBAdaptor>).

This is used to retrieve mappers between any two coordinate systems
whose makeup is described by the assembly table.  Currently one step
(explicit) and two step (implicit) pairwise mapping is supported.  In
one-step mapping an explicit relationship between the coordinate systems
is defined in the assembly table.  In two-step 'chained' mapping no
explicit mapping is present but the coordinate systems must share a
common mapping to an intermediate coordinate system.

=head1 METHODS

=cut

package Bio::EnsEMBL::DBSQL::AssemblyMapperAdaptor;
use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::AssemblyMapper;
use Bio::EnsEMBL::ChainedAssemblyMapper;
use Bio::EnsEMBL::TopLevelAssemblyMapper;

use Bio::EnsEMBL::Utils::Cache; #CPAN LRU cache
use Bio::EnsEMBL::Utils::Exception qw(throw deprecate warning stack_trace_dump);
#use Bio::EnsEMBL::Utils::Exception qw(deprecate throw);
use Bio::EnsEMBL::Utils::SeqRegionCache;

use integer; #do proper arithmetic bitshifts

@ISA = qw(Bio::EnsEMBL::DBSQL::BaseAdaptor);


my $CHUNKFACTOR = 20;  # 2^20 = approx. 10^6

=head2 new

  Arg [1]    : Bio::EnsEMBL::DBAdaptor $dbadaptor the adaptor for
               the database this assembly mapper is using.
  Example    : my $asma = new Bio::EnsEMBL::AssemblyMapperAdaptor($dbadaptor);
  Description: Creates a new AssemblyMapperAdaptor object
  Returntype : Bio::EnsEMBL::DBSQL::AssemblyMapperAdaptor
  Exceptions : none
  Caller     : Bio::EnsEMBL::DBSQL::DBAdaptor
  Status     : Stable

=cut

sub new {
  my($class, $dbadaptor) = @_;

  my $self = $class->SUPER::new($dbadaptor);

  $self->{'_asm_mapper_cache'} = {};

  # use a shared cache (for this database) that contains info about
  # seq regions
  my $seq_region_cache = $self->db->get_SeqRegionCache();
  $self->{'sr_name_cache'} = $seq_region_cache->{'name_cache'};
  $self->{'sr_id_cache'}   = $seq_region_cache->{'id_cache'};

  return $self;
}



=head2  cache_seq_ids_with_mult_assemblys

  Example    : $self->adaptor->cache_seq_ids_with_mult_assemblys();
  Description: Creates a hash of the component seq region ids that
               map to more than one assembly from the assembly table.
  Retruntype : none
  Exceptions : none
  Caller     : AssemblyMapper, ChainedAssemblyMapper
  Status     : At Risk

=cut

sub cache_seq_ids_with_mult_assemblys{
  my $self = shift;
  my %multis;

  return if (defined($self->{'multi_seq_ids'}));

  $self->{'multi_seq_ids'} = {};

  my $sql = qq(
  SELECT    sra.seq_region_id
  FROM      seq_region_attrib sra,
            attrib_type at,
            seq_region sr,
            coord_system cs
  WHERE     sra.attrib_type_id = at.attrib_type_id
    AND     code = "MultAssem"
    AND     sra.seq_region_id = sr.seq_region_id
    AND     sr.coord_system_id = cs.coord_system_id
    AND     cs.species_id = ?);

  my $sth = $self->prepare($sql);

  $sth->bind_param( 1, $self->species_id(), SQL_INTEGER );

  $sth->execute();

  my ($seq_region_id);

  $sth->bind_columns(\$seq_region_id);

  while($sth->fetch()) {
    $self->{'multi_seq_ids'}->{$seq_region_id} = 1;
  }
  $sth->finish;
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
               $mapper = $asma->fetch_by_CoordSystems($cs1,$cs2);
               $mapper = $asma->fetch_by_CoordSystems($cs2,$cs1);
  Returntype : Bio::EnsEMBL::AssemblyMapper
  Exceptions : wrong argument types
  Caller     : general
  Status     : Stable

=cut

sub fetch_by_CoordSystems {
  my $self = shift;
  my $cs1  = shift;
  my $cs2  = shift;

  if(!ref($cs1) || !$cs1->isa('Bio::EnsEMBL::CoordSystem')) {
    throw("cs1 argument must be a Bio::EnsEMBL::CoordSystem.");
  }
  if(!ref($cs2) || !$cs2->isa('Bio::EnsEMBL::CoordSystem')) {
    throw("cs2 argument must be a Bio::EnsEMBL::CoordSystem.");
  }

#  if($cs1->equals($cs2)) {
#    throw("Cannot create mapper between same coord systems: " .
#          $cs1->name . " " . $cs1->version);
#  }

  if($cs1->is_top_level()) {
    return Bio::EnsEMBL::TopLevelAssemblyMapper->new($self, $cs1, $cs2);
  }
  if($cs2->is_top_level()) {
    return Bio::EnsEMBL::TopLevelAssemblyMapper->new($self, $cs2, $cs1);
  }

  my $csa = $self->db->get_CoordSystemAdaptor();

  #retrieve the shortest possible mapping path between these systems
  my @mapping_path = @{$csa->get_mapping_path($cs1,$cs2)};

  if(!@mapping_path) {

    # It is perfectly fine not to have a mapping. No warning needed really
    # Just check the return code!!

#    warning(
#      "There is no mapping defined between these coord systems:\n" .
#      $cs1->name() . " " . $cs1->version() . " and " . $cs2->name() . " " .
#      $cs2->version()
#    );
    return undef;
  }

  my $key = join(':', map({defined($_)?$_->dbID():"-"} @mapping_path));

  my $asm_mapper = $self->{'_asm_mapper_cache'}->{$key};

  return $asm_mapper if($asm_mapper);

  if(@mapping_path == 1) {
    throw("Incorrect mapping path defined in meta table. " .
	  "0 step mapping encountered between:\n" .
	  $cs1->name() . " " . $cs1->version() . " and " . $cs2->name . " " .
	  $cs2->version());
  }

  if(@mapping_path == 2) {
    #1 step regular mapping
    $asm_mapper = Bio::EnsEMBL::AssemblyMapper->new($self, @mapping_path);

#   If you want multiple pieces on two seqRegions to map to each other
#   you need to make an assembly.mapping entry that is seperated with a #
#   instead of an |.

    $self->{'_asm_mapper_cache'}->{$key} = $asm_mapper;
    return $asm_mapper;
  }

  if(@mapping_path == 3) {
   #two step chained mapping
    $asm_mapper = Bio::EnsEMBL::ChainedAssemblyMapper->new($self,@mapping_path);
    #in multi-step mapping it is possible get requests with the
    #coordinate system ordering reversed since both mappings directions
    #cache on both orderings just in case
    #e.g.   chr <-> contig <-> clone   and   clone <-> contig <-> chr

    $self->{'_asm_mapper_cache'}->{$key} = $asm_mapper;
    $key = join(':', map({defined($_)?$_->dbID():"-"} reverse(@mapping_path)));
    $self->{'_asm_mapper_cache'}->{$key} = $asm_mapper;
    return $asm_mapper;
  }

  throw("Only 1 and 2 step coordinate system mapping is currently\n" .
	"supported.  Mapping between " .
          $cs1->name() . " " . $cs1->version() . " and " .
          $cs2->name() . " " . $cs2->version() .
          " requires ". (scalar(@mapping_path)-1) . " steps.");
}



=head2 register_assembled

  Arg [1]    : Bio::EnsEMBL::AssemblyMapper $asm_mapper
               A valid AssemblyMapper object
  Arg [2]    : integer $asm_seq_region
               The dbID of the seq_region to be registered
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
  Status     : Stable

=cut


sub register_assembled {
  my $self = shift;
  my $asm_mapper = shift;
  my $asm_seq_region = shift;
  my $asm_start      = shift;
  my $asm_end        = shift;

  if(!ref($asm_mapper) || !$asm_mapper->isa('Bio::EnsEMBL::AssemblyMapper')) {
    throw("Bio::EnsEMBL::AssemblyMapper argument expected");
  }

  throw("asm_seq_region argument expected") if(!defined($asm_seq_region));
  throw("asm_start argument expected") if(!defined($asm_start));
  throw("asm_end argument expected") if(!defined($asm_end));

  my $asm_cs_id = $asm_mapper->assembled_CoordSystem->dbID();
  my $cmp_cs_id = $asm_mapper->component_CoordSystem->dbID();

  #split up the region to be registered into fixed chunks
  #this allows us to keep track of regions that have already been
  #registered and also works under the assumption that if a small region
  #is requested it is likely that other requests will be made in the
  #vicinity (the minimum size registered the chunksize (2^chunkfactor)

  my @chunk_regions;
  #determine span of chunks
  #bitwise shift right is fast and easy integer division

  my($start_chunk, $end_chunk);

  $start_chunk = $asm_start >> $CHUNKFACTOR;
  $end_chunk   = $asm_end   >> $CHUNKFACTOR;

  # inserts have start = end + 1, on boundary condition start_chunk
  # could be less than end chunk
  if($asm_start == $asm_end + 1) {
    ($start_chunk, $end_chunk) = ($end_chunk, $start_chunk);
  }

  #find regions of continuous unregistered chunks
  my $i;
  my ($begin_chunk_region,$end_chunk_region);
  for ($i = $start_chunk; $i <= $end_chunk; $i++) {
    if($asm_mapper->have_registered_assembled($asm_seq_region, $i)) {
      if(defined($begin_chunk_region)) {
        #this is the end of an unregistered region.
        my $region = [($begin_chunk_region   << $CHUNKFACTOR),
                      (($end_chunk_region+1)     << $CHUNKFACTOR)-1];
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
    my $region = [($begin_chunk_region << $CHUNKFACTOR),
                  (($end_chunk_region+1)   << $CHUNKFACTOR) -1];
    push @chunk_regions, $region;
  }

  return if(!@chunk_regions);

  # keep the Mapper to a reasonable size
  if( $asm_mapper->size() > $asm_mapper->max_pair_count() ) {
    $asm_mapper->flush();
    #we now have to go and register the entire requested region since we
    #just flushed everything

    @chunk_regions = ( [ ( $start_chunk << $CHUNKFACTOR)
                         , (($end_chunk+1) << $CHUNKFACTOR)-1 ] );

    for( my $i = $start_chunk; $i <= $end_chunk; $i++ ) {
      $asm_mapper->register_assembled( $asm_seq_region, $i );
    }
  }

#  my $asm_seq_region_id =
#    $self->_seq_region_name_to_id($asm_seq_region,$asm_cs_id);

  # Retrieve the description of how the assembled region is made from
  # component regions for each of the continuous blocks of unregistered,
  # chunked regions

  my $q = qq{
      SELECT
         asm.cmp_start,
         asm.cmp_end,
         asm.cmp_seq_region_id,
         sr.name,
         sr.length,
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
    $sth->bind_param(1,$asm_seq_region,SQL_INTEGER);
    $sth->bind_param(2,$region_start,SQL_INTEGER);
    $sth->bind_param(3,$region_end,SQL_INTEGER);
    $sth->bind_param(4,$cmp_cs_id,SQL_INTEGER);

    $sth->execute();

    my($cmp_start, $cmp_end, $cmp_seq_region_id, $cmp_seq_region, $ori,
      $asm_start, $asm_end, $cmp_seq_region_length);

    $sth->bind_columns(\$cmp_start, \$cmp_end, \$cmp_seq_region_id,
                       \$cmp_seq_region, \$cmp_seq_region_length, \$ori,
                       \$asm_start, \$asm_end);

    #
    # Load the unregistered regions of the mapper
    #
    while($sth->fetch()) {
      next if($asm_mapper->have_registered_component($cmp_seq_region_id)
               and !defined($self->{'multi_seq_ids'}->{$cmp_seq_region_id}));
      $asm_mapper->register_component($cmp_seq_region_id);
      $asm_mapper->mapper->add_map_coordinates(
                 $asm_seq_region, $asm_start, $asm_end,
                 $ori,
                 $cmp_seq_region_id, $cmp_start, $cmp_end);

      my $arr = [ $cmp_seq_region_id, $cmp_seq_region,
                  $cmp_cs_id, $cmp_seq_region_length ];

      $self->{'sr_name_cache'}->{"$cmp_seq_region:$cmp_cs_id"} = $arr;
      $self->{'sr_id_cache'}->{"$cmp_seq_region_id"} = $arr;
    }
  }

  $sth->finish();
}



sub _seq_region_name_to_id {
  my $self    = shift;
  my $sr_name = shift;
  my $cs_id   = shift;

  if(!defined($sr_name) or
     !defined($cs_id)){
      throw('seq_region_name and coord_system_id args are required');
  }

  my $arr = $self->{'sr_name_cache'}->{"$sr_name:$cs_id"};
  if( $arr ) {
    return $arr->[0];
  }

  # Get the seq_region_id via the name.  This would be quicker if we just
  # used internal ids instead but stored but then we lose the ability
  # the transform accross databases with different internal ids

  my $sth = $self->prepare("SELECT seq_region_id, length " .
			   "FROM   seq_region " .
			   "WHERE  name = ? AND coord_system_id = ?");

  $sth->bind_param(1,$sr_name,SQL_VARCHAR);
  $sth->bind_param(2,$cs_id,SQL_INTEGER);
  $sth->execute();

  my @row = $sth->fetchrow_array();
  unless ( @row ) {
    throw("No-existent seq_region [$sr_name] in coord system $cs_id");
  }
  my @more = $sth->fetchrow_array();
  if ( @more ) {
    throw("Ambiguous seq_region [$sr_name] in coord system $cs_id");
  }

  my ($sr_id, $sr_length) = @row;
  $sth->finish();

  $arr = [ $sr_id, $sr_name, $cs_id, $sr_length ];

  $self->{'sr_name_cache'}->{"$sr_name:$cs_id"} = $arr;
  $self->{'sr_id_cache'}->{"$sr_id"} = $arr;

  return $sr_id;
}

sub _seq_region_id_to_name {
  my $self    = shift;
  my $sr_id = shift;

  ($sr_id) || throw('seq_region_is required');

  my $arr = $self->{'sr_id_cache'}->{"$sr_id"};
  if( $arr ) {
    return $arr->[1];
  }

  # Get the seq_region name via the id.  This would be quicker if we just
  # used internal ids instead but stored but then we lose the ability
  # the transform accross databases with different internal ids

  my $sth = $self->prepare("SELECT name, length ,coord_system_id " .
			   "FROM   seq_region " .
			   "WHERE  seq_region_id = ? ");

  $sth->bind_param(1,$sr_id,SQL_INTEGER);
  $sth->execute();

  my @row = $sth->fetchrow_array();
  if(!$sth->rows() == 1) {
    throw("non-existant seq_region [$sr_id]");
  }

  my ($sr_name, $sr_length, $cs_id) = @row;
  $sth->finish();

  $arr = [ $sr_id, $sr_name, $cs_id, $sr_length ];

  $self->{'sr_name_cache'}->{"$sr_name:$cs_id"} = $arr;
  $self->{'sr_id_cache'}->{"$sr_id"} = $arr;

  return $sr_name;
}


=head2 register_component

  Arg [1]    : Bio::EnsEMBL::AssemblyMapper $asm_mapper
               A valid AssemblyMapper object
  Arg [2]    : integer $cmp_seq_region
               The dbID of the seq_region to be registered
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
  Status     : Stable

=cut

sub register_component {
  my $self = shift;
  my $asm_mapper = shift;
  my $cmp_seq_region = shift;

  if(!ref($asm_mapper) || !$asm_mapper->isa('Bio::EnsEMBL::AssemblyMapper')) {
    throw("Bio::EnsEMBL::AssemblyMapper argument expected");
  }

  if(!defined($cmp_seq_region)) {
    throw("cmp_seq_region argument expected");
  }

  my $cmp_cs_id = $asm_mapper->component_CoordSystem()->dbID();
  my $asm_cs_id = $asm_mapper->assembled_CoordSystem()->dbID();

  #do nothing if this region is already registered or special case
  return if($asm_mapper->have_registered_component($cmp_seq_region)
  and !defined($self->{'multi_seq_ids'}->{$cmp_seq_region}));

#  my $cmp_seq_region_id =
#    $self->_seq_region_name_to_id($cmp_seq_region, $cmp_cs_id);

  # Determine what part of the assembled region this component region makes up

  my $q = qq{
      SELECT
         asm.asm_start,
         asm.asm_end,
         asm.asm_seq_region_id,
         sr.name,
         sr.length
      FROM
         assembly asm, seq_region sr
      WHERE
         asm.cmp_seq_region_id = ? AND
         asm.asm_seq_region_id = sr.seq_region_id AND
         sr.coord_system_id = ?
   };

  my $sth = $self->prepare($q);
  $sth->bind_param(1,$cmp_seq_region,SQL_INTEGER);
  $sth->bind_param(2,$asm_cs_id,SQL_INTEGER);
  $sth->execute();

  my @rows = $sth->fetchrow_array();

  if($sth->rows() == 0) {
    #this component is not used in the assembled part i.e. gap
    $asm_mapper->register_component($cmp_seq_region);
    $sth->finish();
    return;
  }

  #we do not currently support components mapping to multiple assembled
  # make sure that you've got the correct mapping in the meta-table :
  #   chromosome:EquCab2#contig ( use'#' for multiple mappings )
  #   chromosome:EquCab2|contig ( use '|' delimiter for 1-1 mappings )
  #
  my @more = $sth->fetchrow_array();
  if($sth->rows() != 1) {
    $sth->finish();
    throw("Multiple assembled regions for single " .
          "component region cmp_seq_region_id=[$cmp_seq_region]\n".
          "Remember that multiple mappings use the \#-operaator".
          " in the meta-table (i.e. chromosome:EquCab2\#contig\n");
  }

  my ($asm_start, $asm_end, $asm_seq_region_id,
      $asm_seq_region, $asm_seq_region_length) = @rows;

  my $arr = [ $asm_seq_region_id, $asm_seq_region,
              $asm_cs_id, $asm_seq_region_length ];

  $self->{'sr_name_cache'}->{"$asm_seq_region:$asm_cs_id"} = $arr;
  $self->{'sr_id_cache'}->{"$asm_seq_region_id"} = $arr;

  $sth->finish();

  # Register the corresponding assembled region. This allows a us to
  # register things in assembled chunks which allows us to:
  # (1) Keep track of what assembled regions are registered
  # (2) Use locality of reference (if they want something in same general
  #     region it will already be registered).

  $self->register_assembled($asm_mapper,$asm_seq_region_id,$asm_start,$asm_end);
}



=head2 register_chained

  Arg [1]    : Bio::EnsEMBL::ChainedAssemblyMapper $casm_mapper
               The chained assembly mapper to register regions on
  Arg [2]    : string $from ('first' or 'last')
               The direction we are registering from, and the name of the
               internal mapper.
  Arg [3]    : string $seq_region_name
               The name of the seqregion we are registering on
  Arg [4]    : listref $ranges
               A list  of ranges to register (in [$start,$end] tuples).
  Arg [5]    : (optional) $to_slice
               Only register those on this Slice.
  Description: Registers a set of ranges on a chained assembly mapper.
               This function is at the heart of the chained mapping process.
               It retrieves information from the assembly table and
               dynamically constructs the mappings between two coordinate
               systems which are 2 mapping steps apart. It does this by using
               two internal mappers to load up a third mapper which is
               actually used by the ChainedAssemblyMapper to perform the
               mapping.

               This method must be called before any mapping is
               attempted on regions of interest, otherwise only gaps will
               be returned.  Note that the ChainedAssemblyMapper automatically
               calls this method when the need arises.
  Returntype : none
  Exceptions : throw if the seq_region to be registered does not exist
               or if it associated with multiple assembled pieces (bad data
               in assembly table)

               throw if the mapping between the coordinate systems cannot
               be performed in two steps, which means there is an internal
               error in the data in the meta table or in the code that creates
               the mapping paths.
  Caller     : Bio::EnsEMBL::AssemblyMapper
  Status     : Stable

=cut

sub register_chained {
  my $self = shift;
  my $casm_mapper = shift;
  my $from = shift;
  my $seq_region_id = shift;
  my $ranges = shift;
  my $to_slice = shift;

  my $to_seq_region_id;
  if(defined($to_slice)){
   if($casm_mapper->first_CoordSystem()->equals($casm_mapper->last_CoordSystem())){
      return $self->_register_chained_special($casm_mapper, $from, $seq_region_id, $ranges, $to_slice);
    }
    $to_seq_region_id = $to_slice->get_seq_region_id();
    if(!defined($to_seq_region_id)){
      die "Could not get seq_region_id for to_slice".$to_slice->seq_region_name."\n";
    }
  }

  my ($start_name, $start_mid_mapper, $start_cs, $start_registry,
      $end_name, $end_mid_mapper, $end_cs, $end_registry);

  if($from eq 'first') {
    $start_name       = 'first';
    $start_mid_mapper = $casm_mapper->first_middle_mapper();
    $start_cs         = $casm_mapper->first_CoordSystem();
    $start_registry   = $casm_mapper->first_registry();
    $end_mid_mapper   = $casm_mapper->last_middle_mapper();
    $end_cs           = $casm_mapper->last_CoordSystem();
    $end_registry     = $casm_mapper->last_registry();
    $end_name         = 'last';
  } elsif($from eq 'last') {
    $start_name       = 'last';
    $start_mid_mapper = $casm_mapper->last_middle_mapper();
    $start_cs         = $casm_mapper->last_CoordSystem();
    $start_registry   = $casm_mapper->last_registry();
    $end_mid_mapper   = $casm_mapper->first_middle_mapper();
    $end_cs           = $casm_mapper->first_CoordSystem();
    $end_registry     = $casm_mapper->first_registry();
    $end_name         = 'first';
  } else {
    throw("Invalid from argument: [$from], must be 'first' or 'last'");
  }

  my $combined_mapper = $casm_mapper->first_last_mapper();
  my $mid_cs          = $casm_mapper->middle_CoordSystem();
  my $mid_name        = 'middle';
  my $csa             = $self->db->get_CoordSystemAdaptor();

  # Check for the simple case where the ChainedMapper is short
  if( ! defined $mid_cs ) {
      $start_mid_mapper = $combined_mapper;
  }


  ##############
  # obtain the first half of the mappings and load them into the start mapper
  #

  #ascertain which is component and which is actually assembled coord system
  my @path;

  # check for the simple case, where the ChainedMapper is short
  if( defined $mid_cs ) {
      @path = @{$csa->get_mapping_path($start_cs, $mid_cs)};
  } else {
      @path = @{$csa->get_mapping_path( $start_cs, $end_cs )};
  }

  if(@path != 2 && defined( $path[1] )) {
    my $path = join(',', map({$_->name .' '. $_->version} @path));
    my $len  = scalar(@path) - 1;
    throw("Unexpected mapping path between start and intermediate " .
	  "coord systems (". $start_cs->name . " " . $start_cs->version .
	  " and " . $mid_cs->name . " " . $mid_cs->version . ")." .
	  "\nExpected path length 1, got $len. " .
	  "(path=$path)");
  }

  my $sth;
  my ($asm_cs,$cmp_cs);
  $asm_cs = $path[0];
  $cmp_cs = $path[-1];

  #the SQL varies depending on whether we are coming from assembled or
  #component coordinate system

my $asm2cmp = (<<ASMCMP);
   SELECT
         asm.cmp_start,
         asm.cmp_end,
         asm.cmp_seq_region_id,
         sr.name,
         sr.length,
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
ASMCMP


my $cmp2asm = (<<CMPASM);
   SELECT
         asm.asm_start,
         asm.asm_end,
         asm.asm_seq_region_id,
         sr.name,
         sr.length,
         asm.ori,
         asm.cmp_start,
         asm.cmp_end
      FROM
         assembly asm, seq_region sr
      WHERE
         asm.cmp_seq_region_id = ? AND
         ? <= asm.cmp_end AND
         ? >= asm.cmp_start AND
         asm.asm_seq_region_id = sr.seq_region_id AND
	 sr.coord_system_id = ?
CMPASM

  my $asm2cmp_sth;
  my $cmp2asm_sth;
  if(defined($to_slice)){
    my $to_cs = $to_slice->coord_system;
    if($asm_cs->equals($to_cs)){
      $asm2cmp_sth = $self->prepare($asm2cmp);
      $cmp2asm_sth = $self->prepare($cmp2asm." AND asm.asm_seq_region_id = $to_seq_region_id");
    }
    elsif($cmp_cs->equals($to_cs)){
      $asm2cmp_sth = $self->prepare($asm2cmp." AND asm.cmp_seq_region_id = $to_seq_region_id");
      $cmp2asm_sth = $self->prepare($cmp2asm);
    }
    else{
      $asm2cmp_sth = $self->prepare($asm2cmp);
      $cmp2asm_sth = $self->prepare($cmp2asm);
    }
  }	
  else{
    $asm2cmp_sth = $self->prepare($asm2cmp);
    $cmp2asm_sth = $self->prepare($cmp2asm);
  }



  $sth = ($asm_cs->equals($start_cs)) ? $asm2cmp_sth : $cmp2asm_sth;

  my $mid_cs_id;

  # check for the simple case where the ChainedMapper is short
  if( defined $mid_cs ) {
      $mid_cs_id = $mid_cs->dbID();
  } else {
      $mid_cs_id = $end_cs->dbID();
  }

  my @mid_ranges;
  my @start_ranges;

  #need to perform the query for each unregistered range
  foreach my $range (@$ranges) {
    my ($start, $end) = @$range;
    $sth->bind_param(1,$seq_region_id,SQL_INTEGER);
    $sth->bind_param(2,$start,SQL_INTEGER);
    $sth->bind_param(3,$end,SQL_INTEGER);
    $sth->bind_param(4,$mid_cs_id,SQL_INTEGER);
    $sth->execute();

    #load the start <-> mid mapper with the results and record the mid cs
    #ranges we just added to the mapper

    my ($mid_start, $mid_end, $mid_seq_region_id, $mid_seq_region, $mid_length,
        $ori, $start_start, $start_end);

    $sth->bind_columns(\$mid_start, \$mid_end, \$mid_seq_region_id,
		       \$mid_seq_region, \$mid_length, \$ori, \$start_start,
		       \$start_end);

    while($sth->fetch()) {

      if( defined $mid_cs ) {
        $start_mid_mapper->add_map_coordinates
          (
           $seq_region_id,$start_start, $start_end, $ori,
           $mid_seq_region_id, $mid_start, $mid_end
          );
      } else {
        if( $from eq "first" ) {
          $combined_mapper->add_map_coordinates
            (
             $seq_region_id,$start_start, $start_end, $ori,
             $mid_seq_region_id, $mid_start, $mid_end
            );
        } else {
          $combined_mapper->add_map_coordinates
            (
             $mid_seq_region_id, $mid_start, $mid_end, $ori,
             $seq_region_id,$start_start, $start_end
            );
        }
      }

      #update sr_name cache
      my $arr = [ $mid_seq_region_id, $mid_seq_region,
                  $mid_cs_id, $mid_length ];

      $self->{'sr_name_cache'}->{"$mid_seq_region:$mid_cs_id"} = $arr;
      $self->{'sr_id_cache'}->{"$mid_seq_region_id"} = $arr;

      push @mid_ranges,[$mid_seq_region_id,$mid_seq_region,
                        $mid_start,$mid_end];
      push @start_ranges, [ $seq_region_id, $start_start, $start_end ];

      #the region that we actually register may actually be larger or smaller
      #than the region that we wanted to register.
      #register the intersection of the region so we do not end up doing
      #extra work later

      if($start_start < $start || $start_end > $end) {
        $start_registry->check_and_register($seq_region_id,$start_start,
                                            $start_end);
      }
    }
    $sth->finish();
  }

  # in the one step case, we load the mid ranges in the
  # last_registry and we are done
  if( ! defined $mid_cs ) {
    for my $range ( @mid_ranges ) {
      $end_registry->check_and_register( $range->[0], $range->[2],
                                         $range->[3] );
    }

    # and thats it for the simple case ...
    return;
  }


  ###########
  # now the second half of the mapping
  # perform another query and load the mid <-> end mapper using the mid cs
  # ranges
  #

  #ascertain which is component and which is actually assembled coord system
  @path = @{$csa->get_mapping_path($mid_cs, $end_cs)};
  if(@path == 2 || ( @path == 3 && !defined $path[1])) {

    $asm_cs = $path[0];
    $cmp_cs = $path[-1];
  } else {
    my $path = join(',', map({$_->name .' '. $_->version} @path));
    my $len = scalar(@path)-1;
    throw("Unexpected mapping path between intermediate and last" .
	  "coord systems (". $mid_cs->name . " " . $mid_cs->version .
	  " and " . $end_cs->name . " " . $end_cs->version . ")." .
	  "\nExpected path length 1, got $len. " .
	  "(path=$path)");
  }

  if(defined($to_slice)){
    my $to_cs = $to_slice->coord_system;
    if($asm_cs->equals($to_cs)){
      $asm2cmp_sth = $self->prepare($asm2cmp);
      $cmp2asm_sth = $self->prepare($cmp2asm." AND asm.asm_seq_region_id = $to_seq_region_id");
    }
    elsif($cmp_cs->equals($to_cs)){
      $asm2cmp_sth = $self->prepare($asm2cmp." AND asm.cmp_seq_region_id = $to_seq_region_id");
      $cmp2asm_sth = $self->prepare($cmp2asm);
    }
    else{
      $asm2cmp_sth = $self->prepare($asm2cmp);
      $cmp2asm_sth = $self->prepare($cmp2asm);
    }
  }

  $sth = ($asm_cs->equals($mid_cs)) ? $asm2cmp_sth : $cmp2asm_sth;

  my $end_cs_id = $end_cs->dbID();
  foreach my $mid_range (@mid_ranges) {
    my ($mid_seq_region_id, $mid_seq_region,$start, $end) = @$mid_range;
    $sth->bind_param(1,$mid_seq_region_id,SQL_INTEGER);
    $sth->bind_param(2,$start,SQL_INTEGER);
    $sth->bind_param(3,$end,SQL_INTEGER);
    $sth->bind_param(4,$end_cs_id,SQL_INTEGER);
    $sth->execute();

    #load the end <-> mid mapper with the results and record the mid cs
    #ranges we just added to the mapper

    my ($end_start, $end_end, $end_seq_region_id, $end_seq_region, $end_length,
        $ori, $mid_start, $mid_end);

    $sth->bind_columns(\$end_start, \$end_end, \$end_seq_region_id,
		       \$end_seq_region, \$end_length, \$ori, \$mid_start,
		       \$mid_end);

    while($sth->fetch()) {
      $end_mid_mapper->add_map_coordinates
        (
         $end_seq_region_id, $end_start, $end_end, $ori,
         $mid_seq_region_id, $mid_start, $mid_end
        );

      #update sr_name cache
      my $arr = [ $end_seq_region_id,$end_seq_region,$end_cs_id,$end_length ];

      $self->{'sr_name_cache'}->{"$end_seq_region:$end_cs_id"} = $arr;
      $self->{'sr_id_cache'}->{"$end_seq_region_id"} = $arr;

      #register this region on the end coord system
      $end_registry->check_and_register($end_seq_region_id, $end_start, $end_end);
    }
    $sth->finish();
  }

  #########
  # Now that both halves are loaded
  # Do stepwise mapping using both of the loaded mappers to load
  # the final start <-> end mapper
  #

  _build_combined_mapper(\@start_ranges, $start_mid_mapper, $end_mid_mapper,
                         $combined_mapper, $start_name);
  #all done!
  return;
}


=head2 _register_chained_special

  Arg [1]    : Bio::EnsEMBL::ChainedAssemblyMapper $casm_mapper
               The chained assembly mapper to register regions on
  Arg [2]    : string $from ('first' or 'last')
               The direction we are registering from, and the name of the
               internal mapper.
  Arg [3]    : string $seq_region_name
               The name of the seqregion we are registering on
  Arg [4]    : listref $ranges
               A list  of ranges to register (in [$start,$end] tuples).
  Arg [5]    : (optional) $to_slice
               Only register those on this Slice.
  Description: Registers a set of ranges on a chained assembly mapper.
               This function is at the heart of the chained mapping process.
               It retrieves information from the assembly table and
               dynamically constructs the mappings between two coordinate
               systems which are 2 mapping steps apart. It does this by using
               two internal mappers to load up a third mapper which is
               actually used by the ChainedAssemblyMapper to perform the
               mapping.

               This method must be called before any mapping is
               attempted on regions of interest, otherwise only gaps will
               be returned.  Note that the ChainedAssemblyMapper automatically
               calls this method when the need arises.
  Returntype : none
  Exceptions : throw if the seq_region to be registered does not exist
               or if it associated with multiple assembled pieces (bad data
               in assembly table)

               throw if the mapping between the coordinate systems cannot
               be performed in two steps, which means there is an internal
               error in the data in the meta table or in the code that creates
               the mapping paths.
  Caller     : Bio::EnsEMBL::AssemblyMapper
  Status     : Stable

=cut

sub _register_chained_special {
  my $self = shift;
  my $casm_mapper = shift;
  my $from = shift;
  my $seq_region_id = shift;
  my $ranges = shift;
  my $to_slice = shift;
  my $found = 0;

  my $sth = $self->prepare("SELECT
		      asm.cmp_start,
		      asm.cmp_end,
		      asm.cmp_seq_region_id,
		      sr.name,
		      sr.length,
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
		      sr.coord_system_id = ? AND
		      asm.cmp_seq_region_id = ?");


  my ($start_name, $start_mid_mapper, $start_cs, $start_registry,
      $end_name, $end_mid_mapper, $end_cs, $end_registry);

  if($from eq 'first') {
    $start_name       = 'first';
    $start_mid_mapper = $casm_mapper->first_middle_mapper();
    $start_cs         = $casm_mapper->first_CoordSystem();
    $start_registry   = $casm_mapper->first_registry();
    $end_mid_mapper   = $casm_mapper->last_middle_mapper();
    $end_cs           = $casm_mapper->last_CoordSystem();
    $end_registry     = $casm_mapper->last_registry();
    $end_name         = 'last';
  } elsif($from eq 'last') {
    $start_name       = 'last';
    $start_mid_mapper = $casm_mapper->last_middle_mapper();
    $start_cs         = $casm_mapper->last_CoordSystem();
    $start_registry   = $casm_mapper->last_registry();
    $end_mid_mapper   = $casm_mapper->first_middle_mapper();
    $end_cs           = $casm_mapper->first_CoordSystem();
    $end_registry     = $casm_mapper->first_registry();
    $end_name         = 'first';
  } else {
    throw("Invalid from argument: [$from], must be 'first' or 'last'");
  }

  my $combined_mapper = $casm_mapper->first_last_mapper();
  my $mid_cs          = $casm_mapper->middle_CoordSystem();
  my $mid_name        = 'middle';
  my $csa             = $self->db->get_CoordSystemAdaptor();

  # Check for the simple case where the ChainedMapper is short
  if( ! defined $mid_cs ) {
      $start_mid_mapper = $combined_mapper;
  }


  my @path;
  if( defined $mid_cs ) {
      @path = @{$csa->get_mapping_path($start_cs, $mid_cs)};
  } else {
      @path = @{$csa->get_mapping_path( $start_cs, $end_cs )};
  }
  if( ! defined $mid_cs ) {
      $start_mid_mapper = $combined_mapper;
  }

  if(@path != 2 && defined( $path[1] )) {
    my $path = join(',', map({$_->name .' '. $_->version} @path));
    my $len  = scalar(@path) - 1;
    throw("Unexpected mapping path between start and intermediate " .
	  "coord systems (". $start_cs->name . " " . $start_cs->version .
	  " and " . $mid_cs->name . " " . $mid_cs->version . ")." .
	  "\nExpected path length 1, got $len. " .
	  "(path=$path)");
  }

  my ($asm_cs,$cmp_cs);
  $asm_cs = $path[0];
  $cmp_cs = $path[-1];

  $combined_mapper = $casm_mapper->first_last_mapper();
  $mid_cs          = $casm_mapper->middle_CoordSystem();
  $mid_name        = 'middle';
  $csa             = $self->db->get_CoordSystemAdaptor();

  my $mid_cs_id;

  # Check for the simple case where the ChainedMapper is short
  if ( !defined $mid_cs ) {
    $start_mid_mapper = $combined_mapper;
  } else {
    $mid_cs_id = $mid_cs->dbID();
  }

  my @mid_ranges;
  my @start_ranges;

  my $to_cs = $to_slice->coord_system;
  foreach my $direction (1, 0){
    my $id1;
    my $id2;
    if($direction){
      $id1 = $seq_region_id;
      $id2 = $to_slice->get_seq_region_id();
    }
    else{
      $id2 = $seq_region_id;
      $id1 = $to_slice->get_seq_region_id();
    }

    foreach my $range (@$ranges) {
      my ($start, $end) = @$range;
      $sth->bind_param(1,$id1,SQL_INTEGER);
      $sth->bind_param(2,$start,SQL_INTEGER);
      $sth->bind_param(3,$end,SQL_INTEGER);
      $sth->bind_param(4,$to_cs->dbID,SQL_INTEGER);
      $sth->bind_param(5,$id2,SQL_INTEGER);
      $sth->execute();

      my ($mid_start, $mid_end, $mid_seq_region_id, $mid_seq_region, $mid_length,
	  $ori, $start_start, $start_end);

      $sth->bind_columns(\$mid_start, \$mid_end, \$mid_seq_region_id,
			 \$mid_seq_region, \$mid_length, \$ori, \$start_start,
			 \$start_end);

      while($sth->fetch()) {
	$found = 1;

	if( defined $mid_cs ) {
	  $start_mid_mapper->add_map_coordinates
	    (
	     $id1,$start_start, $start_end, $ori,
	     $mid_seq_region_id, $mid_start, $mid_end
	    );
	} else {
	  if( $from eq "first") {
	    if($direction){
	      $combined_mapper->add_map_coordinates
		(
		 $id1,$start_start, $start_end, $ori,
		 $mid_seq_region_id, $mid_start, $mid_end
		);
	    }
	    else{
	      $combined_mapper->add_map_coordinates
		(
		 $mid_seq_region_id, $mid_start, $mid_end, $ori,
		 $id1,$start_start, $start_end
		);
	    }	
	  } else {
	    if($direction){
	      $combined_mapper->add_map_coordinates
		(
		 $mid_seq_region_id, $mid_start, $mid_end, $ori,
		 $id1,$start_start, $start_end
		);
	    }
	    else{
	      $combined_mapper->add_map_coordinates
		(
		 $id1,$start_start, $start_end, $ori,
		 $mid_seq_region_id, $mid_start, $mid_end
		);
	    }
	  }
	}
	
	#update sr_name cache
	my $arr = [ $mid_seq_region_id, $mid_seq_region,
		    $mid_cs_id, $mid_length ];
	
	$self->{'sr_name_cache'}->{"$mid_seq_region:$mid_cs_id"} = $arr;
	$self->{'sr_id_cache'}->{"$mid_seq_region_id"} = $arr;
	
	push @mid_ranges,[$mid_seq_region_id,$mid_seq_region,
			  $mid_start,$mid_end];
	push @start_ranges, [ $id1, $start_start, $start_end ];
	
	#the region that we actually register may actually be larger or smaller
	#than the region that we wanted to register.
	#register the intersection of the region so we do not end up doing
	#extra work later
	
	if($start_start < $start || $start_end > $end) {
	  $start_registry->check_and_register($id1,$start_start,
					      $start_end);
	}
      }
    $sth->finish();
    }
    if($found){
      if( ! defined $mid_cs ) {
	for my $range ( @mid_ranges ) {
	  $end_registry->check_and_register( $range->[0], $range->[2],
					     $range->[3] );
	}
	
	# and thats it for the simple case ...
	return;
      }
    }
  }
}


=head2 register_all

  Arg [1]    : Bio::EnsEMBL::AssemblyMapper $mapper
  Example    : $mapper = $asm_mapper_adaptor->fetch_by_CoordSystems($cs1,$cs2);

               # make cache large enough to hold all of the mappings
               $mapper->max_pair_count(10e6);
               $asm_mapper_adaptor->register_all($mapper);

               # perform mappings as normal
               $mapper->map($slice->seq_region_name(), $sr_start, $sr_end,
                            $sr_strand, $cs1);
               ...
  Description: This function registers the entire set of mappings between
               two coordinate systems in an assembly mapper.
               This will use a lot of memory but will be much more efficient
               when doing a lot of mapping which is spread over the entire
               genome.
  Returntype : none
  Exceptions : none
  Caller     : specialised prograhsm
  Status     : Stable

=cut

sub register_all {
  my $self = shift;
  my $mapper = shift;

  my $asm_cs_id = $mapper->assembled_CoordSystem()->dbID();
  my $cmp_cs_id = $mapper->component_CoordSystem()->dbID();

  # retrieve every relevant assembled/component pair from the assembly table

  my $q = qq{
      SELECT
         asm.cmp_start,
         asm.cmp_end,
         asm.cmp_seq_region_id,
         cmp_sr.name,
         cmp_sr.length,
         asm.ori,
         asm.asm_start,
         asm.asm_end,
         asm.asm_seq_region_id,
         asm_sr.name,
         asm_sr.length
      FROM
         assembly asm, seq_region asm_sr, seq_region cmp_sr
      WHERE
         asm.cmp_seq_region_id = cmp_sr.seq_region_id AND
         asm.asm_seq_region_id = asm_sr.seq_region_id AND
         cmp_sr.coord_system_id = ? AND
         asm_sr.coord_system_id = ?
   };

  my $sth = $self->prepare($q);

  $sth->bind_param(1,$cmp_cs_id,SQL_INTEGER);
  $sth->bind_param(2,$asm_cs_id,SQL_INTEGER);
  $sth->execute();

  # load the mapper with the assembly information

  my ($cmp_start, $cmp_end, $cmp_seq_region_id, $cmp_seq_region, $cmp_length,
      $ori,
      $asm_start, $asm_end, $asm_seq_region_id, $asm_seq_region, $asm_length);

  $sth->bind_columns(\$cmp_start, \$cmp_end, \$cmp_seq_region_id,
                     \$cmp_seq_region, \$cmp_length, \$ori,
                     \$asm_start, \$asm_end, \$asm_seq_region_id,
                     \$asm_seq_region, \$asm_length);

  my %asm_registered;

  while($sth->fetch()) {
    $mapper->register_component($cmp_seq_region_id);
    $mapper->mapper->add_map_coordinates(
                 $asm_seq_region_id, $asm_start, $asm_end,
                 $ori,
                 $cmp_seq_region_id, $cmp_start, $cmp_end);

      my $arr = [$cmp_seq_region_id, $cmp_seq_region, $cmp_cs_id, $cmp_length];

      $self->{'sr_name_cache'}->{"$cmp_seq_region:$cmp_cs_id"} = $arr;
      $self->{'sr_id_cache'}->{"$cmp_seq_region_id"} = $arr;

    # only register each asm seq_region once since it requires some work
    if(!$asm_registered{$asm_seq_region_id}) {
      $asm_registered{$asm_seq_region_id} = 1;

      # register all chunks from start of seq region to end
      my $end_chunk = $asm_length >> $CHUNKFACTOR;
      for(my $i = 0; $i <= $end_chunk; $i++) {
        $mapper->register_assembled($asm_seq_region_id, $i);
      }

      $arr = [$asm_seq_region_id, $asm_seq_region, $asm_cs_id, $asm_length];

      $self->{'sr_name_cache'}->{"$asm_seq_region:$asm_cs_id"} = $arr;
      $self->{'sr_id_cache'}->{"$asm_seq_region_id"} = $arr;
    }
  }

  $sth->finish();

  return;
}



=head2 register_all_chained

  Arg [1]    : Bio::EnsEMBL::ChainedAssemblyMapper $casm_mapper
  Example    : $mapper = $asm_mapper_adaptor->fetch_by_CoordSystems($cs1,$cs2);

               # make the cache large enough to hold all of the mappings
               $mapper->max_pair_count(10e6);
               # load all of the mapping data
               $asm_mapper_adaptor->register_all_chained($mapper);

               # perform mappings as normal
               $mapper->map($slice->seq_region_name(), $sr_start, $sr_end,
                            $sr_strand, $cs1);
               ...
  Description: This function registers the entire set of mappings between
               two coordinate systems in a chained mapper.  This will use a lot
               of memory but will be much more efficient when doing a lot of
               mapping which is spread over the entire genome.
  Returntype : none
  Exceptions : throw if mapper is between coord systems with unexpected
               mapping paths
  Caller     : specialised programs doing a lot of genome-wide mapping
  Status     : Stable

=cut

sub register_all_chained {
  my $self = shift;
  my $casm_mapper = shift;

  my $first_cs = $casm_mapper->first_CoordSystem();
  my $mid_cs   = $casm_mapper->middle_CoordSystem();
  my $last_cs  = $casm_mapper->last_CoordSystem();

  my $start_mid_mapper = $casm_mapper->first_middle_mapper();
  my $end_mid_mapper   = $casm_mapper->last_middle_mapper();
  my $combined_mapper  = $casm_mapper->first_last_mapper();

  my @ranges;

  my $sth = $self->prepare(
     'SELECT
         asm.cmp_start,
         asm.cmp_end,
         asm.cmp_seq_region_id,
         sr_cmp.name,
         sr_cmp.length,
         asm.ori,
         asm.asm_start,
         asm.asm_end,
         asm.asm_seq_region_id,
         sr_asm.name,
         sr_asm.length
      FROM
         assembly asm, seq_region sr_asm, seq_region sr_cmp
      WHERE
         sr_asm.seq_region_id = asm.asm_seq_region_id AND
         sr_cmp.seq_region_id = asm.cmp_seq_region_id AND
         sr_asm.coord_system_id = ? AND
         sr_cmp.coord_system_id = ?');

  my $csa = $self->db()->get_CoordSystemAdaptor();

  my @path;

  if ( !defined $mid_cs ) {
    @path = @{ $csa->get_mapping_path( $first_cs, $last_cs ) };
    if ( !defined( $path[1] ) ) {
      splice( @path, 1, 1 );
    }
  } else {
    @path = @{ $csa->get_mapping_path( $first_cs, $mid_cs ) };
    # fix for when we have something like supercontig#contig#chromosome
    if ( !defined( $path[1] ) ) {
      splice( @path, 1, 1 );
    }
  }

  if ( @path != 2 ) {
    my $path =
      join( ',', map( { $_->name . ' ' . $_->version } @path ) );
    my $len = scalar(@path) - 1;
    throw(   "Unexpected mapping path between start and intermediate "
           . "coord systems ("
           . $first_cs->name . " "
           . $first_cs->version . " and "
           . $mid_cs->name . " "
           . $mid_cs->version . ")."
           . "\nExpected path length 1, got $len. "
           . "(path=$path)" );
  }

  my ($asm_cs,$cmp_cs) = @path;

  $sth->{mysql_use_result} = 1;
  $sth->bind_param(1,$asm_cs->dbID,SQL_INTEGER);
  $sth->bind_param(2,$cmp_cs->dbID,SQL_INTEGER);
  $sth->execute();


  my ($mid_start, $mid_end, $mid_seq_region_id, $mid_seq_region, $mid_length,
      $ori, $start_start, $start_end, $start_seq_region_id, $start_seq_region,
      $start_length);

  if($asm_cs->equals($first_cs)) {
    $sth->bind_columns(\$mid_start, \$mid_end, \$mid_seq_region_id,
                       \$mid_seq_region, \$mid_length, \$ori, \$start_start,
                       \$start_end, \$start_seq_region_id, \$start_seq_region,
                       \$start_length);
  } else {
    $sth->bind_columns(\$start_start, \$start_end, \$start_seq_region_id,
                       \$start_seq_region, \$start_length, \$ori,
                       \$mid_start, \$mid_end, \$mid_seq_region_id,
                       \$mid_seq_region, \$mid_length);

  }

  my ( $mid_cs_id, $start_cs_id, $registry, $mapper );
  if ( !defined $mid_cs ) {
    $mid_cs_id   = $last_cs->dbID();
    $start_cs_id = $first_cs->dbID();
    $mapper      = $combined_mapper;
  } else {
    $mid_cs_id   = $mid_cs->dbID();
    $start_cs_id = $first_cs->dbID();
    $mapper      = $start_mid_mapper;
  }

  $registry =  $casm_mapper->first_registry();

  while($sth->fetch()) {
    $mapper->add_map_coordinates
      (
       $start_seq_region_id, $start_start, $start_end, $ori,
       $mid_seq_region_id, $mid_start, $mid_end
      );
    push( @ranges, [$start_seq_region_id, $start_start, $start_end ] );

    $registry->check_and_register( $start_seq_region_id, 1, $start_length );
    if( ! defined $mid_cs ) {
      $casm_mapper->last_registry()->check_and_register
	( $mid_seq_region_id, $mid_start, $mid_end );
    }

    my $arr = [ $mid_seq_region_id, $mid_seq_region,
                $mid_cs_id, $mid_length ];

    $self->{'sr_name_cache'}->{"$mid_seq_region:$mid_cs_id"} = $arr;
    $self->{'sr_id_cache'}->{"$mid_seq_region_id"} = $arr;

    $arr = [ $start_seq_region_id, $start_seq_region,
             $start_cs_id, $start_length ];

    $self->{'sr_name_cache'}->{"$start_seq_region:$start_cs_id"} = $arr;
    $self->{'sr_id_cache'}->{"$start_seq_region_id"} = $arr;
  }

  if( ! defined $mid_cs ) {
    # thats it for the simple case
    return;
  }


  @path = @{ $csa->get_mapping_path( $last_cs, $mid_cs ) };
  if ( defined($mid_cs) ) {
    if ( !defined( $path[1] ) ) {
      splice( @path, 1, 1 );
    }
  }

  if ( @path != 2 ) {
    my $path =
      join( ',', map( { $_->name . ' ' . $_->version } @path ) );
    my $len = scalar(@path) - 1;
    throw(   "Unexpected mapping path between intermediate and last "
           . "coord systems ("
           . $last_cs->name . " "
           . $last_cs->version . " and "
           . $mid_cs->name . " "
           . $mid_cs->version . ")."
           . "\nExpected path length 1, got $len. "
           . "(path=$path)" );
  }

  ($asm_cs,$cmp_cs) = @path;

  $sth->bind_param(1,$asm_cs->dbID,SQL_INTEGER);
  $sth->bind_param(2,$cmp_cs->dbID,SQL_INTEGER);
  $sth->execute();


  my ($end_start, $end_end, $end_seq_region_id, $end_seq_region,
      $end_length);

  if($asm_cs->equals($mid_cs)) {
    $sth->bind_columns(\$end_start, \$end_end, \$end_seq_region_id,
                       \$end_seq_region, \$end_length, \$ori,
                       \$mid_start, \$mid_end, \$mid_seq_region_id,
                       \$mid_seq_region, \$mid_length);
  } else {
    $sth->bind_columns(\$mid_start, \$mid_end, \$mid_seq_region_id,
                       \$mid_seq_region, \$mid_length, \$ori, \$end_start,
                       \$end_end, \$end_seq_region_id, \$end_seq_region,
                       \$end_length);
  }

  my $end_cs_id = $last_cs->dbID();
  $registry = $casm_mapper->last_registry();

  while($sth->fetch()) {
    $end_mid_mapper->add_map_coordinates
      (
       $end_seq_region_id, $end_start, $end_end, $ori,
       $mid_seq_region_id, $mid_start, $mid_end
      );

    $registry->check_and_register( $end_seq_region_id, 1, $end_length );

    my $arr = [ $end_seq_region_id, $end_seq_region,
             $end_cs_id, $end_length ];
    $self->{'sr_name_cache'}->{"$end_seq_region:$end_cs_id"} = $arr;
    $self->{'sr_id_cache'}->{"$end_seq_region_id"} = $arr;
  }

  _build_combined_mapper( \@ranges, $start_mid_mapper, $end_mid_mapper,
                          $combined_mapper, "first" );

  return;
}



# after both halves of a chained mapper are loaded
# this function maps all ranges in $ranges and loads the
# results into the combined mapper
sub _build_combined_mapper {
  my $ranges = shift;
  my $start_mid_mapper = shift;
  my $end_mid_mapper = shift;
  my $combined_mapper = shift;
  my $start_name = shift;

  my $mid_name = "middle";

  foreach my $range (@$ranges) {
    my ( $seq_region_id, $start, $end) = @$range;

    my $sum = 0;

    my @initial_coords = $start_mid_mapper->map_coordinates($seq_region_id,
                                                            $start,$end,1,
                                                            $start_name);

    foreach my $icoord (@initial_coords) {
      #skip gaps
      if($icoord->isa('Bio::EnsEMBL::Mapper::Gap')) {
        $sum += $icoord->length();
        next;
      }


      #feed the results of the first mapping into the second mapper
      my @final_coords =
        $end_mid_mapper->map_coordinates($icoord->id, $icoord->start,
                                         $icoord->end,
                                         $icoord->strand, $mid_name);


      foreach my $fcoord (@final_coords) {
        #load up the final mapper

        if($fcoord->isa('Bio::EnsEMBL::Mapper::Coordinate')) {
          my $total_start = $start + $sum;
          my $total_end   = $total_start + $fcoord->length - 1;
          my $ori = $fcoord->strand();

          if($start_name eq 'first') { # add coords in consistant order
            $combined_mapper->add_map_coordinates(
                             $seq_region_id, $total_start, $total_end, $ori,
                             $fcoord->id(), $fcoord->start(), $fcoord->end());
          } else {
            $combined_mapper->add_map_coordinates(
                        $fcoord->id(), $fcoord->start(), $fcoord->end(),$ori,
                        $seq_region_id, $total_start, $total_end);
          }

        }
        $sum += $fcoord->length();
      }
    }
  }
  #all done!
}


=head2 seq_regions_to_ids

  Arg [1]    : Bio::EnsEMBL::CoordSystem $coord_system
  Arg [2]    : listref of strings $seq_regions
  Example    : my @ids = @{$asma->seq_regions_to_ids($coord_sys, \@seq_regs)};
  Description: Converts a list of seq_region names to internal identifiers
               using the internal cache that has accumulated while registering
               regions for AssemblyMappers. If any requested regions are
               not  found in the cache an attempt is made to retrieve them
               from the database.
  Returntype : listref of ints
  Exceptions : throw if a non-existant seqregion is provided
  Caller     : general
  Status     : Stable

=cut

sub seq_regions_to_ids {
  my $self = shift;
  my $coord_system = shift;
  my $seq_regions = shift;

  my $cs_id = $coord_system->dbID();

  my @out;

  foreach my $sr (@$seq_regions) {
    my $arr = $self->{'sr_name_cache'}->{"$sr:$cs_id"};
    if( $arr ) {
      push( @out, $arr->[0] );
    } else {
      push @out, $self->_seq_region_name_to_id($sr,$cs_id);
    }
  }

  return \@out;
}


=head2 seq_ids_to_regions

  Arg [1]    : listref of   seq_region ids
  Example    : my @ids = @{$asma->ids_to_seq_regions(\@seq_ids)};
  Description: Converts a list of seq_region ids to seq region names
               using the internal cache that has accumulated while registering
               regions for AssemblyMappers. If any requested regions are
               not  found in the cache an attempt is made to retrieve them
               from the database.
  Returntype : listref of strings
  Exceptions : throw if a non-existant seq_region_id is provided
  Caller     : general
  Status     : Stable

=cut

sub seq_ids_to_regions {
  my $self = shift;
  my $seq_region_ids = shift;

  my @out;

  foreach my $sr (@$seq_region_ids) {
    my $arr = $self->{'sr_id_cache'}->{"$sr"};
    if( $arr ) {
      push( @out, $arr->[1] );
    } else {
      push @out, $self->_seq_region_id_to_name($sr);
    }
  }

  return \@out;
}

=head2 delete_cache

 Description: Delete all the caches for the mappings/seq_regions
 Returntype : none
 Exceptions : none
 Caller     : General
 Status     : At risk

=cut

sub delete_cache{
  my ($self) = @_;

  %{$self->{'sr_name_cache'}}     = ();
  %{$self->{'sr_id_cache'}}       = ();

  foreach my $key (keys %{$self->{'_asm_mapper_cache'}}){
    $self->{'_asm_mapper_cache'}->{$key}->flush();
  }
  %{$self->{'_asm_mapper_cache'}} = ();
  return;
}


1;
