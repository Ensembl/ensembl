
#
# Ensembl module for Bio::EnsEMBL::DBSQL::AssemblyMapperAdaptor
#
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
makeup is described by the assembly table.  Currently one step (explicit) and
two step (implicit) pairwise mapping is supported.  In one-step mapping
an explicit relationship between the coordinate systems is defined in the
assembly table.  In two-step 'chained' mapping no explicit mapping is present 
but the coordinate systems must share a common mapping to an intermediate 
coordinate system.

=head1 CONTACT

Post general queries to B<ensembl-dev@ebi.ac.uk>

=head1 METHODS

=cut


package Bio::EnsEMBL::DBSQL::AssemblyMapperAdaptor;
use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::AssemblyMapper;
use Bio::EnsEMBL::ChainedAssemblyMapper;

use Bio::EnsEMBL::Utils::Cache; #CPAN LRU cache
use Bio::EnsEMBL::Utils::Exception qw(deprecate throw);

use integer; #do proper arithmetic bitshifts

@ISA = qw(Bio::EnsEMBL::DBSQL::BaseAdaptor);


my $CHUNKFACTOR = 20;  # 2^20 = approx. 10^6
# if the mapper is bigger than that its flushed before registering new stuff:
my $MAX_PAIR_COUNT = 1000; 

#number of seq regions to remember ids fo
my $SEQ_REGION_CACHE_SIZE = 2500;

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

  my %cache;
  tie(%cache, 'Bio::EnsEMBL::Utils::Cache', $SEQ_REGION_CACHE_SIZE);
  $self->{'_sr_id_cache'} = \%cache;

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

  if(!ref($cs1) || !$cs1->isa('Bio::EnsEMBL::CoordSystem')) {
    throw("cs1 argument must be a Bio::EnsEMBL::CoordSystem.");
  }
  if(!ref($cs2) || !$cs2->isa('Bio::EnsEMBL::CoordSystem')) {
    throw("cs2 argument must be a Bio::EnsEMBL::CoordSystem.");
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
          $cs1->name() . " " . $cs1->version() . " and " . $cs2->name() . " " .
          $cs2->version());
  }

  my $key = join(':', map({$_->dbID()} @mapping_path));

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
#    $asm_mapper = Bio::EnsEMBL::AssemblyMapper->new($self, @mapping_path);
    $asm_mapper = Bio::EnsEMBL::ChainedAssemblyMapper->new( $self, $mapping_path[0], undef, $mapping_path[1] );  
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
    $key = join(':', map({$_->dbID()} reverse(@mapping_path)));
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
  #determine span of chunks
  #bitwise shift right is fast and easy integer division
 

  my($start_chunk, $end_chunk);

  $start_chunk = $asm_start >> $CHUNKFACTOR;
  $end_chunk   = $asm_end   >> $CHUNKFACTOR;

  #find regions of continuous unregistered chunks
  my $i;
  my ($begin_chunk_region,$end_chunk_region);
  for ($i = $start_chunk; $i <= $end_chunk; $i++) {
    if($asm_mapper->have_registered_assembled($asm_seq_region, $i)) {
      if(defined($begin_chunk_region)) {
        #this is the end of an unregistered region.
        my $region = [($begin_chunk_region   << $CHUNKFACTOR),
                      ($end_chunk_region     << $CHUNKFACTOR)-1];
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
                  ($end_chunk_region   << $CHUNKFACTOR) -1];
    push @chunk_regions, $region;
  }

  return if(!@chunk_regions);

  # keep the Mapper to a reasonable size
  if( $asm_mapper->size() > $MAX_PAIR_COUNT ) {
    $asm_mapper->flush();
    #we now have to go and register the entire requested region since we 
    #just flushed everything
    
    @chunk_regions = ( [ ( $start_chunk << $CHUNKFACTOR)
                         , (($end_chunk+1) << $CHUNKFACTOR)-1 ] );

    for( my $i = $start_chunk; $i <= $end_chunk; $i++ ) {
      $asm_mapper->register_assembled( $asm_seq_region, $i );
    }
  }

  my $asm_seq_region_id =
    $self->_seq_region_name_to_id($asm_seq_region,$asm_cs_id);

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
    $sth->execute($asm_seq_region_id, $region_start, $region_end, $cmp_cs_id);

    my($cmp_start, $cmp_end, $cmp_seq_region_id, $cmp_seq_region, $ori,
      $asm_start, $asm_end);

    $sth->bind_columns(\$cmp_start, \$cmp_end, \$cmp_seq_region_id,
                       \$cmp_seq_region, \$ori, \$asm_start, \$asm_end);

    #
    # Load the unregistered regions of the mapper
    #
    while($sth->fetch()) {
      next if($asm_mapper->have_registered_component($cmp_seq_region));
      $asm_mapper->register_component($cmp_seq_region);
      $asm_mapper->mapper->add_map_coordinates(
                 $asm_seq_region, $asm_start, $asm_end,
                 $ori,
                 $cmp_seq_region, $cmp_start, $cmp_end);
      $self->{'_sr_id_cache'}->{"$cmp_seq_region:$cmp_cs_id"} =
        $cmp_seq_region_id;
    }
  }

  $sth->finish();
}



sub _seq_region_name_to_id {
  my $self    = shift;
  my $sr_name = shift;
  my $cs_id   = shift;

  ($sr_name && $cs_id) || throw('seq_region_name and coord_system_id args ' .
				'are required');

  my $sr_id = $self->{'_sr_id_cache'}->{"$sr_name:$cs_id"};

  return $sr_id if($sr_id);

  # Get the seq_region_id via the name.  This would be quicker if we just
  # used internal ids instead but stored but then we lose the ability
  # the transform accross databases with different internal ids

  my $sth = $self->prepare("SELECT seq_region_id " .
			   "FROM   seq_region " .
			   "WHERE  name = ? AND coord_system_id = ?");

  $sth->execute($sr_name, $cs_id);

  if(!$sth->rows() == 1) {
    throw("Ambiguous or non-existant seq_region [$sr_name] " .
	  "in coord system $cs_id");
  }

  ($sr_id) = $sth->fetchrow_array();
  $sth->finish();

  $self->{'_sr_id_cache'}->{"$sr_name:$cs_id"} = $sr_id;

  return $sr_id;
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
  my $asm_cs_id = $asm_mapper->assembled_CoordSystem()->dbID();

  #do nothing if this region is already registered
  return if($asm_mapper->have_registered_component($cmp_seq_region));

  my $cmp_seq_region_id =
    $self->_seq_region_name_to_id($cmp_seq_region, $cmp_cs_id);

  # Determine what part of the assembled region this component region makes up

  my $q = qq{
      SELECT
         asm.asm_start,
         asm.asm_end,
         asm.asm_seq_region_id,
         sr.name
      FROM
         assembly asm, seq_region sr
      WHERE
         asm.cmp_seq_region_id = ? AND
         asm.asm_seq_region_id = sr.seq_region_id AND
         sr.coord_system_id = ?
   };

  my $sth = $self->prepare($q);
  $sth->execute($cmp_seq_region_id, $asm_cs_id);

  if($sth->rows() == 0) {
    #this component is not used in the assembled part i.e. gap
    $asm_mapper->register_component($cmp_seq_region_id);
    return;
  }

  #we do not currently support components mapping to multiple assembled
  if($sth->rows() != 1) {
    throw("Multiple assembled regions for single " .
          "component region cmp_seq_region_id=[$cmp_seq_region_id]");
  }

  my ($asm_start, $asm_end, $asm_seq_region_id, $asm_seq_region) =
    $sth->fetchrow_array();

  $self->{'_sr_id_cache'}->{"$asm_seq_region:$asm_cs_id"} = $asm_seq_region_id;

  $sth->finish();

  # Register the corresponding assembled region. This allows a us to
  # register things in assembled chunks which allows us to:
  # (1) Keep track of what assembled regions are registered
  # (2) Use locality of reference (if they want something in same general
  #     region it will already be registered).

  $self->register_assembled($asm_mapper,$asm_seq_region,$asm_start,$asm_end);
}



=head register_chained

  Arg [1]    : Bio::EnsEMBL::ChainedAssemblyMapper $casm_mapper
               The chained assembly mapper to register regions on
  Arg [2]    : string $from ('first' or 'last')
               The direction we are registering from, and the name of the
               internal mapper.
  Arg [3]    : string $seq_region_name
               The name of the seqregion we are registering on
  Arg [4]    : listref $ranges
               A list  of ranges to register (in [$start,$end] tuples).

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


=cut

sub register_chained {
  my $self = shift;
  my $casm_mapper = shift;
  my $from = shift;
  my $seq_region_name = shift;
  my $ranges = shift;

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
  my $mid_cs     = $casm_mapper->middle_CoordSystem();
  my $mid_name   = 'middle';
  my $csa = $self->db->get_CoordSystemAdaptor();

  # Check for the simple case where the ChainedMapper is short
  if( ! defined $mid_cs ) {
      $start_mid_mapper = $combined_mapper;
  }

  #the SQL varies depending on whether we are coming from assembled or
  #component coordinate system
  #print STDERR "ASM SQL:";
  my $asm2cmp_sth = $self->prepare(
     'SELECT
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
	 sr.coord_system_id = ?');

  #print STDERR "CMP SQL:";
  my $cmp2asm_sth = $self->prepare(
      'SELECT
         asm.asm_start,
         asm.asm_end,
         asm.asm_seq_region_id,
         sr.name,
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
	 sr.coord_system_id = ?');

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

  if(@path != 2) {
    my $path = join(',', map({$_->name .' '. $_->version} @path));
    my $len  = scalar(@path) - 1;
    throw("Unexpected mapping path between start and intermediate " .
	  "coord systems (". $start_cs->name . " " . $start_cs->version .
	  " and " . $mid_cs->name . " " . $mid_cs->version . ")." .
	  "\nExpected path length 1, got $len. " .
	  "(path=$path)");
  }

  my $sth;
  my ($asm_cs,$cmp_cs) = @path;
  $sth = ($asm_cs->equals($start_cs)) ? $asm2cmp_sth : $cmp2asm_sth;

  my $seq_region_id = $self->_seq_region_name_to_id($seq_region_name,
						    $start_cs->dbID());
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
    $sth->execute($seq_region_id, $start, $end, $mid_cs_id);

    #load the start <-> mid mapper with the results and record the mid cs
    #ranges we just added to the mapper

    my ($mid_start, $mid_end, $mid_seq_region_id, $mid_seq_region,
        $ori, $start_start, $start_end);

    $sth->bind_columns(\$mid_start, \$mid_end, \$mid_seq_region_id,
		       \$mid_seq_region, \$ori, \$start_start,
		       \$start_end);

    while($sth->fetch()) {
      if( defined $mid_cs ) {
        $start_mid_mapper->add_map_coordinates
          (
           $seq_region_name,$start_start, $start_end, $ori,
           $mid_seq_region, $mid_start, $mid_end
          );
      } else {
        if( $from eq "first" ) {
          $combined_mapper->add_map_coordinates
            (
             $seq_region_name,$start_start, $start_end, $ori,
             $mid_seq_region, $mid_start, $mid_end
            );
        } else {
          $combined_mapper->add_map_coordinates
            (
             $mid_seq_region, $mid_start, $mid_end, $ori,
             $seq_region_name,$start_start, $start_end
            );
        }
      }

      #update sr_name cache
      $self->{'_sr_id_cache'}->{"$mid_seq_region:$mid_cs_id"} =
        $mid_seq_region_id;

      push @mid_ranges,[$mid_seq_region_id,$mid_seq_region,
                        $mid_start,$mid_end];
      push @start_ranges, [ $start_start, $start_end ];

      #the region that we actually register may actually be larger or smaller
      #than the region that we wanted to register.
      #register the intersection of the region so we do not end up doing 
      #extra work later

      if($start_start < $start || $start_end > $end) {
        $start_registry->check_and_register($seq_region_name,$start_start,
                                            $start_end);
      }
    }
  }

  # in the one step case, we load the mid ranges in the
  # last_registry and we are done
  if( ! defined $mid_cs ) {
    for my $range ( @mid_ranges ) {
      $end_registry->check_and_register( $range->[1], $range->[2],
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
  if(@path != 2) {
    my $path = join(',', map({$_->name .' '. $_->version} @path));
    my $len = scalar(@path)-1;
    throw("Unexpected mapping path between intermediate and last" .
	  "coord systems (". $mid_cs->name . " " . $mid_cs->version .
	  " and " . $end_cs->name . " " . $end_cs->version . ")." .
	  "\nExpected path length 1, got $len. " .
	  "(path=$path)");
  }

  ($asm_cs,$cmp_cs) = @path;
  $sth = ($asm_cs->equals($mid_cs)) ? $asm2cmp_sth : $cmp2asm_sth;

  my $end_cs_id = $end_cs->dbID();
  foreach my $mid_range (@mid_ranges) {
    my ($mid_seq_region_id, $mid_seq_region,$start, $end) = @$mid_range;
    $sth->execute($mid_seq_region_id, $start, $end, $end_cs_id);
    #print STDERR "bind vals=($mid_seq_region_id, $start, $end, $mid_cs_id)\n";

    #load the end <-> mid mapper with the results and record the mid cs
    #ranges we just added to the mapper

    my ($end_start, $end_end, $end_seq_region_id, $end_seq_region,
        $ori, $mid_start, $mid_end);

    $sth->bind_columns(\$end_start, \$end_end, \$end_seq_region_id,
		       \$end_seq_region, \$ori, \$mid_start,
		       \$mid_end);

    while($sth->fetch()) {
      #print STDERR "Adding to end<->mid mapper:\n" .
      #      "$end_seq_region:$end_start-$end_end<->$mid_seq_region:" .
      #      "$mid_start-$mid_end($ori)\n";
      $end_mid_mapper->add_map_coordinates
        (
         $end_seq_region, $end_start, $end_end, $ori,
         $mid_seq_region, $mid_start, $mid_end
        );

      #update sr_name cache
      $self->{'_sr_id_cache'}->{"$end_seq_region:$end_cs_id"} =
        $end_seq_region_id;

      #register this region on the end coord system
      $end_registry->check_and_register($end_seq_region, $end_start, $end_end);
    }
  }

  #########
  # Now that both halves are loaded
  # Do stepwise mapping using both of the loaded mappers to load
  # the final start <-> end mapper
  #

  foreach my $range (@start_ranges) {
    my ($start, $end) = @$range;

    my $sum = 0;

    my @initial_coords = $start_mid_mapper->map_coordinates($seq_region_name,
                                                            $start,$end,1,
                                                            $start_name);

    foreach my $icoord (@initial_coords) {
      #skip gaps
      if($icoord->isa('Bio::EnsEMBL::Mapper::Gap')) {
        $sum += $icoord->length();
        next;
      }

      #print STDERR "icoord: id=".$icoord->id." start=".$icoord->start." end=".
      #             $icoord->end."\n";

      #feed the results of the first mapping into the second mapper
      my @final_coords =
        $end_mid_mapper->map_coordinates($icoord->id, $icoord->start,
                                         $icoord->end,
                                         $icoord->strand, $mid_name);

      my $istrand = $icoord->strand();
      foreach my $fcoord (@final_coords) {
        #load up the final mapper
        if($fcoord->isa('Bio::EnsEMBL::Mapper::Coordinate')) {
          my $total_start = $start + $sum;
          my $total_end   = $total_start + $fcoord->length - 1;
          my $ori = $fcoord->strand();

          if($from eq 'first') { #the ordering we add coords must be consistant
            $combined_mapper->add_map_coordinates(
                             $seq_region_name, $total_start, $total_end, $ori,
                             $fcoord->id(), $fcoord->start(), $fcoord->end());
          } else {
            $combined_mapper->add_map_coordinates(
                        $fcoord->id(), $fcoord->start(), $fcoord->end(),$ori,
                        $seq_region_name, $total_start, $total_end);
          }

          #print STDERR "  fcoord: id=".$fcoord->id." start=".
          #  $fcoord->start." end=".$fcoord->end."\n";
          #print STDERR "Loading combined mapper with : " ,
          #  "$seq_region_name:$total_start-$total_end, ($ori) <-> "
          #   .$fcoord->id.":".$fcoord->start."-".$fcoord->end."\n";
        }
        $sum += $fcoord->length();
      }
    }
  }
  #all done!
}


=head2 deleteObj

  Arg [1]    : none
  Example    : none
  Description: Cleans up this objects references to other objects so that
               proper garbage collection can occur
  Returntype : none
  Exceptions : none
  Caller     : Bio::EnsEMBL::DBConnection

=cut

sub deleteObj {
  my $self = shift;

  delete $self->{'_asm_mapper_cache'};
  $self->SUPER::deleteObj();
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

=cut

sub seq_regions_to_ids {
  my $self = shift;
  my $coord_system = shift;
  my $seq_regions = shift;

  my $cs_id = $coord_system->dbID();

  my @out;

  foreach my $sr (@$seq_regions) {
    my $id = $self->{'_sr_id_cache'}->{"$sr:$cs_id"};
    $id = $self->_seq_region_name_to_id($sr,$cs_id) if(!$id);
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

  return $self->fetch_by_CoordSystems($cs1,$cs2);
}



1;
