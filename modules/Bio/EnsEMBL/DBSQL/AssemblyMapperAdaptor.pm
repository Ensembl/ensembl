
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

  my $asm_adptr = $dba->get_AssemblyMapperAdaptor;

  my $chr_ctg_mapper = $asm_adptr->fetch_by_coord_systems('chromosome',
                                                          'contig',
                                                          'NCBI33');

  my $asm_asm_mapper = $asm_adptr->fetch_by_coord_systems('chromosome',
                                                          'chromosome',
                                                          'NCBI33',
                                                          'NCBI34');

  my $ctg_clone_mapper = $asm_adptr->fetch_by_coord_systems('contig',
                                                            'clone');

  $asm_adptr->register($chr_ctg_mapper, $chr_slice);


=head1 DESCRIPTION

Adaptor for handling Assembly mappers.  This is a
I<Singleton> class.  ie: There is only one per
database (C<DBAdaptor>).

This is used to retrieve mappers between any two coordinate systems whose
makeup is described by the assembly table.  Currently only pairwise mapping
is supported (i.e there must be an explicit relationship between the coordinate
systems in the assembly table), but in the future 'chained' mapping between
coordinate systems may be possible.

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

  Arg [1]    : Bio::EnsEMBL::DBAdaptor $dbadaptor the adaptor for the database
               this assembly mapper is using.
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

  return $self;
}


=head2 fetch_by_coord_systems

  Arg [1]    : string $cs1
               The name of one of the 1st coord_system to use for mapping
  Arg [2]    : string $cs2
               The name of the 2nd coord_system to use for mapping
  Arg [3]    : (optional) string $cs1_version
               The version of the 2st coord_system to use for mapping
  Arg [4]    : (optional) string $cs2_version
               The version of the 2nd coord_system to use for mapping
  Description: Retrieves an Assembly mapper for two coordinate systems whose
               relationship is described in the assembly table.

               Versions of the coordinate systems may optionally be provided.
               If they are not provided the default version of each system
               will be used instead.

               The ordering of the coodinate systems is arbitrary.  The
               following two statements are equivalent:
               $mapper = $asm_mapper_adaptor->fetch_by_coord_systems('contig',
                                                                     'clone');
               $mapper = $asm_mapper_adaptor->fetch_by_coord_systems('clone',
                                                                     'contig');

  Returntype : Bio::EnsEMBL::AssemblyMapper
  Exceptions : none
  Caller     : Bio::EnsEMBL::DBSQL::AssemblyMapperAdaptor

=cut

sub fetch_by_coord_systems {
  my($self, $cs1, $cs2, $cs1_version, $cs2version) = @_;

  throw('Coord_system1 argument is required') if(!$cs1);
  throw('Coord_system2 argument is required') if(!$cs1);

  my $csa = $self->db()->get_CoordSystemAdaptor();

  my($cs1_id,$cs1,$cs1_version) = $csa->fetch_by_name($cs1,$cs1_version);
  my($cs2_id,$cs2,$cs2_version) = $csa->fetch_by_name($cs2,$cs2_version);

  if($cs1_id == $cs2_id) {
    throw("Cannot create mapper between same coord system and version:\n" .
          "  [$cs1|$cs1_version] and [$cs1|cscs2_version]");
  }

  #retrieve the shortest possible mapping path between these systems
  my @mapping_path = @{$csa->get_mapping_path($cs1_id,$cs2_id)};

  if(!@mapping_path) {
    throw("There is no mapping defined between these coord systems:\n" .
          "  [$cs1|$version] and [$cs2|$cs2_version]");
  }

  my $key = join(':', @mapping_path);

  my $asm_mapper = $self->{'_asm_mapper_cache'}->{$key};

  return $asm_mapper if($asm_mapper);

  if(@mapping_path > 2) {
    throw("Only explicit (1 step) coordinate system mapping is currently\n" .
          "supported.  Mapping between [$cs1|$cs1_version] and " .
          "[$cs2|$cs2_version]\n requires ". (scalar(@mapping_path)-1) .
          " steps.");
  }

  $asm_mapper = Bio::EnsEMBL::AssemblyMapper->new(@mapping_path);

  $self->{'_asm_mapper_cache'}->{$key} = $asm_mapper;

  return $asm_mapper;
}



=head2 fetch_by_type

  Arg [1]    : char $type 
  Description: Fetches a Bio::EnsEMBL::AssemblyMapper object from the adaptor 
               for a particular assembly (golden path) type (e.g. NCBI_xx). 
               The result is cached for a particular assembly type. 
  Returntype : Bio::EnsEMBL::AssemblyMapper 
  Exceptions : none
  Caller     : Bio::EnsEMBL::DBSQL::AssemblyMapperAdaptor 

=cut

sub fetch_by_type{
  my ($self,$type) = @_;

  deprecated('Use fetch_by_coord_systems instead');

  #assume that what the user wanted was a mapper between the sequence coord
  #level and the top coord level

  my $csa = $self->db()->get_CoordSystemAdaptor();

  my ($tl_id,$tl,$tl_version) = $csa->fetch_top_level($type);
  my ($sl_id,$sl,$sl_version) = $csa->fetch_sequence_level();

  return $self->fetch_by_coord_systems($tl,$sl,$tl_version,$sl_version);
}




sub register_Slice {
  my $self = shift;
  my $asm_mapper = shift;
  my $slice = shift;

  if(!$asm_mapper || !ref($asm_mapper) ||
     !$asm_mapper->isa('Bio::EnsEMBL::AssemblyMapper')) {
    throw('AssemblyMapper argument required');
  }

  if(!$slice || !ref($slice) || !$slice->isa('Bio::EnsEMBL::Slice')) {
    throw('Slice argument required');
  }

  my $cmp_cs_id = $asm_mapper->component_coord_system_id();
  my $asm_cs_id = $asm_mapper->assembled_coord_system_id();

  my $csa = $self->db()->get_CoordSystemAdaptor();

  my ($slice_cs_id) =
    $csa->fetch_by_name($slice->coord_system, $slice->version);

  # Get the seq_region id of the slice.  This is probably better stored
  # right on the slice when it is retrieved.
  my $sth = $self->prepare(
    "SELECT seq_region_id" .
    "FROM   seq_region" .
    "WHERE  name = ? AND coord_system_id = ?");

  $sth->execute($slice->seq_region_name, $slice_cs_id);

  if($sth->rows() != 1) {
    throw("Seq_region name=[$name],coord_system_id=[$cs_id] not found");
  }

  my ($seq_region_id) = $sth->fetchrow_array();

  $sth->finish();

  my $start = $slice->start();
  my $end   = $slice->end();

  if($cmp_cs_id == $cmp_cs_id) {
    $self->register_component($asm_mapper, $seq_region_id);
  } elsif($asm_cs_id == $asm_cs_id) {
    $self->register_assembled($asm_mapper, $seq_region_id, $start, $end);
  } else {
    my $name = $slice->name();
    throw("Slice [$slice] is not in one of the coordinate systems ".
          "defined in the provided assembly mapper");
  }
}


sub register_assembled {
  my $self = shift;
  my $asm_mapper = shift;
  my $asm_seq_region_id = shift;
  my $asm_start         = shift;
  my $asm_end           = shift;


  #
  # TBD: somehow chunk assembled regions up and register the chunks
  #

  # Retrieve the description of how the assembled region is made from
  # component regions

  my $q = qq{
      SELECT
         asm.cmp_start,
         asm.cmp_end,
         asm.cmp_seq_region_id,
         asm.ori,
         asm.asm_start,
         asm.asm_end
      FROM
         assembly asm
      WHERE
         asm.asm_seq_region_id = ? AND
         ? <= asm.asm_end AND
         ? >= asm.asm_start
   };

  my $sth = $self->prepare($q);
  $sth->execute($asm_seq_region_id, $asm_start, $asm_end);

  my($cmp_start, $cmp_end, $cmp_seq_region_id, $ori);
  $sth->bind_columns(\$cmp_start, \$cmp_end, \$cmp_seq_region_id, \$ori);

  #
  # Load the unregistered regions of the mapper
  #
  while($sth->fetch()) {
    next if($asm_mapper->have_registered_component($cmp_seq_region_id));

    $asm_mapper->register_component($cmp_seq_region_id);
    $asm_mapper->mapper->add_map_coordinates(
					       $cmp_seq_region_id, $cmp_start, $cmp_end,
                 $ori,
					       $asm_seq_region_id, $asm_start, $asm_end);
  }

  $sth->finish();
}


sub register_component {
  my $self = shift;
  my $asm_mapper = shift;
  my $cmp_seq_region_id = shift;

  #do nothing if this region is already registered
  return if($asm_mapper->have_registered_component($cmp_seq_region_id));

  # Determine what part of the assembled region this component region makes up

  my $q = qq{
      SELECT
         asm.asm_start,
         asm.asm_end,
         asm.asm_seq_region_id
      FROM
         assembly asm
      WHERE
         asm.cmp_seq_region_id = ?
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

  my ($asm_start, $asm_end, $asm_seq_region_id)= $sth->fetchrow_array();

  $sth->finish();

  # Register the corresponding assembled region. This allows a us to
  # register things in assembled chunks which allows us to:
  # (1) Keep track of what assembled regions are registered
  # (2) Use locality of reference (if they want something in same general
  #     region it will already be registered).

  $self->register_assembled($asm_mapper,$asm_seq_region_id,$start,$end);
}


=head2 register_region

  Arg [1]    : Bio::EnsEMBL::AssemblyMapper $assmapper
	       A valid AssemblyMapper object
  Arg [2]    : char $type
	       golden path type (e.g. NCBI_xx)
  Arg [3]    : char $chr_name
	       chromosome name (e.g. '2', 'X')
  Arg [4]    : int $start
	       chromosome start coordinate
  Arg [5]    : int $end
	       chromosome end coordinate
  Description: Declares a chromosomal region to the AssemblyMapper.
               This extracts the relevant data from the assembly
               table and stores it in a Bio::EnsEMBL::Mapper.
               It therefore must be called before any mapping is
               attempted on that region. Otherwise only gaps will
               be returned!
  Returntype : none
  Exceptions : thrown if $assmapper arg is not a Bio::EnsEMBL::AssemblyMapper
  Caller     : Bio::EnsEMBL::DBSQL::AssemblyMapperAdaptor

=cut

sub register_region{
  my ($self, $assmapper, $type, $chr_name, $start, $end) = @_;

  $self->throw("$assmapper is not a Bio::EnsEMBL::AssemblyMapper")
    unless $assmapper->isa("Bio::EnsEMBL::AssemblyMapper");

  my $chr = $self->db->get_ChromosomeAdaptor()->fetch_by_chr_name( $chr_name );
  my $chr_id = $chr->dbID();
  my $max_assembly_contig = $self->db()->get_MetaContainer->get_max_assembly_contig();

  my $select = qq{
      select
         ass.contig_start,
         ass.contig_end,
         ass.contig_id,
         ass.contig_ori,
         ass.chr_start,
         ass.chr_end
      from
         assembly ass
      where
         ass.chromosome_id = $chr_id and
         ? <= ass.chr_end  and
         ? >= ass.chr_start  and
	 ? <= ass.chr_start and
         ass.type = '$type'
   };

  my $sth = $self->prepare($select);
   
  $sth->execute( $start, $end, $start-$max_assembly_contig );
   
  while( my $arrayref = $sth->fetchrow_arrayref ) {
    my ($contig_start, $contig_end, $contig_id, $contig_ori,
        $chr_start, $chr_end) = @$arrayref;
    if( $assmapper->_have_registered_contig($contig_id) == 0 ) {
      $assmapper->_register_contig($contig_id);
      $assmapper->_mapper->add_map_coordinates(
					       $contig_id, $contig_start, 
					       $contig_end, $contig_ori,
					       $chr_name, $chr_start, $chr_end
					      );
    }
  } 
}


=head2 register_contig

  Arg  [1]   : Bio::EnsEMBL::AssemblyMapper $assmapper
	       A valid AssemblyMapper object
  Arg  [2]   : char $type
	       golden path type (e.g. NCBI_xx)
  Arg  [3]   : int $contig_id
	       RawContig internal ID

  Description: Declares a chromosomal region to the AssemblyMapper.
               This extracts the relevant data from the assembly
               table and stores it in a Bio::EnsEMBL::Mapper.
               It therefore must be called before any mapping is
               attempted on that region. Otherwise only gaps will be returned!
  Returntype : 1 if the contig is present in assembly, 
               0 if the contig is absent
  Exceptions : none
  Caller     : Bio::EnsEMBL::DBSQL::AssemblyMapperAdaptor

=cut

sub register_contig {
   my ($self, $assmapper, $type, $contig_id ) = @_;

   my $sth = $self->prepare(qq{
      select
	c.name,
	a.chr_start,
	a.chr_end
      from
	assembly a,
	chromosome c
      where
	 contig_id = '$contig_id' and
	 type = '$type' and
	 c.chromosome_id = a.chromosome_id
   });

   $sth->execute();

   my @ctg_list;

   while (my $row = $sth->fetchrow_arrayref) {
       push @ctg_list, $row;
   }

   if (@ctg_list == 0) {
     #$self->warn("Not found contig $contig_id");
     return ();
   }

   if (@ctg_list > 1) {
     $self->warn("Contig $contig_id is ambiguous in assembly type $type");
     return ();
   }

   my $chr_name = $ctg_list[0]->[0];
   my $chr_start = $ctg_list[0]->[1];
   my $chr_end = $ctg_list[0]->[2];

   return ($chr_name, $chr_start, $chr_end);
}




1;
