=head1 LICENSE

  Copyright (c) 1999-2012 The European Bioinformatics Institute and
  Genome Research Limited.  All rights reserved.

  This software is distributed under a modified Apache license.
  For license details, please see

    http://www.ensembl.org/info/about/code_licence.html

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <dev@ensembl.org>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk@ensembl.org>.

=cut

=head1 NAME

Bio::EnsEMBL::DBSQL::DBAdaptor

=head1 SYNOPSIS

  $db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
    -user   => 'root',
    -dbname => 'pog',
    -host   => 'caldy',
    -driver => 'mysql'
  );

  $gene_adaptor = $db->get_GeneAdaptor();

  $gene = $gene_adaptor->fetch_by_stable_id($stable_id);

  $slice =
    $db->get_SliceAdaptor()->fetch_by_chr_start_end( 'X', 1, 10000 );

=head1 DESCRIPTION

Formerly this class provided database connectivity and a means
to retrieve object adaptors.  This class is now provided for
convenience and backwards compatibility, and delegates its connection
responsibilities to the DBConnection class (no longer inherited from)
and its object adaptor retrieval to the static Bio::EnsEMBL::Registry.

Please use Bio::EnsEMBL::Registry to retrieve object adaptors.

=head1 METHODS

=cut

package Bio::EnsEMBL::DBSQL::DBAdaptor;

use strict;

use Bio::EnsEMBL::DBSQL::DBConnection;
use Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor;
use Bio::EnsEMBL::Utils::SeqRegionCache;
use Bio::EnsEMBL::Utils::Exception qw(throw warning deprecate);
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Utils::ConfigRegistry;

my $reg = "Bio::EnsEMBL::Registry";

=head2 new

  Arg [-DNADB]: (optional) Bio::EnsEMBL::DBSQL::DBAdaptor DNADB 
               All sequence, assembly, contig information etc, will
               be retrieved from this database instead.

  Arg [-NO_CACHE]: (optional) int 1
               This option will turn off caching for slice features,
               so, every time a set of features is retrieved,
               they will come from the database instead of the
               cache.  This option is only recommended for advanced
               users, specially if you need to store and retrieve
               features.  It might reduce performance when querying
               the database if not used properly.  If in doubt, do
               not use it or ask in the developer mailing list.

  Arg [..]   : Other args are passed to superclass
               Bio::EnsEMBL::DBSQL::DBConnection

  Example    : $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
                -user   => 'root',
                -dbname => 'pog',
                -host   => 'caldy',
                -driver => 'mysql'
              );

              $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
                -species => 'Homo_sapiens',
                -group   => 'core',
                -user    => 'root',
                -dbname  => 'pog',
                -host    => 'caldy',
                -driver  => 'mysql'
              );

              $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
                -species         => 'staphylococcus_aureus',
                -group           => 'core',
                -user            => 'root',
                -dbname          => 'staphylococcus_collection_1_52_1a',
                -multispecies_db => 1,
                -host            => 'caldy',
                -driver          => 'mysql'
              );

  Description: Constructor for DBAdaptor.
  Returntype : Bio::EnsEMBL::DBSQL::DBAdaptor
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub new {
  my ( $class, @args ) = @_;

  my $self = bless {}, $class;

  my ( $is_multispecies, $species, $species_id, $group, $con, $dnadb,
    $no_cache, $dbname )
    = rearrange( [
      'MULTISPECIES_DB', 'SPECIES', 'SPECIES_ID', 'GROUP',
      'DBCONN',          'DNADB',   'NO_CACHE',   'DBNAME'
    ],
    @args
    );

  if ( defined($con) ) { $self->dbc($con) }
  else {
    if(! defined $dbname) {
      throw "-DBNAME is a required parameter when creating a DBAdaptor";
    }
    $self->dbc( new Bio::EnsEMBL::DBSQL::DBConnection(@args) );
  }

  if ( defined($species) ) { $self->species($species) }
  if ( defined($group) )   { $self->group($group) }

 
  $self = Bio::EnsEMBL::Utils::ConfigRegistry::gen_load($self);

#  if(!defined($species) ){
#     $reg->find_and_add_aliases($self);
#  }

  $self->species_id( $species_id || 1 );

  $self->is_multispecies( defined($is_multispecies)
                          && $is_multispecies == 1 );

  if ( defined($dnadb) )    { $self->dnadb($dnadb) }
  if ( defined($no_cache) ) { $self->no_cache($no_cache) }

  return $self;
} ## end sub new

=head2 clear_caches

  Example			: $dba->clear_caches();
  Description	: Loops through all linked adaptors and clears their 
                caches if C<clear_cache()> is implemented. Not all caches
                are cleared & the DBAdaptor instance should be removed from
                the registry to clear these remaining essential caches. 
  Returntype 	: None
  Exceptions 	: None

=cut

sub clear_caches {
  my ($self) = @_;
  my $adaptors = Bio::EnsEMBL::Registry->get_all_adaptors(
    $self->species(), $self->group());
  foreach my $adaptor (@{$adaptors}) {
    if($adaptor->can('clear_cache')) {
      $adaptor->clear_cache();
    }
  }
  return;
}

=head2 dbc

  Arg[1]    : (optional) Bio::EnsEMBL::DBSQL::DBConnection

  Example    : $dbc = $dba->dbc();
  Description: Getter/Setter for DBConnection.
  Returntype : Bio::EnsEMBL::DBSQL::DBConnection
  Exceptions : throws if argument not a Bio::EnsEMBL::DBSQL::DBConnection
  Caller     : general
  Status     : Stable

=cut

sub dbc{
  my $self  = shift;
  
  if(@_){
    my $arg = shift;
    if(defined($arg)){
      if(!$arg->isa('Bio::EnsEMBL::DBSQL::DBConnection')){
	throw("$arg is no a DBConnection\n");
      }
    }
    $self->{_dbc} = $arg;
  }
  return $self->{_dbc};
}



=head2 add_db_adaptor

  Arg [1]    : string $name
               the name of the database to attach to this database
  Arg [2]    : Bio::EnsEMBL::DBSQL::DBConnection
               the db adaptor to attach to this database
  Example    : $db->add_db_adaptor('lite', $lite_db_adaptor);
  Description: Attaches another database instance to this database so 
               that it can be used in instances where it is required.
  Returntype : none
  Exceptions : none
  Caller     : EnsWeb
  Status     : At Risk
             : may get deprecated, please use add_db from the registry instead

=cut

sub add_db_adaptor {
  my ($self, $name, $adaptor) = @_;

  unless($name && $adaptor && ref $adaptor) {
    throw('adaptor and name arguments are required');
  }

  Bio::EnsEMBL::Registry->add_db($self, $name, $adaptor);

}


=head2 remove_db_adaptor

  Arg [1]    : string $name
               the name of the database to detach from this database.
  Example    : $lite_db = $db->remove_db_adaptor('lite');
  Description: Detaches a database instance from this database and returns
               it.
  Returntype : none
  Exceptions : none
  Caller     : ?
  Status     : At Risk
             : mey get deprecated, use remove_db instead from the Registry

=cut

sub remove_db_adaptor {
  my ($self, $name) = @_;

  return Bio::EnsEMBL::Registry->remove_db($self, $name);
}


=head2 get_all_db_adaptors

  Arg [1]    : none
  Example    : @attached_dbs = values %{$db->get_all_db_adaptors()};
  Description: returns all of the attached databases as 
               a hash reference of key/value pairs where the keys are
               database names and the values are the attached databases  
  Returntype : hash reference with Bio::EnsEMBL::DBSQL::DBConnection values
  Exceptions : none
  Caller     : Bio::EnsEMBL::DBSQL::ProxyAdaptor
  Status     : At Risk
             : may get deprecated soon
             : please use  Bio::EnsEMBL::Registry->get_all_db_adaptors

=cut

sub get_all_db_adaptors {
  my ($self) = @_;
  return Bio::EnsEMBL::Registry->get_all_db_adaptors($self);
}



=head2 get_db_adaptor

  Arg [1]    : string $name
               the name of the attached database to retrieve
  Example    : $lite_db = $db->get_db_adaptor('lite');
  Description: returns an attached db adaptor of name $name or undef if
               no such attached database exists
  Returntype : Bio::EnsEMBL::DBSQL::DBConnection
  Exceptions : none
  Caller     : ?
  Status     : At Risk
             : may get deprecated soon
             : please use  Bio::EnsEMBL::Registry->get_db_adaptors

=cut

sub get_db_adaptor {
  my ($self, $name) = @_;

  return Bio::EnsEMBL::Registry->get_db($self, $name);
}

=head2 get_available_adaptors

  Example    : my %pairs = %{$dba->get_available_adaptors()};
  Description: gets a hash of the available adaptors
  ReturnType : reference to a hash
  Exceptions : none
  Caller     : Bio::EnsEMBL::Utils::ConfigRegistry
  Status     : Stable

=cut 

sub get_available_adaptors {
  my %pairs = (
    # Firstly those that just have an adaptor named after there object
    # in the main DBSQL directory.
    map( { $_ => "Bio::EnsEMBL::DBSQL::${_}Adaptor" } qw(
        Analysis                 ArchiveStableId      Attribute
        AssemblyExceptionFeature AssemblyMapper       CoordSystem
        CompressedSequence       DBEntry              DnaAlignFeature
        DensityFeature           DensityType          Exon
        Gene                     KaryotypeBand        MiscSet
        MiscFeature              PredictionTranscript PredictionExon
        ProteinFeature           ProteinAlignFeature  RepeatConsensus
        RepeatFeature            Sequence             SeqRegionSynonym  SimpleFeature
        Slice                    SupportingFeature    Transcript
        TranscriptSupportingFeature Translation       UnmappedObject
        UnconventionalTranscriptAssociation           AssemblySlice
        SplicingEvent            SplicingEventFeature SplicingTranscriptPair
        Operon 			 OperonTranscript
        DataFile                 Assembly
        IntronSupportingEvidence
        ) ),
    # Those whose adaptors are in Map::DBSQL
    map( { $_ => "Bio::EnsEMBL::Map::DBSQL::${_}Adaptor" } qw(
        Marker MarkerFeature QtlFeature Qtl Ditag DitagFeature
        ) ),
    # Finally the exceptions... those that have non-standard mapping
    # between object / adaptor ....
    # 'Blast'                => 'Bio::EnsEMBL::External::BlastAdaptor',
    'MetaCoordContainer' => 'Bio::EnsEMBL::DBSQL::MetaCoordContainer',
    'MetaContainer'      => 'Bio::EnsEMBL::DBSQL::MetaContainer',
    'SNP'                => 'Bio::EnsEMBL::DBSQL::ProxySNPAdaptor',
  );

  return ( \%pairs );
} ## end sub get_available_adaptors

###########################################################
#
# Support for DAS
#
###########################################################

=head2 add_DASFeatureFactory

  Arg [1]    : Bio::EnsEMBL::ExternalFeatureFactory $value 
  Example    : none
  Description: Attaches a DAS Feature Factory to this method.  
               ExternalFeatureFactory objects are not really used right now.
               They may be reintroduced or taken out completely.  The fate
               of this function is unknown (although it is presently needed).
  Returntype : none
  Exceptions : none
  Caller     : EnsWeb
  Status     : At Risk
             : with the new web code this may not be needed/supported

=cut

sub add_DASFeatureFactory{
 
 my ($self,$value) = @_;
  
  push(@{$self->{'_das_ff'}},$value);
}


sub remove_all_DASFeatureFactories {
  $_[0]->{'_das_ff'} = [];
}
=head2 _each_DASFeatureFactory

  Args       : none
  Example    : none
  Description: Not sure if this is used, or if it should be removed.  It 
               does not seem to be used at the moment
  Returntype : Bio::EnsEMBL::ExternalFeatureFactory
  Exceptions : none
  Caller     : ??
  Status     : At Risk
             : with the new web code this may not be needed/supported

=cut

sub _each_DASFeatureFactory{
   my ($self) = @_;

   return @{$self->{'_das_ff'}||[]}
}


################################################################## 
# 
# SUPPORT FOR EXTERNAL FEATURE FACTORIES 
# 
##################################################################



=head2 add_ExternalFeatureAdaptor

  Arg [1]    : Bio::EnsEMBL::External::ExternalFeatureAdaptor
  Example    : $db_adaptor->add_ExternalFeatureAdaptor($xfa);
  Description: Adds an external feature adaptor to this database adaptor.
               Adding the external adaptor in this way allows external
               features to be obtained from Slices and from RawContigs.

               The external feature adaptor which is passed to this method
               will have its db attribuite set to this DBAdaptor object via 
               the db accessor method. 

               ExternalFeatureAdaptors passed to this method are stored 
               internally in a hash keyed on the string returned by the 
               ExternalFeatureAdaptors track_name method.
               
               If the track name method is not implemented then the 
               a default key named 'External features' is assigned.  In the
               event of duplicate key names, a number is appended to the
               key name, and incremented for each subsequent adaptor with the
               same track name.  For example, if no track_names are specified 
               then the the external feature adaptors will be stored under the
               keys 'External features', 'External features2' 
               'External features3' etc.
  Returntype : none
  Exceptions : none
  Caller     : general
  
=cut

sub add_ExternalFeatureAdaptor {
  my ($self, $adaptor) = @_;

  unless($adaptor && ref $adaptor && 
	 $adaptor->isa('Bio::EnsEMBL::External::ExternalFeatureAdaptor')) {
     throw("[$adaptor] is not a " .
           "Bio::EnsEMBL::External::ExternalFeatureAdaptor");
  }

  unless(exists $self->{'_xf_adaptors'}) {
    $self->{'_xf_adaptors'} = {};
  }

  my $track_name = $adaptor->{'_track_name'};
  if(!$track_name) {
    $track_name = $adaptor->track_name();
  }

  #use a generic track name if one hasn't been defined
  unless(defined $track_name) {
    $track_name = "External features";
  }

  #if the track name exists add numbers to the end until a free name is found
  if(exists $self->{'_xf_adaptors'}->{"$track_name"}) {
    my $num = 2;
    $num++ while(exists $self->{'_xf_adaptors'}->{"$track_name$num"});
    $self->{'_xf_adaptors'}->{"$track_name$num"} = $adaptor;
  } else {
    $self->{'_xf_adaptors'}->{"$track_name"} = $adaptor;
  }

  $adaptor->ensembl_db($self);
}



=head2 get_ExternalFeatureAdaptors

  Arg [1]    : none
  Example    : @xfas = values %{$db_adaptor->get_ExternalFeatureAdaptors}; 
  Description: Retrieves all of the ExternalFeatureAdaptors which have been
               added to this DBAdaptor.  The ExternalFeatureAdaptors are 
               returned in a reference to a hash keyed on the track names
               of the external adaptors
  Returntype : Reference to a hash of ExternalFeatureAdaptors keyed on 
               their track names.
  Exceptions : none
  Caller     : general

=cut

sub get_ExternalFeatureAdaptors {
  my $self = shift;

  return $self->{'_xf_adaptors'};
}


=head2 add_ExternalFeatureFactory

  Arg [1]    : Bio::EnsEMBL::DB::ExternalFeatureFactoryI $value
  Example    : $db_adaptor->add_ExternalFeatureFactory
  Description: It is recommended that add_ExternalFeatureAdaptor be used 
               instead.  See documentation for 
               Bio::EnsEMBL::External::ExternalFeatureAdaptor

               Adds an external feature factory to the core database
               so that features from external sources can be displayed in 
               ensembl. This method is still available mainly for legacy
               support for external EnsEMBL installations.
  Returntype : none
  Exceptions : none
  Caller     : external

=cut

sub add_ExternalFeatureFactory{
   my ($self,$value) = @_;

   $self->add_ExternalFeatureAdaptor($value);
}

#
# OVERWRITABLE STANDARD ADAPTORS
#

=head2 get_adaptor

  Arg [1]    : Canonical data type for which an adaptor is required.
  Example    : $db_adaptor->get_adaptor("Protein")
  Description: Gets an adaptor object for a standard data type.
  Returntype : Adaptor Object of arbitrary type or undef
  Exceptions : none
  Caller     : external
  Status     : Medium Risk
             : please use the Registry method, as at some time this
             : may no longer be supprted.
 
=cut

sub get_adaptor {
  my ($self, $canonical_name, @other_args) = @_;

  return $reg->get_adaptor($self->species(),$self->group(),$canonical_name);
}



=head2 set_adaptor

  Arg [1]    : Canonical data type for new adaptor.
  Arg [2]    : Object defining the adaptor for arg1.
  Example    : $aa = Bio::EnsEMBL::DBSQL::GeneAdaptor->new($db_adaptor);
             : $db_adaptor->set_adaptor("Gene", $ga)
  Description: Stores the object which represents the adaptor for the
               arg1 data type.
  Returntype : none
  Exceptions : none
  Caller     : external
  Status     : Medium Risk
             : please use the Registry method, as at some time this
             : may no longer be supprted.
 
=cut

sub set_adaptor {
  my ($self, $canonical_name, $module) = @_;

  $reg->add_adaptor($self->species(),$self->group(),$canonical_name,$module);

  return $module;
}


#
# GENERIC FEATURE ADAPTORS
#

=head2 get_GenericFeatureAdaptors

  Arg [1]    : List of names of feature adaptors to get. If no
               adaptor names are given, all the defined adaptors are returned.
  Example    : $db->get_GenericFeature("SomeFeature", "SomeOtherFeature")
  Description: Returns a hash containing the named feature adaptors (or
               all feature adaptors).
  Returntype : reference to a Hash containing the named
               feature adaptors (or all feature adaptors).
  Exceptions : If any of the the named generic feature adaptors do not exist.
  Caller     : external

=cut

sub get_GenericFeatureAdaptors {

  my ($self, @names) = @_;

  my %adaptors = ();

  if (!@names) {
    %adaptors = %{$self->{'generic_feature_adaptors'}};
  } else {
    foreach my $name (@names) {
      if (!exists($self->{'generic_feature_adaptors'}->{$name})) {
        throw("No generic feature adaptor has been defined for $name" );
      }


      $adaptors{$name} = $self->{'generic_feature_adaptors'}->{$name};
    }
  }

  return \%adaptors;
}


=head2 add_GenericFeatureAdaptor

  Arg [1]    : The name of the feature.
  Arg [2]    : Adaptor object for a generic feature.
  Example    : $db->add_GenericFeatureAdaptor("SomeFeature",
                              "Bio::EnsEMBL::DBSQL::SomeFeatureAdaptor")
  Description: Stores the object which represents the adaptor for the
               named feature type.
  Returntype : none
  Exceptions :
  Caller     : external

=cut

sub add_GenericFeatureAdaptor {
  my ($self, $name, $adaptor_obj) = @_;
	
  # check that $adaptor is an object that subclasses BaseFeatureAdaptor	
  if (!$adaptor_obj->isa("Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor")) {
    throw("$name is a " . ref($adaptor_obj) . "which is not a " .
          "subclass of Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor" );
  }

  $self->{'generic_feature_adaptors'}->{$name} = $adaptor_obj;
}

=head2 species

  Arg [1]    : (optional) string $arg
               The new value of the species used by this DBAdaptor. 
  Example    : $species = $dba->species()
  Description: Getter/Setter for the species of to use for 
               this connection.  There is currently no point in setting 
               this value after the connection has already been established 
               by the constructor.
  Returntype : string
  Exceptions : none
  Caller     : new
  Status     : Stable

=cut

sub species {
  my ( $self, $arg ) = @_;

  if ( defined($arg) ) {
    $self->{_species} = $arg;
  }

  $self->{_species};
}

=head2 all_species

  Args       : NONE
  Example    : @all_species = @{$dba->all_species()};
  Description: Returns the names of all species contained in the
               database to which this DBAdaptor is connected.
  Returntype : array reference
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub all_species {
  my ($self) = @_;

  if ( !$self->is_multispecies() ) { return [ $self->species() ] }

  if ( exists( $self->{'_all_species'} ) ) {
    return $self->{'_all_species'};
  }

  my $statement =
      "SELECT meta_value "
    . "FROM meta "
    . "WHERE meta_key = 'species.db_name'";

  my $sth = $self->dbc()->db_handle()->prepare($statement);

  $sth->execute();

  my $species;
  $sth->bind_columns( \$species );

  my @all_species;
  while ( $sth->fetch() ) { push( @all_species, $species ) }

  $self->{'_all_species'} = \@all_species;

  return $self->{'_all_species'};
} ## end sub all_species


=head2 is_multispecies

  Arg [1]    : (optional) boolean $arg
  Example    : if ($dba->is_multispecies()) { }
  Description: Getter/Setter for the is_multispecies boolean of
               to use for this connection.  There is currently no
               point in setting this value after the connection has
               already been established by the constructor.
  Returntype : boolean
  Exceptions : none
  Caller     : new
  Status     : Stable

=cut

sub is_multispecies {
  my ( $self, $arg ) = @_;

  if ( defined($arg) ) {
    $self->{_is_multispecies} = $arg;
  }

  return $self->{_is_multispecies};
}


=head2 species_id

  Arg [1]    : (optional) string $arg
               The new value of the species_id used by this DBAdaptor
               when dealing with multi-species databases.
  Example    : $species_id = $dba->species_id()
  Description: Getter/Setter for the species_id of to use for this
               connection.  There is currently no point in setting
               this value after the connection has already been
               established by the constructor.
  Returntype : string
  Exceptions : none
  Caller     : new
  Status     : Stable

=cut

sub species_id {
  my ( $self, $arg ) = @_;

  if ( defined($arg) ) {
    $self->{_species_id} = $arg;
  }

  return $self->{_species_id};
}


=head2 no_cache

  Arg [1]    : (optional) int $arg
               The new value of the no cache attribute used by this DBAdaptor. 
  Example    : $no_cache = $dba->no_cache();
  Description: Getter/Setter for the no_cache to use for 
               this connection.  There is currently no point in setting 
               this value after the connection has already been established 
               by the constructor.
  Returntype : int
  Exceptions : none
  Caller     : new
  Status     : Stable

=cut

sub no_cache {
  my ($self, $arg ) = @_;

  if ( defined $arg ){
      if ($arg != 1 && $arg != 0){
	  throw("$arg is not allowed for this attribute. Only value 1|0 is allowed");
      }
      $self->{_no_cache} = $arg;
  }
  $self->{_no_cache};
}


=head2 group

  Arg [1]    : (optional) string $arg
               The new value of the group used by this DBAdaptor. 
  Example    : $group = $dba->group()
  Description: Getter/Setter for the group of to use for 
               this connection.  There is currently no point in setting 
               this value after the connection has already been established 
               by the constructor.
  Returntype : string
  Exceptions : none
  Caller     : new
  Status     : Stable

=cut

sub group {
  my ($self, $arg ) = @_;
  ( defined $arg ) &&
    ( $self->{_group} = $arg );
  $self->{_group};
}

=head2 get_SeqRegionCache

  Arg [1]    : none
  Example    : my $srcache = $dba->get_SeqRegionCache();
  Description: Retrieves a seq_region cache for this database
  Returntype : Bio::EnsEMBL::Utils::SeqRegionCache
  Exceptions : none
  Caller     : SliceAdaptor, AssemblyMapperAdaptor
  Status     : Stable

=cut

sub get_SeqRegionCache {
  my $self = shift;

  # use the cache from the database where seq_regions are stored
  if($self != $self->dnadb()) {
    return $self->dnadb()->get_SeqRegionCache();
  }

  if(!$self->{'seq_region_cache'}) {
    $self->{'seq_region_cache'} = Bio::EnsEMBL::Utils::SeqRegionCache->new();
  }

  return $self->{'seq_region_cache'};
}



#convenient method to retrieve the schema_build version for the database being used

sub _get_schema_build{
  my ($self) = @_;

  #avoided using dnadb by default to avoid obfuscation of behaviour
  
  my @dbname = split/_/, $self->dbc->dbname();

  #warn "dbname is $schema_build";

  my $schema_build = pop @dbname;
  $schema_build = pop(@dbname).'_'.$schema_build;


  return $schema_build;
}


=head2 dnadb

 Title   : dnadb
 Usage   : my $dnadb = $db->dnadb();
 Function: returns the database adaptor where the dna lives
           Useful if you only want to keep one copy of the dna
           on disk but have other databases with genes and features in
 Returns : dna database adaptor
 Args    : Bio::EnsEMBL::DBSQL::BaseAdaptor
 Status  : Medium Risk.
         : Use the Registry method add_DNAAdaptor/get_DNAAdaptor instead

=cut

sub dnadb {
  my $self = shift;

  if(@_) {
    my $arg = shift;
    $reg->add_DNAAdaptor($self->species(),$self->group(),$arg->species(),$arg->group());
  }

#  return $self->{'dnadb'} || $self;
  return $reg->get_DNAAdaptor($self->species(),$self->group()) || $self;
}


use vars '$AUTOLOAD';

sub AUTOLOAD {
  my ( $self, @args ) = @_;

  my $type;
  if ( $AUTOLOAD =~ /^.*::get_(\w+)Adaptor$/ ) {
    $type = $1;
  } elsif ( $AUTOLOAD =~ /^.*::get_(\w+)$/ ) {
    $type = $1;
  } else {
    throw( sprintf( "Could not work out type for %s\n", $AUTOLOAD ) );
  }
  
  my $ret = $reg->get_adaptor( $self->species(), $self->group(), $type );

  return $ret if $ret;
  
  warning( sprintf(
    "Could not find %s adaptor in the registry for %s %s\n",
    $type, $self->species(), $self->group() ) );

  throw( sprintf( 
    "Could not get adaptor %s for %s %s\n",
    $type, $self->species(), $self->group() ) );

} ## end sub AUTOLOAD

sub DESTROY { }    # required due to AUTOLOAD


#########################
# sub DEPRECATED METHODS
#########################
=head2 db
  
  Description: DEPRECATED 
  
=cut

sub db{
  my ($self, $arg ) = @_;
 deprecate("db Should no longer be called from the DBAdaptor. DBConnection should now be used OR preferably the object adaptor itself\n");
 return $self->dbc($arg);
}


sub source { deprecate('Do not use - this method does nothing'); }


=head2 assembly_type

  Description: DEPRECATED - Use CoordSystemAdaptor to obtain default coordinate
               system instead.

=cut

sub assembly_type{
  my $self = shift;

  deprecate('Use CoordSystemAdaptor $csa->fetch_all->[0]->version() instead');

  my $csa = $self->get_CoordSystemAdaptor();
  my ($cs) = @{$csa->fetch_all()};
  return ($cs) ? $cs->version() : undef;
}



=head2 list_supported_assemblies

  Description: DEPRECATED - Use CoordSystemAdaptor to obtain list of top-level
               coordinate systems instead

=cut

sub list_supported_assemblies {
  my($self) = @_;
  deprecate('Use CoordSystemAdaptor::fetch_all instead');

  my $csa = $self->get_CoordSystemAdaptor();
  my %versions;
  foreach my $cs (@{$csa->fetch_all()}) {
    $versions{$cs->version()} = 1;
  }

  return keys %versions;
}


sub prepare{
  my ($self, @args) = @_;

 deprecate("prepare Should no longer be called from the DBAdaptor. DBConnection should now be used OR preferably the object adaptor itself\n");
  $self->dbc->prepare(@args);
}

sub dbname{
  my ($self, @args) = @_;

 deprecate("dbname Should no longer be called from the DBAdaptor. DBConnection should now be used OR preferably the object adaptor itself\n");
  $self->dbc->dbname(@args);
} 

sub disconnect_when_inactive{
  my ($self, @args) = @_;

 deprecate("disconnect_when_inactive Should no longer be called from the DBAdaptor. DBConnection should now be used OR preferably the object adaptor itself\n");
  $self->dbc->disconnect_when_inactive(@args);
}

sub reconnect_when_lost{
  my ($self, @args) = @_;

 deprecate("reconnect_when_lost Should no longer be called from the DBAdaptor. DBConnection should now be used OR preferably the object adaptor itself\n");
  $self->dbc->reconnect_when_lost(@args);
}


sub host{
  my ($self, @args) = @_;

 deprecate("host Should no longer be called from the DBAdaptor. DBConnection should now be used OR preferably the object adaptor itself\n");
  $self->dbc->host(@args);
}
sub username{
  my ($self, @args) = @_;

 deprecate("username Should no longer be called from the DBAdaptor. DBConnection should now be used OR preferably the object adaptor itself\n");
  $self->dbc->username(@args);
}
sub password{
  my ($self, @args) = @_;

 deprecate("password Should no longer be called from the DBAdaptor. DBConnection should now be used OR preferably the object adaptor itself\n");
  $self->dbc->password(@args);
}
sub driver{
  my ($self, @args) = @_;

 deprecate("driver Should no longer be called from the DBAdaptor. DBConnection should now be used OR preferably the object adaptor itself\n");
  $self->dbc->driver(@args);
}
sub port{
  my ($self, @args) = @_;

 deprecate("port Should no longer be called from the DBAdaptor. DBConnection should now be used OR preferably the object adaptor itself\n");
  $self->dbc->port(@args);
}

sub db_handle{
  my ($self, @args) = @_;


 deprecate("db_handle Should no longer be called from the DBAdaptor. DBConnection should now be used OR preferably the object adaptor itself\n");
  $self->dbc->db_handle(@args);
}


1;
