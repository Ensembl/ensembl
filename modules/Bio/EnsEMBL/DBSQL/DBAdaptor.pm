
=head1 NAME - Bio::EnsEMBL::DBSQL::DBAdaptor

=head1 SYNOPSIS

    $db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
        -user   => 'root',
        -dbname => 'pog',
        -host   => 'caldy',
        -driver => 'mysql'
        );

    $gene_adaptor = $db->get_GeneAdaptor();

    $gene = $gene_adaptor()->fetch_by_stable_id($stable_id);

    $slice = $db->get_SliceAdaptor()->fetch_by_chr_start_end('X', 1, 10000);

=head1 DESCRIPTION

This is the primary interface to an EnsEMBL database. It maintains an active
connection to the database and allows for the retrieval of ObjectAdaptors,
via a set of get_XxxxAdaptor methods (where Xxxx is the type of adaptor).

ObjectAdaptors can then be used to obtain objects and actual information
from the database.


=head1 CONTACT

Post questions to the EnsEMBL development list <ensembl-dev@ebi.ac.uk>

=head1 METHODS

The rest of the documentation details each of the object methods. Internal
methods are usually preceded with a _

=cut

package Bio::EnsEMBL::DBSQL::DBAdaptor;

use Bio::EnsEMBL::Registry;
my $reg = "Bio::EnsEMBL::Registry";

use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::DBSQL::DBConnection;
use Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor;
use Bio::EnsEMBL::Utils::SeqRegionCache;
use Bio::EnsEMBL::Utils::Exception qw(throw warning deprecate);
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Utils::ConfigRegistry;

=head2 new

  Arg [-DNADB]: (optional) Bio::EnsEMBL::DBSQL::DBAdaptor DNADB 
               All sequence, assembly, contig information etc, will be
               retrieved from this database instead.
  Arg [..]   : Other args are passed to superclass
               Bio::EnsEMBL::DBSQL::DBConnection
  Example    : $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
						    -user   => 'root',
						    -dbname => 'pog',
						    -host   => 'caldy',
						    -driver => 'mysql' );
  Exmaple2   : $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
                                                    -species => 'Homo_sapiens',
                                                    -group   => 'core');
  Description: Constructor for DBAdaptor.
  Returntype : Bio::EnsEMBL::DBSQL::DBAdaptor
  Exceptions : none
  Caller     : general

=cut

sub new {
  my($class, @args) = @_;

  #call superclass constructor
  my $self ={};
  bless $self,$class;

  my ($species, $group, $dbname) =
    rearrange([qw(SPECIES GROUP DBNAME)], @args);

  my ($spec,$gro) = $reg->check_if_already_there(@args);
  if($spec){
#    print STDERR "FOUND $spec $group for $dbname\n";
    $self = $reg->get_DBAdaptor($spec,$gro);
    $self->species($species);
    $self->group($group);
    return $self;
  }
  my $config_sub;
  if(defined($species)){ # NEW style usage of registry
    $self = $reg->get_DBAdaptor($species,$group);
    $self->species($species);
    $self->group($group);
  }
  else{                  # OLD style, so mimic by adding to registry
    my $free=0;
    $species= "DEFAULT";
    if($class->isa('Bio::EnsEMBL::Compara::DBSQL::DBAdaptor')){
      $group = "compara";
      $config_sub =  \&Bio::EnsEMBL::Utils::ConfigRegistry::load_compara;
    }
    elsif($class->isa('Bio::EnsEMBL::Lite::DBAdaptor')){
      $group = 'lite';
      $config_sub =  \&Bio::EnsEMBL::Utils::ConfigRegistry::load_lite;
    }
    elsif($class->isa('Bio::EnsEMBL::External::BlastAdaptor')){
      $group = 'blast';
      $config_sub =  \&Bio::EnsEMBL::Utils::ConfigRegistry::load_blast;
    }
    elsif($class->isa('Bio::EnsEMBL::ExternalData::SNPSQL::DBAdaptor')){
      $group = "SNP";
      $config_sub =  \&Bio::EnsEMBL::Utils::ConfigRegistry::load_SNP;
    }
    elsif($class->isa('Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor')){
      $group = "pipeline";
      $config_sub =  \&Bio::EnsEMBL::Utils::ConfigRegistry::load_pipeline;
    }
    elsif($class->isa('Bio::EnsEMBL::Hive::DBSQL::DBAdaptor')){
      $group = "hive";
      $config_sub =  \&Bio::EnsEMBL::Utils::ConfigRegistry::load_hive;
    }
    elsif($class->isa('Bio::EnsEMBL::DBSQL::DBAdaptor')){
      $group = "core";
      $config_sub =  \&Bio::EnsEMBL::Utils::ConfigRegistry::load_core;
    }
    else{
      throw("Unknown DBAdaptor type $class\n");
    }
    $reg->add_alias($species,$species); #set needed self alias

    my $i = 1;
    while(!$free){
      $reg->add_alias($species.$i,$species.$i); #set needed self alias
      if(!defined($reg->get_DBAdaptor($species.$i, $group))){
	$free =1;
      }
      else{
	$i++;
      }
    }
    $species .= $i;
    push (@args, '-species');
    push (@args, $species);


    &{$config_sub}($class,@args);

    $self = $reg->get_DBAdaptor($species,$group);

  }
  my ( $dnadb ) = rearrange([qw(DNADB)],@args);

  if(defined $dnadb) {
    $self->dnadb($dnadb);
  }
  return $self;
}

=head2 new_fast

  Arg [-CON]: Bio::EnsEMBL::DBSQL::DBConnection

  Exmaple    : $db = new Bio::EnsEMBL::DBSQL::DBAdaptor( -con => $dbc);
  Description: Constructor for DBAdaptor.
  Returntype : Bio::EnsEMBL::DBSQL::DBAdaptor
  Exceptions : none
  Caller     : general

=cut

sub new_fast{
  my($class, @args) = @_;
  my ( $con ) = rearrange([qw(CON)],@args);

  #call superclass constructor
  my $self ={};
  bless $self,$class;
  unless($con && ref $con &&
	 $con->isa('Bio::EnsEMBL::DBSQL::DBConnection')) {
    throw("$con passed is not of type Bio::EnsEMBL::DBSQL::DBConnection");
  }
  $self->db($con);
  $self->species($con->species());
  $self->group($con->group());

  return $self;
}


sub new_merged{
  my($class, $species, @args) = @_;

  my $self ={};
  bless $self,$class;

  $self->species($species);
  $self->group("_MERGED_");

  return $self;
}



=head2 db

  Arg[1]    : (optional) Bio::EnsEMBL::DBSQL::DBConnection

  Exmaple    : $dbc = $dba->db();
  Description: Getter/Setter for DBConnection.
  Returntype : Bio::EnsEMBL::DBSQL::DBConnection
  Exceptions : none
  Caller     : general

=cut

sub db{
  my ($self, $arg ) = @_;
  ( defined $arg ) &&
    ( $self->{_db} = $arg );
  $self->{_db};
}


sub DESTROY{
  my ($self)= @_;

  $self->{'_db'} = undef;
}

=head2 deleteObj

  Arg [1]    : none
  Example    : none
  Description: Cleans up circular reference loops so proper garbage collection
               can occur.
  Returntype : none
  Exceptions : none
  Caller     : DBAdaptorContainer::DESTROY

=cut

sub deleteObj {
  my $self = shift;
  #print "called deleteObj on DBAdaptor\n";

  #clean up external feature adaptor references
  if(exists $self->{'_xf_adaptors'}) {
    foreach my $key (keys %{$self->{'_xf_adaptors'}}) {
      delete $self->{'_xf_adaptors'}->{$key};
    }
  }

  if(exists $self->{'generic_feature_adaptors'}) {
    foreach my $name (keys %{$self->{'generic_feature_adaptors'}}) {
      my $adaptor = $self->{'generic_feature_adaptors'}->{$name};
      if(ref($adaptor) && $adaptor->can('deleteObj')) {
	$adaptor->deleteObj();
      }

      delete $self->{'generic_feature_adaptors'}->{$name};
    }

    delete $self->{'generic_feature_adaptors'};
  }

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

=cut

sub add_db_adaptor {
  my ($self, $name, $adaptor) = @_;

  unless($name && $adaptor && ref $adaptor) {
    throw('adaptor and name arguments are required');
  }

#  print STDERR "ADDING ".$adaptor->db->dbname." to ".$self->db->dbname."  as $name\n";
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

=cut

sub get_all_db_adaptors {
  my ($self) = @_;

  my %ret = %{Bio::EnsEMBL::Registry->get_all_db_adaptors($self)};

#  foreach my $key (%ret){
#    print $key."\n";
#  }
#  print %ret."\n";
  return \%ret;
#  unless(defined $self->{'_db_adaptors'}) {
#    return {};
#  }

#  return $self->{'_db_adaptors'};
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

=cut

sub get_db_adaptor {
  my ($self, $name) = @_;

#  print STDERR "looking for ".$self->db->dbname." and $name\n";
  return Bio::EnsEMBL::Registry->get_db($self, $name);
}

sub get_SNPAdaptor {
  my ($self)  = @_;
 
  my $lite = $self->get_db_adaptor('lite'); #### use register directly here

  my $primary_adaptor;

  if($lite) {
    $primary_adaptor = $lite->get_SNPAdaptor();
  } else {
    my $snp = $self->get_db_adaptor('SNP');
  
    unless($snp) {
      warn("No lite or SNP database, cannot get snp adaptor");
      return undef;
    }

    $primary_adaptor = $snp->get_SNPAdaptor();
    $primary_adaptor->ensembl_db( $self );
  }

  #return a proxy adaptor which can use the lite or the core database
  my $ret = $self->get_adaptor("ProxySNP");
  $ret->set_primary($primary_adaptor);
  return $ret;
#  return $self->get_adaptor("ProxySNP", $primary_adaptor);
}


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

=cut

sub add_DASFeatureFactory{
 
 my ($self,$value) = @_;
  
  push(@{$self->{'_das_ff'}},$value);
}


=head2 _each_DASFeatureFactory

  Args       : none
  Example    : none
  Description: Not sure if this is used, or if it should be removed.  It 
               does not seem to be used at the moment
  Returntype : Bio::EnsEMBL::ExternalFeatureFactory
  Exceptions : none
  Caller     : ??

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

=cut

sub get_adaptor {
  my ($self, $canonical_name, @other_args) = @_;

  return $reg->get_adaptor($self->db->species(),$self->db->group(),$canonical_name);
}



sub prepare{
  my ($self, @args) = @_;

 deprecate("prepare Should no longer be called from the DBAdaptor. DBConnection should now be used OR preferably the object adaptor itself\n");
  $self->db->prepare(@args);
}

sub dbname{
  my ($self, @args) = @_;

 deprecate("dbname Should no longer be called from the DBAdaptor. DBConnection should now be used OR preferably the object adaptor itself\n");
  $self->db->dbname(@args);
}

sub disconnect_when_inactive{
  my ($self, @args) = @_;

 deprecate("disconnect_when_inactive Should no longer be called from the DBAdaptor. DBConnection should now be used OR preferably the object adaptor itself\n");
  $self->db->disconnect_when_inactive(@args);
}

sub host{
  my ($self, @args) = @_;

 deprecate("host Should no longer be called from the DBAdaptor. DBConnection should now be used OR preferably the object adaptor itself\n");
  $self->db->host(@args);
}
sub username{
  my ($self, @args) = @_;

 deprecate("username Should no longer be called from the DBAdaptor. DBConnection should now be used OR preferably the object adaptor itself\n");
  $self->db->username(@args);
}
sub password{
  my ($self, @args) = @_;

 deprecate("password Should no longer be called from the DBAdaptor. DBConnection should now be used OR preferably the object adaptor itself\n");
  $self->db->password(@args);
}
sub driver{
  my ($self, @args) = @_;

 deprecate("driver Should no longer be called from the DBAdaptor. DBConnection should now be used OR preferably the object adaptor itself\n");
  $self->db->driver(@args);
}
sub port{
  my ($self, @args) = @_;

 deprecate("port Should no longer be called from the DBAdaptor. DBConnection should now be used OR preferably the object adaptor itself\n");
  $self->db->port(@args);
}

sub db_handle{
  my ($self, @args) = @_;


 deprecate("db_handle Should no longer be called from the DBAdaptor. DBConnection should now be used OR preferably the object adaptor itself\n");
  $self->db->db_handle(@args);
}

=head2 set_adaptor

  Arg [1]    : Canonical data type for new adaptor.
	Arg [2]    : Object defining the adaptor for arg1.
  Example    : $pa = Bio::EnsEMBL::DBSQL::ProteinAdaptor->new($db_adaptor);
             : $db_adaptor->set_adaptor("Protein", $pa)
  Description: Stores the object which represents the adaptor for the
               arg1 data type.
  Returntype : none
  Exceptions : none
  Caller     : external

=cut

sub set_adaptor {
  my ($self, $canonical_name, $module) = @_;

  my $adaptor = $self->_get_adaptor($module);

  $reg->add_adaptor($self->species(), $self->group(), $canonical_name, $adaptor);

  return $adaptor;
}

=head2 _get_adaptor

  Arg [1]    : string $module
               the fully qualified of the adaptor module to be retrieved
  Arg [2..n] : (optional) arbitrary list @args
               list of arguments to be passed to adaptors constructor
  Example    : $adaptor = $self->_get_adaptor("full::adaptor::name");
  Description: PROTECTED Used by subclasses to obtain adaptor objects
               for this database connection using the fully qualified
               module name of the adaptor. If the adaptor has not been 
               retrieved before it is created, otherwise it is retreived
               from the adaptor cache.
  Returntype : Adaptor Object of arbitrary type
  Exceptions : thrown if $module can not be instantiated
  Caller     : Bio::EnsEMBL::DBAdaptor

=cut

sub _get_adaptor {
  my( $self, $module) = @_;

  my( $adaptor);

  eval "require $module";

  if($@) {
    warning("$module cannot be found.\nException $@\n");
    return undef;
  }

  $adaptor = "$module"->new($self);

  return $adaptor;
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

=cut

sub species {
  my ($self, $arg ) = @_;
  ( defined $arg ) &&
    ( $self->{_species} = $arg );
  $self->{_species};
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


=head2 dnadb

 Title   : dnadb
 Usage   : my $dnadb = $db->dnadb();
 Function: returns the database adaptor where the dna lives
           Useful if you only want to keep one copy of the dna
           on disk but have other databases with genes and features in
 Returns : dna database adaptor
 Args    : Bio::EnsEMBL::DBSQL::BaseAdaptor

=cut

sub dnadb {
  my $self = shift;

  if(@_) {
    my $arg = shift;
#    print STDERR "ADDING DNADB ".$self->species()." ".$self->group()." to use ";
#    print STDERR $arg->db->species()." ".$arg->db->group()." instead for dna\n";
    $reg->add_DNAAdaptor($self->species(),$self->group(),$arg);
  }

#  return $self->{'dnadb'} || $self;
  return $reg->get_DNAAdaptor($self->species(),$self->group()) || $self;
}


use vars '$AUTOLOAD';

sub AUTOLOAD {
  my ($self,@args) = @_;

  my $type;
  if($AUTOLOAD =~ /^.*::get_(\w+)Adaptor$/){ 
    $type = $1;
  }
  elsif($AUTOLOAD =~ /^.*::get_(\w+)$/){ 
    $type = $1;
  }
  else{
    throw("Could not work out type for $AUTOLOAD \n");
  }
#  print STDERR "AUTO ".$self->species()."\t".$self->group()."\t$type\n";
  my $ret = undef;
  if($self->group() eq "_MERGED_"){
    $ret = $reg->get_MergedAdaptor($self->species(),$type);
  }
  else{
    $ret = $reg->get_adaptor($self->species(),$self->group(),$type);
  }
  if($ret){
    return $ret;
  }
  else{
    warning("Could not find $type adaptor in the registry for ".$self->species." ".$self->group."\n");
    return $ret;
  }
  die("No such method: $AUTOLOAD\n");
}

#########################
# sub DEPRECATED METHODS
#########################



=head2 get_MapFragAdaptor

  Description: MapFragAdaptor is deprecated. Use MiscFeatureAdaptor instead.

=cut

sub get_MapFragAdaptor {
  my $self = shift;

  return $self->get_adaptor( "MapFrag" );
}

=head2 get_CloneAdaptor

  Description: CloneAdaptor is deprecated.  Use SliceAdaptor instead.

=cut

sub get_CloneAdaptor {
  my( $self ) = @_;

  return $self->dnadb->get_adaptor("Clone");
}

=head2 get_ChromosomeAdaptor

  Description: ChromosomeAdaptor is deprecated.  Use SliceAdaptor instead.

=cut

sub get_ChromosomeAdaptor {
    my( $self ) = @_;

    return 
      $self->dnadb->get_adaptor("Chromosome");
}

=head2 get_RawContigAdaptor

  Description: RawContigAdaptor is deprecated. Use SliceAdaptor instead.

=cut

sub get_RawContigAdaptor {
    my( $self ) = @_;

    return $self->dnadb->get_adaptor("RawContig");
}


=head2 get_ProteinAdaptor

  Description: ProteinAdaptor is deprecated. Use TranslationAdaptor instead

=cut

sub get_ProteinAdaptor {
    my $self  = shift;
    deprecate("The ProteinAdaptor is deprecated. Use the TranslationAdaptor " .
              "instead of the ProteinAdaptor and Translation instead of " .
              "Protein.");
    return $self->get_adaptor("Protein");
}


sub source {
  deprecate('Do not use - this method does nothing');
}


=head2 assembly_type

  Description: DEPRECATED - Use CoordSystemAdaptor to obtain default coordinate
               system instead.

=cut

sub assembly_type{
  my $self = shift;

  deprecate('Use CoordSystemAdaptor::fetch_top_level instead');

  my $csa = $self->get_CoordSystemAdaptor();

  #get the default top-level coord system
  my ($dbID,$name,$version) = $csa->fetch_top_level();

  return $version;
}



=head2 list_supported_assemblies

  Description: DEPRECATED - Use CoordSystemAdaptor to obtain list of top-level
               coordinate systems instead

=cut

sub list_supported_assemblies {
  my($self) = @_;
  deprecate('Use CoordSystemAdaptor::fetch_all_top_level instead');

  my $csa = $self->get_CoordSystemAdaptor();
  my @versions;
  foreach my $cs (@{$csa->fetch_all_top_level()}) {
    my ($dbID, $name, $version) = @$csa;
    push @versions, $version;
  }
  return @versions;
}


1;
