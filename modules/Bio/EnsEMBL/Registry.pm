#
# Ensembl module for Registry
#
# Copyright EMBL/EBI
##
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::Registry

=head1 SYNOPSIS

$gene_adaptor = Bio::EnsEMBL::Registry->get_adaptor("Homo Sapiens","core","Gene"))


=head1 DESCRIPTION

All Adaptors are stored/registered using this module. This module should then
be used to get the adaptors needed.

The registry can be loaded from a configuration file. If the enviroment
variable ENSEMBL_REGISTRY is set then the file pointed to by it is executed.
If not set then if the file ~/.ensembl_init exists then this will be executed.
An example of such configuration file can be found in 
ensembl/modules/Bio/EnsEMBL/Utils/ensembl_init.example
For the Web server ENSEMBL_REGISTRY should be set in SiteDefs.pm.

The four types of registrys are for db adaptors, dba adaptors, dna adaptors
and the standard type.

=head2 db

These are registrys for backwards compatibillity and enable the subroutines
to add other adaptors to connections.

e.g. get_all_db_adaptors, get_db_adaptor, add_db_adaptor, remove_db_adaptor
are the old DBAdaptor subroutines which are now redirected to the Registry.

So if before we had
   my $sfa = $self->adaptor()->db()->get_db_adaptor('blast');

We now want to change this to
   my $sfa =  Bio::EnsEMBL::Registry->get_db($self->adaptor->db,'blast');

OR preferably if the blast adaptor was set up in configure
   my $sfa = Bio::EnsEMBL::Registry->get_adaptor("Human","core","blast");


=head2 DBA

These are the stores for the DBAdaptors

The Registry will create all the DBConnections needed now if you set up the
configuration correctly. So instead of the old commands like

my $db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(....)
my $exon_adaptor = $db->get_ExonAdaptor;

we should now have just

my  $exon_adaptor = Bio::EnsEMBL::Registry->get_adaptor("Human","core","Exon");


=head2 DNA

This is an internal Registry and allows the configuration of a dnadb. 
An example here is to set the est database to get it's dna data from the core database.

## set the est db to use the core for getting dna data.
#Bio::EnsEMBL::Utils::ConfigRegistry->
#                         dnadb_add("Homo Sapiens","core","Homo Sapiens","est");


=head2 adaptors

This is the registry for all the general types of adaptors like GeneAdaptor, ExonAdaptor, 
Slice Adaptor etc.

These are accessed by the get_adaptor subroutine i.e.

my  $exon_adaptor = Bio::EnsEMBL::Registry->get_adaptor("Human","core","Exon");

=head1 CONTACT

Post questions to the Ensembl developer list: <ensembl-dev@ebi.ac.uk>


=head1 METHODS

=cut


package Bio::EnsEMBL::Registry;

use strict;

use Bio::EnsEMBL::DBSQL::MergedAdaptor;
use Bio::EnsEMBL::Utils::Exception qw( deprecate throw warning );
use Bio::EnsEMBL::Utils::Argument qw(rearrange);

use vars qw(%registry_register);

$registry_register{'_WARN'} = 0;



=head2 load_all
 Will load the registry with the configuration file which is obtained from
 the first in the following and in that order.

  1) if an argument is passed to this method this is used as the conf file.
  2) If the enviroment variable ENSEMBL_REGISTRY is set this is used.
  3) If the file .ensembl_init exist in the home directory it is used

  Arg [1]    : (optional) string $arg file to load the registry from
  Example    : Bio::EnsEMBL::Registry->load_all();
  Returntype : none
  Exceptions : none


=cut
 
sub load_all{
  my $class = shift;
  my $web_reg = shift;

  #$registry_register{'_WARN'} = 0; # default report overwriting
  if(!defined($registry_register{'seen'})){
    $registry_register{'seen'}=1;
    if(defined($web_reg)){
#      print STDERR  "Loading conf from site defs file ".$web_reg."\n";
      if(-e $web_reg){
	unless (my $return = do $web_reg ){
	  throw "Error in Configuration\n $!\n";
	}
	# other wise it gets done again by the web initialisation stuff
	delete $INC{$web_reg}; 
      }
    }
    elsif(defined($ENV{ENSEMBL_REGISTRY}) and -e $ENV{ENSEMBL_REGISTRY}){
#      print STDERR  "Loading conf from ".$ENV{ENSEMBL_REGISTRY}."\n";
      unless (my $return = do $ENV{ENSEMBL_REGISTRY}){
	throw "Error in Configuration\n $!\n";
      }
    }
    elsif(-e $ENV{HOME}."/.ensembl_init"){
      do($ENV{HOME}."/.ensembl_init");
    }
#    else{
#      print STDERR "NO default configuration to load\n";
#    }
  }
#  else{
#    print STDERR "Already configured???\n";
#  }
}

=head2 check_if_already_there

  Arg [DBNAME] : string
                 The name of the database to check for.
  Arg [HOST] : (optional) string
               The domain name of the database host to check for
  Arg [PORT] : int
               The port to check for when connecting to the database
  Arg [DRIVER] : (optional) string
                 The type of database driver to check for

  Description: Check to see if the database is already stored.
  Returntype : 0 if not found else the species and group.
  Exceptions : none
  

=cut

  
#sub check_if_already_there{
#  my ($class) = shift;
#
#  my ($dbname,$host,$driver,$port,$species, $group ) =
#    rearrange([qw(DBNAME HOST DRIVER PORT SPECIES GROUP)], @_);
#
#  if(defined($registry_register{'_DBA'})){
#    foreach my $db (@{$registry_register{'_DBA'}}){
#      my $dbc= $db->db();
#      if($dbc->host() eq $host and $dbc->dbname() eq $dbname
#	 and $dbc->driver() eq $driver and $dbc->port() eq $port){
#	if(defined($species) and defined($group) and 
#	   $db->species eq $species and $db->group eq $group){
#	  return ($db->species(),$db->group());
#	}
#      }
#    }
#  }
#  return 0;
#}


#
# add ons.
#

=head2 add_db

  Arg [1]    : db to add adaptor to.
  Arg [2]    : name of the name to add the adaptor to in the registry.
  Arg [3]    : The adaptor to be added to the registry.
  Example    : Bio::EnsEMBL::Registry->add_db($db, "lite", $dba);
  Returntype : none
  Exceptions : none

=cut

sub add_db{
  my ($class, $db, $name, $adap) = @_;


  if($db->species() ne $adap->species){
    $registry_register{$db->species()}{$db->group()}{'_special'}{$name} = $adap;
  }
}

=head2 remove_db

  Arg [1]    : db to remove adaptor from.
  Arg [2]    : name to remove the adaptor from in the registry.
  Example    : my $db = Bio::EnsEMBL::Registry->remove_db($db, "lite");
  Returntype : adaptor
  Exceptions : none

=cut

sub remove_db{
  my ($class, $db, $name) = @_;

  my $ret = $registry_register{$db->species()}{$db->group()}{'_special'}{$name};
  $registry_register{$db->species()}{$db->group()}{'_special'}{$name} = undef;

  return $ret;
}

=head2 get_db

  Arg [1]    : db to get adaptor from.
  Arg [2]    : name to get the adaptor for in the registry.
  Example    : my $db = Bio::EnsEMBL::Registry->get_db("Human", "core", "lite");
  Returntype : adaptor
  Exceptions : none

=cut

sub get_db{
  my ($class, $db, $name) = @_;

  my $ret = Bio::EnsEMBL::Registry->get_DBAdaptor($db->species,$name);

  if(defined($ret)){
    return $ret;
  }
  return $registry_register{$db->species()}{$db->group()}{'_special'}{$name};
}

=head2 get_all_db_adaptors

  Arg [1]    : db to get all the adaptor from.
  Example    : my $db = Bio::EnsEMBL::Registry->get_all_db_adaptors($db);
  Returntype : adaptor
  Exceptions : none

=cut

sub get_all_db_adaptors{
  my ($class,$db) = @_;
  my %ret=();

 foreach my $key (keys %{$registry_register{$class->get_alias($db->species())}{$db->group()}{'_special'}}){
   $ret{$key} = $registry_register{$class->get_alias($db->species())}{$db->group()}{'_special'}{$key};
 }

  return \%ret;
}


#
# DBAdaptors
#

=head2 add_DBAdaptor

  Arg [1]    : name of the species to add the adaptor to in the registry.
  Arg [2]    : name of the group to add the adaptor to in the registry.
  Arg [3]    : The DBAaptor to be added to the registry.
  Example    : Bio::EnsEMBL::Registry->add_DBAdaptor("Human", "core", $dba);
  Returntype : none
  Exceptions : none

=cut

sub add_DBAdaptor{
  my ($class, $species, $group, $adap) = @_;

  $species = $class->get_alias($species);

  $registry_register{$species}{$group}{'_DB'} = $adap;

  if(!defined($registry_register{'_DBA'})){
    my @list =();
    push(@list,$adap);
    $registry_register{'_DBA'}= \@list;
  }
  else{
    push(@{$registry_register{'_DBA'}},$adap);
  }

}



=head2 get_DBAdaptor

  Arg [1]    : name of the species to get the adaptor for in the registry.
  Arg [2]    : name of the group to get the adaptor for in the registry.
  Example    : $dba = Bio::EnsEMBL::Registry->get_DBAdaptor("Human", "core");
  Returntype : DBAdaptor
  Exceptions : none

=cut

sub get_DBAdaptor{
  my ($class, $species, $group) = @_;

  $species = $class->get_alias($species);

  if(defined($group)){ # group defined so return standard DB Adaptor
    return  $registry_register{$species}{$group}{'_DB'};
  }
  else{ #return a merged db adaptor
    return  new_merged Bio::EnsEMBL::DBSQL::DBAdaptor($species);
  }
}

=head2 get_all_DBAdaptors

  Example    : @dba = @{Bio::EnsEMBL::Registry->get_all_DBAdaptors();
  Returntype : list of DBAdaptors
  Exceptions : none

=cut

sub get_all_DBAdaptors{
  my ($class)=@_;

  return @{$registry_register{'_DBA'}};
}

#
# DNA Adaptors
#

=head2 add_DNAAdaptor

  Arg [1]    : name of the species to add the adaptor to in the registry.
  Arg [2]    : name of the group to add the adaptor to in the registry.
  Arg [3]    : The adaptor to be added to the registry as a DNA adaptor.
  Example    : Bio::EnsEMBL::Registry->add_DNAAdaptor("Human", "core", $dnaAdap);
  Returntype : none
  Exceptions : none

=cut

sub add_DNAAdaptor{
  my ($class, $species, $group, $adap) = @_;

  $species = $class->get_alias($species);

  $registry_register{$species}{$group}{'_DNA'} = $adap;

}

=head2 get_DNAAdaptor

  Arg [1]    : name of the species to get the adaptor for in the registry.
  Arg [2]    : name of the group to get the adaptor for in the registry.
  Example    : $dnaAdap = Bio::EnsEMBL::Registry->get_DNAAdaptor("Human", "core");
  Returntype : adaptor
  Exceptions : none

=cut

sub get_DNAAdaptor{
  my ($class, $species, $group) = @_;

  $species = $class->get_alias($species);
  return  $registry_register{$species}{$group}{'_DNA'};
}

#
# General Adaptors
#

=head2 add_adaptor

  Arg [1]    : name of the species to add the adaptor to in the registry.
  Arg [2]    : name of the group to add the adaptor to in the registry.
  Arg [3]    : name of the type to add the adaptor to in the registry.
  Arg [4]    : The DBAaptor to be added to the registry.
  Arg [5]    : (optional) if set okay to overwrite.
  Example    : Bio::EnsEMBL::Registry->add_adaptor("Human", "core", "Gene", $adap);
  Returntype : none
  Exceptions : none


=cut

sub add_adaptor{
  my ($class,$species,$group,$type,$adap, $reset)= @_;

  $species = $class->get_alias($species);

#
# Becouse the adaptors are not stored initially only there class paths when
# the adaptors are obtained we need to store these instead.
# It is not necessarily an error if the registry is overwritten without
# the reset set but it is an indication that we are overwriting a database
# which should be a warning for now
#

  if(defined($reset)){ # JUST REST THE HASH VLAUE NO MORE PROCESSING NEEDED
    $registry_register{$species}{$group}{$type} = $adap;
    return;
  }
  if(defined($registry_register{$species}{$group}{$type})){ 
    print STDERR ("Overwriting Adaptor in Registry for $species $group $type\n");
    $registry_register{$species}{$group}{$type} = $adap;
   return;
  }
  $registry_register{$species}{$group}{$type} = $adap;

  
  if(!defined ($registry_register{$species}{'list'})){
    my @list =();
    push(@list,$adap);
    $registry_register{$species}{'list'}= \@list;
  }
  else{
    push(@{$registry_register{$species}{'list'}},$adap);
  }

#  print STDERR "REGADD  $species \t $group \t $type\t to the registry\n";

  if(!defined ($registry_register{$type}{$species})){
    my @list =();
    push(@list,$adap);
    $registry_register{$type}{$species}= \@list;
  }
  else{
    push(@{$registry_register{$type}{$species}},$adap);
  }

}


=head2 set_get_via_dnadb_if_set

  set the flag so that for this type of adaptor the data is obtained
  from the dna source and not centrally i.e. estgenes where the sequence
  data is held in the core.

  Arg [1]    : name of the species to set flag for.
  Arg [2]    : name of the type to set flag for. (i.e. Sequence)
  Example    : Bio::EnsEMBL::Registry->set_get_via_dnadb_if_set("Human","Sequence");
  Returntype : none
  Exceptions : none
  
  

=cut

sub set_get_via_dnadb_if_set{
  my ($class,$species,$type) = @_;

  $registry_register{$class->get_alias($species)}{$type}{'DNADB'} = 1;
}

=head2 get_adaptor

  Arg [1]    : name of the species to add the adaptor to in the registry.
  Arg [2]    : name of the group to add the adaptor to in the registry.
  Arg [3]    : name of the type to add the adaptor to in the registry.
  Example    : $adap = Bio::EnsEMBL::Registry->get_adaptor("Human", "core", "Gene");
  Returntype : adaptor
  Exceptions : none

=cut

sub get_adaptor{
  my ($class,$species,$group,$type)= @_;

  $species = $class->get_alias($species);
 
  #throw in a check to see if we should get the dnadb one and not the normal
  if(defined($registry_register{$species}{$type}{'DNADB'}) && $class->get_DNAAdaptor($species,$group)){
    my $dna = $class->get_DNAAdaptor($species,$group);
    $species = $dna->species();
    $group = $dna->group();
  }

  my $ret = $registry_register{$species}{$group}{$type};
  if(!defined($ret)){
    return undef;
#    throw("COULD NOT FIND ADAPTOR species=$species\tgroup=$group\ttype=$type\n");
  }
  if(!ref($ret)){ # not instantiated yet
    my $dba = $registry_register{$species}{$group}{'_DB'};
    my $module = $ret;
    eval "require $module";

    if($@) {
      warning("$module cannot be found.\nException $@\n");
      return undef;
    }
    my $adap = "$module"->new($dba);
    Bio::EnsEMBL::Registry->add_adaptor($species, $group, $type, $adap, "reset");
    $ret = $adap;
  }

  return $ret;
}

=head2 get_all_adaptors

  Arg [1]    : name of the species to get the adaptors for.
  Example    : @adaps = @{Bio::EnsEMBL::Registry->get_all_adaptors()};
  Returntype : list of adaptors
  Exceptions : none

=cut

sub get_all_adaptors{
  my ($class,$species)= @_;

  $species = get_alias($species);
  return $registry_register{$species}{'list'};
}


=head2 get_MergedAdaptor

  Arg [1]    : name of the species to get the adaptors for in the registry.
  Arg [2]    : name of the type to get the adaptors for in the registry.
  Example    : $merged = Bio::EnsEMBL::Registry->get_MergedAdaptor("Mouse","Gene");
  Returntype : Bio::EnsEMBL::DBSQL::MergedAdaptor
  Exceptions : none

=cut

sub get_MergedAdaptor{
  my ($class,$species,$type)=@_;

  $species = $class->get_alias($species);
  my $ret = new Bio::EnsEMBL::DBSQL::MergedAdaptor();
  $ret->add_list(@{$registry_register{$type}{$species}});

  return $ret;
}

=head2 add_alias

  Arg [1]    : name of the species to add alias for
  Arg [2]    : name of the alias
  Example    : Bio::EnsEMBL::Registry->add_alias("Homo Sapiens","Human");
  Description: add alternative name for the species.
  Returntype : none
  Exceptions : none

=cut

sub add_alias{
  my ($class, $species,$key) = @_;

  $registry_register{'_ALIAS'}{$key} = $species;
}

=head2 get_alias

  Arg [1]    : name of the possible alias to get species for
  Arg [2]    : if set will not throw if not found.
  Example    : Bio::EnsEMBL::Registry->get_alias("Human");
  Description: get proper species name.
  Returntype : species name
  Exceptions : if not found and second argument

=cut

sub get_alias{
  my ($class, $key, $no_throw) = @_;

  if(!defined($registry_register{'_ALIAS'}{$key})){
    return $key;
  }
  return $registry_register{'_ALIAS'}{$key};
}

#
# Web specific routines
#


=head2 load_registry_with_web_adaptors
  Will load the registry with all the Adaptors used in the Web server.
  Providing Sitedefs and SpeciesDefs can be found on PERL5LIB path.

  Example    : Bio::EnsEMBL::Registry->load_registry_with_web_adaptors();
  Returntype : none
  Exceptions : Will die if Sitedefs or SpeciesDefs is not found on the
               PERL5LIB path.

=cut

sub load_registy_with_web_adaptors{
  my $class = shift;


  eval{ require SiteDefs };
  if ($@){ die "Can't use SiteDefs.pm - $@\n"; }
    SiteDefs->import(qw(:ALL));

  eval{ require SpeciesDefs };
  if ($@){ die "Can't use SpeciesDefs.pm - $@\n"; }
  my $conf = new SpeciesDefs();

}

=head2 set_default_track
  Sets a flag to say that that this species/group are a default track and do not
  need to be added as another web track.

  Arg [1]    : name of the species to get the adaptors for in the registry.
  Arg [2]    : name of the type to get the adaptors for in the registry.
  Example    : $merged = Bio::EnsEMBL::Registry->set_default_track("Human","core");
  Returntype : none
  Exceptions : none

=cut

sub set_default_track{
  my ($class, $species, $group) = @_;  

  $registry_register{'def_track'}{$species}{$group} = 1;
}

=head2 default_track
  Check flag to see if this is a default track

  Arg [1]    : name of the species to get the adaptors for in the registry.
  Arg [2]    : name of the type to get the adaptors for in the registry.
  Example    : $merged = Bio::EnsEMBL::Registry->set_default_track("Human","core");
  Returntype : int 
  Exceptions : none

=cut

sub default_track{
  my ($class, $species, $group) = @_;  

  if(defined($registry_register{'def_track'}{$species}{$group})){
    return 1;
  }
  
  return 0;
}


=head2 add_new_tracks
  Will add new gene tracks to the configuration of the WEB server if they are
  not of the type default and the configuration already has genes in the display.

  Arg [1]    : hash of the default configuration of the web page
  Returntype : none
  Exceptions : none
  Called by  : UserConfig.pm
  
=cut

sub add_new_tracks{
  my($class, $conf) = @_;

  my $start = 0;
  my $reg = $class;
  my $species_reg = $reg->get_alias($conf->{'species'},"nothrow");
  my $view =  $conf->{'type'};
  if(defined($species_reg)){
    my $config = $conf->{'general'}->{$view};
    foreach my $dba ($reg->get_all_DBAdaptors()){
      #      if($dba->species eq $species_reg and !$reg->default_track($dba->species,$dba->group)){
      if(!$reg->default_track($dba->species,$dba->group)){
	if($start == 0){
	  if(exists($config->{'vega_transcript_lite'}) and defined($config->{'vega_transcript_lite'}->{'pos'})){
	    $start = $config->{'vega_transcript_lite'}->{'pos'};
	  }
	  elsif(exists($config->{'transcript_lite'}) and defined($config->{'transcript_lite'}->{'pos'})){
	    $start = $config->{'transcript_lite'}->{'pos'};
	  }
	  else{ # no transcripts on this view so do not add track here
	    #    print STDERR "no transcript options on this display \n";
	    next;
	  }
	  $start ++;
	}
	else{
	  if(exists($config->{'vega_transcript_lite'}) and defined($config->{'vega_transcript_lite'}->{'pos'})){
	    $start++;
	  }
	  elsif(exists($config->{'transcript_lite'}) and defined($config->{'transcript_lite'}->{'pos'})){
	    $start++;
	  }
	  else{ # no transcripts on this view so do not add track here
	    #    print STDERR "no transcript options on this display \n";
	    next;
	  }
	}
	$reg->_add_new_track( $conf->{'general'}->{$view}, $dba, $start);
      }
    }
  }
}


=head2 _add_new_track
 
  Arg [1]    : hash of the configuration  
  Arg [2]    : dbadaptor to use to get track info from
  Arg [3]    : start index to place track in the right place
  Returntype : none
  Exceptions : none
  Called by  : add_new_tracks

=cut

sub _add_new_track{
  my ($class, $config, $dba, $start ) = @_;


  my $KEY = $dba->group();

  $config->{$KEY} ={
		    'on'    => "on",
		    'compact' => 'yes',
		    'pos'     => $start,
		    'str'     => 'b',
		    'src'     => 'all', # 'ens' or 'all'
		    'available'=> 'species '.$dba->species(),
		    'dba'  => $dba,		   
#		    'colours' => {$config->{'_colourmap'}->colourSet( 'vega_gene' )},
		    'glyphset' => 'generic_transcript',
		    'dep'      => 6,
		   };

  push @{ $config->{'_artefacts'} }, $KEY;
  push @{ $config->{'_settings'}->{'features'}}, [$KEY => $KEY] ;
  
# [ 'transcript_lite'      => "Ensembl Trans."  ],  
  return;
}


1;
