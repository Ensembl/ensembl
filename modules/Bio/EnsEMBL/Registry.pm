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

These are the old DBAdaptors but since the DBAdaptor module will be deleted in a
future version these are here only to aid backwards compatibility.

The Registry will create all the DBConnections needed now if you set up the
configuration correctly. So instead of the old commands like

my $db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(....)
my $exon_adaptor = $db->get_ExonAdaptor;

we should now have just

my  $exon_adaptor = Bio::EnsEMBL::Registry->get_adaptor("Human","core","Exon");


=head2 DNA

This is an internal Registry and allows the configuration of a dnadb. An example here is to
set the est database to get it's dna data from the core database.

## set the est db to use the core for getting dna data.
#Bio::EnsEMBL::Utils::ConfigRegistry->dnadb_add("Homo Sapiens","core","Homo Sapiens","est");

A better way though is to use the Merged adaptors to get the features on the core for all databases.


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
use Bio::EnsEMBL::Utils::ConfigRegistry;
#use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::DBSQL::MergedAdaptor;
#use Bio::EnsEMBL::DBSQL::DBMergedAdaptor;
use Bio::EnsEMBL::Utils::Exception qw( deprecate throw warning );
use Bio::EnsEMBL::Utils::Argument qw(rearrange);


our %register = ();


$register{'_WARN'} = 0; # default report overwriting

if(defined($ENV{ENSEMBL_REGISTRY}) and -e $ENV{ENSEMBL_REGISTRY}){
#  print "Loading conf from ".$ENV{ENSEMBL_REGISTRY}."\n";
  unless (my $return = do $ENV{ENSEMBL_REGISTRY}){
    throw "Error in Configuration\n $!\n";
  }
}
elsif(-e $ENV{HOME}."/.ensembl_init"){
#  print "Loading conf from ".$ENV{HOME}."/.ensembl_init\n";
  do($ENV{HOME}."/.ensembl_init");
}
else{
#  print "NO default configuration to load\n";
}

=head2 warn_on_duplicates

  Arg [1]    : (optional) string $arg
  Example    : Bio::EnsEMBL::Registry->warn_on_duplicates(1);
  Description: Getter / Setter for the production of error 
               messages on overwriting values.
  Returntype : string
  Exceptions : none

=cut

sub warn_on_duplicates{
  my ($class) = shift;

  $register{'_WARN'} = shift if(@_);
  return $register{'_WARN'};
}

sub check_if_already_there{
  my ($class) = shift;

  my ($dbname,$host,$driver,$port ) =
    rearrange([qw(DBNAME HOST DRIVER PORT )], @_);

  if(defined($register{'_DBA'})){
    foreach my $db (@{$register{'_DBA'}}){
      my $dbc= $db->db();
      if($dbc->host() eq $host and $dbc->dbname() eq $dbname
	 and $dbc->driver() eq $driver and $dbc->port() eq $port){
	return ($db->species(),$db->group());
      }
    }
  }
  return 0;
}

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

#  print STDERR "ADDING# ".$adap->dbname." to ".$db->db->dbname."  as $name\n";
  $register{$db->species()}{$db->group()}{'_special'}{$name} = $adap;

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

  my $ret = $register{$db->species()}{$db->group()}{'_special'}{$name};
  $register{$db->species()}{$db->group()}{'_special'}{$name} = undef;

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

  return $register{$db->species()}{$db->group()}{'_special'}{$name};
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

 foreach my $key (keys %{$register{$db->species()}{$db->group()}{'_special'}}){
   $ret{$key} = $register{$db->species()}{$db->group()}{'_special'}{$key};
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

  $species = Bio::EnsEMBL::Utils::ConfigRegistry->get_alias($species);
  if(defined($register{$species}{$group}{'_DB'}) && warn_on_duplicates()){
    warning("Overwriting DBAdaptor in Registry for $species $group $adap\n");
  }

  $register{$species}{$group}{'_DB'} = $adap;

  if(!defined($register{'_DBA'})){
    my @list =();
    push(@list,$adap);
    $register{'_DBA'}= \@list;
  }
  else{
    push(@{$register{'_DBA'}},$adap);
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

  $species = Bio::EnsEMBL::Utils::ConfigRegistry->get_alias($species);

  if(defined($group)){ # group defined so return standard DB Adaptor
    return  $register{$species}{$group}{'_DB'};
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

  return @{$register{'_DBA'}};
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

  $species = Bio::EnsEMBL::Utils::ConfigRegistry->get_alias($species);
  if(defined($register{$species}{$group}{'_DNA'}) && warn_on_duplicates()){
    warning("Overwriting DNAAdaptor in Registry for $species $group $adap\n");
  }

  $register{$species}{$group}{'_DNA'} = $adap;

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

  $species = Bio::EnsEMBL::Utils::ConfigRegistry->get_alias($species);
  return  $register{$species}{$group}{'_DNA'};
}

#
# General Adaptors
#

=head2 add_adaptor

  Arg [1]    : name of the species to add the adaptor to in the registry.
  Arg [2]    : name of the group to add the adaptor to in the registry.
  Arg [3]    : name of the type to add the adaptor to in the registry.
  Arg [4]    : The DBAaptor to be added to the registry.
  Example    : Bio::EnsEMBL::Registry->add_adaptor("Human", "core", "Gene", $adap);
  Returntype : none
  Exceptions : none

=cut

sub add_adaptor{
  my ($class,$species,$group,$type,$adap)= @_;

  $species = Bio::EnsEMBL::Utils::ConfigRegistry->get_alias($species);

  if(defined($register{$species}{$group}{$type}) && warn_on_duplicates()){
    warning("Overwriting Adaptor in Registry for $species $group $type\n");
  }
  $register{$species}{$group}{$type} = $adap;


  if(!defined ($register{$species}{'list'})){
    my @list =();
    push(@list,$adap);
    $register{$species}{'list'}= \@list;
  }
  else{
    push(@{$register{$species}{'list'}},$adap);
  }

  if(!defined ($register{$type}{$species})){
    my @list =();
    push(@list,$adap);
    $register{$type}{$species}= \@list;
  }
  else{
    push(@{$register{$type}{$species}},$adap);
  }

}


sub set_get_via_dnadb_if_set{
  my ($class,$species,$type) = @_;

  $register{$species}{$type}{'DNADB'} = 1;
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

  $species = Bio::EnsEMBL::Utils::ConfigRegistry->get_alias($species);
 
  #throw in a check to see if we should get the dnadb one and not the normal
  if(defined($register{$species}{$type}{'DNADB'}) && $class->get_DNAAdaptor($species,$group)){
    my $dna = $class->get_DNAAdaptor($species,$group);
    $species = $dna->species();
    $group = $dna->group();
  }

  my $ret = $register{$species}{$group}{$type};
  if(!defined($ret)){
    throw("COULD NOT FIND ADAPTOR species=$species\tgroup=$group\ttype=$type\n");
    print STDERR caller();
    print STDERR "\nfin\n";;
  }
  return $ret;
}

=head2 get_all_adaptors

  Arg [1]    : name of the species to get the adaptors for.
  Example    : @adaps = @{Bio::EnsEMBL::Registry->get_all_adaptor()};
  Returntype : list of adaptors
  Exceptions : none

=cut

sub get_all_adaptors{
  my ($class,$species)= @_;

  $species = get_alias($species);
  return $register{$species}{'list'};
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

  $species = Bio::EnsEMBL::Utils::ConfigRegistry->get_alias($species);
  my $ret = new Bio::EnsEMBL::DBSQL::MergedAdaptor();
  $ret->add_list(@{$register{$type}{$species}});

  return $ret;
}

sub add_alias{
  my ($class, $species,$key) = @_;

  $register{'_ALIAS'}{$key} = $species;
}

sub get_alias{
  my ($class, $key) = @_;

  if(!defined($register{'_ALIAS'}{$key})){
    throw("Unknown species $key has it been mistyped??\n");
  }
  return $register{'_ALIAS'}{$key};
}

1;
