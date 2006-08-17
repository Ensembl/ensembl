#
# Ensembl module for Registry
#
# Copyright EMBL/EBI
##
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::Utils::ConfigRegistry;

=head1 SYNOPSIS


Bio::EnsEMBL::Utils::ConfigRegistry->load_core($dba );
 

=head1 DESCRIPTION

The ConfigRegistry will "Register" a set of adaptors based on the type of
database that is being loaded.

=head1 CONTACT

Post questions to the Ensembl developer list: <ensembl-dev@ebi.ac.uk>


=head1 METHODS

=cut


use strict;
use warnings;

package Bio::EnsEMBL::Utils::ConfigRegistry;

use Bio::EnsEMBL::Registry;
my $reg = "Bio::EnsEMBL::Registry";
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::DBSQL::DBConnection;
use Bio::EnsEMBL::DBSQL::DBAdaptor;

use Bio::EnsEMBL::Utils::Exception qw(warning throw  deprecate stack_trace_dump);

=head2 load_core, load_otherfeatures, load_vega, load_compara, load_pipeline, load_SNP, load_lite
  Arg [1]    : DBAdaptor with DBConnection alredy attached
  Returntype : DBAdaptor;
  Exceptions : none

=cut

#
# 1) core. no need to add dnadb
# 2) not core add dnadb
# 3) 
#


sub gen_load{
  my ($dba) = @_;
  my $config_sub;

# at some point we hope to set the group in the DBadaptor
# hence this long check etc should be simpler

  if($dba->isa('Bio::EnsEMBL::Compara::DBSQL::DBAdaptor')){
    if(!defined($dba->group())){
      $dba->group('compara');
    }
    $config_sub =  \&Bio::EnsEMBL::Utils::ConfigRegistry::load_compara; 
  }
  elsif($dba->isa('Bio::EnsEMBL::Lite::DBAdaptor')){
    if(!defined($dba->group())){
      $dba->group('lite');
    }
    $config_sub =  \&Bio::EnsEMBL::Utils::ConfigRegistry::load_lite;
  }
  elsif($dba->isa('Bio::EnsEMBL::External::BlastAdaptor')){
    if(!defined($dba->group())){
      $dba->group('blast');
    }
    $config_sub =  \&Bio::EnsEMBL::Utils::ConfigRegistry::load_blast;
  }
  elsif($dba->isa('Bio::EnsEMBL::ExternalData::SNPSQL::DBAdaptor')){
    if(!defined($dba->group())){
      $dba->group('SNP');
    } 
    $config_sub =  \&Bio::EnsEMBL::Utils::ConfigRegistry::load_SNP;
  }
  elsif($dba->isa('Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor')){
    if(!defined($dba->group())){
      $dba->group('pipeline');
    }
    $config_sub =  \&Bio::EnsEMBL::Utils::ConfigRegistry::load_pipeline;
  }
  elsif($dba->isa('Bio::EnsEMBL::Hive::DBSQL::DBAdaptor')){
    if(!defined($dba->group())){
      $dba->group('hive');
    }
    $config_sub =  \&Bio::EnsEMBL::Utils::ConfigRegistry::load_hive;
  }
  elsif($dba->isa('Bio::EnsEMBL::ExternalData::Haplotype::DBAdaptor')){
    if(!defined($dba->group())){
      $dba->group('haplotype');
    }
    $config_sub =  \&Bio::EnsEMBL::Utils::ConfigRegistry::load_haplotype;    
  }
  elsif($dba->isa('Bio::EnsEMBL::Variation::DBSQL::DBAdaptor')){
    if(!defined($dba->group())){
      $dba->group('variation');
    }
    $config_sub =  \&Bio::EnsEMBL::Utils::ConfigRegistry::load_variation;
  }
  elsif($dba->isa('Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor')){
    if(!defined($dba->group())){
      $dba->group('funcgen');
  }
    $config_sub =  \&Bio::EnsEMBL::Utils::ConfigRegistry::load_funcgen;
  }
  elsif($dba->isa('Bio::EnsEMBL::DBSQL::DBAdaptor')){
    #vega uses the core DBAdaptor so test if vega is in the dbname
    if(!defined($dba->group())){
      $dba->group('core');
    }
    if($dba->group eq "estgene"){
      $config_sub =  \&Bio::EnsEMBL::Utils::ConfigRegistry::load_estgene;
     }
    elsif($dba->group eq "otherfeatures"){
      $config_sub =  \&Bio::EnsEMBL::Utils::ConfigRegistry::load_otherfeatures;
     }
    elsif($dba->group eq "vega"){
      $config_sub =  \&Bio::EnsEMBL::Utils::ConfigRegistry::load_vega;
     }
    else{
      $config_sub =  \&Bio::EnsEMBL::Utils::ConfigRegistry::load_core;
    }
  }
  else{
    # none standard DBA adaptor 
    if(!defined($dba->group())){
      $dba->group('none_standard');
    }
    $config_sub =  \&Bio::EnsEMBL::Utils::ConfigRegistry::load_and_attach_dnadb_to_core;
    #    throw("Unknown DBAdaptor type $dba\n");
  }


  # return if the connection and species, group are the same


  if(defined($dba->species)){
    my $db_reg = $reg->get_DBAdaptor($dba->species,$dba->group);
    if(defined($db_reg)){
      if($dba->dbc->equals($db_reg->dbc)){
	return $db_reg;
      }
      else{
	warn "WARN: Species and group same for two seperate databases\nModify species name for one of these\n";
	$dba->species(find_unique_species($dba->species,$dba->group));
      }
    }
  }
  else{  # no species 
    my @db_reg = @{$reg->get_all_DBAdaptors_by_connection($dba->dbc)};
    foreach my $db_adaptor (@db_reg){
      if($db_adaptor->group eq $dba->group){ # found same db connection and group
	return $db_adaptor;
      }
    }
    $dba->species(find_unique_species("DEFAULT",$dba->group));      
    if($dba->species ne "DEFAULT"){
      warn "WARN: For multiple species use species attribute in DBAdaptor->new\n" 
    }
  }


  Bio::EnsEMBL::Registry->add_DBAdaptor($dba->species(), $dba->group(), $dba);

  #call the loading subroutine. (add the adaptors to the DBAdaptor)
  &{$config_sub}($dba);
   
  return $dba;
}



sub find_unique_species{
  my ($species, $group) = @_;

  $reg->add_alias($species,$species);

  my $i = 0;
  my $free =0;
  while(!$free){
    if($i == 0){
      if(!defined($reg->get_DBAdaptor($species, $group))){
	$free =1;
	$i ="";
      }
      else{
	$i = 1;
      }
    }
    else{
      $reg->add_alias($species.$i,$species.$i); #set needed self alias
      if(!defined($reg->get_DBAdaptor($species.$i, $group))){
	$free =1;
      }
      else{
	$i++;
      }
    }
  }
  
  $species .= $i;
  return ($species);
}



sub load_adaptors{
  my ($dba) = @_;

  my %pairs = %{$dba->get_available_adaptors()};
  
  foreach my $key (keys %pairs){
    Bio::EnsEMBL::Registry->add_adaptor($dba->species, $dba->group, $key, $pairs{$key});
  }

}

sub load_and_attach_dnadb_to_core{
  my ($dba) = @_;

  load_adaptors($dba);
  
  $reg->add_DNAAdaptor($dba->species,$dba->group,$dba->species,"core"); 
}


sub load_core{
  my ($dba) = @_;

  load_adaptors($dba);
}

sub load_compara{
  load_adaptors(@_);
}

sub load_hive{
  load_adaptors(@_);
}

sub load_pipeline{
  load_adaptors(@_);
}

sub load_lite{
  load_adaptors(@_);
}

sub load_SNP{
  load_adaptors(@_);
}

sub load_variation{
  load_and_attach_dnadb_to_core(@_);
}

sub load_funcgen{
  load_and_attach_dnadb_to_core(@_);
}

sub load_haplotype{
  load_adaptors(@_);
}


# these that need to attach to the core to get the sequense data

sub load_estgene{
  load_and_attach_dnadb_to_core(@_);
}

sub load_otherfeatures{
  load_and_attach_dnadb_to_core(@_);
}

sub load_vega{
  load_and_attach_dnadb_to_core(@_);
}



sub add_alias{
  my ($class, @args) = @_;
  my ($species, $aliases) = rearrange([qw(SPECIES ALIAS)],@args);

 #make sure it exists itself
  Bio::EnsEMBL::Registry->add_alias($species,$species);

  if($aliases){
    foreach my $ali (@$aliases){
      Bio::EnsEMBL::Registry->add_alias($species,$ali);
    }
  }

}
#
# overwrite/load new types. Done this way to enable no changes to CVS for
# external users. External users should add there own "GROUPS" in the file
# User_defined_load.
#

eval{ require Bio::EnsEMBL::Utils::User_defined_load };

1;
