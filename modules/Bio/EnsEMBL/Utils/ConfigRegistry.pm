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

=head2 load_core, load_estgene, load_vega, load_compara, load_pipeline, load_SNP, load_lite
  Arg [1]    : DBAdaptor with DBConnection alredy attached
  Returntype : DBAdaptor;
  Exceptions : none

=cut

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
    $config_sub =  \&Bio::EnsEMBL::Utils::ConfigRegistry::load_variation;
  }
  elsif($dba->isa('Bio::EnsEMBL::DBSQL::DBAdaptor')){
    #vega uses the core DBAdaptor so test if vega is in the dbname
    if(!defined($dba->group())){
      $dba->group('core');
    }
    if($dba->group eq "estgene"){
      $config_sub =  \&Bio::EnsEMBL::Utils::ConfigRegistry::load_estgene;
     }
    elsif($dba->group eq "est"){
      $config_sub =  \&Bio::EnsEMBL::Utils::ConfigRegistry::load_est;
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
    $config_sub =  \&Bio::EnsEMBL::Utils::ConfigRegistry::load_and_attach_dnadb_to_core;
    #    throw("Unknown DBAdaptor type $dba\n");
  }


  # return if the connection and species, group are the same
  my $db_reg;
  if($reg->get_alias($dba->species,"no throw")){
    $db_reg = $reg->get_DBAdaptor($dba->species,$dba->group);
  }
  if(defined($db_reg)){
    if($dba->dbc->equals($db_reg->dbc)){
      return $db_reg;
    }
    else{ # diff conn details
      $dba->species(find_unique_species($dba->species,$dba->group));
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
  
  $reg->add_DNAAdaptor($dba->species,$dba->group,"core"); 
}

# those that do not need to attach to core:-#

sub load_core{
  load_adaptors(@_);
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
  load_adaptors(@_);
}


sub load_haplotype{
  load_adaptors(@_);
}


# these that need to attach to the core to get the sequense data

sub load_estgene{
  load_and_attach_dnadb_to_core(@_);
}

sub load_est{
  load_and_attach_dnadb_to_core(@_);
}

sub load_vega{
  load_and_attach_dnadb_to_core(@_);
}



#special cases;-




sub load_blast{
    my ($dba) = @_;
}


sub add_blast_link{
  my ($class, $species, $group) = @_;
  my $dba =undef;
  my $blast = undef;
  print STDERR " add_blast_link called\n";

  if($dba = Bio::EnsEMBL::Registry->get_DBAdaptor($species, $group)){
  }
  else{
    throw("Cannot find $species $group in the registry\n");
  }

  if($blast = Bio::EnsEMBL::Registry->get_adaptor("NONE", "blast", "Blast")){
    $dba->add_db_adaptor("Blast",$blast);
  }
  else{
    throw("Sorry no Blast database has been set up to link to.\n");
  }
}


sub load_go{

 my ($class, $dba) = @_;

 # shouldnt go into the registry ... sorry
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
#if ($@){ print STDERR  "No user defined loads\n"; }

1;
