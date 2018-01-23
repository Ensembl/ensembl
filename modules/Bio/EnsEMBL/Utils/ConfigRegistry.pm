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

Bio::EnsEMBL::Utils::ConfigRegistry;

=head1 SYNOPSIS


  Bio::EnsEMBL::Utils::ConfigRegistry->load_core($dba);


=head1 DESCRIPTION

The ConfigRegistry will "Register" a set of adaptors based on the type
of database that is being loaded.

=head1 METHODS

=cut

package Bio::EnsEMBL::Utils::ConfigRegistry;

use strict;
use warnings;

use Bio::EnsEMBL::Registry;
my $reg = "Bio::EnsEMBL::Registry";
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::DBSQL::DBConnection;
use Bio::EnsEMBL::DBSQL::DBAdaptor;

use Bio::EnsEMBL::Utils::Exception qw(warning throw  deprecate stack_trace_dump);



sub gen_load {
  my ($dba) = @_;
  my $pre_hook;
  my $config_sub;

  # At some point we hope to set the group in the DBadaptor, hence this
  # long check etc. should be simpler.

  if ( $dba->isa('Bio::EnsEMBL::Compara::DBSQL::DBAdaptor') ) {
    if ( !defined( $dba->group() ) ) {
      $dba->group('compara');
    }
    $config_sub = \&Bio::EnsEMBL::Utils::ConfigRegistry::load_compara;
  } elsif ( $dba->isa('Bio::EnsEMBL::Lite::DBAdaptor') ) {
    if ( !defined( $dba->group() ) ) {
      $dba->group('lite');
    }
    $config_sub = \&Bio::EnsEMBL::Utils::ConfigRegistry::load_lite;
  } elsif ( $dba->isa('Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor') ) {
    if ( !defined( $dba->group() ) ) {
      $dba->group('pipeline');
    }
    $config_sub = \&Bio::EnsEMBL::Utils::ConfigRegistry::load_pipeline;
  } elsif ( $dba->isa('Bio::EnsEMBL::Hive::DBSQL::DBAdaptor') ) {
    if ( !defined( $dba->group() ) ) {
      $dba->group('hive');
    }
    $config_sub = \&Bio::EnsEMBL::Utils::ConfigRegistry::load_hive;
  } elsif ( $dba->isa('Bio::EnsEMBL::Variation::DBSQL::DBAdaptor') ) {
    if ( !defined( $dba->group() ) ) {
      $dba->group('variation');
    }
    $config_sub = \&Bio::EnsEMBL::Utils::ConfigRegistry::load_variation;
  } elsif ( $dba->isa('Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor') ) {
    if ( !defined( $dba->group() ) ) {
      $dba->group('funcgen');
    }
    $pre_hook = \&Bio::EnsEMBL::Utils::ConfigRegistry::pre_funcgen_hook;
    $config_sub = \&Bio::EnsEMBL::Utils::ConfigRegistry::load_funcgen;
  } elsif ( $dba->isa('Bio::EnsEMBL::DBSQL::OntologyDBAdaptor') || $dba->isa('Bio::Ensembl::DBSQL::OntologyTermAdaptor') ) {
    if ( !defined( $dba->group() ) ) {
      $dba->group('ontology');
    }
    $config_sub = \&Bio::EnsEMBL::Utils::ConfigRegistry::load_ontology;
  } elsif ( $dba->isa('Bio::EnsEMBL::DBSQL::DBAdaptor') ) {
    #vega uses the core DBAdaptor so test if vega is in the dbname
    if ( !defined( $dba->group() ) ) {
      $dba->group('core');
    }

    if ( $dba->group eq "estgene" ) {
      $config_sub = \&Bio::EnsEMBL::Utils::ConfigRegistry::load_estgene;
    } elsif ( $dba->group eq "otherfeatures" ) {
      $config_sub =
        \&Bio::EnsEMBL::Utils::ConfigRegistry::load_otherfeatures;
    } elsif ( $dba->group eq "rnaseq" ) {
      $config_sub =
        \&Bio::EnsEMBL::Utils::ConfigRegistry::load_rnaseq;
    } elsif ( $dba->group eq 'vega' || $dba->group eq 'vega_update' ) {
      $config_sub = \&Bio::EnsEMBL::Utils::ConfigRegistry::load_vega;
    } else {
      $config_sub = \&Bio::EnsEMBL::Utils::ConfigRegistry::load_core;
    }

  } else {
    # none standard DBA adaptor
    if ( !defined( $dba->group() ) ) {
      $dba->group('none_standard');
    }
    $config_sub =
      \&Bio::EnsEMBL::Utils::ConfigRegistry::load_and_attach_dnadb_to_core;
    #    throw("Unknown DBAdaptor type $dba\n");
  }
  
  #Run the pre-hook if one was defined
  $pre_hook->($dba) if $pre_hook;

  # return if the connection and species, group are the same

  if ( defined( $dba->species ) ) {
    my $db_reg = $reg->get_DBAdaptor( $dba->species, $dba->group, 1 );
    if ( defined($db_reg) ) {
      if ( $dba->dbc->equals( $db_reg->dbc ) ) { return $db_reg }
      else {
        my $msg =
          sprintf( 'WARN: Species (%s) and group (%s) '
            . 'same for two seperate databases',
          $dba->species(), $dba->group() );

        warn "${msg}\nModify species name for one of these\n";
        $dba->species(
          find_unique_species( $dba->species, $dba->group ) );
      }
    }
  } else {    # no species
    
    my @db_reg =
      @{ $reg->get_all_DBAdaptors_by_connection( $dba->dbc ) };

    foreach my $db_adaptor (@db_reg) {
      if ( $db_adaptor->group eq $dba->group ) {
        # found same db connection and group
        return $db_adaptor;
      }
    }

    $dba->species( find_unique_species( "DEFAULT", $dba->group ) );
    if ( $dba->species ne "DEFAULT" ) {
      warn "WARN: For multiple species "
        . "use species attribute in DBAdaptor->new()\n";
    }
  }

  Bio::EnsEMBL::Registry->add_DBAdaptor( $dba->species(), $dba->group(),
    $dba );

  #call the loading subroutine. (add the adaptors to the DBAdaptor)
  &{$config_sub}($dba);

  return $dba;
} ## end sub gen_load



sub find_unique_species {
  my ( $species, $group ) = @_;

  $reg->add_alias( $species, $species );

  my $i    = 0;
  my $free = 0;

  while ( !$free ) {
    if ( $i == 0 ) {
      if ( !defined( $reg->get_DBAdaptor( $species, $group ) ) ) {
        $free = 1;
        $i    = "";
      } else {
        $i = 1;
      }
    } else {
      # set needed self alias
      $reg->add_alias( $species . $i, $species . $i );

      if ( !defined( $reg->get_DBAdaptor( $species . $i, $group ) ) ) {
        $free = 1;
      } else {
        $i++;
      }
    }
  }

  $species .= $i;
  return ($species);
} ## end sub find_unique_species



sub load_adaptors {
  my ($dba) = @_;

  my %pairs = %{ $dba->get_available_adaptors() };

  while ( my ( $key, $value ) = each(%pairs) ) {
    Bio::EnsEMBL::Registry->add_adaptor( $dba->species(), $dba->group(),
      $key, $value );
  }
}

sub load_and_attach_dnadb_to_core {
  my ($dba) = @_;

  load_adaptors($dba);
  $reg->add_DNAAdaptor( $dba->species(), $dba->group(), $dba->species(),
    'core' );
}


=head2 load_core
  Arg [1]    : DBAdaptor with DBConnection already attached
  Returntype : DBAdaptor
  Exceptions : none
=cut
sub load_core      { load_adaptors(@_) }


#
# 1) core. no need to add dnadb
# 2) not core add dnadb
# 3) 
#

=head2 load_compara
  Arg [1]    : DBAdaptor with DBConnection already attached
  Returntype : DBAdaptor
  Exceptions : none
=cut
sub load_compara   { load_adaptors(@_) }

=head2 load_hive
  Arg [1]    : DBAdaptor with DBConnection already attached
  Returntype : DBAdaptor
  Exceptions : none
=cut
sub load_hive      { load_adaptors(@_) }

=head2 load_pipeline
  Arg [1]    : DBAdaptor with DBConnection already attached
  Returntype : DBAdaptor
  Exceptions : none
=cut
sub load_pipeline  { load_adaptors(@_) }

=head2 load_SNP
  Arg [1]    : DBAdaptor with DBConnection already attached
  Returntype : DBAdaptor
  Exceptions : none
=cut
sub load_SNP       { load_adaptors(@_) }

sub load_haplotype { load_adaptors(@_) }

sub load_ontology  { load_adaptors(@_) }


# these that need to attach to the core to get the sequence data

sub load_estgene       { load_and_attach_dnadb_to_core(@_) }

sub load_variation { load_and_attach_dnadb_to_core(@_) }

sub load_funcgen   { load_and_attach_dnadb_to_core(@_) }

=head2 load_otherfeatures
  Arg [1]    : DBAdaptor with DBConnection alredy attached
  Returntype : DBAdaptor
  Exceptions : none

=cut
sub load_otherfeatures { load_and_attach_dnadb_to_core(@_) }

sub load_rnaseq { load_and_attach_dnadb_to_core(@_) }

=head2 load_vega
  Arg [1]    : DBAdaptor with DBConnection already attached
  Returntype : DBAdaptor
  Exceptions : none
=cut
sub load_vega { load_and_attach_dnadb_to_core(@_) }


sub add_alias {
  my ( $class, @args ) = @_;

  my ( $species, $aliases ) = rearrange( [qw(SPECIES ALIAS)], @args );

  # Make sure it exists itself
  Bio::EnsEMBL::Registry->add_alias( $species, $species );

  if ( defined($aliases) ) {
    foreach my $ali (@$aliases) {
      Bio::EnsEMBL::Registry->add_alias( $species, $ali );
    }
  }
}

# WARNING:  "CONVENIENCE METHOD" for retriving the species name when one was 
#           not set. Regulation DB requirement
sub pre_funcgen_hook {
  my ($dba) = @_;
  if(! $dba->species() ) {
    warn "Setting name";
    my $name = $dba->dbc()->sql_helper()->execute_single_result(
      -SQL => 'select meta_value from meta where meta_key =?',
      -PARAMS => ['species.production_name'],
    );
    $dba->dbc()->disconnect_if_idle();
    $dba->species($name);
  }
  return;
}

#
# overwrite/load new types. Done this way to enable no changes to CVS for
# external users. External users should add there own "GROUPS" in the file
# User_defined_load.
#

eval{ require Bio::EnsEMBL::Utils::User_defined_load };

1;
