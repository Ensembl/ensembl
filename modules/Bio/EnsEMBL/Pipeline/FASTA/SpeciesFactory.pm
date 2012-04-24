=pod

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

=head1 NAME

Bio::EnsEMBL::Pipeline::FASTA::SpeciesFactory

=head1 DESCRIPTION

A module which generates dump jobs for each species it finds in the Ensembl
Registry. The type of dump performed is controlled by the I<sequence_type_list>
parameter. The species we run the code on can be controlled by specifying
the I<species> parameter or by reducing the number of DBAdaptors loaded into
the registry. 

Allowed parameters are:

=over 8

=item sequence_type_list -  The type of dump to perform. Should be an array and 
                            can contain I<dna>, I<cdna> and I<ncrna>. Defaults
                            to all of these.

=item species - Can be an array of species to perform dumps for or a single
                species name. If specified only jobs will be created for
                those species. Defaults to nothing so all species are processed

item db_types - Specify the types of database to dump. Defaults to core and
                should be an array.

=back

The code flows to two outputs. Please take note if you are reusing this module

=over 8

=item 2 - Perform DNA dumps

=item 3 - Perform Gene dumps

=back

Multiple types of DB can be specifed with the I<db_types> method call but
be aware that this is flowed as 1 job per species for all types.

=cut

package Bio::EnsEMBL::Pipeline::FASTA::SpeciesFactory;

use strict;
use warnings;

use base qw/Bio::EnsEMBL::Pipeline::FASTA::Base/;

use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Utils::Scalar qw/wrap_array/;

sub param_defaults {
  my ($self) = @_;
  return {
    sequence_type_list => [qw/dna cdna ncrna/],
    db_types => [qw/core/],
    species => []
  };
}

sub fetch_input {
  my ($self) = @_;
  
  my $core_dbas = $self->get_DBAdaptors();
  $self->info('Found %d core DBAdaptor(s) to process', scalar(@{$core_dbas}));
  $self->param('dbas', $core_dbas);
  
  my %species_lookup = 
    map { $_ => 1 } 
    map { Bio::EnsEMBL::Registry->get_alias($_)  } 
    @{$self->param('species')};
  $self->param('species_lookup', \%species_lookup);
  
  my %sequence_types = map { $_ => 1 } @{ $self->param('sequence_type_list') };
  $self->param('sequence_types', \%sequence_types);
  
  return;
}
  
sub run {
  my ($self) = @_;
  my @dna;
  my @genes;
  my @species;
  foreach my $dba (@{$self->param('dbas')}) {
    if(!$self->process_dba($dba)) {
      $self->fine('Skipping %s', $dba->species());
      next;
    }
    
    my $dna_flow = $self->dna_flow($dba);
    if($dna_flow) {
      push(@dna, [$self->input_id($dba, 'dna'), $dna_flow]);
    }
    
    my $genes_flow = $self->genes_flow($dba);
    if($genes_flow) {
      push(@genes, [$self->input_id($dba, 'genes'), $genes_flow]);
    }
    
    push(@species, [ { species => $dba->species() }, 5 ]);
  }
  $self->param('dna', \@dna);
  $self->param('genes', \@genes);
  $self->param('species', \@species);
  return;
}

sub write_output {
  my ($self) = @_;
  $self->do_flow('dna');
  $self->do_flow('genes');
  $self->do_flow('species');
  return;
}

sub get_DBAdaptors {
  my ($self) = @_;
  return Bio::EnsEMBL::Registry->get_all_DBAdaptors(-GROUP => 'core');
}

sub do_flow {
  my ($self, $key) = @_;
  my $targets = $self->param($key);
  foreach my $entry (@{$targets}) {
    my ($input_id, $flow) = @{$entry};
    $self->fine('Flowing %s to %d for %s', $input_id->{species}, $flow, $key);
    $self->dataflow_output_id($input_id, $flow);
  }
  return;
}

sub process_dba {
  my ($self, $dba) = @_;
  
  #Reject if DB was ancestral sequences
  return 0 if $dba->species() =~ /ancestral/i;
  
  #If species is defined then make sure we only allow those species through
  if(@{$self->param('species')}) {
    my $lookup = $self->param('species_lookup');
    my $name = $dba->species();
    my $aliases = Bio::EnsEMBL::Registry->get_all_aliases($name);
    push(@{$aliases}, $name);
    my $found = 0;
    foreach my $alias (@{$aliases}) {
      if($lookup->{$alias}) {
        $found = 1;
        last;
      }
    }
    return $found;
  }
  
  #Otherwise just accept
  return 1;
}

# return 0 if we do not want to do any flowing otherwise return 2

sub dna_flow {
  my ($self, $dba) = @_;
  return 0 unless $self->param('sequence_types')->{dna};
  return 2;
}

# return 0 if we do not want to do any flowing otherwise return 3

sub genes_flow {
  my ($self, $dba) = @_;
  my $types = $self->param('sequence_types');
  return 0 if ! $types->{cdna} && ! $types->{ncrna};
  return 3;
}

sub input_id {
  my ($self, $dba, $type) = @_;
  my $mc = $dba->get_MetaContainer();
  my $input_id = {
    db_types => $self->db_types($dba),
    species => $mc->get_production_name(),
  };
  if($type eq 'dna') {
    $input_id->{sequence_type_list} = ['dna']; 
  }
  else {
    my $types = $self->param('sequence_types');
    my @types;
    push(@types, 'cdna') if $types->{cdna};
    push(@types, 'ncrna') if $types->{ncrna};
    $input_id->{sequence_type_list} = \@types;
  }
  return $input_id;
}

sub db_types {
  my ($self, $dba) = @_;
  return $self->param('db_types');
}

1;
