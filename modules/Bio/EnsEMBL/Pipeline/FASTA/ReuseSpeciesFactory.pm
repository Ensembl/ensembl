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

Bio::EnsEMBL::Pipeline::FASTA::ReuseSpeciesFactory

=head1 DESCRIPTION

An extension of the SpeciesFactory code. This uses the ensembl production
database to decide if a species has had an update to its DNA or not. An update
means any change to the assembly or repeat masking.

Allowed parameters are:

=over 8

=item release - Needed to query production with

=item ftp_dir - If not specified then we cannot reuse

=item force_species - Specify species we want to redump even though 
                      our queries of production could say otherwise

=back

The registry should also have a DBAdaptor for the production schema 
registered under the species B<multi> and the group B<production>. 

The code adds an additional flow output:

=over 8

=item 4 - Perform DNA reuse

=back

=cut

package Bio::EnsEMBL::Pipeline::FASTA::ReuseSpeciesFactory;

use strict;
use warnings;

use base qw/Bio::EnsEMBL::Pipeline::FASTA::SpeciesFactory/;

use Bio::EnsEMBL::Registry;
use File::Spec;

sub param_defaults {
  my ($self) = @_;
  my $p = {
    %{$self->SUPER::param_defaults()},
    
    force_species => [],
  };
  return $p;
}

sub fetch_input {
  my ($self) = @_;
  $self->SUPER::fetch_input();
  $self->throw("Need to define a release parameter") unless $self->param('release');
  
  #Calculate quick lookup
  my %force_species_lookup;
  foreach my $s (@{$self->param('force_species')}) {
    my $new = Bio::EnsEMBL::Registry->get_alias($s);
    $force_species_lookup{$new} = 1;
  }
  $self->param('force_species_lookup', \%force_species_lookup);
  
  return;
}

sub dna_flow {
  my ($self, $dba) = @_;
  my $parent_flow = $self->SUPER::dna_flow($dba);
  return 0 if ! $parent_flow; #return if parent decided to skip
  if(! $self->param('ftp_dir')) {
    $self->fine('No ftp_dir param defined so will flow %s to %d', $dba->species(), $parent_flow);
    return $parent_flow;
  }
  my $requires_new_dna = $self->requires_new_dna($dba);
  return $parent_flow if $requires_new_dna;
  return 4; #nominated flow for copying
}

sub requires_new_dna {
  my ($self, $dba) = @_;
  
  if($self->force_run($dba)) {
    $self->fine('Automatically forcing DNA rerun for %s', $dba->species());
    return 1;
  }
  
  if(!$self->has_source_dir($dba)) {
    $self->fine('Source directory is missing; forcing DNA rerun for %s', $dba->species());
    return 1;
  }
  
  my $sql = <<'SQL';
select count(*)
from changelog c
join changelog_species cs using (changelog_id)
join species s using (species_id)
where c.release_id = ?
and (c.assembly = ? or c.repeat_masking = ?)
and c.team = ?
and c.status = ?
and production_name = ?
SQL
  my $production_name  = $dba->get_MetaContainer()->get_production_name();
  $dba->dbc()->disconnect_if_idle();
  my $release = $self->param('release');
  my $params = [ $release, 'Y', 'Y', 'Genebuild', 'handed_over', $production_name ];
  my $prod_dba = $self->get_production_DBAdaptor();
  my $result = $prod_dba->dbc()->sql_helper()->execute_single_result(-SQL => $sql, -PARAMS => $params);
  $prod_dba->dbc()->disconnect_if_idle();
  return $result;
}

sub force_run {
  my ($self, $dba) = @_;
  my $new = Bio::EnsEMBL::Registry->get_alias($dba->species());
  return ($self->param('force_species_lookup')->{$new}) ? 1 : 0;
}

sub has_source_dir {
  my ($self, $dba) = @_;
  my $dir = $self->old_path($dba->species());
  $dba->dbc()->disconnect_if_idle();
  return (-d $dir) ? 1 : 0;
}

sub get_production_DBAdaptor {
  my ($self) = @_;
  return Bio::EnsEMBL::Registry->get_DBAdaptor('multi', 'production');
}

1;
