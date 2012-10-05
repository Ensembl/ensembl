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

Bio::EnsEMBL::Pipeline::Production::ClassSpeciesFactory

=head1 DESCRIPTION

An extension of the SpeciesFactory code. This uses the ensembl production
database to decide if 
- there has been a change to the species
- there is a variation database associated

Allowed parameters are:

=over 8

=item release - Needed to query production with

=back

The registry should also have a DBAdaptor for the production schema 
registered under the species B<multi> and the group B<production>. 

The code adds an additional flow output:

=over 8

=item 4 - Perform DNA reuse

=back

=cut

package Bio::EnsEMBL::Pipeline::Production::ClassSpeciesFactory;

use strict;
use warnings;

use base qw/Bio::EnsEMBL::Pipeline::SpeciesFactory/;

use Bio::EnsEMBL::Registry;
use File::Spec;



sub run {
  my ($self) = @_;
  my @dbs;
  foreach my $dba (@{$self->param('dbas')}) {
    if(!$self->process_dba($dba)) {
      $self->fine('Skipping %s', $dba->species());
      next;
    }

    my $all = $self->production_flow($dba, 'all');
    if ($self->param('run_all')) {
      $all = 2;
    }
    if($all) {
      push(@dbs, [$self->input_id($dba), $all]);
      my $vega = $self->production_flow($dba, 'vega');
      if ($vega) {
        push(@dbs, [$self->input_id($dba), $vega]);
      }
      my $variation = $self->production_flow($dba, 'variation');
      if ($variation) {
        push(@dbs, [$self->input_id($dba), $variation]);
      }
      my $karyotype = $self->production_flow($dba, 'karyotype');
      if ($karyotype) {
        push(@dbs, [$self->input_id($dba), $karyotype]);
      }
    }
  }
  $self->param('dbs', \@dbs);
  return;
}


sub input_id {
  my ($self, $dba) = @_;
  my $mc = $dba->get_MetaContainer();
  $dba->dbc()->disconnect_if_idle();
  my $input_id = {
    species => $mc->get_production_name(),
  };
  return $input_id;
}


sub has_karyotype {
  my ($self, $dba) = @_;
  my $helper = $dba->dbc()->sql_helper();
  my $sql = q{
    SELECT count(*)
    FROM seq_region_attrib sa, attrib_type at, seq_region s, coord_system cs
    WHERE s.seq_region_id = sa.seq_region_id
    AND cs.coord_system_id = s.coord_system_id
    AND at.attrib_type_id = sa.attrib_type_id
    AND cs.species_id = ?
    AND at.code = 'karyotype_rank' };
  my $count = $helper->execute_single_result(-SQL => $sql, -PARAMS => [$dba->species_id()]);
  $dba->dbc()->disconnect_if_idle();
  return $count;
}


sub has_vega {
  my ($self, $dba) = @_;
  my $production_name  = $dba->get_MetaContainer()->get_production_name();
  my $sql = q{
     SELECT count(*)
     FROM db d, species s
     WHERE db_type = 'vega'
     AND d.is_current = 1
     AND s.species_id = d.species_id
     AND db_name = ?
     AND db_release = ? };
  my $prod_dba = $self->get_production_DBAdaptor();
  my @params = ($production_name, $self->param('release'));
  my $result = $prod_dba->dbc()->sql_helper()->execute_single_result(-SQL => $sql, -PARAMS => [@params]);
  $prod_dba->dbc()->disconnect_if_idle();
  return $result;
}


sub production_flow {
  my ($self, $dba, $class) = @_;
  if($self->is_run($dba, $class)) {
    if ($class =~ 'vega') {
      return 5;
    }
    if ($class =~ 'variation') {
      return 4;
    }
    if ($class =~ 'karyotype') {
      return 3;
    }
    if ($class =~ 'all') {
      return 2;
    }
  }
}


sub is_run {
  my ($self, $dba, $class) = @_;
  my $production_name  = $dba->get_MetaContainer()->get_production_name();

  if ($class =~ 'karyotype') {
    return $self->has_karyotype($dba);
  }

  if ($class =~ 'vega') {
    return $self->has_vega($dba);
  }
  
  my $sql = <<'SQL';
     SELECT count(*)
     FROM   db_list dl, db d
     WHERE  dl.db_id = d.db_id and db_type = 'core' and is_current = 1 
     AND full_db_name like ?
     AND    species_id IN (
     SELECT species_id 
     FROM   changelog c, changelog_species cs 
     WHERE  c.changelog_id = cs.changelog_id 
     AND    release_id = ?
     AND    status not in ('cancelled', 'postponed') 
     AND    (gene_set = 'Y' OR assembly = 'Y' OR repeat_masking = 'Y' OR variation_pos_changed = 'Y'))
SQL

  my @params = ("$production_name%", $self->param('release'));

  if ($class =~ 'variation') {
    $sql .= <<'SQL';
     AND    species_id IN (
     SELECT distinct species_id 
     FROM   db 
     WHERE  db_release = ? AND db_type = 'variation')
SQL
    push (@params, $self->param('release'));
  }

  $dba->dbc()->disconnect_if_idle();
  my $prod_dba = $self->get_production_DBAdaptor();
  my $result = $prod_dba->dbc()->sql_helper()->execute_single_result(-SQL => $sql, -PARAMS => [@params]);
  $prod_dba->dbc()->disconnect_if_idle();
  return $result;
}


sub write_output {
  my ($self) = @_;
  $self->do_flow('dbs');
  return;
}



1;
