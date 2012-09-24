
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

Bio::EnsEMBL::Pipeline::Production::EGSpeciesFactory

=head1 DESCRIPTION

An extension of the ClassSpeciesFactory code, for use with
EnsemblGenomes, which uses the production database differently
and thus needs a simpler 'is_run' function.

=cut

package Bio::EnsEMBL::Pipeline::Production::EGSpeciesFactory;

use strict;
use warnings;

use base qw/Bio::EnsEMBL::Pipeline::Production::ClassSpeciesFactory/;

sub is_run {
	my ( $self, $dba, $class ) = @_;
	my $production_name = $dba->get_MetaContainer()->get_production_name();

	if ( $class =~ 'karyotype' ) {
		return $self->has_karyotype($dba);
	}
	$dba->dbc()->disconnect_if_idle();
	return 1;
}

sub process_dba {
	my ( $self, $dba ) = @_;
	my $result = $self->SUPER::process_dba($dba);
	if ( $result == 1 && @{ $self->param('division') } ) {
		$result = 0;
		for my $division (@{$self->param('division')}) {
			if($dba->get_MetaContainer()->get_division() eq $division) {
				$result = 1;
				last;
			}
		}
		$dba->dbc()->disconnect_if_idle();
	}
	return $result;
}
1;
