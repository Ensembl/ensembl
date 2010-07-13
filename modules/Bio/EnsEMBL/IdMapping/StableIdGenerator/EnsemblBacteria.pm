=head1 LICENSE

  Copyright (c) 1999-2010 The European Bioinformatics Institute and
  Genome Research Limited.  All rights reserved.

  This software is distributed under a modified Apache license.
  For license details, please see

    http://www.ensembl.org/info/about/code_licence.html

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <ensembl-dev@ebi.ac.uk>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk@ensembl.org>.

=cut

package Bio::EnsEMBL::IdMapping::StableIdGenerator::EnsemblBacteria;
use strict;
use warnings;
no warnings 'uninitialized';
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use base qw(Bio::EnsEMBL::IdMapping::StableIdGenerator::EnsemblGeneric);

#
# new generator for EnsemblBacteria
# 1. generates EB style IDs
# 2. updates base ID to allow for multi-species DBs to be correctly incremented
#
sub initial_stable_id {
	my $self = shift;
	my $type = shift;
	my $base = $self->get_base();
	# retrieve last good stable ID from hash
	my $init_stable_id = $self->{stable_id_list}{$type};
	if ( !$init_stable_id ) {
		$self->logger->debug(
			"Finding new init_stable_id as base for new $type stable IDs.\n");

		# use stable ID from configuration if set
		if ( $init_stable_id =
			$self->conf->param("starting_${type}_stable_id") )
		{
			$self->logger->debug(
"Using pre-configured $init_stable_id as base for new $type stable IDs.\n"
			);
			return $init_stable_id;
		}
		my $s_dba = $self->cache->get_DBAdaptor('source');
		my $s_dbh = $s_dba->dbc->db_handle;

		# look in the ${type}_stable_id table first
		my $sql =
qq(SELECT MAX(stable_id) FROM ${type}_stable_id where stable_id like '${base}%');
		print $sql;
		$init_stable_id = $self->fetch_value_from_db( $s_dbh, $sql );

		# also look in gene_archive to make sure there are no larger Ids there
		unless ( $type eq 'exon' ) {
			$sql = qq(SELECT MAX(${type}_stable_id) FROM gene_archive);
			my $archived_stable_id = $self->fetch_value_from_db( $s_dbh, $sql );
			if (    $archived_stable_id
				and $self->is_valid($archived_stable_id)
				and ( $archived_stable_id gt $init_stable_id ) )
			{
				$init_stable_id = $archived_stable_id;
			}
		}
		$self->{stable_id_list}{$type} = $init_stable_id;
	} else {
		$self->logger->debug(
"Using preexisting initial $init_stable_id as base for new $type stable IDs.\n"
		);
	}
	if ($init_stable_id) {

	 # since $init_stable_id now is the highest existing stable Id for this
	 # object type, we need to increment it to find the first one we want to use
	 # for new assignments
		$init_stable_id = $self->increment_stable_id( $init_stable_id, $type );
		$self->logger->debug(
			"Using $init_stable_id as base for new $type stable IDs.\n");
	} else {
		$self->logger->warning(
			"Can't find highest ${type}_stable_id in source db.\n");
		my $pref =
		  $self->cache->get_DBAdaptor('target')->get_MetaContainer()
		  ->list_value_by_key('species.stable_id_prefix');
		if ($pref) {
			$init_stable_id =
			  $pref->[0] . substr( uc($type), 0, 1 ) . '00000000000';
			$self->logger->debug(
				"Using $init_stable_id as base for new $type stable IDs.\n");
		}
	}
	return $init_stable_id;
}

sub increment_stable_id {
	my $self      = shift;
	my $stable_id = shift;
	my $type      = shift;
	unless ( $self->is_valid($stable_id) ) {
		throw("Unknown or missing stable ID: $stable_id.");
	}
	my $base = $self->get_base();
	$stable_id =~ /$base([A-Z]{1,4})(\d{11})/;
	my $number = $2;
	my $new_stable_id = $base . $1 . ( ++$number );
	$self->{stable_id_list}{$type} = $new_stable_id;
	return $new_stable_id;
}

=head2 is_valid

  Arg[1]      : String $stable_id - the stable Id to check
  Example     : unless ($generator->is_valid($stable_id)) {
                  die "Invalid stable Id: $stable_id.\n";
                }
  Description : Tests a stable Id to be valid (according to the Ensembl stable
                Id format definition).
  Return type : Boolean - TRUE if valid, FALSE otherwise
  Exceptions  : none
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub is_valid {
  my ( $self, $stable_id ) = @_;

  my $base = $self->get_base();

  return ( $stable_id
           and ( $stable_id =~ /$base([A-z]{1,4})(\d{11})/ ) );
}

sub get_base { return 'EB' }

1;
