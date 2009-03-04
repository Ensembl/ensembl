=head1 LICENSE

  Copyright (c) 1999-2009 The European Bioinformatics Institute and
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

=head1 NAME

Bio::EnsEMBL::DBSQL::OligoArrayAdaptor - A database adaptor for fetching and
storing OligoArray objects.

=head1 SYNOPSIS

  my $oaa = $db->get_OligoArrayAdaptor();

  my $array  = $oaa->fetch_by_name('Array-1');
  my @arrays = @{ $oaa->fetch_all() };

=head1 DESCRIPTION

The OligoArrayAdaptor is a database adaptor for storing and retrieving
OligoArray objects.

=head1 METHODS

=cut

package Bio::EnsEMBL::DBSQL::OligoArrayAdaptor;

use strict;
use warnings;

use Bio::EnsEMBL::Utils::Exception qw( warning );
use Bio::EnsEMBL::OligoArray;
use Bio::EnsEMBL::DBSQL::BaseAdaptor;

use vars qw(@ISA);
@ISA = qw(Bio::EnsEMBL::DBSQL::BaseAdaptor);


=head2 fetch_by_name

  Arg [1]    : string - name of an array
  Example    : my $array = $oaa->fetch_by_name('Array-1');
  Description: Retrieves a named OligoArray object from the database.
  Returntype : Bio::EnsEMBL::OligoArray
  Exceptions : None
  Caller     : General
  Status     : Medium Risk

=cut

sub fetch_by_name {
    my $self = shift;
    my $name = shift;
    
    $self->bind_param_generic_fetch($name,SQL_VARCHAR);

    my $result = $self->generic_fetch("oa.name = ?");
	
    if (scalar @$result > 1) {
		warning("Array $name is not unique in the database, but only one result has been returned");
    } 

    return $result->[0];
}

=head2 fetch_all_by_type

  Arg [1]    : List of strings - type(s) (e.g. AFFY or OLIGO)
  Example    : my @arrays = @{$oaa->fetch_all_by_type('OLIGO')};
  Description: Fetch all arrays of a particular type.
  Returntype : Listref of Bio::EnsEMBL::OligoArray objects
  Exceptions : Throws if no type is provided
  Caller     : General
  Status     : Medium Risk

=cut

sub fetch_all_by_type {
	my ($self, @types) = @_;
	
	throw('Need type as parameter') if !@types;
	
	my $constraint;
	if (scalar @types == 1) {
		$constraint = qq( oa.type = ? );
		$self->bind_param_generic_fetch($types[0],SQL_VARCHAR);
	} else {
		$constraint = join q(','), @types;
		$constraint = qq( oa.type IN ('$constraint') );
	}

	return $self->generic_fetch($constraint);
}

=head2 fetch_attributes

  Arg [1]    : Bio::EnsEMBL::OligoArray - array to fetch attributes for
  Example    : None
  Description: This function is solely intended to lazy load attributes into
               empty OligoArray objects. You should not need to call this.
  Returntype : None
  Exceptions : None
  Caller     : Bio::EnsEMBL::OligoArray getters
  Status     : Medium Risk

=cut

sub fetch_attributes {
    my $self = shift;
    my $array = shift;

    my $tmp_array = $self->fetch_by_dbID( $array->dbID() );
    %$array = %$tmp_array;
}

=head2 _tables

  Args       : None
  Example    : None
  Description: PROTECTED implementation of superclass abstract method.
               Returns the names and aliases of the tables to use for queries.
  Returntype : List of listrefs of strings
  Exceptions : None
  Caller     : Internal
  Status     : Medium Risk

=cut

sub _tables {
	my $self = shift;
	
	return ['oligo_array', 'oa'];
}

=head2 _columns

  Args       : None
  Example    : None
  Description: PROTECTED implementation of superclass abstract method.
               Returns a list of columns to use for queries.
  Returntype : List of strings
  Exceptions : None
  Caller     : Internal
  Status     : Medium Risk

=cut

sub _columns {
	my $self = shift;
	
	return qw( oa.oligo_array_id oa.parent_array_id oa.probe_setsize oa.name oa.type );
}

=head2 _objs_from_sth

  Arg [1]    : DBI statement handle object
  Example    : None
  Description: PROTECTED implementation of superclass abstract method.
               Creates OligoArray objects from an executed DBI statement
			   handle.
  Returntype : Listref of Bio::EnsEMBL::OligoArray objects
  Exceptions : None
  Caller     : Internal
  Status     : Medium Risk

=cut

sub _objs_from_sth {
	my ($self, $sth) = @_;
	
	my (@result, $array_id, $parent_id, $setsize, $name, $type);
	
	$sth->bind_columns( \$array_id, \$parent_id, \$setsize, \$name, \$type );
	
	while ( $sth->fetch() ) {
		my $array = Bio::EnsEMBL::OligoArray->new(
			-dbID    => $array_id,
			-adaptor => $self,
			-name    => $name,
			-setsize => $setsize,
			-type    => $type,
	  	);
		push @result, $array;
		if ($parent_id) {
			my $parent_array = Bio::EnsEMBL::OligoArray->new(
				-dbID    => $parent_id,
				-adaptor => $self,
			);
			$array->superset($parent_array);
		}
	}
	return \@result;
}

=head2 store

  Args       : List of Bio::EnsEMBL::OligoArray objects
  Example    : $oaa->store($array1, $array2, $array3);
  Description: Stores given OligoArray objects in the database. Should only be
               called once per array because no checks are made for duplicates.
			   Sets dbID and adaptor on the objects that it stores.
  Returntype : None
  Exceptions : None
  Caller     : General
  Status     : Medium Risk

=cut

sub store {
    my $self = shift;
    my @args = @_;
    
    foreach my $array (@args) {
		if ( !$array->isa('Bio::EnsEMBL::OligoArray') ) {
			warning('Can only store OligoArray objects');
			next;
		}
		
		# Has array already been stored?
		next if ( $array->dbID() && $array->adaptor() == $self );

		my $superset = $array->superset();
		if ( defined $superset && !$superset->dbID() ) {
			$self->store($superset);
		}

		my $sth = $self->prepare("
			INSERT INTO oligo_array
			(name, probe_setsize, type, parent_array_id)
			VALUES (?, ?, ?, ?)
		");
		$sth->bind_param(1, $array->name(),    SQL_VARCHAR);
		$sth->bind_param(2, $array->setsize(), SQL_INTEGER);
		$sth->bind_param(3, $array->type(),    SQL_VARCHAR);
		if (defined $superset) {
			$sth->bind_param(4, $superset->dbID(), SQL_INTEGER);
		} else {
			$sth->bind_param(4, undef);
		}
		$sth->execute();
		my $dbID = $sth->{'mysql_insertid'};
		$array->dbID($dbID);
		$array->adaptor($self);
	}
}

=head2 list_dbIDs

  Args       : None
  Example    : my @array_ids = @{$oaa->list_dbIDs()};
  Description: Gets an array of internal IDs for all OligoArray objects in the
               current database.
  Returntype : List of ints
  Exceptions : None
  Caller     : ?
  Status     : Medium Risk

=cut

sub list_dbIDs {
    my ($self) = @_;
	
    return $self->_list_dbIDs('oligo_array');
}

1;

