#
# Ensembl module for Bio::EnsEMBL::DBSQL::AffyArrayAdaptor
#
# You may distribute this module under the same terms as Perl itself

=head1 NAME

Bio::EnsEMBL::DBSQL::AffyArrayAdaptor - A database adaptor for fetching and
storing AffyArray objects.

=head1 SYNOPSIS

my $aaa = $db->get_AffyArrayAdaptor();

my $array = $aaa->fetch_by_name('Affy-1');
my @arrays = @{$aaa->fetch_all()};

=head1 DESCRIPTION

The AffyArrayAdaptor is a database adaptor for storing and retrieving
AffyArray objects.

=head1 AUTHOR

This module was originally written by Arne Stabenau, but was changed to be a
subclass of OligoArrayAdaptor by Ian Sealy.

This module is part of the Ensembl project: http://www.ensembl.org/

=head1 CONTACT

Post comments or questions to the Ensembl development list: ensembl-dev@ebi.ac.uk

=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::DBSQL::AffyArrayAdaptor;

use Bio::EnsEMBL::AffyArray;
use Bio::EnsEMBL::DBSQL::OligoArrayAdaptor;

use vars qw(@ISA);
@ISA = qw(Bio::EnsEMBL::DBSQL::OligoArrayAdaptor);


=head2 _default_where_clause

  Args       : None
  Example    : None
  Description: PROTECTED implementation of superclass abstract method.
               Ensures this adaptor only returns Affy arrays.
  Returntype : string
  Exceptions : None
  Caller     : Internal
  Status     : Medium Risk

=cut

sub _default_where_clause {
  my $self = shift;

  return "oa.type='AFFY'";
}

=head2 _objs_from_sth

  Arg [1]    : DBI statement handle object
  Example    : None
  Description: PROTECTED implementation of superclass abstract method.
               Creates AffyArray objects from an executed DBI statement
			   handle.
  Returntype : Listref of Bio::EnsEMBL::AffyArray objects
  Exceptions : None
  Caller     : Internal
  Status     : Medium Risk

=cut

sub _objs_from_sth {
	my ($self, $sth) = @_;
	
	my (@result, $array_id, $parent_id, $setsize, $name, $type);
	
	$sth->bind_columns( \$array_id, \$parent_id, \$setsize, \$name, \$type );
	
	while ( $sth->fetch() ) {
		my $array = Bio::EnsEMBL::AffyArray->new(
			-dbID    => $array_id,
			-adaptor => $self,
			-name    => $name,
			-setsize => $setsize,
	  	);
		push @result, $array;
		if ($parent_id) {
			my $parent_array = Bio::EnsEMBL::AffyArray->new(
				-dbID    => $parent_id,
				-adaptor => $self,
			);
			$array->superset($parent_array);
		}
	}
	return \@result;
}

=head2 list_dbIDs

  Args       : None
  Example    : my @array_ids = @{$aaa->list_dbIDs()};
  Description: Gets an array of internal IDs for all AffyArray objects in the
               current database.
  Returntype : List of ints
  Exceptions : None
  Caller     : ?
  Status     : Medium Risk

=cut

sub list_dbIDs {
	my ($self) = @_;
	
	#return $self->_list_dbIDs('oligo_array');
	# Can't use _list_dbIDs because only want OligoArray objects of type AFFY
	
	my @out;
	my $sql = "SELECT oligo_array_id  FROM oligo_array WHERE type='AFFY'";
	my $sth = $self->prepare($sql);
	$sth->execute;
	
	while (my ($id) = $sth->fetchrow() ) {
		push @out, $id;
	}
	
	$sth->finish;
	
	return \@out;
}

1;

