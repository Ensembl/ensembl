#
# EnsEMBL module for Bio::EnsEMBL::DBSQL::AssemblyExceptionFeatureAdaptor
#
# Cared for by Ewan Birney <birney@ebi.ac.uk>
#
# Copyright EMBL/EBI
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::DBSQL::AssemblyExceptionFeatureAdaptor 

=head1 SYNOPSIS

my $assembly_exception_feature_adaptor = $database_adaptor->get_AssemblyExceptionFeatureAdaptor();
@assembly_exception_features = $assembly_exception_feature_adaptor->fetch_by_Slice($slice);

=head1 DESCRIPTION

Assembly Exception Feature Adaptor - database access for assembly exception features

=head1 AUTHOR - Glenn Proctor

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

package Bio::EnsEMBL::DBSQL::AssemblyExceptionFeatureAdaptor;
use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor;
use Bio::EnsEMBL::AssemblyExceptionFeature;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);

@ISA = qw(Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor);

=head2 generic_fetch

  Arg [1]    : (optional) string $constraint
               An SQL query constraint (i.e. part of the WHERE clause)
  Arg [2]    : (optional) Bio::EnsEMBL::AssemblyMapper $mapper
               A mapper object used to remap features
               as they are retrieved from the database
  Arg [3]    : (optional) Bio::EnsEMBL::Slice $slice
               A slice that features should be remapped to
  Arg [4]    : (optional) boolean $keep_all
               Set to 1 if all features, even ones entirely off slice,
               should be kept
  Example    : $fts = $a->generic_fetch('contig_id in (1234, 1235)', 'Swall');
  Description: Overrides the default generic fetch for this object. Will
               ignore slice constraints.
  Returntype : listref of Bio::EnsEMBL::SeqFeature in contig coordinates
  Exceptions : none
  Caller     : BaseFeatureAdaptor

=cut

sub generic_fetch {
  my ($self, $constraint, $mapper, $slice, $keep_all) = @_;

  my $columns = join(', ', $self->_columns());

  my $db = $self->db();

  my $sql = "SELECT $columns FROM assembly_exception ae";

  # if there is a slice constraint we remove the constraints ...
  unless( $constraint =~ /seq_region_id/ || $constraint eq "" ) {
    $sql .= "WHERE $constraint";
  }

  my $sth = $db->prepare($sql);

  $sth->execute;

  my $res = $self->_objs_from_sth($sth, $mapper, $slice, $keep_all);

  return $res;
}

=head2 _tablename

  Arg [1]    : none
  Example    : none
  Description: PROTECTED implementation of superclass abstract method
               returns the names, aliases of the tables to use for queries
  Returntype : list of listrefs of strings
  Exceptions : none
  Caller     : internal

=cut

sub _tables {
  my $self = shift;

  return ['assembly_exception', 'ae'];
}


=head2 _columns

  Arg [1]    : none
  Example    : none
  Description: PROTECTED implementation of superclass abstract method
               returns a list of columns to use for queries
  Returntype : list of strings
  Exceptions : none
  Caller     : internal

=cut

sub _columns {
  my $self = shift;

  return qw( ae.seq_region_id ae.seq_region_start ae.seq_region_end ae.exc_type ae.exc_seq_region_id ae.exc_seq_region_start ae.exc_seq_region_end );
}


=head2 _objs_from_sth

  Arg [1]    : hash reference $hashref
  Example    : none
  Description: PROTECTED implementation of superclass abstract method.
               creates AssemblyExceptionFeatures from an executed DBI statement handle.
  Returntype : list reference to Bio::EnsEMBL::AssemblyExceptionFeature objects
  Exceptions : none
  Caller     : internal

=cut

sub _objs_from_sth {

  my ($self, $sth, $mapper, $dest_slice) = @_;

  my @features;

  my $slice_adaptor = $self->db()->get_SliceAdaptor();

  my ($seq_region_id, $seq_region_start, $seq_region_end, $type,
      $exc_seq_region_id, $exc_seq_region_start, $exc_seq_region_end);

  $sth->bind_columns(\$seq_region_id, \$seq_region_start, \$seq_region_end, \$type,
		     \$exc_seq_region_id, \$exc_seq_region_start, \$exc_seq_region_end);

  while ($sth->fetch()) {

    my $sr_slice  = $slice_adaptor->fetch_by_seq_region_id($seq_region_id);
    my $exc_slice = $slice_adaptor->fetch_by_seq_region_id($exc_seq_region_id);

    # each row creates TWO features, each of which has alternate_slice pointing to the "other" one

    push @features,
      Bio::EnsEMBL::AssemblyExceptionFeature->new('-start'           => $seq_region_start,
						  '-end'             => $seq_region_end,
						  '-strand'          => 1,
						  '-adaptor'         => $self,
						  '-slice'           => $sr_slice,
						  '-alternate_slice' => $exc_slice,
						  '-type'            => $type);

    push @features, 
      Bio::EnsEMBL::AssemblyExceptionFeature->new('-start'           => $exc_seq_region_start,
						  '-end'             => $exc_seq_region_end,
						  '-strand'          => 1,
						  '-adaptor'         => $self,
						  '-slice'           => $exc_slice,
						  '-alternate_slice' => $sr_slice,
						  '-type'            => $type);

  }

  return \@features;
}


1;
