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

=cut

package Bio::EnsEMBL::DBSQL::MetaCoordContainer;

use strict;
use warnings;

use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::Utils::Exception;

@ISA = qw(Bio::EnsEMBL::DBSQL::BaseAdaptor);

sub new {
  my $class = shift;

  my $self = $class->SUPER::new(@_);

  #
  # Retrieve the list of the coordinate systems that features are stored
  # in and cache them.
  #

  my @coord_systems =
    @{ $self->db()->dnadb()->get_CoordSystemAdaptor->fetch_all() };

  my @cs_ids;
  foreach my $cs (@coord_systems) { push( @cs_ids, $cs->dbID() ) }

  my $sth = $self->prepare(
              'SELECT mc.table_name, mc.coord_system_id, mc.max_length '
                . 'FROM meta_coord mc '
                . 'WHERE mc.coord_system_id in ('
                . join( ',', @cs_ids )
                . ')' );

  $sth->execute();

  while ( my ( $table_name, $cs_id, $max_length ) =
          $sth->fetchrow_array() )
  {
    $table_name = lc($table_name);

    $self->{'_feature_cache'}->{$table_name} ||= [];

    push( @{ $self->{'_feature_cache'}->{$table_name} }, $cs_id );
    $self->{'_max_len_cache'}->{$cs_id}->{$table_name} = $max_length;
  }
  $sth->finish();

  return $self;
} ## end sub new




=head2 fetch_all_CoordSystems_by_feature_type

  Arg [1]    : string $table - the name of the table to retrieve coord systems
               for.  E.g. 'gene', 'exon', 'dna_align_feature'
  Example    : @css = @{$mcc->fetch_all_CoordSystems_by_feature_type('gene')};
  Description: This retrieves the list of coordinate systems that features
               in a particular table are stored.  It is used internally by
               the API to perform queries to these tables and to ensure that
               features are only stored in appropriate coordinate systems.
  Returntype : listref of Bio::EnsEMBL::CoordSystem objects
  Exceptions : throw if name argument not provided
  Caller     : BaseFeatureAdaptor
  Status     : Stable

=cut

sub fetch_all_CoordSystems_by_feature_type {
  my $self = shift;
  my $table = lc(shift); #case insensitive matching

  throw('Name argument is required') unless $table;

  if(!$self->{'_feature_cache'}->{$table}) {
    return [];
  }

  my @cs_ids = @{$self->{'_feature_cache'}->{$table}};
  my @coord_systems;

  my $csa = $self->db->get_CoordSystemAdaptor();

  foreach my $cs_id (@cs_ids) {
    my $cs = $csa->fetch_by_dbID($cs_id);

    if(!$cs) {
      throw("meta_coord table refers to non-existant coord_system $cs_id");
    }

    push @coord_systems, $cs;
  }

  return \@coord_systems;
}



=head2 fetch_max_length_by_CoordSystem_feature_type

  Arg [1]    : Bio::EnsEMBL::CoordSystem $cs
  Arg [2]    : string $table
  Example    : $max_len = $mcc->fetch_max_length_by_CoordSystem_feature_type($cs,'gene');
  Description: Returns the maximum length of features of a given type in
               a given coordinate system.
  Returntype : int or undef
  Exceptions : throw on incorrect argument
  Caller     : BaseFeatureAdaptor
  Status     : Stable

=cut


sub fetch_max_length_by_CoordSystem_feature_type {
  my $self = shift;
  my $cs = shift;
  my $table = shift;

  if(!ref($cs) || !$cs->isa('Bio::EnsEMBL::CoordSystem')) {
    throw('Bio::EnsEMBL::CoordSystem argument expected');
  }

  throw("Table name argument is required") unless $table;

  return $self->{'_max_len_cache'}->{$cs->dbID()}->{lc($table)};
}



=head2 add_feature_type

  Arg [1]    : Bio::EnsEMBL::CoordSystem $cs
               The coordinate system to associate with a feature table
  Arg [2]    : string $table - the name of the table in which features of
               a given coordinate system will be stored in
  Arg [3]    : int $length
               This length is used to update the max_length in the database
               and the internal cache. 
  Example    : $csa->add_feature_table($chr_coord_system, 'gene');
  Description: This function tells the coordinate system adaptor that
               features from a specified table will be stored in a certain
               coordinate system.  If this information is not already stored
               in the database it will be added.
  Returntype : none
  Exceptions : none
  Caller     : BaseFeatureAdaptor
  Status     : Stable

=cut


sub add_feature_type {
  my $self = shift;
  my $cs   = shift;
  my $table = lc(shift);
  my $length = shift;
  if(!ref($cs) || !$cs->isa('Bio::EnsEMBL::CoordSystem')) {
    throw('CoordSystem argument is required.');
  }

  if(!$table) {
    throw('Table argument is required.');
  }

  my $cs_ids = $self->{'_feature_cache'}->{$table} || [];

  my ($exists) = grep {$cs->dbID() == $_} @$cs_ids;
  if( $exists ) {
    if( !$self->{'_max_len_cache'}->{$cs->dbID()}->{$table} ||
        $self->{'_max_len_cache'}->{$cs->dbID()}->{$table} < $length ) {
      my $sth = $self->prepare('UPDATE meta_coord ' .
                               "SET max_length = $length " .
                               'WHERE coord_system_id = ? ' .
                               'AND table_name = ? '.
                               "AND (max_length<$length ".
                               "OR max_length is null)");
      $sth->execute( $cs->dbID(), $table );
      $self->{'_max_len_cache'}->{$cs->dbID()}->{$table} = $length;
    }
    return;
  }

  #store the new tablename -> coord system relationship in the db
  #ignore failures b/c during the pipeline multiple processes may try
  #to update this table and only the first will be successful
  my $sth = $self->prepare('INSERT IGNORE INTO meta_coord ' .
                              'SET coord_system_id = ?, ' .
                                  'table_name = ?, ' .
			   'max_length = ? ' 
			  );

  $sth->execute($cs->dbID, $table, $length );

  #update the internal cache
  $self->{'_feature_cache'}->{$table} ||= [];
  push @{$self->{'_feature_cache'}->{$table}}, $cs->dbID();
  $self->{'_max_len_cache'}->{$cs->dbID()}->{$table} = $length;

  return;
}


1;
