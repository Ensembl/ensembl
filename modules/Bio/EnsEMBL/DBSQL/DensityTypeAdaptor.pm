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

Bio::EnsEMBL::DBSQL::DensityTypeAdaptor

=head1 SYNOPSIS

  my $density_type_adaptor =
    $registry->get_adaptor( 'Human', 'Core', 'DensityType' );

  my @density_types = @{ $density_type_adaptor->fetch_all() };

  my $dt = $density_type_adaptor->fetch_by_dbID(12);

=head1 DESCRIPTION

DensityTypeAdaptor - Performs database interaction for DensityType objects.

=head1 METHODS

=cut

package Bio::EnsEMBL::DBSQL::DensityTypeAdaptor;
use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::DBSQL::BaseAdaptor;

use Bio::EnsEMBL::DensityType;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);

@ISA = qw(Bio::EnsEMBL::DBSQL::BaseAdaptor);



=head2 new

  Arg [1]    : see superclass (Bio::EnsEMBL::DBSQL::BaseAdaptor) arguments
  Example    : #use this instead of the constructor directly:
               my $dta = $db_adaptor->get_DensityTypeAdaptor();
  Description: Constructor. Creates a new DensityTypeAdaptor
  Returntype : Bio::EnsEMBL::DBSQL::DensityTypeAdaptor
  Exceptions : none
  Caller     : DBAdaptor
  Status     : Stable

=cut

sub new {
  my $class = shift;

  my $self = $class->SUPER::new(@_);

  $self->{'dbID_cache'} = {};

  return $self;
}


# _tables
#  Arg [1]    : none
#  Description: PROTECTED implementation of superclass abstract method.
#               Returns the names, aliases of the tables to use for queries.
#  Returntype : list of listrefs of strings
#  Exceptions : none
#  Caller     : internal
#  Status     : Stable

sub _tables {
  return (['density_type', 'dt']);
}


# _columns
#  Arg [1]    : none
#  Example    : none
#  Description: PROTECTED implementation of superclass abstract method.
#               Returns a list of columns to use for queries.
#  Returntype : list of strings
#  Exceptions : none
#  Caller     : internal
#  Status     : Stable

sub _columns {
  return ('dt.density_type_id', 'dt.analysis_id', 'dt.block_size', 'dt.region_features', 'dt.value_type');
}




=head2 fetch_all

  Arg [1]    : none
  Example    : my @density_types = @{$density_type_adaptor->fetch_all};
  Description: Retrieves every density type in the database.
               NOTE:  In a multi-species database, this method will
               return all the entries, not just the ones associated with
               the current species.
  Returntype : reference to list of Bio::EnsEMBL::DensityType objects
  Exceptions : none
  Caller     : general, new
  Status     : Stable

=cut

sub fetch_all {
  my $self = shift;

  my @out;

  my $sth = $self->prepare("SELECT density_type_id, analysis_id, block_size,".
                           "       value_type, region_features " .
                           "FROM density_type");

  $sth->execute();

  my($dbID, $analysis_id, $blk_size, $vtype, $region_features );
  $sth->bind_columns(\$dbID, \$analysis_id, \$blk_size, \$vtype, \$region_features );

  my $analysis_adaptor = $self->db->get_AnalysisAdaptor();

  while($sth->fetch()) {
    my $analysis = $analysis_adaptor->fetch_by_dbID($analysis_id);


    my $dt = Bio::EnsEMBL::DensityType->new(-ADAPTOR => $self,
                                            -DBID    => $dbID,
                                            -ANALYSIS => $analysis,
                                            -BLOCK_SIZE => $blk_size,
					    -REGION_FEATURES => $region_features,
                                            -VALUE_TYPE => $vtype);

    $self->{'dbID_cache'}->{$dbID} = $dt;

    push @out, $dt;
  }

  return \@out;
}



=head2 fetch_by_dbID

  Arg [1]    : int $dbID
  Example    : my $dt = $density_type_adaptor->fetch_by_dbID($dbID);
  Description: Retrieves a density type object via its internal identifier
  Returntype : Bio::EnsEMBL::DensityType
  Exceptions : throw if dbID argument not defined
  Caller     : general
  Status     : Stable

=cut

sub fetch_by_dbID {
  my $self = shift;
  my $dbID = shift;

  if(!defined($dbID)) {
    throw("dbID argument must be defined");
  }

  if($self->{'dbID_cache'}->{$dbID}) {
    return $self->{'dbID_cache'}->{$dbID};
  }

  # go back to database and refill caches
  $self->fetch_all();

  return $self->{'dbID_cache'}->{$dbID};
}


=head2 fetch_all_by_logic_name

  Arg [1]    : string $logic_name
  Example    : my @dts = @{$dtype_adaptor->fetch_all('repeat_coverage')};
  Description: Retrieves all density types with a given logic name.
               NOTE:  In a multi-species database, this method will
               return all the entries matching the search criteria, not
               just the ones associated with the current species.
  Returntype : reference to list of Bio::EnsEMBL::DensityTypes
  Exceptions : thrown if logic_name argument is not provided
  Caller     : general
  Status     : Stable

=cut

sub fetch_all_by_logic_name {
  my $self = shift;
  my $logic_name = shift;

  if(!defined($logic_name)) {
    throw("logic_name argument is required.");
  }

  my $analysis_adaptor = $self->db()->get_AnalysisAdaptor();
  my $analysis = $analysis_adaptor->fetch_by_logic_name($logic_name);

  return [] if(!$analysis);

  my $sth = $self->prepare("SELECT density_type_id, block_size,".
                           "       value_type, region_features " .
                           "FROM density_type " .
                           "WHERE analysis_id = ?");
  $sth->bind_param(1,$analysis->dbID,SQL_INTEGER);
  $sth->execute();

  my($dbID, $blk_size, $vtype, $region_features );
  $sth->bind_columns(\$dbID, \$blk_size, \$vtype, \$region_features);

  my @out;

  while($sth->fetch()) {

    my $dt = Bio::EnsEMBL::DensityType->new(-ADAPTOR => $self,
                                            -DBID    => $dbID,
                                            -ANALYSIS => $analysis,
                                            -BLOCK_SIZE => $blk_size,
					    -REGION_FEATURES => $region_features,
                                            -VALUE_TYPE => $vtype);

    $self->{'dbID_cache'}->{$dbID} = $dt;

    push @out, $dt;
  }

  return \@out;
}


=head2 store

  Arg [1]    : list of Bio::EnsEMBL::DensityType @dt
               the density types to store in the database
  Example    : $density_type->store(@density_types);
  Description: Stores a list of density type objects in the database
  Returntype : none
  Exceptions : thrown if @dt is not defined
               or if any elements of @dt are not Bio::EnsEMBL::DensityType 
  Caller     : general
  Status     : Stable

=cut

sub store {
  my ($self,@dt) = @_;

  if( scalar(@dt) == 0 ) {
    throw("Must call store with list of Density Types");
  }

  my $insert_ignore = $self->insert_ignore_clause();
  my $sth = $self->prepare
    ("${insert_ignore} INTO density_type (analysis_id,".
                                  "block_size, value_type, region_features ) ". 
    "VALUES (?, ?, ?, ?)");

  my $db = $self->db();
  my $analysis_adaptor = $db->get_AnalysisAdaptor();

 FEATURE: foreach my $dt ( @dt ) {
    if( !ref $dt || !$dt->isa("Bio::EnsEMBL::DensityType") ) {
      throw("Density Type must be an Ensembl DensityType, " .
            "not a [".ref($dt)."]");
    }

    if($dt->is_stored($db)) {
      next FEATURE;
    }

    if(!defined($dt->analysis())) {
      throw("An analysis must be attached to the density type to be stored.");
    }

    #store the analysis if it has not been stored yet
    if(!$dt->analysis->is_stored($db)) {
      $analysis_adaptor->store($dt->analysis());
    }
	
    my $block_size = $dt->block_size();
    $block_size |= 0;
    my $region_features = $dt->region_features();
    $region_features |= 0;

    $sth->bind_param(1,$dt->analysis->dbID,SQL_INTEGER);
    $sth->bind_param(2,$block_size,SQL_INTEGER);
    $sth->bind_param(3,$dt->value_type,SQL_VARCHAR);
    $sth->bind_param(4,$region_features, SQL_VARCHAR);
    my $inserted = $sth->execute();

    my $dbID;

    # $inserted can be 0E0 which is true but equal to 0
    if(!$inserted || $inserted == 0) {
      # insert failed, presumably because was already stored in database

      my @dts=@{$self->fetch_all_by_logic_name($dt->analysis()->logic_name())};
      my ($stored_dt) = grep {$_->block_size() == $dt->block_size()} @dts;
      if(!$stored_dt) {
        throw("Could not retrieve or store DensityType from database.\n" .
              "Incorrect db permissions or missing density_type table?\n");
      }
      $dbID = $stored_dt->dbID();
    } else {
      $dbID = $self->last_insert_id('density_type_id', undef, 'density_type');
    }

    # next two lines are to set the density type as stored
    $dt->dbID($dbID);
    $dt->adaptor($self);

    $self->{'dbID_cache'}->{$dbID} = $dt;
  }
}

1;
