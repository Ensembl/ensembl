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

Bio::EnsEMBL::DBSQL::MiscSetAdaptor - Provides database interaction for
Bio::EnsEMBL::MiscSet objects.

=head1 SYNOPSIS

  my $msa = $registry->get_adaptor( 'Human', 'Core', 'MiscSet' );

  my $misc_set = $msa->fetch_by_dbID(1234);

  $misc_set = $msa->fetch_by_code('clone');

=head1 DESCRIPTION

This class provides database interactivity for MiscSet objects.
MiscSets are used to classify MiscFeatures into groups.

=head1 METHODS

=cut

package Bio::EnsEMBL::DBSQL::MiscSetAdaptor;

use strict;
use warnings;

use Bio::EnsEMBL::MiscSet;
use Bio::EnsEMBL::DBSQL::BaseAdaptor;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);

use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::DBSQL::BaseAdaptor);


=head2 new

  Arg [...]  : Superclass args.  See Bio::EnsEMBL::DBSQL::BaseAdaptor
  Description: Instantiates a Bio::EnsEMBL::DBSQL::MiscSetAdaptor and
               caches the contents of the MiscSet table.
  Returntype : Bio::EnsEMBL::MiscSet
  Exceptions : none
  Caller     : MiscFeatureAdaptor
  Status     : Stable

=cut


sub new {
  my $class = shift;

  my $self = $class->SUPER::new(@_);

  $self->{'_id_cache'} = {};
  $self->{'_code_cache'} = {};

  # cache the entire contents of the misc set table
  # the table is small and it removes the need to repeatedly query the
  # table or join to the table

  $self->fetch_all();

  return $self;
}




=head2 fetch_all

  Arg [1]    : none
  Example    : foreach my $ms (@{$msa->fetch_all()}) {
                 print $ms->code(), ' ', $ms->name(), "\n";
               }
  Description: Retrieves every MiscSet defined in the DB.
               NOTE:  In a multi-species database, this method will
               return all the entries matching the search criteria, not
               just the ones associated with the current species.
  Returntype : listref of Bio::EnsEMBL::MiscSets
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub fetch_all {
  my $self = shift;

  my $sth = $self->prepare
    ('SELECT misc_set_id, code, name, description, max_length FROM misc_set');

  $sth->execute();

  my ($dbID, $code, $name, $desc, $max_len);
  $sth->bind_columns(\$dbID, \$code, \$name, \$desc, \$max_len);

  my @all;

  while($sth->fetch()) {
    my $ms = Bio::EnsEMBL::MiscSet->new
      (-DBID     => $dbID,
       -ADAPTOR  => $self,
       -CODE     => $code,
       -NAME     =>  $name,
       -DESCRIPTION => $desc,
       -LONGEST_FEATURE => $max_len);

    $self->{'_id_cache'}->{$dbID} = $ms;
    $self->{'_code_cache'}->{lc($code)} = $ms;
    push @all, $ms;
  }

  $sth->finish();

  return \@all;
}



=head2 fetch_by_dbID

  Arg [1]    : int $dbID
               The internal identifier of the misc set to retrieve
  Example    : my $ms = $msa->fetch_by_dbID($dbID);
  Description: Retrieves a misc set via its internal identifier
  Returntype : Bio::EnsEMBL::MiscSet
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub fetch_by_dbID {
  my $self = shift;
  my $dbID = shift;

  if(!$self->{'_id_cache'}->{$dbID}) {
    # on a cache miss reread the whole table and reload the cache
    $self->fetch_all();
  }

  return $self->{'_id_cache'}->{$dbID};
}



=head2 fetch_by_code

  Arg [1]    : string $code
               The unique code of the MiscSet to retrieve
  Example    : my $ms = $msa->fetch_by_code('clone');
  Description: Retrieves a MiscSet via its code
  Returntype : Bio::EnsEMBL::MiscSet
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub fetch_by_code {
  my $self = shift;
  my $code = shift;

  if(!$self->{'_code_cache'}->{lc($code)}) {
    # on cache miss, reread whole table and reload cache
    $self->fetch_all();
  }

  return $self->{'_code_cache'}->{lc($code)};
}



=head2 store

  Arg [1]    : list of MiscSets @mist_sets
  Example    : $misc_set_adaptor->store(@misc_sets);
  Description: Stores a list of MiscSets in the database, and sets the
               dbID and adaptor attributes of the stored sets.
  Returntype : none
  Exceptions : throw on incorrect arguments
               warning if a feature is already stored in this database
  Caller     : MiscFeatureAdaptor::store
  Status     : Stable

=cut

sub store {
  my $self = shift;
  my @misc_sets = @_;

  # we use 'insert ignore' so that inserts can occur safely on the farm
  # otherwise 2 processes could try to insert at the same time and one
  # would fail

  my $insert_ignore = $self->insert_ignore_clause();
  my $sth = $self->prepare(
    qq{${insert_ignore} INTO misc_set (
         code,
         name,
         description,
         max_length
      ) VALUES (?, ?, ?, ?)
    });

  my $db = $self->db();

 SET:
  foreach my $ms (@misc_sets) {
    if(!ref($ms) || !$ms->isa('Bio::EnsEMBL::MiscSet')) {
      throw("List of MiscSet arguments expected.");
    }

    if($ms->is_stored($db)) {
      warning("MiscSet [".$ms->dbID."] is already stored in this database.");
      next SET;
    }

    $sth->bind_param(1,$ms->code,SQL_VARCHAR);
    $sth->bind_param(2,$ms->name,SQL_VARCHAR);
    $sth->bind_param(3,$ms->description,SQL_LONGVARCHAR);
    $sth->bind_param(4,$ms->longest_feature,SQL_INTEGER);

    my $num_inserted = $sth->execute();

    my $dbID;

    if($num_inserted == 0) {
      # insert failed because set with this code already exists
      my $sth2 = $self->prepare("SELECT misc_set_id from misc_set " .
                                "WHERE code = ?");
      $sth2->bind_param(1,$ms->code,SQL_VARCHAR);
      $sth2->execute();

      ($dbID) = $sth2->fetchrow_array();

      if($sth2->rows() != 1) {
        throw("Could not retrieve or store MiscSet, code=[".$ms->code."]\n".
              "Wrong database user/permissions?");
      }
    } else {
      $dbID = $self->last_insert_id('misc_set_id', undef, 'misc_set');
    }

    $ms->dbID($dbID);
    $ms->adaptor($self);

    # update the internal caches
    $self->{'_id_cache'}->{$dbID} = $ms;
    $self->{'_code_cache'}->{lc($ms->code())} = $ms;
  }

  return;
}

=head2 update

  Arg [1]    : Bio::EnsEMBL::MiscSet $miscset
  Example    : $adaptor->update($miscset)
  Description: Updates this misc_set in the database
  Returntype : int 1 if update is performed, undef if it is not
  Exceptions : throw if arg is not an misc_set object
  Caller     : ?
  Status     : Stable

=cut

sub update {
  my $self = shift;
  my $m    = shift;

  if (!ref($m) || !$m->isa('Bio::EnsEMBL::MiscSet')) {
    throw("Expected Bio::EnsEMBL::MiscSet argument.");
  }

  if(!$m->is_stored($self->db())) {
    return undef;
  }

  my $sth = $self->prepare("UPDATE misc_set ".
			   "SET code =?, name =?, description = ?, max_length = ? ".
			   "WHERE misc_set_id = ?");

  $sth->bind_param(1,$m->code,SQL_VARCHAR);
  $sth->bind_param(2,$m->name,SQL_VARCHAR);
  $sth->bind_param(3,$m->description,SQL_VARCHAR);
  $sth->bind_param(4,$m->longest_feature,SQL_INTEGER);
  $sth->bind_param(5,$m->dbID,SQL_INTEGER);

  $sth->execute();
  $sth->finish();

 # update the internal caches
  $self->{'_id_cache'}->{$m->dbID} = $m;
  $self->{'_code_cache'}->{lc($m->code())} = $m;
}

1;
