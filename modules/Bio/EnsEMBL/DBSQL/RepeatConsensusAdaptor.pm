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

Bio::EnsEMBL::DBSQL::RepeatConsensusAdaptor

=head1 SYNOPSIS

  $rca = $database_adaptor->get_RepeatConsensusAdaptor();

  $repeat_consensus = $rca->fetch_by_dbID(132);
  $repeat_consensus = $rca->fetch_by_name_class( 'AluSx', 'SINE/Alu' );

  $rca->store( $rc1, $rc2, $rc3 );

=head1 DESCRIPTION

This is an adaptor for the retrieval and storage of RepeatConsensus
objects.

=head1 METHODS

=cut

package Bio::EnsEMBL::DBSQL::RepeatConsensusAdaptor;

use strict;
use warnings;
use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::RepeatConsensus;
use Bio::EnsEMBL::Utils::Exception qw(throw deprecate);

use base qw(Bio::EnsEMBL::DBSQL::BaseAdaptor);

=head2 fetch_all_repeat_types

  Example			: my $types = $rca->fetch_all_repeat_types();
  Description	: Returns the distinct repeat types available from a database
  Returntype 	: Array
  Exceptions 	: -

=cut


sub fetch_all_repeat_types {
  my ($self) = @_;
  return $self->dbc()->sql_helper()->execute_simple(
    -SQL => 'SELECT DISTINCT repeat_type FROM repeat_consensus');
}


=head2 fetch_by_dbID

  Arg [1]    : int $db_id
               The database identifier for the RepeatConsensus to obtain
  Example    : $repeat_consensus = $repeat_consensus_adaptor->fetch_by_dbID(4);
  Description: Obtains a RepeatConsensus object from the database via its
               primary key. 
  Returntype : Bio::EnsEMBL::RepeatConsensus
  Exceptions : none
  Caller     : general, Bio::EnsEMBL::RepeatFeatureAdaptor
  Status     : Stable

=cut

sub fetch_by_dbID {
    my( $self, $db_id ) = @_;

    my ($rc) = @{$self->_generic_fetch("repeat_consensus_id = $db_id")};

    return $rc;
}



=head2 fetch_by_name

  Arg [1]    : string $name
               the name of the repeat consensus to obtain
  Example    : $rc = $repeat_consensus_adaptor->fetch_by_name('AluSx');
  Description: Obtains a repeat consensus from the database via its name
  Returntype : Bio::EnsEMBL::RepeatConsensus
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub fetch_by_name {
    my( $self, $name ) = @_;

    my ($rc) = @{$self->_generic_fetch("repeat_name = '$name'")};

    return $rc;
}


=head2 fetch_by_name_class

  Arg [1]    : string $name
               the name of the repeat consensus to obtain
  Arg [2]    : string $class
               the class of the repeat consensus to obtain
  Example    : $rc = $repeat_consensus_adaptor->
                 fetch_by_name_class('AluSx', 'SINE/Alu');
  Description: Obtains a repeat consensus from the database
               via its name and class
  Returntype : Bio::EnsEMBL::RepeatConsensus
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub fetch_by_name_class {
    my( $self, $name, $class ) = @_;


    my ($rc) = @{$self->_generic_fetch(qq{
      repeat_name  = '$name'
      AND repeat_class = '$class'
    })};

    return $rc;
}


=head2 fetch_all_by_class_seq

  Arg [1]    : string $class
               the class of the repeat consensus to obtain
  Arg [2]    : string $seq
               the sequence of the repeat consensus to obtain
  Example    : $rc = $repeat_consensus_adaptor->
                 fetch_all_by_class_seq('trf', 'ATGGTGTCA');
  Description: Obtains a repeat consensus from the database
               via its class and sequence
  Returntype : listREF of Bio::EnsEMBL::RepeatConsensus
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub fetch_all_by_class_seq {
    my( $self, $class, $seq ) = @_;

    return $self->_generic_fetch(qq{
            repeat_class     = '$class'
        AND repeat_consensus = '$seq'
    });
}


=head2 _generic_fetch

  Arg [1]    : string $where_clause
  Example    : none
  Description: PRIVATE used to create RepeatConsensus features from an 
               SQL constraint
  Returntype : listref of Bio::EnsEMBL::RepeatConsensus objects
  Exceptions : none
  Caller     : internal
  Status     : Stable

=cut

sub _generic_fetch {
    my( $self, $where_clause ) = @_;

    my( $repeat_consensus_id, $repeat_name, $repeat_class,$repeat_length,
        $repeat_consensus, $repeat_type );

    my $sth = $self->prepare(qq{
        SELECT repeat_consensus_id
          , repeat_name
          , repeat_class
          , repeat_type
          , repeat_consensus
        FROM repeat_consensus
        WHERE }. $where_clause);

    $sth->execute;
    $sth->bind_columns(
        \$repeat_consensus_id,
        \$repeat_name,
        \$repeat_class,
        \$repeat_type,
        \$repeat_consensus
        );

    
    my @consensi;
    while ($sth->fetch) {
      if ($repeat_consensus =~ /^(\d+)\(N\)$/) {
	$repeat_length = $1;
      } else {
	$repeat_length = CORE::length($repeat_consensus);
      }

      push @consensi, Bio::EnsEMBL::RepeatConsensus->new
          (-DBID => $repeat_consensus_id,
           -NAME => $repeat_name,
           -REPEAT_CLASS => $repeat_class,
           -REPEAT_TYPE => $repeat_type,
           -LENGTH => $repeat_length,
           -ADAPTOR => $self,
           -REPEAT_CONSENSUS => $repeat_consensus);
    }
    return \@consensi;
}


=head2 store

  Arg [1]    : list of Bio::EnsEMBL::RepeatConsensus @consensi
  Example    : $repeat_consensus_adaptor->store(@consensi);
  Description: stores a list of RepeatConsensus objects in the database
  Returntype : none
  Exceptions : none
  Caller     : ?
  Status     : Stable

=cut

sub store {
  my( $self, @consensi ) = @_;

  my $sth = $self->prepare(q{
    INSERT into repeat_consensus( repeat_consensus_id
          , repeat_name
          , repeat_class
          , repeat_type
          , repeat_consensus )
      VALUES (NULL, ?,?,?,?)
    });

  foreach my $rc (@consensi) {
    my $name  = $rc->name
      or throw("name not set");
    my $class = $rc->repeat_class
      or throw("repeat_class not set");
    my $type  = $rc->repeat_type();
    $type = "" unless defined $type;
    my $seq   = $rc->repeat_consensus
      or throw("repeat_consensus not set");

    $sth->bind_param(1,$name,SQL_VARCHAR);
    $sth->bind_param(2,$class,SQL_VARCHAR);
    $sth->bind_param(3,$type,SQL_VARCHAR);
    $sth->bind_param(4,$seq,SQL_LONGVARCHAR);

    $sth->execute();

    my $db_id = $self->last_insert_id('repeat_consensus_id', undef, 'repeat_consensus')
    or throw("Didn't get an insertid from the INSERT statement");

    $rc->dbID($db_id);
    $rc->adaptor($self);
  }
}

1;

__END__

