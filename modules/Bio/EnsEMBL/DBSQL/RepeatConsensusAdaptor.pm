#
# EnsEMBL module for Bio::EnsEMBL::DBSQL::RepeatFeatureAdaptor
#
# Copyright EMBL/EBI
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::DBSQL::RepeatConsensusAdaptor

=head1 SYNOPSIS

$repeat_consensus_adaptor = $database_adaptor->get_RepeatConsensusAdaptor();

=head1 DESCRIPTION

This is an adaptor for the retrieval and storage of RepeatConsensus objects.

=head1 AUTHOR - James Gilbert

Email jgrg@ebi.ac.uk

=head1 CONTACT

Arne Stabenau - stabenau@ebi.ac.uk
Graham McVicker - mcvicker@ebi.ac.uk
Ewan Birney - birney@ebi.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

package Bio::EnsEMBL::DBSQL::RepeatConsensusAdaptor;

use strict;
use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::RepeatConsensus;
use vars qw(@ISA);

@ISA = ('Bio::EnsEMBL::DBSQL::BaseAdaptor');


=head2 fetch_by_dbID

  Arg [1]    : int $db_id
               The database identifier for the RepeatConsensus to obtain
  Example    : $repeat_consensus = $repeat_consensus_adaptor->fetch_by_dbID(4);
  Description: Obtains a RepeatConsensus object from the database via its
               primary key. 
  Returntype : list of Bio::EnsEMBL::RepeatConsensus
  Exceptions : none
  Caller     : general, Bio::EnsEMBL::RepeatFeatureAdaptor

=cut

sub fetch_by_dbID {
    my( $self, $db_id ) = @_;

    return $self->_generic_fetch("repeat_consensus_id = $db_id");   
}



=head2 fetch_by_name

  Arg [1]    : string $name
               the name of the repeat consensus to obtain
  Example    : $rc = $repeat_consensus_adaptor->fetch_by_name('AluSx');
  Description: Obtains a repeat consensus from the database via its name
  Returntype : list of Bio::EnsEMBL::RepeatConsensus
  Exceptions : none
  Caller     : general

=cut

sub fetch_by_name {
    my( $self, $name ) = @_;

    return $self->_generic_fetch("repeat_name = '$name'");   
}


=head2 _generic_fetch

  Arg [1]    : string $where_clause
  Example    : none
  Description: PRIVATE used to create RepeatConsensus features from an 
               SQL constraint
  Returntype : list of Bio::EnsEMBL::RepeatConsensus
  Exceptions : none
  Caller     : internal

=cut

sub _generic_fetch {
    my( $self, $where_clause ) = @_;
    
    my( $repeat_consensus_id,
        $repeat_name,
        $repeat_class,
        $repeat_length,
        );
    
    my $sth = $self->prepare(qq{
        SELECT repeat_consensus_id
          , repeat_name
          , repeat_class
          , LENGTH(repeat_consensus)
        FROM repeat_consensus
        WHERE }. $where_clause);
    $sth->execute;
    $sth->bind_columns(
        \$repeat_consensus_id,
        \$repeat_name,
        \$repeat_class,
        \$repeat_length,
        );
    
    my( @consensi );
    while ($sth->fetch) {
        my $rc = Bio::EnsEMBL::RepeatConsensus->new;
        $rc->dbID($repeat_consensus_id);
        $rc->name($repeat_name);
        $rc->repeat_class($repeat_class);
        $rc->length($repeat_length);
        $rc->adaptor($self);
        push(@consensi, $rc);
    }
    return @consensi;
}


=head2 fetch_seq_string_for_dbID

  Arg [1]    : int $db_id
  Example    : none
  Description: Should probably be deprecated or renamed since it does not 
               return an object and it is a fetch method - not consistent.
               Retrieves the repeat_consensus string of a feature.  Looking
               at the database this always apears to be 'N' anyway.  
  Returntype : string
  Exceptions : thrown if the RepeatConsensus with $db_id cannot be found
  Caller     : ?

=cut

sub fetch_seq_string_for_dbID {
    my( $self, $db_id ) = @_;
    
    my $sth = $self->prepare(qq{
        SELECT repeat_consensus
        FROM repeat_consensus
        WHERE repeat_consensus_id = $db_id
        });
    $sth->execute;
    
    my ($seq) = $sth->fetchrow  
        or $self->throw("Can't fetch repeat_consensus for repeat_consensus_id = '$db_id'");
    return $seq;
}


=head2 store

  Arg [1]    : list of Bio::EnsEMBL::RepeatConsensus @consensi
  Example    : $repeat_consensus_adaptor->store(@consensi);
  Description: stores a list of RepeatConsensus objects in the database
  Returntype : none
  Exceptions : none
  Caller     : ?

=cut

sub store {
  my( $self, @consensi ) = @_;
  
  my $sth = $self->prepare(q{
    INSERT into repeat_consensus( repeat_consensus_id
          , repeat_name
          , repeat_class
          , repeat_consensus )
      VALUES (NULL, ?,?,?)
    });
    
  foreach my $rc (@consensi) {
    my $name  = $rc->name
      or $self->throw("name not set");
    my $class = $rc->repeat_class
      or $self->throw("repeat_class not set");
    my $seq   = $rc->repeat_consensus
      or $self->throw("repeat_consensus not set");
    
    $sth->execute($name, $class, $seq);
    
    my $db_id = $sth->{'mysql_insertid'}
    or $self->throw("Didn't get an insertid from the INSERT statement");
    
    $rc->dbID($db_id);
    $rc->adaptor($self);
  }
}

1;

__END__

