
### Bio::EnsEMBL::DBSQL::RepeatConsensusAdaptor

package Bio::EnsEMBL::DBSQL::RepeatConsensusAdaptor;

use strict;
use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::RepeatConsensus;
use vars qw(@ISA);

@ISA = ('Bio::EnsEMBL::DBSQL::BaseAdaptor');

sub fetch_by_dbID {
    my( $self, $db_id ) = @_;

    return $self->_generic_fetch(
        "repeat_id = $db_id"
        );   
}

sub fetch_by_name {
    my( $self, $name ) = @_;

    return $self->_generic_fetch(
        "repeat_name = '$name'"
        );   
}

sub _generic_fetch {
    my( $self, $where_clause ) = @_;
    
    my( $repeat_id,
        $repeat_name,
        $repeat_class,
        $repeat_length,
        );
    
    my $sth = $self->prepare(qq{
        SELECT repeat_id
          , repeat_name
          , repeat_class
          , LENGTH(repeat_consensus)
        FROM repeat_consensus
        WHERE }. $where_clause);
    $sth->execute;
    $sth->bind_columns(
        \$repeat_id,
        \$repeat_name,
        \$repeat_class,
        \$repeat_length,
        );
    
    my( @consensi );
    while ($sth->fetch) {
        my $rc = Bio::EnsEMBL::RepeatConsensus->new;
        $rc->dbID($repeat_id);
        $rc->name($repeat_name);
        $rc->repeat_class($repeat_class);
        $rc->length($repeat_length);
        $rc->adaptor($self);
        push(@consensi, $rc);
    }
    return @consensi;
}

sub fetch_seq_string_for_dbID {
    my( $self, $db_id ) = @_;
    
    my $sth = $self->prepare(qq{
        SELECT repeat_consensus
        FROM repeat_consensus
        WHERE repeat_id = $db_id
        });
    $sth->execute;
    
    my ($seq) = $sth->fetchrow  
        or $self->throw("Can't fetch repeat_consensus for repeat_id = '$db_id'");
    return $seq;
}

sub store {
    my( $self, @consensi ) = @_;
    
    my $sth = $self->prepare(q{
        INSERT into repeat_consensus( repeat_id
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

=head1 NAME - Bio::EnsEMBL::DBSQL::RepeatConsensusAdaptor

=head1 AUTHOR

James Gilbert B<email> jgrg@sanger.ac.uk

