
### Bio::EnsEMBL::RepeatConsensus

package Bio::EnsEMBL::RepeatConsensus;

use strict;
use Bio::PrimarySeqI;
use Bio::EnsEMBL::Root;
use vars '@ISA';

@ISA = qw(Bio::EnsEMBL::Root Bio::PrimarySeqI);

sub new {
    my( $pkg ) = @_;
    
    return bless {}, $pkg;
}

sub moltype  { return 'dna' };
sub alphabet { return 'dna' };

sub dbID {
    my( $self, $db_id ) = @_;
    
    if ($db_id) {
        $self->{'_db_id'} = $db_id;
    }
    return $self->{'_db_id'};
}

# Alias primary_id method to dbID
*primary_id = \&dbID;

sub name {
    my( $self, $name ) = @_;
    
    if ($name) {
        $self->{'_name'} = $name;
    }
    return $self->{'_name'};
}

# Alias display_id and accession_number methods to name
*display_id       = \&name;
*accession_number = \&name;

sub length {
    my( $self, $length ) = @_;
    
    if ($length) {
        $self->{'_length'} = $length;
    }
    return $self->{'_length'};
}

sub repeat_class {
    my( $self, $repeat_class ) = @_;
    
    if ($repeat_class) {
        $self->{'_repeat_class'} = $repeat_class;
    }
    return $self->{'_repeat_class'};
}

sub desc {
    my( $self ) = @_;
    
    my $class = $self->repeat_class or return;
    return "class=$class";
}

sub adaptor {
    my( $self, $adaptor ) = @_;
    
    if ($adaptor) {
        $self->{'_adaptor'} = $adaptor;
    }
    return $self->{'_adaptor'};
}

sub repeat_consensus {
    my( $self, $repeat_consensus ) = @_;
    
    if ($repeat_consensus) {
        $self->{'_repeat_consensus'} = $repeat_consensus;
    }
    return $self->{'_repeat_consensus'};
}

sub seq {
    my( $self ) = @_;
    
    # Stored in the object?
    if (my $rc = $self->repeat_consensus) {
        return $rc;
    }
    
    # No, in the database.
    my $db_id = $self->dbID
        or $self->throw("No dbID - is the consensus stored in the database?");
    return $self->adaptor->fetch_seq_string_for_dbID($db_id);
}

1;

__END__

=head1 NAME - Bio::EnsEMBL::RepeatConsensus

=head1 DESCRIPTION

This object represents an entry in the
repeat_consensus table.

It can contain the consensus sequence for a
repeat such as a particular Alu, or "cag" for a
simple triplet repeat.

=head1 AUTHOR

James Gilbert B<email> jgrg@sanger.ac.uk

