
### Bio::EnsEMBL::DBSQL::RepeatFeatureAdaptor

package Bio::EnsEMBL::DBSQL::RepeatFeatureAdaptor;

use strict;
use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::RepeatFeature;
use vars qw(@ISA);

@ISA = ('Bio::EnsEMBL::DBSQL::BaseAdaptor');

sub fetch_for_contig {
    my( $self, $contig ) = @_;
    
    my $contig_id = $contig->dbID;
    my @repeats = $self->_generic_fetch(
        qq{ AND f.contig_id = $contig_id }
        );
    foreach my $r (@repeats) {
        $r->attach_seq($contig);
    }
}

sub _generic_fetch {
    my( $self, $where_clause ) = @_;
    
    my $sth = $self->prepare(qq{
        SELECT f.repeat_feature_id
          , f.contig_id
          , f.contig_start
          , f.contig_end
          , f.contig_strand
          , f.repeat_id
          , f.repeat_start
          , f.repeat_end
          , f.analysis_id
        FROM repeat_feature f
          , repeat r
        WHERE f.repeat_id = r.repeat_id
        }. $where_clause);
    $sth->execute;
    
    my( @repeats );
    while (my @data = $sth->fetchrow_arrayref) {
        my $r = Bio::EnsEMBL::RepeatFeature->new(@data);
        push(@repeats, $r);
    }
    return( @repeats );
}

1;

__END__

=head1 NAME - Bio::EnsEMBL::DBSQL::RepeatFeatureAdaptor

=head1 AUTHOR

James Gilbert B<email> jgrg@sanger.ac.uk

