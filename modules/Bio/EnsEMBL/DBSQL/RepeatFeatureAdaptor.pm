
### Bio::EnsEMBL::DBSQL::RepeatFeatureAdaptor

package Bio::EnsEMBL::DBSQL::RepeatFeatureAdaptor;

use strict;
use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::RepeatFeature;
use vars qw(@ISA);

@ISA = ('Bio::EnsEMBL::DBSQL::BaseAdaptor');

sub fetch_by_RawContig {
    my( $self, $contig ) = @_;

    my @repeats = $self->fetch_by_contig_id($contig->dbID);
    foreach my $r (@repeats) {
        $r->attach_seq($contig);
    }
    return @repeats;
}

sub fetch_by_contig_id {
    my( $self, $contig_id ) = @_;

    return $self->_generic_fetch(
        qq{ contig_id = $contig_id }
        );
}

sub fetch_by_dbID {
    my( $self, $db_id ) = @_;

    my ($rf) = $self->_generic_fetch(
        qq{ repeat_feature_id = $db_id }
        );
    return $rf;
}

sub _generic_fetch {
    my( $self, $where_clause ) = @_;

    my( $repeat_feature_id,
        $contig_id,
        $contig_start,
        $contig_end,
        $contig_strand,
        $repeat_id,
        $repeat_start,
        $repeat_end,
        $analysis_id,
        );

    my $sth = $self->prepare(qq{
        SELECT repeat_feature_id
          , contig_id
          , contig_start
          , contig_end
          , contig_strand
          , repeat_id
          , repeat_start
          , repeat_end
          , analysis_id
        FROM repeat_feature f
        WHERE }. $where_clause);

    $sth->execute;
    $sth->bind_columns(
        \$repeat_feature_id,
        \$contig_id,
        \$contig_start,
        \$contig_end,
        \$contig_strand,
        \$repeat_id,
        \$repeat_start,
        \$repeat_end,
        \$analysis_id,
        );

    my $rca = $self->db->get_RepeatConsensusAdaptor;
    my $aa  = $self->db->get_AnalysisAdaptor;

    my( @repeats, %analysis_cache );
    while ($sth->fetch) {
        # new in RepeatFeature takes no arguments
        my $r = Bio::EnsEMBL::RepeatFeature->new;
        $r->dbID($repeat_feature_id);

        # So RepeatFeature can get its repeat
        $r->repeat_consensus_adaptor($rca);
        $r->repeat_id($repeat_id);

        $r->contig_id( $contig_id     );
        $r->start    ( $contig_start  );
        $r->end      ( $contig_end    );
        $r->strand   ( $contig_strand );
        $r->hstart   ( $repeat_start  );
        $r->hend     ( $repeat_end    );

        my( $ana_obj );
        unless ($ana_obj = $analysis_cache{$analysis_id}) {
            $ana_obj = $aa->fetch_by_dbID($analysis_id)
                or $self->throw("No analysis object for ID '$analysis_id'");
            $analysis_cache{$analysis_id} = $ana_obj;
        }
        $r->analysis($ana_obj);

        push(@repeats, $r);
    }
    return( @repeats );
}

sub store {
    my( $self, $contig_id, @repeats ) = @_;

    my $rca = $self->db->get_RepeatConsensusAdaptor;
    my ($cons, $db_id);

    $self->throw("Can't store repeats without a contig_id (got '$contig_id')")
        unless $contig_id =~ /^\d+$/;

    my $sth = $self->prepare(qq{
        INSERT into repeat_feature( repeat_feature_id
          , contig_id
          , contig_start
          , contig_end
          , contig_strand
          , repeat_id
          , repeat_start
          , repeat_end
	  , score
          , analysis_id )
        VALUES(NULL, ?,?,?,?,?,?,?,?,?)
        });
    foreach my $rf (@repeats) {

	# must have a consensus attached

        $self->throw("Must have a RepeatConsensus attached")
	 unless defined ($cons = $rf->repeat_consensus);

	# for tandem repeats - simply store consensus and repeat
	# one pair per hit. don't need to check consensi stored
	# already. consensus has name and class set to 'trf'

	if ($cons->name eq 'trf') {
	    $rca->store($cons);
	    $rf->repeat_id($cons->dbID);

	} else {

	# for other repeats - need to see if a consensus is stored already

            unless ($rf->repeat_id) {

	        # need to get the consensus seq object for this repeat

                unless ($cons->dbID) {
	            my @match = ($rca->fetch_by_name($cons->name));

		    if (@match > 1) {
		        $self->warn(@match . " consensi for " . $cons->name . "\n");
		    } elsif (@match == 0) {
		    # if we don't match a consensus already stored
		    # create a fake one
		    # set consensus to 'N' as null seq not allowed
		    # FIXME: not happy with this, but ho hum ...
		        $self->warn("Can't find " . $cons->name . "\n");
		        $cons->repeat_consensus("N");
		        $rca->store($cons);
		    }
	            $db_id = ($rca->fetch_by_name($cons->name))[0]->dbID;

	            $cons->dbID($db_id);
	        }
	        $rf->repeat_id($cons->dbID);
	    }
	}

        $sth->execute(
            $contig_id,
            $rf->start,
            $rf->end,
            $rf->strand,
            $rf->repeat_id,
            $rf->hstart,
            $rf->hend,
	    $rf->score,
            $rf->analysis->dbID,
            );
        my $db_id = $sth->{'mysql_insertid'}
            or $self->throw("Didn't get an insertid from the INSERT statement");
        $rf->dbID($db_id);
    }
}

1;

__END__

=head1 NAME - Bio::EnsEMBL::DBSQL::RepeatFeatureAdaptor

=head1 AUTHOR

James Gilbert B<email> jgrg@sanger.ac.uk

