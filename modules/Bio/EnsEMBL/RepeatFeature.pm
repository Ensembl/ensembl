
### Bio::EnsEMBL::RepeatFeature

package Bio::EnsEMBL::RepeatFeature;

use strict;
use Bio::EnsEMBL::SeqFeature;

use vars '@ISA';

@ISA = qw{ Bio::EnsEMBL::SeqFeature Bio::LocationI };


# new() comes from SeqFeature::new()


# Unique to RepeatFeature

sub score {
    my( $self, $score ) = @_;
    
    if (defined $score) {
        $self->{'_score'} = $score;
    }
    return $self->{'_score'};
}

sub repeat_id {
    my( $self, $repeat_id ) = @_;
    
    if ($repeat_id) {
        $self->{'_repeat_id'} = $repeat_id;
    }
    return $self->{'_repeat_id'};
}

sub contig_id {
    my( $self, $contig_id ) = @_;
    
    if ($contig_id) {
        $self->{'_contig_id'} = $contig_id;
    }
    return $self->{'_contig_id'};
}

sub repeat_consensus_adaptor {
    my( $self, $rca ) = @_;
    
    if ($rca) {
        unless (ref($rca) and $rca->isa('Bio::EnsEMBL::DBSQL::RepeatConsensusAdaptor')) {
            $self->throw("Not a 'Bio::EnsEMBL::DBSQL::RepeatConsenusAdpator': $rca")
        }
        $self->{'_repeat_consensus_adaptor'} = $rca;
    }
    return $self->{'_repeat_consensus_adaptor'};
}

sub repeat_consensus {
    my( $self, $con ) = @_;
    
    if (defined $con) {
	$self->throw("$con is not a RepeatConsensus")
	 unless $con->isa("Bio::EnsEMBL::RepeatConsensus");
        $self->{'_repeat_consensus'} = $con;
    }
    return $self->{'_repeat_consensus'};
}

sub get_RepeatConsensus {
    my( $self ) = @_;
    
    my $repeat_id = $self->repeat_id;
    return $self->repeat_consensus_adaptor->fetch_by_dbID($repeat_id);
}

# Bio::EnsEMBL::SeqFeatureI methods

# Should this be dbID?  -- I'm implementing Bio::EnsEMBL::SeqFeatureI
sub id {
    my $self = shift;
    
    $self->warn("Delegating id call to dbID");
    
    $self->dbID(@_);
}

sub dbID {
    my( $self, $db_id ) = @_;
    
    if ($db_id) {
        $self->{'_db_id'} = $db_id;
    }
    return $self->{'_db_id'};
}

sub analysis {
    my( $self, $analysis ) = @_;
    
    if ($analysis) {
        $self->{'_analysis'} = $analysis;
    }
    return $self->{'_analysis'};
}

sub validate {

}

# Should these return anything?
sub       phase { return -1  }
sub   end_phase { return -1  }
sub      hphase { return -1  }
sub  hend_phase { return -1  }
sub  percent_id { return 100 }
sub hpercent_id { return 100 }
sub     e_value { return 0   }
sub     p_value { return 0   }
sub    hp_value { return 0   }

# LocationI methods

sub start {
    my( $self, $start ) = @_;
    
    if ($start) {
        $self->{'_start'} = $start;
    }
    return $self->{'_start'};
}

sub end {
    my( $self, $end ) = @_;
    
    if ($end) {
        $self->{'_end'} = $end;
    }
    return $self->{'_end'};
}

sub strand {
    my( $self, $strand ) = @_;
    
    if ($strand) {
        $self->{'_strand'} = $strand;
    }
    return $self->{'_strand'};
}

sub hstart {
    my( $self, $hstart ) = @_;
    
    if ($hstart) {
        $self->{'_hstart'} = $hstart;
    }
    return $self->{'_hstart'};
}

sub hend {
    my( $self, $hend ) = @_;
    
    if ($hend) {
        $self->{'_hend'} = $hend;
    }
    return $self->{'_hend'};
}

sub hstrand { return 1 }

sub location_type { return 'EXACT' }

sub min_start { return };
sub min_end   { return };
sub max_start { return };
sub max_end   { return };

sub start_pos_type { return 'EXACT' }
sub   end_pos_type { return 'EXACT' }

sub to_FTString {
    my( $self ) = @_;
    
    my $start  = $self->start or $self->throw("start not set");
    my $end    = $self->end   or $self->throw("end not set");
    my $strand = $self->strand;
    $self->throw("strand not set") unless defined($strand);
    
    my $loc = "$start..$end";
    if ($strand == -1) {
        $loc = "complement($loc)";
    }
}



1;

__END__

=head1 NAME - Bio::EnsEMBL::RepeatFeature

=head1 AUTHOR

James Gilbert B<email> jgrg@sanger.ac.uk

