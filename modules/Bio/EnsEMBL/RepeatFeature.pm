
### Bio::EnsEMBL::RepeatFeature

package Bio::EnsEMBL::RepeatFeature;

use strict;
use Bio::EnsEMBL::SeqFeature;
use Bio::LocationI;

use vars '@ISA';

@ISA = qw{ Bio::EnsEMBL::SeqFeature Bio::LocationI };


#ultra fast hacky constructor for rapid feature creation
sub new_fast {
  my ($class, $hashref) = @_;

  return bless $hashref, $class;
}


sub repeat_consensus_id {
    my( $self, $repeat_consensus_id ) = @_;
    
    if ($repeat_consensus_id) {
        $self->{'_repeat_consensus_id'} = $repeat_consensus_id;
    }
    return $self->{'_repeat_consensus_id'};
}


sub adaptor {
  my ($self, $adaptor) = @_;

  if(defined $adaptor) {
    $self->{'_adaptor'} = $adaptor;
  }

  return $self->{'_adaptor'};
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
    
    my $rc_id = $self->repeat_consensus_id;
    return $self->adaptor->db->get_RepeatConsensusAdaptor
      ->fetch_by_dbID($rc_id);
}


sub dbID {
    my( $self, $db_id ) = @_;
    
    if (defined $db_id) {
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

