
### Bio::EnsEMBL::RepeatFeature

package Bio::EnsEMBL::RepeatFeature;

use strict;
use Bio::EnsEMBL::SeqFeature;

use vars '@ISA';

@ISA = qw{ Bio::EnsEMBL::SeqFeature };


=head2 new_fast

  Arg [1]    : hash reference $hashref
  Example    : none
  Description: This is an ultra fast constructor which requires knowledge of
               the objects internals to be used.  It is only used by 
               RepeatFeatureAdaptors (when thousands of repeats need to be
               quickly created).  The SeqFeature superclass constructor 'new'
               should be used in most instances.
  Returntype : Bio::EnsEMBL::RepeatFeature
  Exceptions : none
  Caller     : RepeatFeatureAdaptors

=cut

sub new_fast {
  my ($class, $hashref) = @_;
  
  return bless $hashref, $class;
}



=head2 adaptor

  Arg [1]    : (optional) Bio::EnsEMBL::DBSQL::RepeatFeatureAdaptor $adaptor
  Example    : $adaptor = $repeat->adaptor;
  Description: The adaptor which performs database requests for this object.
               This should be set when the object is stored in the database or
               retrieved from the database.
  Returntype : Bio::EnsEMBL::DBSQL::RepeatFeatureAdaptor
  Exceptions : none
  Caller     : general

=cut

sub adaptor {
  my ($self, $adaptor) = @_;
  
  if(defined $adaptor) {
    $self->{'_adaptor'} = $adaptor;
  }
  
  return $self->{'_adaptor'};
}



=head2 repeat_consensus

  Arg [1]    : (optional) Bio::EnsEMBL::RepeatConsensus
  Example    : $repeat_consensus = $repeat->repeat_consensus;
  Description: Getter/Setter for the repeat consensus of this repeat
  Returntype : Bio::EnsEMBL::RepeatConsensus
  Exceptions : none
  Caller     : general

=cut

sub repeat_consensus {
  my( $self, $con ) = @_;
    
  if (defined $con) {
    $self->throw("$con is not a RepeatConsensus")
      unless $con->isa("Bio::EnsEMBL::RepeatConsensus");
    $self->{'_repeat_consensus'} = $con;
  }
  return $self->{'_repeat_consensus'};
}



=head2 dbID

  Arg [1]    : (optional) $db_id
  Example    : $dbID = $repeat->dbID
  Description: getter/setter for this objects internal database identifier
  Returntype : int
  Exceptions : none
  Caller     : general

=cut

sub dbID {
    my( $self, $db_id ) = @_;
    
    if (defined $db_id) {
        $self->{'_db_id'} = $db_id;
    }
    return $self->{'_db_id'};
}


=head2 analysis

  Arg [1]    : (optional) Bio::EnsEMBL::Analysis
  Example    : $analysis = $repeat_feat->analysis;
  Description: Getter/Setter for the analysis that was used to generate this
               object
  Returntype : Bio::EnsEMBL::Analysis
  Exceptions : none
  Caller     : general

=cut

sub analysis {
    my( $self, $analysis ) = @_;
    
    if ($analysis) {
        $self->{'_analysis'} = $analysis;
    }
    return $self->{'_analysis'};
}


=head2 hstart

  Arg [1]    : (optional) int $hstart
  Example    : $hit_start = $repeat->hstart;
  Description: Getter/Setter for the start bp of this repeat match on the 
               consensus sequence.
  Returntype : int
  Exceptions : none 
  Caller     : general

=cut

sub hstart {
    my( $self, $hstart ) = @_;
    
    if ($hstart) {
        $self->{'_hstart'} = $hstart;
    }
    return $self->{'_hstart'};
}



=head2 hend

  Arg [1]    : (optional) int $hend
  Example    : $hit_end = $repeat->hend;
  Description: Getter/Setter for the end bp of this repeat match on the 
               consensus sequence.
  Returntype : int
  Exceptions : none
  Caller     : general

=cut

sub hend {
    my( $self, $hend ) = @_;
    
    if ($hend) {
        $self->{'_hend'} = $hend;
    }
    return $self->{'_hend'};
}



=head2 hstrand

  Arg [1]    : none 
  Example    : none
  Description: always returns 1. method exists for consistancy with other 
               features.
  Returntype : int
  Exceptions : none
  Caller     : 

=cut

sub hstrand { 
  return 1; 
}


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



=head2 repeat_consensus_id

  Arg [1]    : none
  Example    : none
  Description: DEPRECATED use repeat_consensus->dbID instead
  Returntype : none
  Exceptions : none
  Caller     : none

=cut


sub repeat_consensus_id {
  my( $self, $repeat_consensus_id ) = @_;

    my ($f,$p, $l) = caller;
    $self->warn("repeat_consensus_id is deprecated, use repeat_consensus ".
		"instead. caller = $f, $p, $l"); 
  
  if ($repeat_consensus_id) {
    $self->{'_repeat_consensus_id'} = $repeat_consensus_id;
  }
  return $self->{'_repeat_consensus_id'};
}




=head2 get_RepeatConsensus

  Arg [1]    : none
  Example    : none
  Description: DEPRECATED use repeat_consensus instead
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub get_RepeatConsensus {
    my( $self ) = @_;

    my ($f,$p, $l) = caller;
    $self->warn("get_RepeatConsensus is deprecated, use repeat_consensus ".
		"instead. caller = $f, $p, $l"); 
    
    my $rc_id = $self->repeat_consensus_id;
    return $self->adaptor->db->get_RepeatConsensusAdaptor
      ->fetch_by_dbID($rc_id);
}



1;

__END__

=head1 NAME - Bio::EnsEMBL::RepeatFeature

=head1 AUTHOR

James Gilbert B<email> jgrg@sanger.ac.uk

