### Bio::EnsEMBL::Lite::RepeatFeatureAdaptor
#
# Copyright EMBL-EBI 2001
#
# Author: Graham McVicker
# 
# Date : 12.08.2002
#

=head1 NAME

Bio::EnsEMBL::Lite::RepeatFeatureAdaptor

=head1 SYNOPSIS

A repeat feature adaptor for the lite database which is designed to
obtain lite-weight repeat feature objects from the denormalized, 
locationally-indexed Lite database.  

=head1 CONTACT

Post general questions to ensembl-dev mailing list: <ensembl-dev@ebi.ac.uk>

=head1 APPENDIX

=cut

use strict;


package Bio::EnsEMBL::Lite::RepeatFeatureAdaptor;

use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::RepeatFeature;
use Bio::EnsEMBL::RepeatConsensus;


use vars '@ISA';

@ISA = qw(Bio::EnsEMBL::DBSQL::BaseAdaptor);

my $MAX_REPEAT_LENGTH = 100000;


=head2 fetch_by_Slice_and_type

  Arg [1]    : Bio::EnsEMBL::Slice $slice
  Arg [2]    : string $type
  Example    : @rfs = $rf_adaptor->fetch_by_Slice_and_type($slice, 'Alu'); 
  Description: Obtains a list of repeat features of a given repeat type on
               the region defined by slice
  Returntype : listreference of Bio::EnsEMBL::RepeatFeatures
  Exceptions : none
  Caller     : contigview

=cut

sub fetch_all_by_Slice_and_type {
  my ($self, $slice, $type) = @_;
 
  my $sth = $self->prepare(
	"SELECT r.id, r.hid,  r.chr_name, r.chr_start, r.chr_end, r.chr_strand
         FROM   repeat r
         WHERE  r.chr_name = ? AND r.chr_start <= ? AND r.chr_start >= ? 
                AND r.chr_end >= ? AND r.type = ?" );
 
  $sth->execute( $slice->chr_name(), $slice->chr_end(), 
		 $slice->chr_start() - $MAX_REPEAT_LENGTH, 
		 $slice->chr_start(), $type );

  return $self->_objs_from_sth($sth, $slice);
}


=head2 fetch_all_by_Slice

  Arg [1]    : Bio::EnsEMBL::Slice $slice 
  Example    : $repeats = $lite_repeat_apaptor->fetch_all_by_Slice($slice);
  Description: Retrieves a list of RepeatFeatures from the lite database.
  Returntype : Bio::EnsEMBL::RepeatFeature
  Exceptions : none
  Caller     : ProxyRepeatFeatureAdaptor (via Slice via contigview)

=cut

sub fetch_all_by_Slice {
    my ( $self, $slice ) = @_;
    
    my $sth = $self->prepare(
	"SELECT r.id, r.hid,  r.chr_name, r.chr_start, r.chr_end, r.chr_strand
         FROM   repeat r
         WHERE  r.chr_name = ? AND r.chr_start <= ? AND r.chr_start >= ? 
                AND r.chr_end >= ?" );
 
    $sth->execute( $slice->chr_name(), $slice->chr_end(), 
		   $slice->chr_start() - $MAX_REPEAT_LENGTH, 
		   $slice->chr_start() );
    
    return $self->_objs_from_sth($sth, $slice);
}



=head2 _objs_from_sth

  Arg [1]    : DBI:st $sth 
               An executed DBI statement handle
  Arg [2]    : Bio::EnsEMBL::Slice $slice
               The slice to map the coordinates to
  Example    : $repeats = $self->_objs_from_sth($sth, $slice);
  Description: Creates a list of repeat objects from an executed DBI statement
               handle.
  Returntype : list reference to Bio::EnsEMBL::RepeatFeatures
  Exceptions : none
  Caller     : internal

=cut


sub _objs_from_sth {
  my ($self, $sth, $slice) = @_;

  my @repeats = ();
  my($rc, $core, $repeat_adaptor, $rc_adaptor);
  my $slice_start = $slice->chr_start;
  my $slice_end   = $slice->chr_end;
  my $slice_strand = $slice->strand;
  my %rc_hash;

  if($core = $self->db->get_db_adaptor('core')) {
    $repeat_adaptor = $core->get_RepeatFeatureAdaptor();
    $rc_adaptor = $core->get_RepeatConsensusAdaptor();
  } else {
    $repeat_adaptor = $self;
    $self->warn("Core Database not attached to lite database.  Not able to 
                   Retrieve repeat consensus adaptor for repeat features\n");
  }

  my ($id, $hid, $chr_name, $start, $end, $strand);
  $sth->bind_columns(\$id, \$hid, \$chr_name, \$start, \$end, \$strand);

  my ($feat_start, $feat_end, $feat_strand);
  
  while($sth->fetch) {
    #skip features that are entirely outsied the slice area
    next if ($start > $slice_end || $end < $slice_start);

    #cache repeat_consensi to reduce object construction overhead
    unless($rc = $rc_hash{$hid}) {
      $rc = new Bio::EnsEMBL::RepeatConsensus();
      $rc->name($hid);
      if($rc_adaptor) {
	$rc->adaptor($rc_adaptor);
      }
      $rc_hash{$hid} = $rc;
    }

    #convert chromosomal coordinates to slice coordinates
    if($slice_strand == 1) {
      $feat_start  = $start - $slice_start + 1;
      $feat_end    = $end   - $slice_start   + 1;
      $feat_strand = $strand;
    } else {
      $feat_start  = $slice_end - $end   + 1;
      $feat_end    = $slice_end - $start + 1;
      $feat_strand = $strand * -1;
    }

    #create partially filled repeat object using fast (hacky) constructor 
    push @repeats, Bio::EnsEMBL::RepeatFeature->new_fast(
			   { '_gsf_tag_hash'  => {},
			     '_gsf_sub_array' => [],
			     '_parse_h'       => {},
                             '_gsf_start'         => $feat_start,
                             '_gsf_end'           => $feat_end,
			     '_gsf_strand'        => $feat_strand,
			     '_adaptor'       => $repeat_adaptor,
			     '_repeat_consensus' => $rc,
			     '_db_id'         => $id } );
    }

  return \@repeats; 
}


=head2 fetch_by_Slice_and_type

  Arg [1]    : none
  Example    : none
  Description: DEPRECATED use fetch_all_by_Slice_and_type instead
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub fetch_by_Slice_and_type {
  my ($self, @args) = @_;

  $self->warn("fetch_by_Slice_and_type has been renamed fetch_all_by_Slice_and_type\n" . caller);

  return $self->fetch_all_by_Slice_and_type(@args);
}


=head2 fetch_by_Slice

  Arg [1]    : none
  Example    : none
  Description: DEPRECATED use fetch_all_by_Slice instead
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub fetch_by_Slice {
  my ($self, @args) = @_;

  $self->warn("fetch_by_Slice has been renamed fetch_all_by_Slice\n" . caller);

  return $self->fetch_all_by_Slice(@args);
}



1;

__END__


