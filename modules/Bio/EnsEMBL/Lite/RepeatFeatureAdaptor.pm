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

sub fetch_by_Slice_and_type {
  my ($self, $slice, $type) = @_;
 
  my $sth = $self->prepare(
	"SELECT r.id, r.hid,  r.chr_name, r.chr_start, r.chr_end, r.chr_strand
         FROM   repeat r
         WHERE  r.chr_name = ? AND r.chr_start <= ? AND r.chr_start >= ? 
                AND r.chr_end >= ? AND r.type = ?" );
 
  $sth->execute( $slice->chr_name(), $slice->chr_end(), 
		 $slice->chr_start() - $MAX_REPEAT_LENGTH, 
		 $slice->chr_start(), $type );

  my @repeats = ();

  my ($rc, $core, $repeat_adaptor, $rc_adaptor);
  my @repeats;
  my($rc, $core, $repeat_adaptor, $rc_adaptor);
  my $slice_start = $slice->chr_start + 1;
  my %rc_hash;

  if($core = $self->db()->get_db_adaptor('core')) {
    $repeat_adaptor = $core->get_RepeatFeatureAdaptor();
    $rc_adaptor = $core->get_RepeatConsensusAdaptor();
  } else {
    $repeat_adaptor = $self;
    $self->warn("Core Database not attached to lite database.  Not able to 
                 Retrieve repeat consensus adaptor for repeat features\n");
  }

  my ($id, $hid, $chr_name, $start, $end, $strand);
  $sth->bind_columns(\$id, \$hid, \$chr_name, \$start, \$end, \$strand);

  while($sth->fetch) {
    #cache repeat_consensi to reduce object construction overhead
    unless($rc = $rc_hash{$hid}) {
      $rc = new Bio::EnsEMBL::RepeatConsensus();
      $rc->name($hid);
      if($rc_adaptor) {
	$rc->adaptor($rc_adaptor);
      }
      $rc_hash{$hid} = $rc;
    }

    #create partially filled repeat object using fast (hacky) constructor 
    push @repeats, Bio::EnsEMBL::RepeatFeature->new_fast(
			   { '_gsf_tag_hash'  => {},
			     '_gsf_sub_array' => [],
			     '_parse_h'       => {},
                             '_gsf_start'         => $start - $slice_start,
                             '_gsf_end'           => $end - $slice_start,
			     '_gsf_strand'        => $strand,
			     '_adaptor'       => $repeat_adaptor,
			     '_repeat_consensus' => $rc,
			     '_db_id'         => $id } );
  }

  return \@repeats;
}


sub fetch_by_Slice {
    my ( $self, $slice ) = @_;
    
    my $sth = $self->prepare(
	"SELECT r.id, r.hid,  r.chr_name, r.chr_start, r.chr_end, r.chr_strand
         FROM   repeat r
         WHERE  r.chr_name = ? AND r.chr_start <= ? AND r.chr_start >= ? 
                AND r.chr_end >= ?" );
 
    $sth->execute( $slice->chr_name(), $slice->chr_end(), 
		   $slice->chr_start() - $MAX_REPEAT_LENGTH, 
		   $slice->chr_start() );
    
    my @repeats = ();
    my($rc, $core, $repeat_adaptor, $rc_adaptor);
    my $slice_start = $slice->chr_start + 1;
    my %rc_hash;

    if($core = $self->db()->get_db_adaptor('core')) {
      $repeat_adaptor = $core->get_RepeatFeatureAdaptor();
      $rc_adaptor = $core->get_RepeatConsensusAdaptor();
    } else {
      $repeat_adaptor = $self;
      $self->warn("Core Database not attached to lite database.  Not able to 
                   Retrieve repeat consensus adaptor for repeat features\n");
    }

    my ($id, $hid, $chr_name, $start, $end, $strand);
    $sth->bind_columns(\$id, \$hid, \$chr_name, \$start, \$end, \$strand);

    while($sth->fetch) {
      #cache repeat_consensi to reduce object construction overhead
      unless($rc = $rc_hash{$hid}) {
	$rc = new Bio::EnsEMBL::RepeatConsensus();
	$rc->name($hid);
	if($rc_adaptor) {
	  $rc->adaptor($rc_adaptor);
	}
	$rc_hash{$hid} = $rc;
      }

      #create partially filled repeat object using fast (hacky) constructor 
      push @repeats, Bio::EnsEMBL::RepeatFeature->new_fast(
			   { '_gsf_tag_hash'  => {},
			     '_gsf_sub_array' => [],
			     '_parse_h'       => {},
                             '_gsf_start'         => $start - $slice_start,
                             '_gsf_end'           => $end - $slice_start,
			     '_gsf_strand'        => $strand,
			     '_adaptor'       => $repeat_adaptor,
			     '_repeat_consensus' => $rc,
			     '_db_id'         => $id } );
    }


    return \@repeats; 
}

1;

__END__


