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

  Arne Stabenau: stabenau@ebi.ac.uk
  Ewan Birney  : birney@ebi.ac.uk
  Graham McVicker : mcvicker@ebi.ac.uk

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
    
    my @repeats;
    while( my $row = $sth->fetchrow_arrayref() ) {
      my $id = $row->[0];
      my $hid = $row->[1];
      my $chr_name = $row->[2];
      my $start = $row->[3];
      my $end = $row->[4];
      my $strand = $row->[5];

      #create a partially filled repeat consensus object
      my $rc = new Bio::EnsEMBL::RepeatConsensus;
      my $core = $self->db()->get_db_adaptor('core');

      #create a partially filled repeat object
      my $r = new Bio::EnsEMBL::RepeatFeature();
            $rc->name($hid);
      $r->repeat_id($id);

      $r->start($start - $slice->chr_start() + 1);
      $r->end($end - $slice->chr_start() + 1);
      $r->strand($strand);

      if($core) {
	$rc->adaptor($core->get_RepeatConsensusAdaptor());
	$rc->name($hid);
	$rc->dbID($id);
	$r->repeat_consensus($rc);

	#set the adaptor to be the proxy repeat feature adaptor
	$r->adaptor($core->get_RepeatFeatureAdaptor());
      } else {
	$self->warn("Core Database not attached to lite database.  Not able to 
                     Retrieve repeat consensi for repeat features\n");
	$r->adaptor($self);
      }

      
      push @repeats, $r;
    }

    return @repeats;
}

1;

__END__


