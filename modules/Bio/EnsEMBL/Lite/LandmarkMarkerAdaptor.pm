#
# EnsEMBL module for Bio::EnsEMBL::DBSQL::LandmarkMarkerAdaptor
#
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::DBSQL::LandmarkMarkerAdaptor - DESCRIPTION of Object

=head1 SYNOPSIS

Give standard usage here

=head1 DESCRIPTION

Simple Feature Adaptor - database access for simple features 

=head1 AUTHOR - Ewan Birney

Email birney@ebi.ac.uk

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::Lite::LandmarkMarkerAdaptor;
use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::SimpleFeature;
use Bio::EnsEMBL::MarkerFeature;


@ISA = qw(Bio::EnsEMBL::DBSQL::BaseAdaptor);

my $MAX_FEATURE_LENGTH = 5000;



sub fetch_by_dbID {
  my ($self,$id) = @_;
  
  $self->warn( "No dbID in landmark_features" );
  return;
}



sub fetch_by_contig_id_constraint{
}

sub fetch_by_contig_id{
  my($self, $cid, $logic_name) = @_;

}


=head2 fetch_by_Slice

  Arg  1    : Bio::EnsEMBL::Slice $slice
              the area you want to retrieve features from
  Function  : retrieve Landmark Markers from core database
  Returntype: list Bio::EnsEMBL::MarkerFeature
  Exceptions: none
  Caller    : Bio::EnsEMBL::Slice

=cut



sub fetch_by_Slice {
  my( $self, $slice ) = @_;
 
  my $sth = $self->prepare( "
    SELECT marker, name, chr_start, chr_end, chr_strand, chr_name
      FROM landmark_marker
     WHERE chr_name = ?
       AND chr_start > ? 
       AND chr_start < ?
       AND chr_end > ? " );

  $sth->execute( $slice->chr_name(), $slice->chr_start()-$MAX_FEATURE_LENGTH, 
		 $slice->chr_end(), $slice->chr_start() );

  my @result;
  while( my $hr = $sth->fetchrow_hashref() ) {
    my $start = $hr->{'chr_start'} - $slice->chr_start() + 1;
    my $end = $hr->{'chr_end'} - $slice->chr_start() + 1;
    my $strand = $hr->{'chr_strand'} * $slice->strand();
    
    my $mfeature = Bio::EnsEMBL::MarkerFeature->new();
    
    $mfeature->start( $start );
    $mfeature->end( $end );
    $mfeature->strand( $strand );
    $mfeature->attach_seq( $slice );
    $mfeature->display_label( $hr->{'name'} );
    $mfeature->marker_name( $hr->{'marker'} );

    push( @result, $mfeature );
  }
  return @result;
}




sub store{
   my ($self,$contig_id,@sf) = @_;
   $self->warn( "Cannot store MarkerFeature" );
   return;
}


1;
