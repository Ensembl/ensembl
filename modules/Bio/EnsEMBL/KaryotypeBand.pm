#
# Ensembl module for Bio::EnsEMBL::DBSQL::KaryotypeBand
#
# Cared for by James Stalker <jws@sanger.ac.uk>
#
# Copyright James Stalker
#
# You may distribute this module under the same terms as perl itself

=head1 NAME

Bio::EnsEMBL::DBSQL::KaryotypeBand

=head1 SYNOPSIS

use Bio::EnsEMBL::KaryotypeBand;

# create and populate a karyotype band
$kb = Bio::EnsEMBL::KaryotyeBand;
$kb->name('q31');
$kb->chr_name('1');
$kb->start(1);
$kb->end(1_000_000);
$kb->stain('gpos50');


=head1 DESCRIPTION

KaryotypeBand objects encapsulate data pertaining to a single karyotype band.
Access these objects through a Bio::EnsEMBL::DBSQL::KaryotypeBandAdaptor

=head1 AUTHOR

James Stalker

This modules is part of the Ensembl project http://www.ensembl.org

=head1 CONTACT

Post questions to the EnsEMBL developer mailing list <ensembl-dev@ebi.ac.uk>  

=cut

package Bio::EnsEMBL::KaryotypeBand;

use strict;
use vars qw(@ISA);
use Bio::EnsEMBL::Root;
@ISA = qw(Bio::EnsEMBL::Root);



=head2 new

  Arg [1]    : none
  Example    : $kb = Bio::EnsEMBL::KaryotypeBand->new;
  Description: Constructor.  Creates a new KaryotypeBand object
  Returntype : Bio::EnsEMBL::KarytotypeBand
  Exceptions : none
  Caller     : Bio::EnsEMBL::KaryotypeBandAdaptor

=cut

sub new {
    my ($class) = @_;

    my $self = {};
    bless $self,$class;

    return $self;
}



=head2 name

  Arg [1]    : (optional) string $value
  Example    : my $band_name = $band->name(); 
  Description: Getter/Setter for the name of this band
  Returntype : string
  Exceptions : none
  Caller     : general

=cut

sub name{
   my $self = shift;
   if( @_ ) {
      my $value = shift;
      $self->{'name'} = $value;
    }
    return $self->{'name'};
}



=head2 chr_name

  Arg [1]    : (optional) string $chr_name
               Name of the chromosome this Karyotype band is on 
  Example    : $chr_name = $band->chr_name; 
  Description: Getter/Setter for the name of the chromosome this band is on 
  Returntype : string
               Name of the chromosom this band is on
  Exceptions : none
  Caller     : general

=cut

sub chr_name {
   my $self = shift;
   if( @_ ) {
      my $value = shift;
      $self->{'_chr_name'} = $value;
    }
    return $self->{'_chr_name'};
}


=head2 start

  Arg [1]    : (optional) int $newvalue
  Example    : my $band_start = $band->start(); 
  Description: get/set for the band start (e.g. 10000) in absolute basepairs  
  Returntype : int 
  Exceptions : none
  Caller     : general

=cut

sub start{
   my $self = shift;
   if( @_ ) {
      my $value = shift;
      $self->{'start'} = $value;
    }
    return $self->{'start'};
}



=head2 end

  Arg [1]    : (optional) int $value
  Example    : $band_end = $band->end;
  Description: get/set for the band end (e.g. 10000) in absolute basepairs  
  Returntype : int
  Exceptions : none
  Caller     : general

=cut

sub end{
   my $self = shift;
   if( @_ ) {
      my $value = shift;
      $self->{'end'} = $value;
    }
    return $self->{'end'};
}



=head2 stain

  Arg [1]    : (optional) string $value
  Example    : my $band_stain = $band->stain(); 
  Description: get/set for the band stain (e.g. 'gpos50') 
  Returntype : string 
  Exceptions : none
  Caller     : general

=cut

sub stain{
   my $self = shift;
   if( @_ ) {
      my $value = shift;
      $self->{'stain'} = $value;
    }
    return $self->{'stain'};
}


1;
