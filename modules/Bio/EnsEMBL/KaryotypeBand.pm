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

Give standard usage here

=head1 DESCRIPTION

KaryotypeBand objects encapsulate data pertaining to a single karyotype band.
Access these objects through a Bio::EnsEMBL::DBSQL::KaryotypeBandAdaptor

=head1 AUTHOR

James Stalker

This modules is part of the Ensembl project http://www.ensembl.org

=head1 CONTACT

Email jws@sanger.ac.uk

=cut

package Bio::EnsEMBL::KaryotypeBand;

use strict;
use vars qw(@ISA);
use Bio::EnsEMBL::Root;
@ISA = qw(Bio::EnsEMBL::Root);


sub new {
    my ($class) = @_;

    my $self = {};
    bless $self,$class;

    return $self;
}

=head2 name

 Title   : name
 Usage   : my $band_name = $band->name();
           $band->name($newvalue);
 Function: get/set for the band name (e.g.'p34.1') 
 Returns : band name
 Args    : newvalue (optional)

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

  Arg  1    : (optional) string $chr_name
              Name of the chromosome this Karyotype band is on
  Function  : Getter/Setter for the name of the chromosome this band is on
  Returntype: string 
              the name of the chromosome this band is on 
  Exceptions: none
  Caller    : general

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

 Title   : start
 Usage   : my $band_start = $band->start();
           $band->start($newvalue);
 Function: get/set for the band start (e.g. 10000) in absolute basepairs 
 Returns : band start
 Args    : newvalue (optional)

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

 Title   : end
 Usage   : my $band_end = $band->end();
           $band->end($newvalue);
 Function: get/set for the band end (e.g. 10000) in absolute basepairs 
 Returns : band end
 Args    : newvalue (optional)

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

 Title   : stain
 Usage   : my $band_stain = $band->stain();
           $band->stain($newvalue);
 Function: get/set for the band stain (e.g. 'gpos50')
 Returns : band stain
 Args    : newvalue (optional)

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
