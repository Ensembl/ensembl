#
# Ensembl module for Bio::EnsEMBL::DBSQL::KaryotypeBandAdaptor
#
# Cared for by James Stalker <jws@sanger.ac.uk>
#
# Copyright James Stalker
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::DBSQL::KaryotypeBandAdaptor

=head1 SYNOPSIS

$kary_adaptor = $db_adaptor->get_KaryotypeBandAdaptor();
foreach $band ( $kary_adaptor->fetch_by_Slice($slice) ) {
  #do something with band
}

=head1 DESCRIPTION

Database adaptor to provide access to KaryotypeBand objects

=head1 AUTHOR

James Stalker

This modules is part of the Ensembl project http://www.ensembl.org

=head1 CONTACT

Email jws@sanger.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal 
methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::DBSQL::KaryotypeBandAdaptor;
use vars qw(@ISA);
use Bio::EnsEMBL::KaryotypeBand;
use strict;

# Object preamble - inherits from Bio::EnsEMBL::BaseAdaptor

use Bio::EnsEMBL::DBSQL::BaseAdaptor;

@ISA = qw(Bio::EnsEMBL::DBSQL::BaseAdaptor);


# inherit new from BaseAdaptor


=head2 fetch_all_by_Slice

  Arg [1]    : Bio::EnsEMBL::Slice $slice
               The slice object covering the region to retrieve bands from 
  Example    : @bands = @{$karyotype_band_adaptor->fetch_all_by_Slice($slice)};
  Description: Fetches karyotype band object from the database for the
               region given by the slice.
  Returntype : listref of Bio::EnsEMBL::KaryotypeBand objects in 
               Slice coordinates
  Exceptions : none
  Caller     : Bio::EnsEMBL::Slice::get_KarytypeBands 

=cut

sub fetch_all_by_Slice {
  my ($self,$slice) = @_;
  
  my $start = $slice->chr_start();
  my $end = $slice->chr_end();
  my $chr_name = $slice->chr_name();

  my $chr_adaptor = $self->db()->get_ChromosomeAdaptor();
  my $chr_id = $chr_adaptor->fetch_by_chr_name($chr_name)->dbID();

  my $sth = $self->prepare("	SELECT	chr_start,
					chr_end,
					band,
					stain
				FROM	karyotype 
				WHERE	chromosome_id = $chr_id
				AND	$start <= chr_end 
				AND	$end > chr_start 
                                order   by chr_start
			     ");

  $sth->execute;
  my @bands = ();
  my ($chr_start,$chr_end,$band,$stain) = ();
  
  while (($chr_start,$chr_end,$band,$stain) = $sth->fetchrow_array()){
    last unless defined $band;
    my $band_obj = Bio::EnsEMBL::KaryotypeBand->new();
    $band_obj->name($band);
    $band_obj->chr_name($chr_name);

    #convert to slice coordinates
    $band_obj->start($chr_start - $start + 1);
    $band_obj->end($chr_end - $start + 1);

    $band_obj->stain($stain);
    push (@bands, $band_obj);
  }

  return \@bands;
}


=head2 fetch_all_by_chr_name

  Arg [1]    : string $chr_name
               Name of the chromosome from which to retrieve band objects 
  Example    : @bands=@{$karyotype_band_adaptor->fetch_all_by_chr_name('X')}; 
  Description: Fetches all the karyotype band objects from the database for the
               given chromosome. 
  Returntype : listref of Bio::EnsEMBL::KaryotypeBand in chromosomal 
               (assembly) coordinates 
  Exceptions : none 
  Caller     : general 

=cut

sub fetch_all_by_chr_name {
    my ($self,$chr_name) = @_;

    my $chr_adaptor = $self->db()->get_ChromosomeAdaptor();
    my $chr_id = $chr_adaptor->fetch_by_chr_name($chr_name)->dbID();

    my $sth = $self->prepare(
        "SELECT	chr_start, chr_end, stain, band
         FROM karyotype 
         WHERE chromosome_id = ?
         ORDER BY chr_start"
    );
    $sth->execute( $chr_id );
    
    my ($chr_start,$chr_end,$stain, $band) = ();
    my @out = ();

    while(($chr_start, $chr_end, $stain, $band) = $sth->fetchrow_array()) {
      last unless defined $chr_start;
      
      my $band_obj = Bio::EnsEMBL::KaryotypeBand->new();
      $band_obj->name($band);
      $band_obj->chr_name($chr_name);
      $band_obj->start($chr_start);
      $band_obj->end($chr_end);
      $band_obj->stain($stain);
      
      push @out, $band_obj;
    }
    
    return \@out;
}



=head2 fetch_all_by_chr_band

  Arg  [1]   : string $chr_name
               Name of the chromosome from which to retrieve the band
  Arg  [2]   : string $band
               The name of the band to retrieve from the specified chromosome
  Example    : $band = $kary_adaptor->fetch_all_by_chr_band('4', 'q23'); 
  Description: Fetches the karyotype band object from the database 
               for the given chromosome and band name.
  Returntype : Bio::EnsEMBL::KaryotypeBand in 
               chromosomal (assembly) coordinates
  Exceptions : none
  Caller     : general

=cut

sub fetch_all_by_chr_band {
    my ($self,$chr_name, $band) = @_;

    $self->throw("Need band name") unless defined $band;

    my $chr_adaptor = $self->db()->get_ChromosomeAdaptor();
   
    my $chr_id = $chr_adaptor->fetch_by_chr_name($chr_name)->dbID();

    my $sth = $self->prepare(
        "SELECT	chr_start, chr_end, stain
         FROM karyotype 
         WHERE band = ? AND chromosome_id = ?");

    $sth->execute( $band, $chr_id );

    my ($chr, $chr_start,$chr_end,$stain) = $sth->fetchrow_array;

    return undef unless defined $chr_start;

    my $band_obj = Bio::EnsEMBL::KaryotypeBand->new();
    $band_obj->name($band);
    $band_obj->chr_name($chr_name);
    $band_obj->start($chr_start);
    $band_obj->end($chr_end);
    $band_obj->stain($stain);

    return $band_obj;
}


=head2 fetch_by_chr_band

  Arg [1]    : none
  Example    : none
  Description: DEPRECATED use fetch_all_by_chr_band instead
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub fetch_by_chr_band {
  my ($self, @args) = @_;

  $self->warn("fetch_by_chr_band has been renamed fetch_all_by_chr_band\n" . caller);

  return $self->fetch_all_by_chr_band(@args);
}



=head2 fetch_by_chr_name

  Arg [1]    : none
  Example    : none
  Description: DEPRECATED use fetch_all_by_chr_name instead
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub fetch_by_chr_name {
  my ($self, @args) = @_;

  $self->warn("fetch_by_chr_name has been renamed fetch_all_by_chr_name\n" . caller);

  return $self->fetch_all_by_chr_name(@args);
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
