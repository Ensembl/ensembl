=head1 LICENSE

  Copyright (c) 1999-2012 The European Bioinformatics Institute and
  Genome Research Limited.  All rights reserved.

  This software is distributed under a modified Apache license.
  For license details, please see

    http://www.ensembl.org/info/about/code_licence.html

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <dev@ensembl.org>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk@ensembl.org>.

=cut

=head1 NAME

Bio::EnsEMBL::DBSQL::KaryotypeBandAdaptor

=head1 SYNOPSIS

  $kary_adaptor = $db_adaptor->get_KaryotypeBandAdaptor();

  foreach $band ( @{ $kary_adaptor->fetch_all_by_Slice($slice) } ) {
    # do something with band
  }

  $band = $kary_adaptor->fetch_by_dbID($id);

  my @bands = @{ $kary_adaptor->fetch_all_by_chr_name('X') };

  my $band = $kary_adaptor->fetch_by_chr_band( '4', 'q23' );

=head1 DESCRIPTION

Database adaptor to provide access to KaryotypeBand objects

=head1 METHODS

=cut

package Bio::EnsEMBL::DBSQL::KaryotypeBandAdaptor;

use strict;

use vars qw(@ISA);

use Bio::EnsEMBL::KaryotypeBand;
use Bio::EnsEMBL::Utils::Exception qw(throw warning deprecate);
use Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor;

@ISA = qw(Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor);

#_tables
#
#  Arg [1]    : none
#  Example    : none
#  Description: PROTECTED Implementation of abstract superclass method to
#               provide the name of the tables to query
#  Returntype : string
#  Exceptions : none
#  Caller     : internal


sub _tables {
  my $self = shift;

  return (['karyotype','k'])
}


#_columns

#  Arg [1]    : none
#  Example    : none
#  Description: PROTECTED Implementation of abstract superclass method to 
#               provide the name of the columns to query 
#  Returntype : list of strings
#  Exceptions : none
#  Caller     : internal

sub _columns {
  my $self = shift;

  #warning _objs_from_sth implementation depends on ordering
  return qw (
       k.karyotype_id
       k.seq_region_id
       k.seq_region_start
       k.seq_region_end
       k.band
       k.stain );
}


sub _objs_from_sth {
  my ($self, $sth) = @_;
  my $db = $self->db();
  my $slice_adaptor = $db->get_SliceAdaptor();

  my @features;
  my %slice_cache;

  my($karyotype_id,$seq_region_id,$seq_region_start,$seq_region_end,
     $band,$stain);

  $sth->bind_columns(\$karyotype_id, \$seq_region_id, \$seq_region_start,
                     \$seq_region_end, \$band, \$stain);

  while ( $sth->fetch() ) {
    #need to get the internal_seq_region, if present
    $seq_region_id = $self->get_seq_region_id_internal($seq_region_id);

    my $slice = $slice_cache{$seq_region_id} ||=
      $slice_adaptor->fetch_by_seq_region_id($seq_region_id);

    push( @features,
          $self->_create_feature( 'Bio::EnsEMBL::KaryotypeBand', {
                                    -START   => $seq_region_start,
                                    -END     => $seq_region_end,
                                    -SLICE   => $slice,
                                    -ADAPTOR => $self,
                                    -DBID    => $karyotype_id,
                                    -NAME    => $band,
                                    -STAIN   => $stain
                                  } ) );

  }

  return \@features;
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
  Status     : Stable

=cut

sub fetch_all_by_chr_name {
    my ($self,$chr_name) = @_;
    
    throw('Chromosome name argument expected') if(!$chr_name);

    my $slice =
      $self->db->get_SliceAdaptor->fetch_by_region(undef, $chr_name);
    unless ($slice){
        warning("Cannot retrieve chromosome $chr_name");
        return;
    }
    return $self->fetch_all_by_Slice($slice);
}



sub fetch_all_by_chr_band {
  my ($self, $chr_name, $band) = @_;

  throw('Chromosome name argument expected') if(!$chr_name);
  throw('Band argument expected') if(!$band);

  my $slice = $self->db->get_SliceAdaptor->fetch_by_region(undef,
                                                           $chr_name);

  my $constraint = "k.band like '$band%'";
  return $self->fetch_all_by_Slice_constraint($slice,$constraint);
}


=head2 fetch_by_chr_band

  Arg  [1]   : string $chr_name
               Name of the chromosome from which to retrieve the band
  Arg  [2]   : string $band
               The name of the band to retrieve from the specified chromosome
  Example    : @bands = @{$kary_adaptor->fetch_all_by_chr_band('4', 'q23')};
  Description: Fetches the karyotype band object from the database
               for the given chromosome and band name.  If no such band
               exists, undef is returned instead.  This function uses fuzzy
               matching of the band name. For example the bands 'q23.1' and
               'q23.4' could be matched by fetch_all_by_chr_band('20', 'q23');
  Returntype : Bio::EnsEMBL::KaryotypeBand in chromosomal coordinates.
  Exceptions : throws if chr or band is missing in arguments
  Caller     : general
  Status     : Stable

=cut

sub fetch_by_chr_band {
  my $self = shift;
  deprecate('Use fetch_all_by_chr_band instead.');

  my ($band) = @{$self->fetch_all_by_chr_band(@_)};
  return $band;
}


=head2 list_dbIDs

  Arg [1]    : none
  Example    : @kary_ids = @{$karyotype_band_adaptor->list_dbIDs()};
  Description: Gets an array of internal ids for all karyotype bands in the
               current db
  Arg[1]     : <optional> int. not 0 for the ids to be sorted by the seq_region.
  Returntype : reference to a list of ints
  Exceptions : none
  Caller     : ?
  Status     : Stable

=cut

sub list_dbIDs {
  my ($self, $ordered) = shift;

  return $self->_list_dbIDs("karyotype",undef, $ordered);
}


1;
