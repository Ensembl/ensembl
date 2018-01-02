=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=cut


=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

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

=head2 store

  Arg [1]    : list of Bio::EnsEMBL::KaryotypeBand @k
               the karyotype bands to store in the database
  Example    : $karyotype_adaptor->store(@karyotype_band);
  Description: Stores a list of karyotype band objects in the database
  Returntype : none
  Exceptions : thrown if @k is not defined, if any of the features do not
               have an attached slice.
               or if any elements of @k are not Bio::EnsEMBL::KaryotypeBand
  Caller     : general
  Status     : Stable

=cut

sub store{
  my ($self,@k) = @_;

  if( scalar(@k) == 0 ) {
    throw("Must call store with list of karyotype bands");
  }

  my $sth = $self->prepare
    ("INSERT INTO karyotype (seq_region_id, seq_region_start, " .
                            "seq_region_end, band, " .
                            "stain) " .
     "VALUES (?,?,?,?,?)");

  my $db = $self->db();
  my $slice_adaptor = $db->get_SliceAdaptor();

 BAND: foreach my $k ( @k ) {

    if( !ref $k || !$k->isa("Bio::EnsEMBL::KaryotypeBand") ) {
      throw("KaryotypeBand must be an Ensembl KaryotypeBand, " .
            "not a [".ref($k)."]");
    }

    if($k->is_stored($db)) {
      warning("KaryotypeBand [".$k->dbID."] is already stored" .
              " in this database.");
      next BAND;
    }

    my $original = $k;
    my $seq_region_id;
    ($k, $seq_region_id) = $self->_pre_store($k);

    $sth->bind_param(1,$seq_region_id,SQL_INTEGER);
    $sth->bind_param(2,$k->start,SQL_INTEGER);
    $sth->bind_param(3,$k->end,SQL_INTEGER);
    $sth->bind_param(4,$k->name,SQL_VARCHAR);
    $sth->bind_param(5,$k->stain,SQL_VARCHAR);

    $sth->execute();

    $original->dbID($self->last_insert_id('karyotype_id', undef, 'karyotype'));
    $original->adaptor($self);
  }
}


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
  my ($self, $ordered) = @_;

  return $self->_list_dbIDs("karyotype",undef, $ordered);
}


1;
