=head1 LICENSE

  Copyright (c) 1999-2013 The European Bioinformatics Institute and
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

=head1 NAME

Bio::EnsEMBL::DBSQL::GenomeContainer - Encapsulates all access to 
genome related information

=head1 SYNOPSIS

  use Bio::EnsEMBL::Registry;

  Bio::EnsEMBL::Registry->load_registry_from_db(
    -host => 'ensembldb.ensembl.org',
    -user => 'anonymous'
  );

  $genome =
    Bio::EnsEMBL::Registry->get_adaptor( "human", "core", "GenomeContainer" );

  my $version = $genome->get_version;

  my $ref_length = $genome->get_ref_length;

  my $coord_systems = $genome->get_coord_systems;



=head1 DESCRIPTION

This module is responsible for fetching and storing genome-wide information.
Genome is an abstract object which contains information linking the species, the assembly and the ensembl annotation.

=head1 METHODS

=cut

package Bio::EnsEMBL::DBSQL::GenomeContainer;

use strict;
use warnings;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::Utils::Exception qw( deprecate throw warning );
use Bio::EnsEMBL::Utils::Scalar qw( assert_ref );

use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::DBSQL::BaseAdaptor);



=head2 new

  Arg [...]  : Superclass args.  See Bio::EnsEMBL::DBSQL::BaseAdaptor
  Description: Instantiates a Bio::EnsEMBL::DBSQL::GenomeContainer
  Returntype : Bio::EnsEMBL::GenomeContainer
  Exceptions : none
  Caller     : DBAdaptor
  Status     : Stable

=cut


sub new {
  my $class = shift;

  my $self = $class->SUPER::new(@_);

  # cache creation could go here
  return $self;
}

=head2 _meta_container

  Arg [1]    : none
  Example    : $meta_container = $genome->_meta_container();
  Description: Internal method to return a MetaContainer object for the genome
  Returntype : Bio::EnsEMBL::DBSQL::MetaContainer
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub _meta_container {
  my $self = shift;
  return $self->db->get_adaptor('MetaContainer');
}


=head2 get_version

  Arg [1]    : (optional) assembly version
  Example    : $version = $genome->get_version();
  Description: Getter/Setter for the assembly version

  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_version {
  my ($self, $version) = @_;
  if (defined $version) {
    $self->{'version'} = $version;
  }
  if (!defined $self->{'version'}) {
    my $csa = $self->db()->get_adaptor('CoordSystem');
    my @cs = @{ $csa->fetch_all() };
    $self->{'version'} = $cs[0]->version();
  }
  return $self->{'version'};
}

=head2 get_accession

  Arg [1]    : (optional) assembly accession
  Example    : $accession = $genome->get_accession();
  Description: Getter/setter for the accession of the assembly currently used

  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_accession {
  my ($self, $accession) = @_;
  if (defined $accession) {
    $self->{'accession'} = $accession;
  }
  if (!defined $self->{'accession'}) {
    $self->{'accession'} = $self->_meta_container->single_value_by_key('assembly.accession');
  }
  return $self->{'accession'};
}

=head2 _get_length

  Arg [1]    : none
  Example    : $length = $genome->_get_length('toplevel');
  Description: Internal method to return the length for a type of slices
  Returntype : integer
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub _get_length {
  my ($self, $cs_name) = @_;
  my $slice_adaptor = $self->db->get_adaptor('Slice');
  my $seqlevel = $slice_adaptor->fetch_all($cs_name);
  my $count;
  foreach my $seq (@$seqlevel) {
    $count += $seq->length();
  }
  return $count;
}



=head2 get_ref_length

  Arg [1]    : (optional) golden path length
  Example    : $ref_length = $genome->get_ref_length();
  Description: Getter/setter for the golden path of the assembly currently used

  Returntype : integer
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_ref_length {
  my ($self, $ref_length) = @_;
  if (defined $ref_length) {
    $self->{'ref_length'} = $ref_length;
  }
  if (!defined $self->{'ref_length'}) {
    $self->{'ref_length'} = $self->_get_length('toplevel');
  }
  return $self->{'ref_length'};
}


=head2 get_total_length

  Arg [1]    : (optional) base pair length
  Example    : $total_length = $genome->get_total_length();
  Description: Getter/setter for the total length (number of base pairs) for the assembly currently used

  Returntype : integer
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_total_length {
  my ($self, $total_length) = @_;
  if (defined $total_length) {
    $self->{'total_length'} = $total_length;
  }
  if (!defined $self->{'total_length'}) {
    $self->{'total_length'} = $self->_get_length('seqlevel');
  }
  return $self->{'total_length'};
}

=head2 get_toplevel

  Arg [1]    : none
  Example    : $toplevel = $genome->get_toplevel();
  Description: Returns the toplevel for the assembly currently used

  Returntype : ListRef of Bio::EnsEMBL::Slice
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_toplevel {
  my $self = @_;
  if (!defined $self->{'toplevel'}) {
    my $sa = $self->db()->get_adaptor('Slice');
    $self->{'toplevel'} = $sa->fetch_all('toplevel');
  }
  return $self->{'toplevel'};
}


=head2 get_karyotype

  Arg [1]    : none
  Example    : $karyotype = $genome->get_karyotype();
  Description: Returns the karyotype for the assembly currently used

  Returntype : ListRef of Bio::EnsEMBL::Slice
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_karyotype {
  my $self = @_;
  if (!defined $self->{'karyotype'}) {
    my $sa = $self->db()->get_adaptor('Slice');
    $self->{'karyotype'} = $sa->fetch_all_karyotype;
  }
  return $self->{'karyotype'};
}

=head2 get_coord_systems

  Arg [1]    : none
  Example    : $coord_systems = $genome->get_coord_systems();
  Description: Returns the coord_systems for the assembly currently used

  Returntype : ListRef of Bio::EnsEMBL::CoordSystem
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_coord_systems {
  my $self = @_;
  if (!defined $self->{'coord_systems'}) {
    my $version = $self->version();
    my $csa = $self->db->get_adaptor('CoordSystem');
    $self->{'coord_systems'} = $csa->fetch_all_by_version($version);
  }
  return $self->{'coord_systems'};
}

=head2 _get_count

  Arg [1]    : none
  Example    : $count = $genome->_get_count('coding_cnt');
  Description: Internal method to return a count for a given attribute code
  Returntype : integer
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub _get_count {
  my ($self, $code) = @_;
  my $aa = $self->db()->get_adaptor('Attribute');
  my $attributes = $aa->fetch_all_by_Object(undef, 'seq_region', $code);
  my $count;
  foreach my $attribute (@$attributes) {
    $count += $attribute->value();
  }
  return $count;
}

=head2 get_coding_count

  Arg [1]    : none
  Example    : $coding_count = $genome->get_coding_count();
  Description: Returns the number of coding genes in the current build

  Returntype : integer
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_coding_count {
  my ($self, $coding_count) = @_;
  if (defined $coding_count) {
    $self->{'coding_count'} = $coding_count;
  }
  if (!defined $self->{'coding_count'}) {
    $self->{'coding_count'} = $self->_get_count('coding_cnt');
  }
  return $self->{'coding_count'};
}


=head2 get_snoncoding_count

  Arg [1]    : none
  Example    : $snoncoding_count = $genome->get_snoncoding_count();
  Description: Returns the number of short non coding genes in the current build

  Returntype : integer
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_snoncoding_count {
  my ($self, $snoncoding_count) = @_;
  if (defined $snoncoding_count) {
    $self->{'snoncoding_count'} = $snoncoding_count;
  }
  if (!defined $self->{'snoncoding_count'}) {
    $self->{'snoncoding_count'} = $self->_get_count('snoncoding_cnt');
  }
  return $self->{'snoncoding_count'};
}

=head2 get_lnoncoding_count

  Arg [1]    : none
  Example    : $lnoncoding_count = $genome->get_lnoncoding_count();
  Description: Returns the number of long non coding genes in the current build

  Returntype : integer
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_lnoncoding_count {
  my ($self, $lnoncoding_count) = @_;
  if (defined $lnoncoding_count) {
    $self->{'lnoncoding_count'} = $lnoncoding_count;
  }
  if (!defined $self->{'lnoncoding_count'}) {
    $self->{'lnoncoding_count'} = $self->_get_count('lnoncoding_cnt');
  }
  return $self->{'lnoncoding_count'};
}

=head2 get_pseudogene_count

  Arg [1]    : none
  Example    : $pseudogene_count = $genome->get_pseudogene_count();
  Description: Returns the number of pseudogenes in the current build

  Returntype : integer
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut


sub get_pseudogene_count {
  my ($self, $pseudogene_count) = @_;
  if (defined $pseudogene_count) {
    $self->{'pseudogene_count'} = $pseudogene_count;
  }
  if (!defined $self->{'pseudogene_count'}) {
    $self->{'pseudogene_count'} = $self->_get_count('pseudogene_cnt');
  }
  return $self->{'pseudogene_count'};
}

=head2 get_alt_coding_count

  Arg [1]    : none
  Example    : $alt_coding_count = $genome->get_alt_coding_count();
  Description: Returns the number of coding genes on alternate sequences

  Returntype : integer
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_alt_coding_count {
  my ($self, $alt_coding_count) = @_;
  if (defined $alt_coding_count) {
    $self->{'alt_coding_count'} = $alt_coding_count;
  }
  if (!defined $self->{'alt_coding_count'}) {
    $self->{'alt_coding_count'} = $self->_get_count('coding_acnt');
  }
  return $self->{'alt_coding_count'};
}


=head2 get_alt_snoncoding_count

  Arg [1]    : none
  Example    : $alt_snoncoding_count = $genome->get_alt_snoncoding_count();
  Description: Returns the number of short non coding genes on alternate sequences

  Returntype : integer
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_alt_snoncoding_count {
  my ($self, $alt_snoncoding_count) = @_;
  if (defined $alt_snoncoding_count) {
    $self->{'alt_snoncoding_count'} = $alt_snoncoding_count;
  }
  if (!defined $self->{'alt_snoncoding_count'}) {
    $self->{'alt_snoncoding_count'} = $self->_get_count('snoncoding_acnt');
  }
  return $self->{'alt_snoncoding_count'};
}

=head2 get_alt_lnoncoding_count

  Arg [1]    : none
  Example    : $alt_lnoncoding_count = $genome->get_alt_lnoncoding_count();
  Description: Returns the number of long non coding genes on alternate sequences

  Returntype : integer
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_alt_lnoncoding_count {
  my ($self, $alt_lnoncoding_count) = @_;
  if (defined $alt_lnoncoding_count) {
    $self->{'alt_lnoncoding_count'} = $alt_lnoncoding_count;
  }
  if (!defined $self->{'alt_lnoncoding_count'}) {
    $self->{'alt_lnoncoding_count'} = $self->_get_count('lnoncoding_acnt');
  }
  return $self->{'alt_lnoncoding_count'};
}


=head2 get_alt_pseudogene_count

  Arg [1]    : none
  Example    : $alt_pseudogene_count = $genome->get_alt_pseudogene_count();
  Description: Returns the number of pseudogenes on alternate sequences

  Returntype : integer
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_alt_pseudogene_count {
  my ($self, $alt_pseudogene_count) = @_;
  if (defined $alt_pseudogene_count) {
    $self->{'alt_pseudogene_count'} = $alt_pseudogene_count;
  }
  if (!defined $self->{'alt_pseudogene_count'}) {
    $self->{'alt_pseudogene_count'} = $self->_get_count('pseudogene_acnt');
  }
  return $self->{'alt_pseudogene_count'};
}

=head2 get_short_variation_count

  Arg [1]    : none
  Example    : $short_variation_count = $genome->get_short_variation_count();
  Description: Returns the number of short variants in the current build

  Returntype : integer
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_short_variation_count {
  my ($self, $short_variation_count) = @_;
  if (defined $short_variation_count) {
    $self->{'short_variation_count'} = $short_variation_count;
  }
  if (!defined $self->{'short_variation_count'}) {
    $self->{'short_variation_count'} = $self->_get_count('short_variation_count');
  }
  return $self->{'short_variation_count'};
}


=head2 get_prediction_count

  Arg [1]    : (optional) logic_name
  Arg [2]    : (optional) prediction_count
  Example    : $prediction_count = $genome->get_prediction_count();
  Description: Return the number of predicted genes in the current build
               Can be restricted to a given analysis

  Returntype : integer
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_prediction_count {
  my ($self, $logic_name, $prediction_count) = @_;
  if (defined $prediction_count) {
    $self->{'prediction_count'} = $prediction_count;
  }
  if (!defined $self->{'prediction_count'}) {
    my $constraint;
    my $pa = $self->db->get_adaptor('PredictionTranscript');
    if ($logic_name) {
      my $analysis_adaptor = $self->get_adaptor('Analysis');
      my $analysis = $analysis_adaptor->fetch_by_logic_name($logic_name);
      $constraint = 'analysis_id = ' . $analysis->dbID;
    }
    $self->{'prediction_count'} = $pa->generic_count($constraint);
  }
  return $self->{'prediction_count'};
}


=head2 get_structural_variation_count

  Arg [1]    : none
  Example    : $structural_variation_count = $genome->get_structural_variation_count();
  Description: Return the number of structural variations in the current build
  Returntype : integer
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_structural_variation_count {
  my $self = @_;
  my $slice_adaptor = $self->db->get_adaptor('Slice');
  my $slices = $slice_adaptor->fetch_all('toplevel');
  my $count;
  foreach my $seq (@$slices) {
    $count += scalar(@{ $seq->get_all_StructuralVariationFeatures() });
  }
  return $count;
}


1;
