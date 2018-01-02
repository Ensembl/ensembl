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

Bio::EnsEMBL::DBSQL::KaryotypeBand

=head1 SYNOPSIS

  use Bio::EnsEMBL::KaryotypeBand;

  # Create and populate a karyotype band (normally done by adaptor)
  $kb = Bio::EnsEMBL::KaryotyeBand(
    -START   => 1,
    -END     => 1_000_000,
    -SLICE   => $chrX_slice,
    -NAME    => 'q31',
    -STAIN   => 'gpos50',
    -ADAPTOR => $db->get_KaryotypeBandAdaptor(),
    -DBID    => 10
  );

  # Can tranform this band into other coord systems, just like other
  # features
  $kb = $kb->transform('supercontig');

  $start      = $kb->start();
  $end        = $kb->end();
  $seq_region = $kb->slice->seq_region_name();

  # Karyotypes have internal ids as well
  $kary_id = $kb->dbID();

=head1 DESCRIPTION

KaryotypeBand objects encapsulate data pertaining to a
single karyotype band.  Access these objects through a
Bio::EnsEMBL::DBSQL::KaryotypeBandAdaptor.

KarytoypeBand inherits from Bio::EnsEMBL::Feature and can be used just
as any other feature can be.

=head1 METHODS

=cut

package Bio::EnsEMBL::KaryotypeBand;

use strict;
use vars qw(@ISA);

use Bio::EnsEMBL::Feature;
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Utils::Exception qw(deprecate warning);

@ISA = qw(Bio::EnsEMBL::Feature);


=head2 new

  Arg [NAME] : string (optional)
               The name of this band
  Arg [STAIN]: string (optional)
               The stain of this band
  Arg [...]  : Arguments passed to superclass constructor.
               See Bio::EnsEMBL::Feature
  Example    : $kb = Bio::EnsEMBL::KaryotypeBand->new(-START => $start,
                                                      -END   => $end,
                                                      -SLICE => $slice,
                                                      -NAME  => 'q11.21',
                                                      -STAIN => 'gneg');
  Description: Constructor.  Creates a new KaryotypeBand object, which can be
               treated as any other feature object. Note that karyotypes
               bands always have strand = 0.
  Returntype : Bio::EnsEMBL::KarytotypeBand
  Exceptions : none
  Caller     : Bio::EnsEMBL::KaryotypeBandAdaptor
  Status     : Stable

=cut

sub new {
  my $class = shift;

  my $self = $class->SUPER::new(@_);

  my ($name, $stain) = rearrange(['NAME','STAIN'],@_);
  $self->{'name'} = $name;
  $self->{'stain'} = $stain;
  $self->{'strand'} = 0;

  return $self;
}


=head2 name

  Arg [1]    : (optional) string $value
  Example    : my $band_name = $band->name(); 
  Description: Getter/Setter for the name of this band
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub name{
  my $self = shift;
  $self->{'name'} = shift if(@_);
  return $self->{'name'};
}



=head2 stain

  Arg [1]    : (optional) string $value
  Example    : my $band_stain = $band->stain();
  Description: get/set for the band stain (e.g. 'gpos50')
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub stain{
  my $self = shift;
  $self->{'stain'} = shift if(@_);
  return $self->{'stain'};
}



=head2 strand

  Arg [1]    : none
	Example    : $strand = $karyotype_feat->strand();
  Description: Overrides the Feature strand method to always return a
               value of 0 for karyotype features (they are unstranded features)
  Returntype : int (always 0)
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub strand {
	my $self = shift;
  return 0;
}


=head2 move

  Arg [1]    : $start - The new end of this band
  Arg [2]    : $end - The new start of this band
  Arg [3]    : $strand - ignored always set to 0
  Example    : $kb->move(1, 10_000);
  Description: Overrides superclass move() method to ensure strand is always 0.
               See Bio::EnsEMBL::Feature::move
  Returntype : none
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub move {
  my ($self, $start, $end, $strand) = @_;

  #maintain a strandedness of 0
  return $self->SUPER::move($start,$end,0);
}


=head2 display_id

  Arg [1]    : none
  Example    : print $kb->display_id();
  Description: This method returns a string that is considered to be
               the 'display' identifier.  For karyotype bands this is the
               name of the karyotype band or '' if no name is defined.
  Returntype : string
  Exceptions : none
  Caller     : web drawing code
  Status     : Stable

=cut

sub display_id {
  my $self = shift;
  return $self->{'name'} || '';
}

=head2 summary_as_hash

  Example       : $karyotype_band_summary = $band->summary_as_hash();
  Description   : Extends Feature::summary_as_hash
                  Retrieves a summary of this Karyotype object.
                    
  Returns       : hashref of arrays of descriptive strings
  Status        : Intended for internal use
=cut

sub summary_as_hash {
  my $self = shift;
  my $summary_ref = $self->SUPER::summary_as_hash;
  $summary_ref->{'stain'} = $self->stain;
  
  return $summary_ref;
}


1;
