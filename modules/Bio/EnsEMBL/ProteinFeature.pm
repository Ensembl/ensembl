=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2020] EMBL-European Bioinformatics Institute

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

Bio::EnsEMBL::ProteinFeature

=head1 SYNOPSIS

  my $feature = Bio::EnsEMBL::ProteinFeature->new(
    -start    => $start,
    -end      => $end,
    -hstart   => $hit_start,
    -hend     => $hit_end,
    -hseqname => $hit_name
  );

=head1 DESCRIPTION

ProteinFeature objects represent domains or other features of interest
on a peptide sequence.

=head1 METHODS

=cut

package Bio::EnsEMBL::ProteinFeature;

use strict;

use Bio::EnsEMBL::FeaturePair;
use Bio::EnsEMBL::BaseAlignFeature;
use Bio::EnsEMBL::Utils::Argument qw(rearrange);

use parent qw(Bio::EnsEMBL::BaseAlignFeature);

=head2 new

  Arg [IDESC]           : (optional) string An interpro description
  Arg [INTERPRO_AC]     : (optional) string An interpro accession
  Arg [TRANSLATION_ID]  : (optional) integer A translation dbID
  Arg [...]             : named arguments to FeaturePair superclass
  Example    :

    $pf =
      Bio::EnsEMBL::ProteinFeature->new( -IDESC       => $idesc,
                                         -INTERPRO_AC => $iac,
                                         @fp_args );

  Description: Instantiates a Bio::EnsEMBL::ProteinFeature
  Returntype : Bio::EnsEMBL::FeaturePair
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut


sub new {
  my ($proto, @args) = @_;

  my $class = ref($proto) || $proto;

  my $self;
  my ($idesc, $ilabel, $interpro_ac, $translation_id, $external_data, $hit_description, $cigar_string, $align_type, $slice) = rearrange(['IDESC', 'ILABEL', 'INTERPRO_AC', 'TRANSLATION_ID', 'EXTERNAL_DATA', 'HDESCRIPTION', 'CIGAR_STRING', 'ALIGN_TYPE', 'SLICE'], @args);

# BaseAlignFeature expects cigar_line or features
  if($cigar_string && $align_type){
    $self = $class->SUPER::new(@args);
  }else{
  #call the grand parent directly
    $self = $class->Bio::EnsEMBL::FeaturePair::new(@args);
  }

  # the strand of protein features is always 0
  $self->{'strand'}         = 0;
  $self->{'idesc'}          = $idesc || '';
  $self->{'ilabel'}         = $ilabel || '';
  $self->{'interpro_ac'}    = $interpro_ac || '';
  $self->{'translation_id'} = $translation_id || '';
  $self->{'external_data'} = $external_data || '';
  $self->{'hit_description'} = $hit_description || '';
  $self->{'cigar_string'} = $cigar_string || '';
  $self->{'align_type'} = $align_type;

  return $self;
}

=head2 strand

  Arg [1]    : Ignored
  Description: Overwrites Bio::EnsEMBL::Feature->strand to not allow
             : the strand to be set.
  Returntype : int
  Status     : Stable

=cut

#do not allow the strand to be set
sub strand {
  my $self = shift;
  return $self->{'strand'};
}

=head2 idesc

  Arg [1]    : (optional) string The interpro description
  Example    : print $protein_feature->idesc();
  Description: Getter/Setter for the interpro description of this protein
               feature.
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub idesc {
  my $self = shift;
  $self->{'idesc'} = shift if (@_);
  return $self->{'idesc'};
}

=head2 ilabel

  Arg [1]    : (optional) string The interpro label
  Example    : print $protein_feature->ilabel();
  Description: Getter/Setter for the interpro label of this protein
               feature.
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub ilabel {
  my $self = shift;
  $self->{'ilabel'} = shift if (@_);
  return $self->{'ilabel'};
}

=head2 interpro_ac

  Arg [1]    : (optional) string The interpro accession
  Example    : print $protein_feature->interpro_ac();
  Description: Getter/Setter for the interpro accession of this protein
               feature.
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub interpro_ac {
  my $self = shift;
  $self->{'interpro_ac'} = shift if (@_);
  return $self->{'interpro_ac'};
}

=head2 translation_id

  Arg [1]    : (optional) integer The dbID of the translation
  Example    : print $protein_feature->translation_id();
  Description: Getter/Setter for the translation dbID of this protein
               feature.
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub translation_id {
  my $self = shift;
  $self->{'translation_id'} = shift if (@_);
  return $self->{'translation_id'};
}

sub external_data {
  my $self = shift;
  $self->{'external_data'} = shift if (@_);
  return $self->{'external_data'};
}


=head2 summary_as_hash

  Example       : $protein_feature_summary = $protein_feature->summary_as_hash();
  Description   : Retrieves a textual summary of this Protein feature.
                  Not inherited from Feature.
  Returns       : hashref of arrays of descriptive strings
  Status        : Intended for internal use
=cut

sub summary_as_hash {
  my $self = shift;
  my %summary;
  $summary{'type'} = $self->analysis->db;
  $summary{'id'} = $self->display_id;
  $summary{'start'} = $self->start;
  $summary{'end'} = $self->end;
  $summary{'interpro'} = $self->interpro_ac;
  $summary{'description'} = $self->idesc;
  $summary{'hit_start'} = $self->hstart;
  $summary{'hit_end'} = $self->hend;
  $summary{'cigar_string'} = $self->cigar_string;
  $summary{'align_type'} = $self->align_type;
  $summary{'hseqname'} = $self->hseqname;
  $summary{'translation_id'} = $self->translation_id;
  
  return \%summary;
}


=head2 alignment_strings

  Arg [1]    : list of string $flags
  Example    : $pf->alignment_strings
  Description: Allows to rebuild the alignment string of both the query and target sequence
               using the sequence from translation object and
               MD Z String for mismatching positions. Regex : [0-9]+(([A-Z]|\^[A-Z]+)[0-9]+)* (Refer:  SAM/BAM specification)
               eg: MD:Z:96^RHKTDSFVGLMGKRALNS0V14
  Returntype : array reference containing 2 strings
               the first corresponds to seq
               the second corresponds to hseq
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut


sub alignment_strings {
  my $self = shift;

  #Translations
  my $transl_adaptor = $self->adaptor->db->get_TranslationAdaptor();
  my $transl_object = $transl_adaptor->fetch_by_dbID($self->translation_id);
  my $seq;
  if(defined $transl_object && $transl_object->isa('Bio::EnsEMBL::Translation')) {
    $seq = $transl_object->transcript()->translate()->seq();
  }

  if ($self->align_type eq 'mdtag') {
    if(defined $seq && defined $self->cigar_string){
      return $self->_mdz_alignment_string($seq,$self->cigar_string);
    }else{
      warn "sequence or cigar_line not found for  " .  $self->translation_id;
    }
  } else {
    throw("alignment_strings method not implemented for " . $self->align_type);
  }
  return;
}


sub transform {
  my $self = shift;

  $self->throw( "ProteinFeature cant be transformed directly as".
    " they are not on EnsEMBL coord system" );
  return;
}


=head2 _hit_unit

  Arg [1]    : none
  Description: PRIVATE implementation of abstract superclass method.  Returns
               1 as the 'unit' used for the hit sequence.
  Returntype : int
  Exceptions : none
  Caller     : Bio::EnsEMBL::BaseAlignFeature
  Status     : Stable


=cut

sub _hit_unit {
  return 3;
}


=head2 _query_unit

  Arg [1]    : none
  Description: PRIVATE implementation of abstract superclass method.  Returns
               3 as the 'unit' used for the query sequence.
  Returntype : int
  Exceptions : none
  Caller     : Bio::EnsEMBL::BaseAlignFeature
  Status     : Stable


=cut

sub _query_unit {
  return 3;
}





1;
