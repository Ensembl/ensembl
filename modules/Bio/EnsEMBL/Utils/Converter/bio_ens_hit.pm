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

=head1 DESCRIPTION

Sequence alignment hits were previously stored within the core database
as ungapped alignments. This imposed 2 major constraints on alignments:

  a) alignments for a single hit record would require multiple rows in
     the database, and
  b) it was not possible to accurately retrieve the exact original alignment.

Therefore, in the new branch sequence alignments are now stored as
ungapped alignments in the cigar line format (where CIGAR stands for
Concise Idiosyncratic Gapped Alignment Report).

In the cigar line format alignments are sotred as follows:

  M: Match
  D: Deletino
  I: Insertion

An example of an alignment for a hypthetical protein match is shown
below:


  Query:   42 PGPAGLP----GSVGLQGPRGLRGPLP-GPLGPPL...
              PG    P    G     GP   R      PLGP
  Sbjct: 1672 PGTP*TPLVPLGPWVPLGPSSPR--LPSGPLGPTD...

protein_align_feature table as the following cigar line:

  7M4D12M2I2MD7M

=cut

package Bio::EnsEMBL::Utils::Converter::bio_ens_hit;

use strict;
use vars qw(@ISA);
use Bio::EnsEMBL::Utils::Converter::bio_ens;
use Bio::EnsEMBL::DnaDnaAlignFeature;
use Bio::EnsEMBL::DnaPepAlignFeature;
use Bio::EnsEMBL::PepDnaAlignFeature;
use Bio::EnsEMBL::ProteinFeature;

@ISA = qw(Bio::EnsEMBL::Utils::Converter::bio_ens);

sub _initialize {
    my ($self, @args) = @_;
    $self->SUPER::_initialize(@args);

    # After super initialized, analysis and contig are ready.
    my $bio_ens_seqFeature_converter = new Bio::EnsEMBL::Utils::Converter(
        -in => 'Bio::SeqFeature::Generic',
        -out => 'Bio::EnsEMBL::SeqFeature',
        -analysis => $self->analysis,
        -contig => $self->contig
    );
    $self->_bio_ens_seqFeature_converter($bio_ens_seqFeature_converter);

}

sub _convert_single {
    my ($self, $input) = @_;
    
    my $in = $self->in;
    my $out = $self->out;
    
    if($in =~ /Bio::Search::Hit::GenericHit/){
        return $self->_convert_single_hit($input);
    }elsif($in =~ /Bio::Search::HSP::GenericHSP/){
        return $self->_convert_single_hsp($input);
    }else{
        $self->throw("[$in]->[$out], not implemented");
    }
}

sub _convert_single_hit {


}

sub _convert_single_hsp {
    my ($self, $hsp) = @_;

    unless(ref($hsp) && $hsp->isa('Bio::Search::HSP::GenericHSP')){
        $self->throw("a GenericHSP object needed");
    }

    my $bio_ens_seqFeature_converter = $self->_bio_ens_seqFeature_converter;
    my $ens_feature1 = $bio_ens_seqFeature_converter->_convert_single(
        $hsp->feature1);
    my $ens_feature2 = $bio_ens_seqFeature_converter->_convert_single(
        $hsp->feature2);

    $ens_feature1->p_value($hsp->evalue);
    $ens_feature1->score($hsp->score);
    $ens_feature1->percent_id($hsp->percent_identity);
    $ens_feature2->p_value($hsp->evalue);
    $ens_feature2->score($hsp->score);
    $ens_feature2->percent_id($hsp->percent_identity);
    
    my $cigar_string = $hsp->cigar_string;
    my @args = (
        -feature1 => $ens_feature1,
        -feature2 => $ens_feature2,
        -cigar_string => $cigar_string
    );

    my $contig = $self->contig;
    # choose the AlignFeature based on the blast program
    my $program = $hsp->algorithm;

    $self->throw("HSP does not have algorithm value") unless(defined($program));
    my $align_feature;
    
    if($program =~ /blastn/i){
        $align_feature = new Bio::EnsEMBL::DnaDnaAlignFeature(@args);
        $align_feature->attach_seq($contig);
    }elsif($program =~ /blastx/i){
        $align_feature = new Bio::EnsEMBL::DnaPepAlignFeature(@args);
        $align_feature->attach_seq($contig);
    }else{
        $self->throw("$program is not supported yet");
    }
    
    return $align_feature;
}

# an internal getter/setter for a converter used for seq feature conversion.

sub _bio_ens_seqFeature_converter {
    my ($self, $arg) = @_;
    if(defined $arg){
        $self->{_bio_ens_seqFeature_converter} = $arg;
    }
    return $self->{_bio_ens_seqFeature_converter};
}

1;
