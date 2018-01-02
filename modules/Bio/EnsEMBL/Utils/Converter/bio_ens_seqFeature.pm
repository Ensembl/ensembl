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

=head1 AUTHOR

Juguang Xiao <juguang@tll.org.sg>

=cut

=head1 NAME

Bio::EnsEMBL::Utils::Converter::bio_ens_seqFeature

=head1 SYNOPISIS

Please read Bio::EnsEMBL::Utils::Converter

=head1 DESCRIPTION

=head1 METHODS

=cut

package Bio::EnsEMBL::Utils::Converter::bio_ens_seqFeature;

use strict;
use vars qw(@ISA);
use Bio::EnsEMBL::SeqFeature;
use Bio::EnsEMBL::SimpleFeature;
use Bio::EnsEMBL::Exon;
use Bio::EnsEMBL::Utils::Converter::bio_ens;
@ISA = qw(Bio::EnsEMBL::Utils::Converter::bio_ens);

sub _convert_single {
    my ($self, $in) = @_;
    
    unless($in && defined($in) && $in->isa('Bio::SeqFeature::Generic')){
        $self->throw("a Bio::SeqFeature::Generic object needed");
    }
    
    my $seqFeature = $in;
    my $seq_id = $seqFeature->seq_id;
    unless(defined($seq_id) && $seq_id){
        $self->warn("No seq_id value. EnsEMBL SeqFeature will validate it");
        $seq_id = 'Unknown';
    }
    
    # Debated issue here. There are p_value and percent_id in EnsEMBL API and DB
    # schema, but not in bioperl. If in bioperl there are tags called p_value or
    # percent_id, then the values are passed, otherwise set the default 1.
    #
    # the problem arise when I try to converter the seqfeature for tmhmm to 
    # EnsEMBL seqFeature.
    # -- Juguang, 11 July '03
    my $score = $in->score || 0;
    my $percent_id;
    if($in->has_tag('percent_id')){
        ($percent_id) = $in->get_tag_values('percent_id');
    }else{
        $percent_id ||= 0;
    }
    my $p_value;
    if($in->has_tag('p_value')){
        ($p_value) = $in->get_tag_values('p_value');
    }elsif($in->has_tag('evalue')){
        ($p_value) = $in->get_tag_values('evalue');
    }else{
        $p_value ||= 1;
    }
    my $ens_seqFeature;
    my %args = (
        -start => $in->start,
        -end => $in->end,
        -strand => $in->strand,
        -score => $score,
        -analysis => $self->analysis,
        -source_tag => $in->source_tag,
        -seqname => $seq_id,
        -percent_id => $percent_id,
        -p_value => $p_value
    );

    my $output_module = $self->out;
    
    if($output_module eq 'Bio::EnsEMBL::SeqFeature'){
        
        $ens_seqFeature = new Bio::EnsEMBL::SeqFeature(%args);
    }elsif($self->out eq 'Bio::EnsEMBL::SimpleFeature'){
        $ens_seqFeature = new Bio::EnsEMBL::SimpleFeature(%args);
        # The field that there is in SimpleFeature, but not in SeqFeature.
        $ens_seqFeature->display_label('__NONE__');
    }elsif($self->out eq 'Bio::EnsEMBL::Exon'){
        $ens_seqFeature = Bio::EnsEMBL::Exon->new_fast(
            $self->contig, $seqFeature->start, $seqFeature->end, 
            $seqFeature->strand);
    }elsif($self->out eq 'Bio::EnsEMBL::ProteinFeature'){
        my $seq_id2 = $self->analysis->logic_name;
        unless(defined $self->translation_id){
            $self->throw('translation_id unset, in ProteinFeature conversion');
        }
        $args{'-seqname'} = $self->translation_id;
        $ens_seqFeature = Bio::EnsEMBL::ProteinFeature->new(
            -feature1 => Bio::EnsEMBL::SeqFeature->new(%args),
            -feature2 => Bio::EnsEMBL::SeqFeature->new(
                -start => 0,
                -end => 0,
                -seqname => $seq_id2
            )
        );
    }else{
        $self->throw("[$output_module] as -out, not supported");
    }
    
    $ens_seqFeature->attach_seq($self->contig);
    return $ens_seqFeature;
}

1;
