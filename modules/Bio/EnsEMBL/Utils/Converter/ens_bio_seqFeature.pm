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

Juguang Xiao <juguang@fugu-sg.org>

=cut

=head1 NAME

Bio::EnsEMBL::Utils::Converter::ens_bio_seqFeature

=head1 SYNOPISIS

=head1 DESCRIPTION

=head1 METHODS

=cut

package Bio::EnsEMBL::Utils::Converter::ens_bio_seqFeature;

use strict;
use vars qw(@ISA);
use Bio::EnsEMBL::Utils::Converter::ens_bio;
@ISA = qw(Bio::EnsEMBL::Utils::Converter::ens_bio);

sub _convert_single {
    my ($self, $in) = @_;
    
    $self->throw("Input not defined") unless($in && defined($in));
    unless(ref($in) && $in->isa('Bio::EnsEMBL::SeqFeature')){
        $self->throw('A Bio::EnsEMBL::SeqFeature object needed');
    }

    my @args = (
        -start => $in->start,
        -end => $in->end,
        -strand => $in->strand,
        -score => $in->score,
        -source_tag => $in->source_tag,
        -seq_id => $in->seqname
    );
    
    my $seqFeature = new Bio::SeqFeature::Generic(@args);
    
    return $seqFeature;
}

1;
