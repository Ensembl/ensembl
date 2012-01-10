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
