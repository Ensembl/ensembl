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

Bio::EnsEMBL::Utils::Converter::ens_bio

=head1 SYNOPISIS

You are not supposed to use this module directly. Please read
Bio::EnsEMBL::Utils::Converter

=head1 DESCRIPTION

This is a helper module to assist Bio::EnsEMBL::Utils::Converter find
which converter instance should be used, based on the -in and -out
parameters.

=head1 METHODS

=cut

package Bio::EnsEMBL::Utils::Converter::ens_bio;

use strict;
use vars qw(@ISA);
use Bio::EnsEMBL::Utils::Converter;
@ISA = qw(Bio::EnsEMBL::Utils::Converter);

=head2 new
  Please see Bio::EnsEMBL::Utils::Converter::new
=cut

sub new {
    my ($caller, @args) = @_;
    my $class = ref($caller) || $caller;

    if($class eq 'Bio::EnsEMBL::Utils::Converter::ens_bio'){
        my %params = @args;
        @params{map{lc $_} keys %params} = values %params;
        my $module = $class->_guess_module($params{-in}, $params{-out});
        return undef unless ($class->_load_module($module));
        return "$module"->new(@args);
    }else{
        my $self = $class->SUPER::new(@args);
#        $self->_initialize(@args);
        return $self;
    }

}

# Unlike bio_ens, ens_bio does not need _initialize method for analysis and 
# contig information.
#

sub _guess_module {
    my ($self, $in, $out) = @_;
    my $tail;
    if($in eq 'Bio::EnsEMBL::SeqFeature'){
        $tail = 'ens_bio_seqFeature';
    }elsif($in eq 'Bio::Ens::EMBL::FeaturePair'){
        $tail = 'ens_bio_featurePair';
    }else{
        $self->throw("[$in] to [$out], not supported");
    }
    return "Bio::EnsEMBL::Utils::Converter::$tail";
}

1;
