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
