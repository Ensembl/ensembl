# Bio::EnsEMBL::Utils::Converter::bio_ens
#
# Created and cared for by Juguang Xiao <juguang@fugu-sg.org>
# Created date: 4/3/2003
# 
# Copyright Juguang Xiao
# 
# You may distribute this module under the same terms as perl itself
#
# POD documentation
#

=head1 NAME

Bio::EnsEMBL::Utils::Converter::bio_ens

=head1 SYNOPISIS

You should not use this module directly. Please check out the 
Bio::EnsEMBL::Utils::Converter module.

=head1 DESCRIPTION


=head1 FEEDBACK

=head2 Mailing Lists

=head2 Reporting Bugs


=head1 AUTHOR Juguang Xiao

Juguang Xiao <juguang@fugu-sg.org>

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

# Let the code begin ...

package Bio::EnsEMBL::Utils::Converter::bio_ens;

use strict;
use vars qw(@ISA);
use Bio::EnsEMBL::Analysis;
use Bio::EnsEMBL::Utils::Converter;
@ISA = qw(Bio::EnsEMBL::Utils::Converter);

=head2 new

Please see Bio::EnsEMBL::Utils::Converter::new

=cut

sub new {
    my ($caller, @args) = @_;
    my $class = ref($caller) || $caller;

    if($class eq 'Bio::EnsEMBL::Utils::Converter::bio_ens'){
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

sub _initialize {
    my ($self, @args) = @_;
    $self->SUPER::_initialize(@args);
    
    my ($analysis, $contig) =
        $self->_rearrange([qw(analysis contig)], @args);

    if(defined($analysis) && $analysis->isa('Bio::Pipeline::Analysis')){
        my $converter_for_analysis = new Bio::EnsEMBL::Utils::Converter(
            -in => 'Bio::Pipeline::Analysis',
            -out => 'Bio::EnsEMBL::Analysis'
        );
        ($analysis) = @{$converter_for_analysis->convert([$analysis])};
    }

    $self->analysis($analysis);
    $self->contig($contig);
}


sub _guess_module {
    my ($self, $in, $out) = @_;
    my $tail;
    if($in eq 'Bio::Search::HSP::GenericHSP'){
        $tail = 'bio_ens_hit';
    }elsif($in eq 'Bio::SeqFeature::Generic'){
        $tail = 'bio_ens_seqFeature';
    }elsif($in eq 'Bio::SeqFeature::FeaturePair'){
        $tail = 'bio_ens_featurePair';
    }elsif($in eq 'Bio::Pipeline::Analysis'){
        $tail = 'bio_ens_analysis';
    }else{
        $self->throw("[$in] to [$out], not supported");
    }
    return "Bio::EnsEMBL::Utils::Converter::$tail";
}


=head2 analysis
  Title   : analysis
  Usage   : $self->analysis
  Function: get and set for analysis
  Return  : L<Bio::EnsEMBL::Analysis>
  Args    : L<Bio::EnsEMBL::Analysis>   
=cut

sub analysis {
    my ($self, $arg) = @_;
    if(defined($arg)){
        $self->throws("A Bio::EnsEMBL::Analysis object expected.") unless(defined $arg);
        $self->{_analysis} = $arg;
    }
    return $self->{_analysis};
}
    


=head2 contig
  Title   : contig
  Usage   : $self->contig
  Function: get and set for contig
  Return  : 
  Args    :    
=cut

sub contig {
    my ($self, $arg) = @_;
    if(defined($arg)){
        $self->{_contig} = $arg;
    }
    return $self->{_contig};
}
    


1;
