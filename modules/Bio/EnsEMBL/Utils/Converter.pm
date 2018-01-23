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

Bio::EnsEMBL::Utils::Converter, a converter factory

=head1 SYNOPSIS

  my $converter = Bio::EnsEMBL::Utils::Converter->new(
    -in  => 'Bio::SeqFeature::Generic',
    -out => 'Bio::EnsEMBL::SimpleFeature'
  );

  my ( $fearture1, $feature2 );
  my $ens_simple_features =
    $converter->convert( [ $feature1, $feature2 ] );
  my @ens_simple_features = @{$ens_simple_features};

=head1 DESCRIPTION

Module to converter the business objects between EnsEMBL and any other
projects, currently BioPerl.

What the ready conversions are,

    Bio::SeqFeature::Generic <-> Bio::EnsEMBL::SeqFeature, Bio::EnsEMBL::SimpleFeature
    Bio::SeqFeature::FeaturePair <-> Bio::EnsEMBL::SeqFeature, Bio::EnsEMBL::RepeatFeature
    Bio::Search::HSP::GenericHSP -> Bio::EnsEMBL::BaseAlignFeature's submodules
    Bio::Tools::Prediction::Gene -> Bio::EnsEMBL::PredictionTranscript
    Bio::Tools::Prediction::Exon -> Bio::EnsEMBL::Exon
    Bio::Pipeline::Analysis -> Bio::EnsEMBL::Analysis

=head1 METHODS

=cut


package Bio::EnsEMBL::Utils::Converter;

use strict;

=head2 new

  Title   : new
  Usage   : 
        my $converter = Bio::EnsEMBL::Utils::Converter->new(
            -in => 'Bio::SeqFeature::Generic',
             -out => 'Bio::EnsEMBL::SimpleFeature'
         );

  Function: constructor for converter object
  Returns : L<Bio::EnsEMBL::Utils::Converter>
  Args    : 
    in - the module name of the input.
    out - the module name of the output.
    analysis - a Bio::EnsEMBL::Analysis object, if converting other objects to EnsEMBL features.
    contig - a Bio::EnsEMBL::RawContig object, if converting other objects to EnsEMBL features.

=cut

sub new {
    my ($caller, @args) = @_;
    my $class = ref($caller) || $caller;
    
    if($class =~ /Bio::EnsEMBL::Utils::Converter::(\S+)/){
        my $self = $class->SUPER::new(@args);
        $self->_initialize(@args);
        return $self;
    }else{
        my %params = @args;
        @params{map {lc $_} keys %params} = values %params;
        my $module = $class->_guess_module($params{-in}, $params{-out});

        return undef unless($class->_load_module($module));
        return "$module"->new(@args);
    }
}

# This would be invoked by sub-module's _initialize.

sub _initialize {
    my ($self, @args) = @_;
    
    my ($in, $out) = $self->_rearrange([qw(IN OUT)], @args);

    $self->in($in);
    $self->out($out);
}

=head2 _guess_module
  
  Usage   : $module = $class->_guess_module(
    'Bio::EnsEMBL::SimpleFeature',
    'Bio::EnsEMBL::Generic'
  );

=cut

sub _guess_module {
    my ($self, $in, $out) = @_;
    if($in =~ /^Bio::EnsEMBL::(\S+)/ and $out =~ /^Bio::EnsEMBL::(\S+)/){
        $self->throw("Cannot convert between EnsEMBL objects.\n[$in] to [$out]");
    }elsif($in =~ /^Bio::EnsEMBL::(\S+)/){
        return 'Bio::EnsEMBL::Utils::Converter::ens_bio';
    }elsif($out =~ /^Bio::EnsEMBL::(\S+)/){
        return 'Bio::EnsEMBL::Utils::Converter::bio_ens';
    }else{
        $self->throw("Cannot convert between non-EnsEMBL objects.\n[$in] to [$out]");
    }
}

=head2 convert
    
    Title   : convert
    Usage   : my $array_ref = $converter->convert(\@input);
    Function: does the actual conversion
    Returns : an array ref of converted objects
    Args    : an array ref of converting objects

=cut

sub convert{
    my ($self, $input) = @_;

    $input || $self->throw("Need a ref of array of input objects to convert");
    
    my $output_module = $self->out;
    $self->throw("Cannot load [$output_module] perl module") 
        unless $self->_load_module($output_module);

    unless(ref($input) eq 'ARRAY'){
        $self->warn("The input is supposed to be an array ref");
        return $self->_convert_single($input);
    }

    my @output = ();
    foreach(@{$input}){
        push(@output, $self->_convert_single($_));
    }

    return \@output;
}

sub _convert_single{
    shift->throw("Not implemented. Please check the instance subclass");
}

foreach my $field (qw(in out)){
    my $slot=__PACKAGE__ ."::$field";
    no strict 'refs'; ## no critic
    *$field=sub{
        my $self=shift;
        $self->{$slot}=shift if @_;
        return $self->{$slot};
    };
}

=head2 _load_module
  
  This method is copied from Bio::Root::Root

=cut

sub _load_module {
    my ($self, $name) = @_;
    my ($module, $load, $m);
    $module = "_<$name.pm";
    return 1 if $main::{$module};

    # untaint operation for safe web-based running (modified after a fix
    # a fix by Lincoln) HL
    if ($name !~ /^([\w:]+)$/) {
        $self->throw("$name is an illegal perl package name");
    }

    $load = "$name.pm";
    my $io = Bio::Root::IO->new();
    # catfile comes from IO
    $load = $io->catfile((split(/::/,$load)));
    eval {
        require $load;
    };
    if ( $@ ) {
        $self->throw("Failed to load module $name. ".$@);
    }
    return 1;
}

1;
