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

Bio::EnsEMBL::Utils::Converter::bio_ens

=head1 SYNOPISIS

You should not use this module directly. Please check out the
Bio::EnsEMBL::Utils::Converter module.

=head1 DESCRIPTION

=head1 METHODS

=cut

package Bio::EnsEMBL::Utils::Converter::bio_ens;

use strict;
use vars qw(@ISA);
use Bio::EnsEMBL::Analysis;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Utils::Converter;
use Scalar::Util qw(weaken);
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
    
    my ($dbadaptor, 
        $dbdriver, $dbhost, $dbport, $dbuser, $dbpass, $dbname,
        $analysis, $analysis_dbid, $analysis_logic_name, 
        $contig, $contig_dbid, $contig_name, 
        $translation_id) =
        
        $self->_rearrange([qw(DBADAPTOR 
            DBDRIVER DBHOST DBPORT DBUSER DBPASS DBNAME
            ANALYSIS ANALYSIS_DBID ANALYSIS_LOGIC_NAME 
            CONTIG CONTIG_DBID CONTIG_NAME
            TRANSLATION_ID)], @args);

    
    if(defined $dbadaptor){
        $self->dbadaptor($dbadaptor);
    }elsif(defined $dbname){
        $self->ensembl_db(@args);
    }else{
        # No db information.
    }
    
    if(defined $analysis){
        $self->analysis($analysis);
        # then ignore the analysis_dbid and analysis_logic_name
    }elsif(defined $analysis_dbid){
        $self->analysis_dbID($analysis_dbid);
    }elsif(defined $analysis_logic_name){
        $self->analysis_logic_name($analysis_logic_name);
    }else{
        # No analysis information offered
    }
    
    if(defined $contig){
        ($contig) = ref($contig) eq 'ARRAY' ? @{$contig} : $contig;
        $self->contig($contig);
    }elsif(defined $contig_dbid){
        $self->contig_dbID($contig_dbid);
    }elsif(defined $contig_name){
        $self->contig_name($contig_name);
    }else{
        # No contig information
    }

    if(defined $translation_id){
        $self->translation_id($translation_id);
    }
}


sub _guess_module {
    my ($self, $in, $out) = @_;
    my $tail;
    if($in eq 'Bio::Search::HSP::GenericHSP'){
        $tail = 'bio_ens_hsp';
    }elsif($in eq 'Bio::SeqFeature::Generic'){
        $tail = 'bio_ens_seqFeature';
    }elsif($in eq 'Bio::SeqFeature::FeaturePair'){
        $tail = 'bio_ens_featurePair';
    }elsif($in eq 'Bio::Pipeline::Analysis'){
        $tail = 'bio_ens_analysis';
    }elsif($in eq 'Bio::Tools::Prediction::Gene'){
        $tail = 'bio_ens_predictionGene';
    }elsif($in eq 'Bio::Tools::Prediction::Exon'){
        $tail = 'bio_ens_predictionExon';
    }elsif($in eq 'Bio::SeqFeature::Gene::GeneStructure'){
        $tail = 'bio_ens_gene';
    }elsif($in eq 'Bio::SeqFeature::Gene::Transcript'){
        $tail = 'bio_ens_transcript';
    }elsif($in eq 'Bio::SeqFeature::Gene::Exon'){
        $tail = 'bio_ens_exon';
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
        # convert the analysis, if it's not Bio::Pipeline::Analysis
        if($arg->isa('Bio::Pipeline::Analysis')){
            my $converter_for_analysis = new Bio::EnsEMBL::Utils::Converter(
                -in => 'Bio::Pipeline::Analysis',
                -out => 'Bio::EnsEMBL::Analysis'
            );
            ($arg) = @{ $converter_for_analysis->convert([$arg]) };
        }

        $self->throws("A Bio::EnsEMBL::Analysis object expected.") 
            unless($arg->isa('Bio::EnsEMBL::Analysis'));
        $self->{_analysis} = $arg;
        $self->{_analysis_dbid} = $arg->dbID;
        $self->{_analysis_logic_name} = $arg->logic_name;
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
        if($arg->isa('Bio::EnsEMBL::RawContig')){
            $self->{_contig_dbid} = $arg->dbID;
            $self->{_contig_name} = $arg->name;
        }elsif($arg->isa('Bio::EnsEMBL::Slice')){
                $self->{_slice_dbid} = $arg->dbID;
        }elsif($arg->isa('Bio::PrimarySeqI')){
            ;
        }else{
            $self->throw("a Bio::EnsEMBL::RawContig needed");
        }
        $self->{_contig} = $arg;
        
    }
    return $self->{_contig};
}
    
=head2 dbadaptor

  Title   : dbadaptor
  Usage   : $self->dbadaptor
  Function: get and set for dbadaptor
  Return  : L<Bio::EnsEMBL::DBSQL::DBAdaptor>
  Args    : L<Bio::EnsEMBL::DBSQL::DBAdaptor>   

=cut

sub dbadaptor {
    my ($self, $arg) = @_;
    if(defined($arg)){
        $self->throws("A Bio::EnsEMBL::DBSQL::DBAdaptor object expected.") unless(defined $arg);
        weaken($self->{_dbadaptor} = $arg);
    }
    return $self->{_dbadaptor};
}

=head2 ensembl_db

  Title   : ensembl_db
  Usage   : 
  Function: 
  Return  :
  Args    :

=cut

sub ensembl_db {
    my ($self, @args) = @_;
    
    my ($dbdriver, $dbhost, $dbport, $dbuser, $dbpass, $dbname) = $self->_rearrange(
        [qw(DBDRIVER DBHOST DBPORT DBUSER DBPASS DBNAME)], @args);

    my $dbadaptor = new Bio::EnsEMBL::DBSQL::DBAdaptor(
        -driver => $dbdriver,
        -host => $dbhost,
        -port => $dbport,
        -user => $dbuser,
        -pass => $dbpass,
        -dbname => $dbname
    );
    $self->dbadaptor($dbadaptor);
}

=head2 analysis_dbID

  Title   : analysis_dbID
  Usage   : 
  Function: 
  Return  :
  Args    :

=cut

sub analysis_dbID {
    my ($self, $arg) = @_;

    if(defined $arg){
        my $analysis;
        eval{
            $analysis = $self->dbadaptor->get_AnalysisAdaptor->fetch_by_dbID($arg);
        };
        $self->throw("Failed during fetching analysis by dbID\n$@") if($@);
        $self->analysis($analysis);
    }
    $self->{_analysis_dbid};    
}


=head2 analysis_logic_name

  Title   : analysis_logic_name
  Usage   : 
  Function: 
  Return  :
  Args    :

=cut

sub analysis_logic_name {
    my ($self, $arg) = @_;
    
    return $self->{_analysis_logic_name} unless(defined $arg);
    my $analysis;
    eval{
        $analysis = 
            $self->dbadaptor->get_AnalysisAdaptor->fetch_by_logic_name($arg);
    };
    $self->throw("Not found analysis with logic name as \[$arg\]\n$@") if($@);

    $self->analysis($analysis);
    return $self->{_analysis_logic_name};
}

=head2 contig_dbID

  Title   : contig_dbID
  Usage   : $self->contig_dbID
  Function: get and set for contig_dbID
  Return  : 
  Args    :    

=cut

sub contig_dbID {
    my ($self, $arg) = @_;
    if(defined($arg)){
        my $contig;
        eval{
            $contig = 
                $self->dbadaptor->get_RawContigAdaptor->fetch_by_dbID($arg);
        };
        $self->throw("Failed during fetching contig by dbID\n$@") if($@);
        $self->contig($contig);
    }
    return $self->{_contig_dbid};
}

=head2 contig_name
  Title   : contig_name
  Usage   : $self->contig_name
  Function: get and set for contig_name
  Return  : 
  Args    :    
=cut

sub contig_name {
    my ($self, $arg) = @_;
    if(defined($arg)){
        my $contig;
        eval{
            $contig = 
                $self->dbadaptor->get_RawContigAdaptor->fetch_by_name($arg);
        };
        $self->throw("Failed during fetching contig by dbID\n$@") if($@);
        $self->contig($contig);
    }
    return $self->{_contig_name};
}

=head2 slice_dbID
  Title   : slice
  Usage   : $self->slice
  Function: get and set for slice
  Return  : L<Bio::EnsEMBL::Slice>
  Args    : L<Bio::EnsEMBL::Slice>   
=cut

sub slice_dbID {
    my ($self, $arg) = @_;
    if(defined($arg)){
        my $slice;
        $self->throw("undefined dbadpator") unless defined $self->dbadpaotr;

        eval{
            my $sliceAdaptor = $self->dbadaptor->get_SliceAdaptor;
            $slice = $sliceAdaptor->fetch_by_dbID($arg);
        };
            
        $self->throw("Failed to fetch slice by dbID\n$@") if($@);
        $self->contig($slice);
    }
}

=head2 slice_chr_start_end
  Title   : slice_chr_start_end
  Usage   : my $slice = $self->slice_chr_start_end($chr, $start, $end);
  Function: get and set for slice_chr_start_end
  Return  : 
  Args    :    
=cut

sub slice_chr_start_end {
    my ($self, $chr, $start, $end) = @_;
    if(defined($chr) && defined($start) && defined($end)){
        my $slice;
        eval{
            my $sliceAdaptor = $self->dbadaptor->get_SliceAdaptor;
            $slice = $sliceAdaptor->fetch_by_chr_start_end($chr, $start, $end);
        };
        $self->throw("Failed to fetch slice by chr start end\n$@") if($@);
        $self->contig($slice);
    }
}
 
sub translation_id {
    my ($self, $arg) = @_;
    return $self->{_translation_id} = $arg if(defined($arg));
    return $self->{_translation_id};
}
1;
