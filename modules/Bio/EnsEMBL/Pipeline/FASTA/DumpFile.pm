=pod

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

=head1 NAME

Bio::EnsEMBL::Pipeline::FASTA::DumpFile

=head1 DESCRIPTION

The main workhorse of the FASTA dumping pipeline. This module has two
functions

=over 8 

=item 1 - Dumping Genomic DNA sequences in a memory efficient manner in unmasked, softmasked & hardmasked formats

=item 2 - Dumping Genes as cDNA, proteins and ncRNA transcripts (abinitio included)

=back

The script is responsible for creating the filenames of these target
files, taking data from the database and the formatting of the FASTA
headers. It is also responsible for the creation of README files pertaining
to the type of dumps produced. The final files are all Gzipped at normal
levels of compression.

B<N.B.> This code will remove any files already found in the target directory
on its first run as it assumes all data will be dumped in the one process. It
is selective of its directory meaning a rerun of DNA dumps will not cause
the protein/cdna files to be removed.

Allowed parameters are:

=over 8

=item species - The species to dump

=item sequence_type_list - The data to dump. I<dna>, I<cdna> and I<ncrna> are allowed

=item release - A required parameter for the version of Ensembl we are dumping for

=item db_types - Array reference of the database groups to use. Defaults to core

=item process_logic_names - Array reference of transcript logic names to only process (only produce dumps for these). Applied before skip_logic_names

=item skip_logic_names - Array reference of transcript logic names to skip over (we do not produce dumps for these)


=item base_path - The base of the dumps

=item dna_chunk_size - Indicates the number of 60bp chunks to retrieve and 
                       process when formatting FASTA files. Normally do not 
                       touch

=item allow_appending - If the same file name is generated we will 
                        append into that file rather than overwriting

=item overwrite_files - If the same file name is generated we will overwrite 
                        the into that file rather than appending

=back

=cut

package Bio::EnsEMBL::Pipeline::FASTA::DumpFile;

use strict;
use warnings;

use base qw(Bio::EnsEMBL::Pipeline::FASTA::Base);

use File::Spec;
use IO::Compress::Gzip qw/gzip $GzipError/;
use IO::File;
use Bio::EnsEMBL::PaddedSlice;
use Bio::EnsEMBL::Utils::BiotypeMapper;
use Bio::EnsEMBL::Utils::Exception qw/throw/;
use Bio::EnsEMBL::Utils::Scalar qw/check_ref/;
use Bio::EnsEMBL::Utils::IO::FASTASerializer;
use Bio::EnsEMBL::Utils::IO qw/work_with_file gz_work_with_file/;

my $DNA_INDEXING_FLOW = 1;
my $PEPTIDE_INDEXING_FLOW = 2;
my $GENE_INDEXING_FLOW = 3;

sub param_defaults {
  my ($self) = @_;
  return {
    #user configurable
    allow_appending => 1,
    overwrite_files => 0,
    
    dna_chunk_size => 17000,
    
    skip_logic_names => [],
    process_logic_names => [],
    
    #DON'T MESS
    #used to track if we need to reopen a file in append mode or not
    generated_files => {}, 
    remove_files_from_dir => {},
    dataflows => []
  };
}

sub fetch_input {
  my ($self) = @_;
  
  my %sequence_types = map { $_ => 1 } @{ $self->param('sequence_type_list') };
  $self->param('sequence_types', \%sequence_types);
  
  my $dba = $self->get_DBAdaptor();
  my $analyses = $dba->get_MetaContainer()->list_value_by_key('repeat.analysis');
  $self->param('analyses', $analyses);
  
  my $types = $self->param('db_types');
  $types = ['core'] unless $types;
  $self->param('db_types', $types);
  
  my %skip_logic_names = map { $_ => 1 } @{$self->param('skip_logic_names')};
  $self->param('skip_logic', \%skip_logic_names);
  $self->param('skip_logic_active', 1) if @{$self->param('skip_logic_names')};
  my %process_logic_names = map { $_ => 1 } @{$self->param('process_logic_names')};
  $self->param('process_logic', \%process_logic_names);
  $self->param('process_logic_active', 1) if @{$self->param('process_logic_names')};
  
  return;
}

sub run {
  my ($self) = @_;
  my $types = $self->param('db_types');
  foreach my $type (@{$types}) {
    my $dba = $self->get_DBAdaptor($type);
    if(! $dba) {
      $self->info("Cannot continue with %s as we cannot find a DBAdaptor", $type);
      next;
    }
    $self->run_type($type);
  }
  return;
}

sub write_output {
  my ($self) = @_;
  my $dataflows = $self->param('dataflows');
  foreach my $flow (@{$dataflows}) {
    $self->dataflow_output_id(@{$flow});
  }
  return;
}

sub run_type {
  my ($self, $type) = @_;

  my $species = $self->param('species');
  my $sequence_types = $self->param('sequence_types');

  # dump file for each type on a per slice basis
  # types are dna,cDNA, peptide, ncRNA

  #Only run if we are told to & the current DBA is the same as the attached DNADB by checking the Stringified ref
  my $dba = $self->get_DBAdaptor($type);
  if ( $sequence_types->{dna} && $dba eq $dba->dnadb() ) {
    $self->info( "Starting dna dump for " . $species );
    $self->_dump_dna($type);
    $self->_create_README('dna');
  }
  
  if ( $sequence_types->{cdna} ) { #includes peptides whether you like it or not
    $self->info( "Starting cdna dump for " . $species );
    my ($transcripts, $peptide) = $self->_dump_transcripts('cdna', $type);
    
    $self->info( "Starting prediction transcript dumps for " . $species );
    my ($pred_transcripts, $pred_proteins) = $self->_dump_prediction_transcripts($type);
    
    $self->_create_README('cdna') if $transcripts || $pred_transcripts;
    $self->_create_README('pep') if $peptide || $pred_proteins;
  }
  if ( $sequence_types->{ncrna} ) {
    $self->info( "Starting ncRNA dump for " . $species );
    my ($ncrna_transcripts) = $self->_dump_transcripts('ncrna', $type);
    $self->_create_README('ncrna') if $ncrna_transcripts;
  }

  $self->cleanup_DBAdaptor($type);
}

# Dump entire sequence, also dump data into chromosome files as appropriate
sub _dump_dna {
  my ($self,$type) = @_;
  
  my @chromosomes;
  my @non_chromosomes;
  my $filter_human = 1;
  foreach my $s (@{$self->get_Slices($type, $filter_human)}) {
    my $chr = $s->is_chromosome();
    push(@chromosomes, $s) if $chr;
    push(@non_chromosomes, $s) if ! $chr;
  }
  
  ############ NON CHROMOSOME WORK
  $self->info('Processing %d non-chromosome(s)', scalar(@non_chromosomes));
  if(@non_chromosomes) {
    my ( $non_specific_file, $non_specific_fh, $other_serializer ) =
      $self->_generate_fasta_serializer( 'dna', 'nonchromosomal' );
    my ( $rm_non_specific_file, $rm_non_specific_fh, $other_rm_serializer ) =
      $self->_generate_fasta_serializer( 'dna_sm', 'nonchromosomal' );
    foreach my $s (@non_chromosomes) {
      $self->_dump_slice($s, $other_serializer, $other_rm_serializer);
    }
    #Quick close of the SM FH to flush all data out to disk; skip gzipping & leave that to the next call
    $self->tidy_file_handle($rm_non_specific_fh, $rm_non_specific_file, 1);
    my ($hard_mask_fh, $hard_mask_file) = $self->_convert_softmask_to_hardmask($rm_non_specific_file, $rm_non_specific_fh);
    
    $self->tidy_file_handle( $non_specific_fh, $non_specific_file );
    $self->tidy_file_handle( $rm_non_specific_fh, $rm_non_specific_file );
    $self->tidy_file_handle( $hard_mask_fh, $hard_mask_file);
    $self->info('Dumped non-chromosomes');
  }

  ############ CHROMOSOME WORK 
  $self->info('Processing %d chromosome(s)', scalar(@chromosomes));
  foreach my $s (@chromosomes) {
    my ( $chromo_file_name, $chromo_fh, $chromo_serializer ) =
      $self->_generate_fasta_serializer( 'dna', 'chromosome',
      $s->seq_region_name(), undef);
    # repeat masked data too
    my ( $rm_chromo_file_name, $rm_chromo_fh, $rm_chromo_serializer ) =
      $self->_generate_fasta_serializer( 'dna_sm', 'chromosome',
      $s->seq_region_name(), undef);
    
    $self->_dump_slice($s, $chromo_serializer, $rm_chromo_serializer);
    
    #Quick close of the SM FH to flush all data out to disk; skip gzipping & leave that to the next call
    $self->tidy_file_handle($rm_chromo_fh, $rm_chromo_file_name, 1);
    my ($chromo_hard_mask_fh, $chromo_hard_mask_file) = $self->_convert_softmask_to_hardmask($rm_chromo_file_name, $rm_chromo_fh);
    
    $self->tidy_file_handle($chromo_fh, $chromo_file_name);
    $self->tidy_file_handle($rm_chromo_fh, $rm_chromo_file_name);
    $self->tidy_file_handle($chromo_hard_mask_fh, $chromo_hard_mask_file);
  }
  $self->info("Dumped chromosomes");
  
  #input_id
  push(@{$self->param('dataflows')}, [{ data_type => 'dna', species => $self->param('species') }, $DNA_INDEXING_FLOW]);
  push(@{$self->param('dataflows')}, [{ data_type => 'dna_sm', species => $self->param('species') }, $DNA_INDEXING_FLOW]);
  push(@{$self->param('dataflows')}, [{ data_type => 'dna_rm', species => $self->param('species') }, $DNA_INDEXING_FLOW]);
  
  return;
}

sub _dump_slice {
  my ($self, $s, $serialiser, $rm_serialiser) = @_;
  
  my $analyses = $self->param('analyses');
  
  my $chr = $s->is_chromosome();
  $self->info('Starting slice - %s:%d-%d', $s->seq_region_name(), $s->start(), $s->end());
  $self->info('    Slice is a chromosome') if $chr;
  $self->info('    Slice is non-chromosomal') if ! $chr;
  
  # Make a padded slice (to automatically pad with N's outside of known regions)
  # and make a repeat-masked slice and then pad that too.
  my $padded_slice = Bio::EnsEMBL::PaddedSlice->new(-SLICE => $s);
  $serialiser->print_Seq($padded_slice);
  
  my $soft_mask = 1;
  my $masked_slice = $s->get_repeatmasked_seq($analyses, $soft_mask);
  my $padded_masked_slice = Bio::EnsEMBL::PaddedSlice->new(-SLICE => $masked_slice);
  $rm_serialiser->print_Seq($padded_masked_slice);  
  
  return;
}

#Assumes we are working with un-compressed files
sub _convert_softmask_to_hardmask {
  my ($self, $soft_mask_file, $soft_mask_fh) = @_;
  if(! -f $soft_mask_file) {
    $self->info('Skipping as the target file %s does not exist. Must have been deleted', $soft_mask_file);
    return;
  }
  my $hard_mask_file = $soft_mask_file;
  $hard_mask_file =~ s/\.dna_sm\./.dna_rm./;
  my $hm_fh = IO::File->new($hard_mask_file, 'w');
  $self->info('Converting soft-masked file %s into hard-masked file %s', $soft_mask_file, $hard_mask_file);
  work_with_file($soft_mask_file, 'r', sub {
    my ($sm_fh) = @_;
    while(my $line = <$sm_fh>) {
      if(index($line, '>') == 0) {
        $line =~ s/dna_sm/dna_rm/;
      }
      else {
        $line =~ tr/[acgtn]/N/;
      }
      print $hm_fh $line;
    };
    return;
  });
  return ($hm_fh, $hard_mask_file);
}

sub _dump_transcripts {
  my ($self, $transcript_type, $type) = @_;
  
  my $has_transcript_data = 0;
  my $has_protein_data = 0;

  my $transcript_level = ($transcript_type ne 'ncrna') ? 'all' : undef;
  my ( $filename, $fh, $transcript_serializer ) =
    $self->_generate_fasta_serializer( $transcript_type, $transcript_level );

  my ( $peptide_filename, $pep_fh, $peptide_serializer );

  # some cDNAs are translated, make a file to receive them.
  if ( $transcript_type eq 'cdna') {
    ( $peptide_filename, $pep_fh, $peptide_serializer ) =
      $self->_generate_fasta_serializer( 'pep', 'all' );
  }

  # work out what biotypes correspond to $transcript_type
  my $biotype_mapper = Bio::EnsEMBL::Utils::BiotypeMapper->new();
  my $biotypes_list  = $biotype_mapper->group_members($transcript_type);

  my $dba          = $self->get_DBAdaptor($type);
  my $gene_adaptor = $dba->get_GeneAdaptor();

  # get all the transcripts that are $transcript_type e.g. cdna, ncrna,
  foreach my $biotype ( @{$biotypes_list} ) {
    my $gene_list = $gene_adaptor->fetch_all_by_biotype($biotype);
    $self->info("Biotype %s has %d gene(s)", $biotype, scalar( @{$gene_list} ));
    while ( my $gene = shift @{$gene_list} ) {
      $self->fine( 'Gene %s', $gene->display_id );
      my $transcript_list = $gene->get_all_Transcripts();
      foreach my $transcript ( @{$transcript_list} ) {
        $self->fine( 'Transcript %s', $transcript->display_id );
        next unless $self->ok_to_process_logic_name($transcript);

        # foreach transcripts of all genes with biotypes classed as cdna
        my $transcript_seq = $transcript->seq();
        $self->_create_display_id($transcript, $transcript_seq, $transcript_type);
        $transcript_serializer->print_Seq($transcript_seq);
        if ($biotype_mapper->member_of_group( $biotype, 'peptide_producing')) {
          my $translation = $transcript->translation();
          if ($translation) {
            my $translation_seq = $transcript->translate();
            $self->_create_display_id($translation, $translation_seq, $transcript_type);
            $peptide_serializer->print_Seq($translation_seq);
            
            $has_protein_data = 1;
          }
        }
        
        $has_transcript_data = 1;
      }
    }
  }

  $self->tidy_file_handle( $fh, $filename );
  if ( $transcript_type eq 'cdna' ) {
    $self->tidy_file_handle( $pep_fh, $peptide_filename );
  }
  
  if($has_protein_data) {
    push(@{$self->param('dataflows')}, [{ file => $self->_final_filename($peptide_filename), species => $self->param('species') }, $PEPTIDE_INDEXING_FLOW]);
  }
  if($has_transcript_data) {
    push(@{$self->param('dataflows')}, [{ file => $self->_final_filename($filename), species => $self->param('species') }, $GENE_INDEXING_FLOW]);
  }
  
  return ($has_transcript_data, $has_protein_data);
}

# Dump prediction transcripts and peptides. All predicted transcripts have translations
sub _dump_prediction_transcripts {
  my ($self, $type) = @_;
  my $dba  = $self->get_DBAdaptor($type);
  
  my $has_transcript_data = 0;
  my $has_protein_data = 0;
  
  my $prediction_transcript_adaptor = $dba->get_PredictionTranscriptAdaptor();
  my $transcript_list = $prediction_transcript_adaptor->fetch_all();
  my $count = scalar(@{$transcript_list});
  $self->info('Found %d prediction transcript(s)', $count);
  if($count) {
    my ( $abinitio_filename, $fh, $abinitio_serializer ) =
      $self->_generate_fasta_serializer( 'cdna', 'abinitio' );
    my ( $abinitio_peptide_filename, $pep_fh, $abinitio_peptide_serializer ) =
      $self->_generate_fasta_serializer( 'pep', 'abinitio' );
    
    while ( my $transcript = shift @{$transcript_list} ) {
      next unless $self->ok_to_process_logic_name($transcript);
      
      $has_transcript_data = 1;
      my $transcript_seq = $transcript->seq();
      $self->_create_display_id( $transcript, $transcript_seq, 'cdna' );
      $abinitio_serializer->print_Seq($transcript_seq);
    
      my $translation_seq = $transcript->translate();
      if ( $transcript->translation() ) {
        $has_protein_data = 1;
        $self->_create_display_id( $transcript, $translation_seq, 'pep' );
        $abinitio_peptide_serializer->print_Seq($translation_seq);
      }
    }
    
    $self->tidy_file_handle( $fh,     $abinitio_filename );
    $self->tidy_file_handle( $pep_fh, $abinitio_peptide_filename );
    
    if($has_protein_data) {
      push(@{$self->param('dataflows')}, [{ file => $self->_final_filename($abinitio_peptide_filename), species => $self->param('species') }, $PEPTIDE_INDEXING_FLOW]);
    }
    if($has_transcript_data) {
      push(@{$self->param('dataflows')}, [{ file => $self->_final_filename($abinitio_filename), species => $self->param('species') }, $GENE_INDEXING_FLOW]);
    }
  }
  

  return ($has_transcript_data, $has_protein_data);
}

# We can optionally skip the Gzip process & just delegate to the super class
# for it's cleanup routines which only work with an open file handle. Therefore 
# only pass it onto the super implementation *if* the handle was open. 
# Also only Gzip if the source file exists (it could have been unlinked from
# an earlier call)

sub tidy_file_handle {
  my ($self, $fh, $path, $no_gzip) = @_;
  if($fh->opened()) {
    my $tidy = $self->SUPER::tidy_file_handle($fh, $path);
    return 1 if $tidy;
  }
  
  return if $no_gzip; #don't gzip if we were told to skip
  return if ! -f $path; #don't gzip if we had no file
  
  my $target = $path.".gz";
  $self->info('Gzipping "%s"', $path);
  my %args;
  if($self->param('generated_files')->{$target}) {
    if($self->param('allow_appending')) {
      $self->info('Going to append to the file %s as we have created two files of the same name in the same session', $target);
      $args{Append} = 1;
    }
    elsif($self->param('overwrite_files')) {
      $self->info('Overwriting the file %s as we have created two files of the same name in the same session', $target);
    }
    else {
      $self->throw("Cannot continue. The file %s has already been created this session. Fail!");
    }
  }
  gzip $path => $target, %args or throw "GZip error compressing $path to $target: $GzipError";
  $self->info('    Removing original file from filesystem');
  unlink $path or throw "Could not delete $path: $!";
  $self->info('    Finished');
  $self->param('generated_files')->{$target} = 1;
  return 0;
}

#We assume a transcript is ok to process unless proven otherwise
sub ok_to_process_logic_name {
  my ($self, $transcript) = @_;
  my $ok = 1;
  my $logic_name = $transcript->analysis()->logic_name();
  if($self->param('process_logic_active')) {
    if(! $self->param('process_logic')->{$logic_name}) {
      $self->fine('Transcript %s has been filtered because logic_name %s is not in the active logic name list', $transcript->stable_id(), $logic_name);
      $ok = 0;
    }
  }
  if($self->param('skip_logic_active')) {
    if($self->param('skip_logic')->{$logic_name}) {
      $self->fine('Transcript %s has been filtered because logic_name %s is in the skip logic name list', $transcript->stable_id(), $logic_name);
      $ok = 0;
    }
  }
  return $ok;
}

#Generates a FASTA serializer but returns the (filename, handle & instance)
sub _generate_fasta_serializer {
  my ( $self, $datatype, $level, $section, $header_formatter ) = @_;
  $header_formatter ||= $self->_custom_header();
  my $chunk = $self->param('dna_chunk_size');
  my $filename = $self->_generate_file_name( $datatype, $level, $section );
  my $fh = IO::File->new($filename, '>') or throw "Cannot open $filename for writing: $!";
  my $ser = Bio::EnsEMBL::Utils::IO::FASTASerializer->new($fh, $header_formatter, $chunk);
  return ( $filename, $fh, $ser );
}

#
# _generate_file_name(data type, level, section )
#                      dna       toplevel  undef
#                      dna       chromosome 6

sub _generate_file_name {
  my ( $self, $data_type, $level, $section ) = @_;    #level & section is optional

  # File name format looks like:
  # <species>.<assembly>.<release>.<sequence type>.<id type>.<id>.fa
  # e.g. Homo_sapiens.GRCh37.64.dna_rm.chromosome.HG905_PATCH.fa
  #      Homo_sapiens.GRCh37.64.dna.chromosome.20.fa
  #      Ciona_savignyi.CSAV2.0.65.dna.toplevel.fa
  my @name_bits;
  push @name_bits, $self->web_name();
  push @name_bits, $self->assembly();
  push @name_bits, $self->param('release');
  push @name_bits, lc($data_type);
  push @name_bits, $level if $level;
  push @name_bits, $section if $section;
  push @name_bits, 'fa';

  my $file_name = join( '.', @name_bits );

  $data_type =~ s/_[rs]m$//;    # remove repeatmask or softmask designation from path component
  my $data_type_dir = $self->fasta_path($data_type);
  $self->_remove_files_from_dir($data_type_dir);
  return File::Spec->catfile( $data_type_dir, $file_name );
}

# Attempts to remove any generated files previously present for the instance
# of the Process
sub _remove_files_from_dir {
  my ($self, $dir) = @_;
  if(! $self->param('remove_files_from_dir')->{$dir}) {
    $self->unlink_all_files($dir);
    $self->param('remove_files_from_dir')->{$dir} = 1;
  }
  return;
}

##Logic used to generate the expected format for a FASTA header
sub _create_display_id {
  my ($self, $object, $seq, $type) = @_;
  
  my $stable_id;
  my $location;
  my $decoded_type;
  my $decoded_status;
  my %attributes;
  
  if(check_ref( $object, 'Bio::EnsEMBL::Transcript')) {
    $attributes{transcript_biotype} = $object->biotype();
    
    #If pred transcript then no gene but type & status are different
    if(check_ref($object, 'Bio::EnsEMBL::PredictionTranscript')) {
      $stable_id = $object->stable_id();
      $location  = $object->feature_Slice()->name();
      $decoded_type = $type;
      $decoded_status = lc($object->analysis()->logic_name());
      if($type eq 'pep') {
        $attributes{transcript} = $stable_id;
      }
    }
    #Must be a real "transcript"
    else {
      $stable_id = $object->stable_id();
      $location  = $object->feature_Slice()->name();
      my $gene = $object->get_Gene();
      $attributes{gene} = $gene->stable_id();
      $attributes{gene_biotype} = $gene->biotype();
      
      #If ncRNA then we set type to the logic name and status to gene's biotype (taken from original script)
      if($type eq 'ncrna') {
        $decoded_type = lc($object->analysis()->logic_name());
        $decoded_status = $gene->biotype();
      }
      elsif($object->biotype() =~ /pseudogene/i && ! $object->translation()) {
        $decoded_type = $type;
        $decoded_status = 'pseudogene';
      }
      #Otherwise use type & object's transcript's status
      else {
        $decoded_type   = $type;
        $decoded_status = lc($object->status());
      }
    }
  }
  #If it's a translation then grab the transcript and gene then set accordingly
  elsif(check_ref($object, 'Bio::EnsEMBL::Translation')) {
    my $transcript = $object->transcript();
    my $gene = $transcript->get_Gene();
    $stable_id = $object->stable_id();
    $location  = $transcript->feature_Slice()->name();
    %attributes = (
      gene => $gene->stable_id(),
      gene_biotype => $gene->biotype(),
      transcript => $transcript->stable_id(),
      transcript_biotype => $transcript->biotype()
    );
    $decoded_type = 'pep';
    $decoded_status = lc($transcript->status());
  }
  else {
    throw sprintf( 'Do not understand how to format a display_id for type "%s"',
      ref($object) );
  }
  
  my $attr_str = join(q{ }, 
    map { $_.':'.$attributes{$_} } 
    grep { exists $attributes{$_} } 
    qw/gene transcript gene_biotype transcript_biotype/);

  my $format = '%s %s:%s %s %s';
  
  my $id = sprintf( $format, $stable_id, $decoded_type, $decoded_status, $location, $attr_str);
  $seq->display_id($id);

  return;
}

sub _custom_header {
  my ($self) = @_;
  return sub {
    my $slice = shift;
    if ( !$slice->isa('Bio::EnsEMBL::Slice') ) {
      return $slice->display_id();
    }

    #RMS means masked data. soft_mask() true means it was softmasked
    my $dna_type = 'dna';
    if($slice->isa('Bio::EnsEMBL::RepeatMaskedSlice')) {
      $dna_type .= ($slice->soft_mask()) ? '_sm' : '_rm';
    }

    my $id        = $slice->seq_region_name;
    my $idtype    = $slice->coord_system->name;
    my $location  = $slice->name;
    my $ref       = $slice->assembly_exception_type();
    my $header    = sprintf('%s %s:%s %s %s', $id, $dna_type, $idtype, $location, $ref);  
    return $header;
  };
}

sub _final_filename {
  my ($self, $filename) = @_;
  return $filename if $filename =~ /\.gz$/;
  return $filename.'.gz';
}

sub assembly_accession {
  my ($self) = @_;
  my $mc = $self->get_DBAdaptor()->get_MetaContainer();
  return $mc->single_value_by_key('assembly.accession');
}

sub assembly_accession_type {
  my ($self) = @_;
  my $mc = $self->get_DBAdaptor()->get_MetaContainer();
  return $mc->single_value_by_key('assembly.web_accession_type');
}

sub _create_README {

  #Text for readme files

  my %text = (
    dna => <<'README',
#######################
Fasta DNA dumps
#######################

-----------
FILE NAMES
------------
The files are consistently named following this pattern:
   <species>.<assembly>.<release>.<sequence type>.<id type>.<id>.fa.gz

<species>:   The systematic name of the species.
<assembly>:  The assembly build name.
<release>:   The release number.
<sequence type>:
 * 'dna' - unmasked genomic DNA sequences.
  * 'dna_rm' - masked genomic DNA.  Interspersed repeats and low
     complexity regions are detected with the RepeatMasker tool and masked
     by replacing repeats with 'N's.
  * 'dna_sm' - soft-masked genomic DNA. All repeats and low complexity regions
    have been replaced with lowercased versions of their nucleic base
<id type> One of the following:
  * 'chromosome'     - The top-level coordinate system in most species in Ensembl
  * 'nonchromosomal' - Contains DNA that has not been assigned a chromosome
  * 'seqlevel'       - This is usually sequence scaffolds, chunks or clones.
     -- 'scaffold'   - Larger sequence contigs from the assembly of shorter
        sequencing reads (often from whole genome shotgun, WGS) which could
        not yet be assembled into chromosomes. Often more genome sequencing
        is needed to narrow gaps and establish a tiling path.
     -- 'chunk' -  While contig sequences can be assembled into large entities,
        they sometimes have to be artificially broken down into smaller entities
        called 'chunks'. This is due to limitations in the annotation
        pipeline and the finite record size imposed by MySQL which stores the
        sequence and annotation information.
     -- 'clone' - In general this is the smallest sequence entity.  It is often
        identical to the sequence of one BAC clone, or sequence region
        of one BAC clone which forms the tiling path.
<id>:     The actual sequence identifier. Depending on the <id type> the <id>
          could represent the name of a chromosome, a scaffold, a contig, a clone ..
          Field is empty for seqlevel files
fa: All files in these directories represent FASTA database files
gz: All files are compacted with GNU Zip for storage efficiency.


EXAMPLES
   The genomic sequence of human chromosome 1:
     Homo_sapiens.GRCh37.57.dna.chromosome.1.fa.gz

   The masked version of the genome sequence on human chromosome 1
   (contains '_rm' or '_sm' in the name):
     Homo_sapiens.GRCh37.57.dna_rm.chromosome.1.fa.gz
     Homo_sapiens.GRCh37.57.dna_sm.chromosome.1.fa.gz

   Non-chromosomal assembly sequences:
   e.g. mitochondrial genome, sequence contigs not yet mapped on chromosomes
     Homo_sapiens.GRCh37.57.dna.nonchromosomal.fa.gz
     Homo_sapiens.GRCh37.57.dna_rm.nonchromosomal.fa.gz
     Homo_sapiens.GRCh37.57.dna_sm.nonchromosomal.fa.gz

---------
TOPLEVEL
---------
These files contains all sequence regions flagged as toplevel in an Ensembl
schema. This includes chromsomes, regions not assembled into chromosomes and
N padded haplotype/patch regions.

EXAMPLES

  Toplevel sequences unmasked:
    Homo_sapiens.GRCh37.57.dna.toplevel.fa.gz
  
  Toplevel soft/hard masked sequences:
    Homo_sapiens.GRCh37.57.dna_sm.toplevel.fa.gz
    Homo_sapiens.GRCh37.57.dna_rm.toplevel.fa.gz

-----------------
PRIMARY ASSEMBLY
-----------------
Primary assembly contains all toplevel sequence regions excluding haplotypes
and patches. This file is best used for performing sequence similarity searches
where patch and haplotype sequences would confuse analysis.   

EXAMPLES

  Primary assembly sequences unmasked:
    Homo_sapiens.GRCh37.57.dna.primary_assembly.fa.gz
  
  Primary assembly soft/hard masked sequences:
    Homo_sapiens.GRCh37.57.dna_sm.primary_assembly.fa.gz
    Homo_sapiens.GRCh37.57.dna_rm.primary_assembly.fa.gz

--------------
SPECIAL CASES
--------------
Some chromosomes have alternate haplotypes which are presented in files with 
the haplotype sequence only:
   Homo_sapiens.GRCh37.56.dna_rm.chromosome.HSCHR6_MHC_QBL.fa.gz
   Homo_sapiens.GRCh37.56.dna_rm.chromosome.HSCHR17_1.fa.gz

All alternative assembly and patch regions have their sequence padded 
with N's to ensure alignment programs can report the correct index
regions

e.g. A patch region with a start position of 1,000,001 will have 1e6 N's added
its start so an alignment program will report coordinates with respect to the
whole chromosome.

Some species have sequenced Y chromosomes and the pseudoautosomal region (PAR)
on the Y is annotated.  By definition the PAR region is identical on the 
X and Y chromosome.  The Y chromosome file contains the Y chromsome 
minus the PAR regions.

README

    pep => <<'README',
####################
Fasta Peptide dumps
####################
These files hold the protein translations of Ensembl gene predictions.

-----------
FILE NAMES
------------
The files are consistently named following this pattern:
   <species>.<assembly>.<release>.<sequence type>.<status>.fa.gz

<species>:       The systematic name of the species.
<assembly>:      The assembly build name.
<release>:       The release number.
<sequence type>: pep for peptide sequences
<status>
  * 'pep.all' - the super-set of all translations resulting from Ensembl known
     or novel gene predictions.
  * 'pep.abinitio' translations resulting from 'ab initio' gene
     prediction algorithms such as SNAP and GENSCAN. In general, all
     'ab initio' predictions are based solely on the genomic sequence and
     not any other experimental evidence. Therefore, not all GENSCAN
     or SNAP predictions represent biologically real proteins.
fa : All files in these directories represent FASTA database files
gz : All files are compacted with GNU Zip for storage efficiency.

EXAMPLES (Note: Most species do not sequences for each different <status>)
 for Human:
    Homo_sapiens.NCBI36.40.pep.all.fa.gz
      contains all known and novel peptides
    Homo_sapiens.NCBI36.40.pep.abinitio.fa.gz
      contains all abinitio predicted peptide

Difference between known and novel
----------------------------------
Protein models that can be mapped to species-specific entries in
Swiss-Prot, RefSeq or SPTrEMBL are referred to in Ensembl as
known genes.  Those that cannot be mapped are called novel
(e.g. genes predicted on the basis of evidence from closely related species).

For models annotated by HAVANA the status is set manually. Models that have 
an HGNC name are referred to as known and the remaining models are referred to
as novel.

-------------------------------
FASTA Sequence Header Lines
------------------------------
The FASTA sequence header lines are designed to be consistent across
all types of Ensembl FASTA sequences.  This gives enough information
for the sequence to be identified outside the context of the FASTA
database file.

General format:

>ID SEQTYPE:STATUS LOCATION GENE TRANSCRIPT

Example of Ensembl Peptide header:

>ENSP00000328693 pep:novel chromosome:NCBI35:1:904515:910768:1 gene:ENSG00000158815:transcript:ENST00000328693 gene_biotype:protein_coding transcript_biotype:protein_coding
 ^               ^   ^     ^                                   ^                    ^                          ^                            ^ 
 ID              |   |  LOCATION                          GENE:stable gene ID       |                       GENE: gene biotype           TRANSCRIPT: transcript biotype
                 | STATUS                                           TRANSCRIPT: stable transcript ID
               SEQTYPE

README

    cdna => <<'README',
##################
Fasta cDNA dumps
#################

These files hold the cDNA sequences corresponding to Ensembl gene predictions.

------------
FILE NAMES
------------
The files are consistently named following this pattern:
<species>.<assembly>.<release>.<sequence type>.<status>.fa.gz

<species>: The systematic name of the species.
<assembly>: The assembly build name.
<release>: The release number.
<sequence type>: cdna for cDNA sequences
<status>
  * 'cdna.all' - the super-set of all transcripts resulting from
     Ensembl known, novel and pseudo gene predictions (see more below).
  * 'cdna.abinitio' - transcripts resulting from 'ab initio' gene prediction
     algorithms such as SNAP and GENSCAN. In general all 'ab initio'
     predictions are solely based on the genomic sequence and do not
     use other experimental evidence. Therefore, not all GENSCAN or SNAP
     cDNA predictions represent biologically real cDNAs.
     Consequently, these predictions should be used with care.

EXAMPLES  (Note: Most species do not sequences for each different <status>)
  for Human:
    Homo_sapiens.NCBI36.40.cdna.all.fa.gz
      cDNA sequences for all transcripts: known, novel and pseudo
    Homo_sapiens.NCBI36.40.cdna.abinitio.fa.gz
      cDNA sequences for 'ab-initio' prediction transcripts.

Difference between known and novel transcripts
-----------------------------------------------
Transcript or protein models that can be mapped to species-specific entries
in Swiss-Prot, RefSeq or SPTrEMBL are referred to as known genes in Ensembl.
Those that cannot be mapped are called novel genes (e.g. genes predicted on
the basis of evidence from closely related species).

For models annotated by HAVANA the status is set manually. Models that have 
an HGNC name are referred to as known and the remaining models are referred to
as novel.

-------------------------------
FASTA Sequence Header Lines
------------------------------
The FASTA sequence header lines are designed to be consistent across
all types of Ensembl FASTA sequences.  This gives enough information
for the sequence to be identified outside the context of the FASTA file.

General format:

>ID SEQTYPE:STATUS LOCATION GENE

Example of an Ensembl cDNA header:

>ENST00000289823 cdna:known chromosome:NCBI35:8:21922367:21927699:1 gene:ENSG00000158815 gene_biotype:protein_coding transcript_biotype:protein_coding
 ^               ^    ^     ^                                       ^                    ^                           ^       
 ID              |    |  LOCATION                         GENE: gene stable ID        GENE: gene biotype        TRANSCRIPT: transcript biotype
                 |  STATUS
              SEQTYPE


README

    ncrna => <<'README',
##################
Fasta RNA dumps
#################

These files hold the transcript sequences corresponding to non-coding RNA genes (ncRNA).

------------
FILE NAMES
------------
The files are consistently named following this pattern:
<species>.<assembly>.<release>.<sequence type>.fa.gz

<species>: The systematic name of the species.
<assembly>: The assembly build name.
<release>: The release number.
<sequence type>: ncrna for non-coding RNA sequences

EXAMPLES
  for Human:
    Homo_sapiens.NCBI36.40.ncrna.fa.gz
      Transcript sequences for all ncRNA gene types.


-------------------------------
FASTA Sequence Header Lines
------------------------------
The FASTA sequence header lines are designed to be consistent across
all types of Ensembl FASTA sequences.  This gives enough information
for the sequence to be identified outside the context of the FASTA file.

General format:

>ENST00000347977 ncrna:miRNA chromosome:NCBI35:1:217347790:217347874:-1 gene:ENSG00000195671 gene_biotype:ncRNA transcript_biotype:ncRNA
   ^             ^     ^     ^                                          ^                    ^                           ^ 
   ID            |     |  LOCATION                            GENE: gene stable ID       GENE: gene biotype           TRANSCRIPT: transcript biotype   
                 |   STATUS
              SEQTYPE


README
  );

  my $warning = <<'README';
#### README ####

IMPORTANT: Please note you can download correlation data tables,
supported by Ensembl, via the highly customisable BioMart and
EnsMart data mining tools. See http://www.ensembl.org/biomart/martview or
http://www.ebi.ac.uk/biomart/ for more information.

README

  my ( $self, $data_type ) = @_;
  my $base_path = $self->fasta_path();
  my $path      = File::Spec->catfile( $base_path, $data_type, 'README' );
  my $accession = $self->assembly_accession();
  my $txt       = $text{$data_type};
  throw "Cannot find README text for type $data_type" unless $txt;

  #Add accession information if it is available
  if($data_type eq 'dna' && $accession) {
    my $type = $self->assembly_accession_type();
    $warning .= <<EXTRA;
The genome assembly represented here corresponds to $type 
$accession

EXTRA
  }
  
  work_with_file($path, 'w', sub {
    my ($fh) = @_;
    print $fh $warning;
    print $fh $txt;
    return;
  });
  return;
}

1;

