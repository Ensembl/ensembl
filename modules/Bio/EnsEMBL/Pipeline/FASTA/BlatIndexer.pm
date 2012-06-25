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

Bio::EnsEMBL::Pipeline::FASTA::BlastIndexer

=head1 DESCRIPTION

Creates 2bit file of the given GZipped file. The resulting index
is created under the parameter location I<base_path> in blat/index. The filename
is prefixed with the port number of the blat server this file should be
run on.

The module also performs filtering of non-reference sequence regions
and can filter the redundant Y chromosome piece for human (as 2bit does
not like repeated sequence region names).

Allowed parameters are:

=over 8

=item file - The file to index

=item program - The location of the faToTwoBit program

=item port_offset - Value to add onto the species_id from the website DB
                    to name the file correctly

=item base_path - The base of the dumps

=item index     - The type of file to index; supported values are empty, 
                  I<dna>, I<dna_sm> or I<dna_rm>. If specified we will look for this
                  string in the filename surrounded by '.' e.g. .dna.

=back

The registry should also have a DBAdaptor for the website schema 
registered under the species B<multi> and the group B<web> for species_id to
Blat port number. 

=cut


package Bio::EnsEMBL::Pipeline::FASTA::BlatIndexer;

use strict;
use warnings;
use base qw/Bio::EnsEMBL::Pipeline::FASTA::Indexer/;

use File::Spec;
use File::stat;
use Bio::EnsEMBL::Utils::IO qw/work_with_file/;
use Bio::EnsEMBL::Utils::Exception qw/throw/;
use Bio::EnsEMBL::Registry;

sub param_defaults {
  my ($self) = @_;
  return {
    program => 'faToTwoBit',
    port_offset => 30000,
    'index' => 'dna', #or dna_rm and dna_sm
  };
}

sub fetch_input {
  my ($self) = @_;
  $self->assert_executable($self->param('program'));
  $self->assert_executable('zcat');
  $self->assert_executable('gunzip');
  return;
}

sub run {
  my ($self) = @_;
  if($self->run_indexing()) {
    $self->SUPER::run();
  }
  return;
}

sub run_indexing {
  my ($self) = @_;
  my $index = $self->param('index');
  if($index) {
    my $file = $self->param('file');
    return (index($file, ".${index}.") > -1) ? 1 : 0;
  }
  return 1;
}

sub index_file {
  my ($self, $file) = @_;
  
  my $target_file = $self->target_file();
  my $cmd = sprintf(q{%s %s %s}, 
    $self->param('program'), $file, $target_file);
  
  $self->info('About to run "%s"', $cmd);
  my $output = `$cmd 2>&1`;
  my $rc = $? >> 8;
  throw "Cannot run program '$cmd'. Return code was ${rc}. Program output was $output" if $rc;
  unlink $file or throw "Cannot remove the file '$file' from the filesystem: $!";
  
  #Check the file size. If it's 16 bytes then reject as that is an empty file for 2bit
  my $filesize = stat($target_file)->size();
  if($filesize <= 16) {
    unlink $file;
    my $msg = sprintf(
      'The file %s produced a 2bit file %d byte(s). Lower than 17 bytes therefore empty 2 bit file',
      $file, $filesize
    );
    $self->throw($msg);
  }
  
  return;
}

sub decompress {
  my ($self) = @_;
  
  #If we have no non-reference seq regions then use normal decompress
  if(! $self->has_non_refs()) {
    return $self->SUPER::decompress();
  }
  
  #Filter for non-refs
  my $source = $self->param('file');
  my $target_dir = $self->target_dir();
  my ($vol, $dir, $file) = File::Spec->splitpath($source);
  $file =~ s/.gz$//;
  my $target = File::Spec->catdir($target_dir, $file);
  
  my $allowed_regions = $self->allowed_regions();
  
  $self->info('Decompressing %s to %s', $source, $target);
  
  open my $src_fh, '-|', "zcat $source" or throw "Cannot decompress $source to $target";
  work_with_file($target, 'w', sub {
    my ($trg_fh) = @_;
    my $transfer = 0;
    while(my $line = <$src_fh>) {
      #HEADER
      if(index($line, '>') == 0) {
        #regex is looking for NNN:NNN:NNN:1:80:1 i.e. the name
        my ($name) = $line =~ />.+\s(.+:.+:.+:\d+:\d+:\d+)/; 
        $transfer = ($allowed_regions->{$name}) ? 1 : 0;
        if($transfer) {
          $self->info('%s was an allowed Slice', $name);
        }
        else {
          $self->info('%s will be skipped; not a reference Slice', $name);
        }
      }
      print $trg_fh $line if $transfer;
    }
  });
  close($src_fh);

  return $target;
}

sub allowed_regions {
  my ($self) = @_;
  my $filter_human = 1;
  my @slices = grep { $_->is_reference() } @{$self->get_Slices('core', $filter_human)};
  my %hash = map { $_->name() => 1 } @slices;
  return \%hash;
}

#Filename like 30061.Homo_sapiens.GRCh37.2bit
sub target_filename {
  my ($self) = @_;
  my $port = $self->blat_port();
  my $name = $self->web_name();
  my $assembly = $self->assembly();
  return join(q{.}, $port, $name, $assembly, '2bit');
}

sub target_file {
  my ($self) = @_;
  my $target_dir = $self->target_dir();
  my $target_filename = $self->target_filename();
  return File::Spec->catfile($target_dir, $target_filename);
  return;
}

sub target_dir {
  my ($self) = @_;
  my $index = $self->param('index');
  my $base_dir  = ($index eq 'dna')     ? 'blat' 
                : ($index eq 'dna_rm')  ? 'blat_rm' 
                : ($index eq 'dna_sm')  ? 'blat_sm' 
                :                         q{}
                ;
  $self->throw(sprintf('Cannot decode a directory from index type "%s"', $index));
  return $self->get_dir($base_dir, $self->param('index'));
}

sub blat_port {
  my ($self) = @_;
  my $dba = Bio::EnsEMBL::Registry->get_DBAdaptor('multi', 'web');
  my $id = $dba->dbc()->sql_helper()->execute_single_result(
    -SQL => 'select species_id from species where name =?',
    -PARAMS => [$self->web_name()]
  );
  return $id + $self->param('port_offset');
}

sub has_non_refs {
  my ($self) = @_;
  my $sql = <<'SQL';
select count(*)
from attrib_type at
join seq_region_attrib sra using (attrib_type_id)
join seq_region sr using (seq_region_id)
join coord_system cs using (coord_system_id)
where cs.species_id =?
and at.code =?
SQL
  my $dba = $self->get_DBAdaptor();
  return $dba->dbc()->sql_helper()->execute_single_result(
    -SQL => $sql, -PARAMS => [$dba->species_id(), 'non_ref']);
}

1;
