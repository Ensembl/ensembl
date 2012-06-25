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

Bio::EnsEMBL::Pipeline::FASTA::ConcatFiles

=head1 DESCRIPTION

Performs a find in the DNA dumps directory for the given species and then
concats files which match a specified name pattern. We only allow
two types of concats; DNA and RM DNA. The concat file is a series
of cat command calls from all other Gzipped FASTA dumps (allowed under
the GZip specification). 

Allowed parameters are:

=over 8

=item release - Needed to build the target path

=item species - Required to indicate which species we are working with

=item data_type - The type of data to work with. Can be I<dna>, I<dn_sm> or I<dna_rm>

=item base_path - The base of the dumps

=back

=cut

package Bio::EnsEMBL::Pipeline::FASTA::ConcatFiles;

use strict;
use warnings;
use base qw/Bio::EnsEMBL::Pipeline::FASTA::Base/;

use File::Spec;
use File::stat;

sub param_defaults {
  my ($self) = @_;
  return {
    dna => {
      regex => qr/.+\.dna\..+\.fa\.gz$/,
    },
    dna_rm => {
      regex => qr/.+\.dna_rm\..+\.fa\.gz$/,
    },
    dna_sm => {
      regex => qr/.+\.dna_sm\..+\.fa\.gz$/,
    },
  };
}

sub fetch_input {
  my ($self) = @_;
  foreach my $key (qw/data_type species release base_path/) {
    $self->throw("Cannot find the required parameter $key") unless $self->param($key);
  }
  return;
}

# sticks ends of files together into one big file.
sub run {
  my ($self) = @_;
  
  my @file_list = @{$self->get_dna_files()};
  my $count = scalar(@file_list);
  my $running_total_size = 0;
  
  if($count) {
    my $target_file = $self->target_file();
    $self->info("Concatting type %s with %d file(s) into %s", $self->param('data_type'), $count, $target_file);
    
    if(-f $target_file) {
      $self->info("Target already exists. Removing");
      unlink $target_file or $self->throw("Could not remove $target_file: $!");
    }
    
    $self->info('Running concat');
    foreach my $file (@file_list) {
      $self->fine('Processing %s', $file);
      $running_total_size += stat($file)->size;
      system("cat $file >> $target_file") 
        and $self->throw( sprintf('Cannot concat %s into %s. RC %d', $file, $target_file, ($?>>8)));
    }

    $self->info("Catted files together");
    
    my $catted_size = stat($target_file)->size;
    
    if($running_total_size != $catted_size) {
      $self->throw(sprintf('The total size of the files catted together should be %d but was in fact %d. Failing as we expect the catted size to be the same', $running_total_size, $catted_size));
    }
    
    $self->param('target_file', $target_file);
  }
  else {
    $self->throw("Cannot continue as we found no files to concat");
  }
  return;
}

sub write_output {
  my ($self) = @_;
  my $file = $self->param('target_file');
  if($file) {
    $self->dataflow_output_id({ file => $file, species => $self->param('species') }, 1);
  }
  return;
}

sub get_dna_files {
  my ($self) = @_;
  my $path = $self->fasta_path('dna');
  my $data_type = $self->param('data_type'); 
  my $regex_hash = $self->param($data_type); 
  if(! $regex_hash ) {
    $self->throw("We do not have an entry for the data_type $data_type in our regex lookup hash. Edit this module");
  }
  my $regex = $regex_hash->{regex};
  my $filter = sub {
    my ($filename) = @_;
    return ($filename =~ $regex && $filename !~ /\.toplevel\./) ? 1 : 0;
  };
  my $files = $self->find_files($path, $filter);
  return [ sort @{$files} ];
}


sub target_file {
  my ($self) = @_;
  # File name format looks like:
  # <species>.<assembly>.<release>.<sequence type>.<id type>.<id>.fa.gz
  # e.g. Homo_sapiens.GRCh37.64.dna_rm.toplevel.fa.gz
  my @name_bits;
  push @name_bits, $self->web_name();
  push @name_bits, $self->assembly();
  push @name_bits, $self->param('release');
  push @name_bits, $self->param('data_type');
  push @name_bits, 'toplevel';
  push @name_bits, 'fa', 'gz';
  my $file_name = join( '.', @name_bits );
  my $dir = $self->fasta_path('dna');
  return File::Spec->catfile( $dir, $file_name );
}

1;
