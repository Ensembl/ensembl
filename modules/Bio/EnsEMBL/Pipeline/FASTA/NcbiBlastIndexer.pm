=pod

=head1 LICENSE

  Copyright (c) 1999-2013 The European Bioinformatics Institute and
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

Bio::EnsEMBL::Pipeline::FASTA::NcbiBlastIndexer

=head1 DESCRIPTION

Creates NCBI Blast indexes of the given GZipped file. The resulting index
is created under the parameter location I<base_path> in ncbi_blast and then in a
directory defined by the type of dump. The type of dump also changes the file
name generated. Genomic dumps have their release number replaced with the
last repeat masked date. 

Allowed parameters are:

=over 8

=item file - The file to index

=item program - The location of the xdformat program

=item molecule - The type of molecule to index. I<dna> and I<pep> are allowed

=item type - Type of index we are creating. I<genomic> and I<genes> are allowed

=item base_path - The base of the dumps

=item release - Required for correct DB naming

=back

=cut

package Bio::EnsEMBL::Pipeline::FASTA::NcbiBlastIndexer;

use strict;
use warnings;
use base qw/Bio::EnsEMBL::Pipeline::FASTA::BlastIndexer/;

sub param_defaults {
  my ($self) = @_;
  return {
    program => 'makeblastdb',
    blast_dir => 'ncbi_blast',
  };
}

sub index_file {
  my ($self, $file) = @_;
  my $molecule_arg = ($self->param('molecule') eq 'dna') ? 'nucl' : 'prot' ;
  my $target_dir = $self->target_dir();
  my $target_file = $self->target_file($file);
  my $db_title = $self->db_title($file);
  
  my $cmd = sprintf(q{cd %s && %s -in %s -out %s -dbtype %s -title %s -input_type fasta}, 
    $target_dir, $self->param('program'), $file, $target_file, $molecule_arg, $db_title);
  $self->run_cmd($cmd);  
  unlink $file or $self->throw("Cannot remove the file '$file' from the filesystem: $!");
  $self->param('index_base', $target_file);
  return;
}

1;
