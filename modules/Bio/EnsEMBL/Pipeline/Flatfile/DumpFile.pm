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

Bio::EnsEMBL::Pipeline::Flatfile::DumpFile

=head1 DESCRIPTION

The main workhorse of the Flatfile dumping pipeline.

The script is responsible for creating the filenames of these target
files, taking data from the database and the formatting of the flat files
headers. The final files are all Gzipped at normal levels of compression.

Allowed parameters are:

=over 8

=item species - The species to dump

=item base_path - The base of the dumps

=item release - The current release we are emitting

=item type - The type of data we are emitting. Should be embl or genbank

=back

=cut

package Bio::EnsEMBL::Pipeline::Flatfile::DumpFile;

use strict;
use warnings;

use base qw(Bio::EnsEMBL::Pipeline::Flatfile::Base);

use Bio::EnsEMBL::Utils::Exception qw/throw/;
use Bio::EnsEMBL::Utils::SeqDumper;
use Bio::EnsEMBL::Utils::IO qw/gz_work_with_file work_with_file/;
use File::Path qw/rmtree/;

sub param_defaults {
  my ($self) = @_;
  return {
    supported_types => {embl => 1, genbank => 1},
  };
}

sub fetch_input {
  my ($self) = @_;
  
  my $type = $self->param('type');
  throw "No type specified" unless $type;
  throw "Unsupported type '$type' specified" unless $self->param('supported_types')->{$type};
  
  throw "Need a species" unless $self->param('species');
  throw "Need a release" unless $self->param('release');
  throw "Need a base_path" unless $self->param('base_path');
  
  return;
}

sub run {
  my ($self) = @_;
  
  my $root = $self->data_path();
  if(-d $root) {
    $self->info('Directory "%s" already exists; removing', $root);
    rmtree($root);
  }
  
  my $type = $self->param('type');
  my $target = "dump_${type}";
  my $seq_dumper = $self->_seq_dumper();
  
  my @chromosomes;
  my @non_chromosomes;
  foreach my $s (@{$self->get_Slices()}) {
    my $chr = $s->is_chromosome();
    push(@chromosomes, $s) if $chr;
    push(@non_chromosomes, $s) if ! $chr;
  }
  
  if(@non_chromosomes) {
    my $path = $self->_generate_file_name('nonchromosomal');
    $self->info('Dumping non-chromosomal data to %s', $path);
    gz_work_with_file($path, 'w', sub {
      my ($fh) = @_;
      foreach my $slice (@non_chromosomes) {
        $self->fine('Dumping non-chromosomal %s', $slice->name());
        $seq_dumper->$target($slice, $fh);
      }
      return;
    });
  }
  else {
    $self->info('Did not find any non-chromosomal data');
  }
  
  foreach my $slice (@chromosomes) {
    $self->fine('Dumping chromosome %s', $slice->name());
    my $path = $self->_generate_file_name($slice->coord_system_name(), $slice->seq_region_name());
    my $args = {};
    if(-f $path) {
      $self->fine('Path "%s" already exists; appending', $path);
      $args->{Append} = 1;
    }
    gz_work_with_file($path, 'w', sub {
      my ($fh) = @_;
      $seq_dumper->$target($slice, $fh);
      return;
    }, $args);
  }
  
  $self->_create_README();
  
  return;
}

sub _seq_dumper {
  my ($self) = @_;
  my $seq_dumper = Bio::EnsEMBL::Utils::SeqDumper->new();
  $seq_dumper->disable_feature_type('similarity');
  $seq_dumper->disable_feature_type('genscan');
  $seq_dumper->disable_feature_type('variation');
  $seq_dumper->disable_feature_type('repeat');
  return $seq_dumper;
}

sub _generate_file_name {
  my ($self, $section, $name) = @_;

  # File name format looks like:
  # <species>.<assembly>.<release>.<section.name|section>.dat.gz
  # e.g. Homo_sapiens.GRCh37.64.chromosome.20.dat.gz
  #      Homo_sapiens.GRCh37.64.nonchromosomal.dat.gz
  my @name_bits;
  push @name_bits, $self->web_name();
  push @name_bits, $self->assembly();
  push @name_bits, $self->param('release');
  push @name_bits, $section if $section;
  push @name_bits, $name if $name;
  push @name_bits, 'dat', 'gz';

  my $file_name = join( '.', @name_bits );
  my $path = $self->data_path();
  return File::Spec->catfile($path, $file_name);
}

sub _create_README {
  my ($self) = @_;
  my $species = $self->scientific_name();
  my $format = uc($self->param('type'));
  
  my $readme = <<README;
#### README ####

IMPORTANT: Please note you can download correlation data tables, 
supported by Ensembl, via the highly customisable BioMart and 
EnsMart data mining tools. See http://www.ensembl.org/biomart/martview or
http://www.ebi.ac.uk/biomart/ for more information.

-----------------------
$format FLATFILE DUMPS
-----------------------
This directory contains $species $format flatfile dumps.  To ease 
downloading of the files, the $format format entries are bundled 
into groups of chromosomes and non-chromosomal regions.  
All files are then compacted with gzip.

Ensembl provides an automatic reannotation of $species genomic data.
These data will be dumped in a number of forms - one of them being 
$format flat files.  As the annotation of this form comes from Ensembl, 
and not the original sequence entry, the two annotations are 
likely to be different.

$format flat file format dumping provides all the confirmed protein coding 
genes known by Ensembl. Considerably more information is stored in Ensembl: 
the flat file just gives a representation which is compatible with 
existing tools.

The main body of the entry gives the same information as is in the main 
$format flat file entry.

    * ID - the $format id
    * AC - the EMBL/GenBank/DDBJ accession number (only the primary 
           accession number used)
    * SV - The accession.version pair which gives the exact reference to 
           a particular sequence
    * CC - comment lines to help you interpret the entry 

Currently the following features are dumped into the feature table of 
the Ensembl entry:

    * Transcripts as CDS entries. Each transcript has the following 
      attributes attached
          o Transcript id - a stable id, which Ensembl will attempt to 
            preserve as sensibly as possible during updates of the data
          o Gene id - indication of the gene that this transcript belongs 
            to. gene ids are stable and preserved as sensibly as possible 
            during updates of the data
          o Translation - the peptide translation of the transcript. 
    * Exons as exon entries. Each exon has the following information
          o Exon id. The exon id is stable and preserved as sensibly 
            as possible during sequence updates
          o start_phase. The phase of the splice site at the 5' end 
            of the exon. Phase 0 means between two codons, phase 1 
            means between the first and the second base of the codon 
            (meaning that there are 2 bases until the reading frame of 
            the exon) and phase 2 means between the second and the third 
            base of the codon (one base until the reading frame starts).
          o end_phase. The phase of the splice site at the 3' end of the 
            exon: same definition as above (though of course, being end_phase, 
            the position relative to the exon's reading frame is different 
            for phase 1 and 2). 

We are considering other information that should be made dumpable. In 
general we would prefer people to use database access over flat file 
access if you want to do something serious with the data. 

README
  
  my $path = File::Spec->catfile($self->data_path(), 'README');
  work_with_file($path, 'w', sub {
    my ($fh) = @_;
    print $fh $readme;
    return;
  });
  return;
}


1;

