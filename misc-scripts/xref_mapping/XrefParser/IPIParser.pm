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

package XrefParser::IPIParser;

use strict;
use warnings;
use Carp;
use File::Basename;

use base qw( XrefParser::BaseParser );

# IPI file format: fasta, e.g.
# >IPI:IPI00000005.1|SWISS-PROT:P01111|TREMBL:Q15104|REFSEQ_NP:NP_002515|ENSEMBL:ENSP00000261444 Tax_Id=9606 Transforming protein N-Ras
# MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAG
# PTRTVDTKQAHELAKSYGIPFIETSAKTRQGVEDAFYTLVREIRQYRMKKLNSSDDGTQG
# CMGLPCVVM

sub run {

  my ($self, $ref_arg) = @_;
  my $source_id    = $ref_arg->{source_id};
  my $species_id   = $ref_arg->{species_id};
  my $files        = $ref_arg->{files};
  my $verbose      = $ref_arg->{verbose};

  if((!defined $source_id) or (!defined $species_id) or (!defined $files) ){
    croak "Need to pass source_id, species_id, files and rel_file as pairs";
  }
  $verbose |=0;

  my $file = @{$files}[0];

  my @xrefs;

  local $/ = "\n>";

  my $ipi_io = $self->get_filehandle($file);

  if ( !defined $ipi_io ) {
    print STDERR "ERROR: Could not open $file\n";
    return 1;    # 1 = error
  }

  my %species_tax_id = %{$self->get_taxonomy_from_species_id($species_id)};

  while ( $_ = $ipi_io->getline() ) {
    my $xref;

    my ($header, $sequence) = $_ =~ /^>?(.+?)\n([^>]*)/s or warn("Can't parse FASTA entry: $_\n");

    # deconstruct header
    my @header = split /\|/, $header;
    my ($ipi) = $header[0] =~ /^IPI:(IPI(\d)+(\.\d+)?)/ or carp ("Can't deduce IPI identifier from " .  $header[0]);
    my ($ipi_ac, $ipi_ver) = $ipi =~ /(IPI\d+)\.(\d+)/;

    my ($tax_id, $description) = $header[-1] =~ /.*Tax_Id=(\d+)\s+(.*)/;

    # note currently we ignore all the other cross-references in the IPI file

    # only interested in species with the taxonomy ID were looking for
    next if ( !defined $tax_id || !defined $species_tax_id{$tax_id});

    # make sequence into one long string
    $sequence =~ s/\n//g;

    # build the xref object and store it
    $xref->{ACCESSION}     = $ipi_ac;
    $xref->{VERSION}       = $ipi_ver;
    $xref->{LABEL}         = $ipi;
    $xref->{DESCRIPTION}   = $description;
    $xref->{SEQUENCE}      = $sequence;
    $xref->{SOURCE_ID}     = $source_id;
    $xref->{SPECIES_ID}    = $species_id;
    $xref->{SEQUENCE_TYPE} = 'peptide';
    $xref->{STATUS}        = 'experimental';
    $xref->{INFO_TYPE}     = "SEQUENCE_MATCH";

    push @xrefs, $xref;

  }

  $ipi_io->close();


  $self->upload_xref_object_graphs(\@xrefs);

  print scalar(@xrefs) . " IPI xrefs succesfully parsed\n" if($verbose);

  return 0; #successful
}

1;
