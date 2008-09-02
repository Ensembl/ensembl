package XrefParser::CeleraParser;

use strict;
use File::Basename;

use base qw( XrefParser::BaseParser );

# Celera database dump for anopheles - FASTA format
#
# >agCP5429,cg_name=agCG43843,ga_name=GA_x9P1GAV56A9,transcript_name=agCT42178
# MNPNSTGSSSAAGSSISTSSLPGIERLIGRENWETWKFAVQTFLELEDLWCAVKPKKNDD
# GSYESVDTAKDRKARAKIILLLEPVNYVHVKEATTAKEVWSKLEKAFDDSGLTRRVGLLH
#
# This is the parser that provides most functionality, subclasses 
# (CeleraProteinParser, CeleraTranscriptParser) just set sequence type)

sub run {

  my $self = shift if (defined(caller(1)));

  my $source_id = shift;
  my $species_id = shift;
  my $files       = shift;
  my $release_file   = shift;
  my $verbose       = shift;

  my $file = @{$files}[0];

  my $celera_gene_source_id = $self->get_source_id_for_source_name('Celera_Gene');

  my @xrefs;

  local $/ = "\n>";

  my $file_io = $self->get_filehandle($file);

  if ( !defined $file_io ) {
    print STDERR "Could not open $file\n";
    return 1;
  }

  while ( $_ = $file_io->getline() ) {
    next if (/^File:/);   # skip header

    my $xref;

    my ($header, $sequence) = $_ =~ /^>?(.+?)\n([^>]*)/s or warn("Can't parse FASTA entry: $_\n");

    # deconstruct header - just use first part
    my ($accession, $cg) = split /,/, $header;

    # make sequence into one long string
    $sequence =~ s/\n//g;

    # build the xref object and store it
    $xref->{ACCESSION}     = $accession;
    $xref->{LABEL}         = $accession;
    $xref->{SEQUENCE}      = $sequence;
    $xref->{SOURCE_ID}     = $source_id;
    $xref->{SPECIES_ID}    = $species_id;
    $xref->{SEQUENCE_TYPE} = $self->get_sequence_type();
    $xref->{STATUS}        = 'experimental';

    # pull cg_name from peptide files as well and create dependent xrefs
    if ($self->get_sequence_type() =~ /peptide/) {
      my ($cg_name) = $cg =~ /cg_name=(.*)/;
      my %dep;
      $dep{SOURCE_NAME} = 'Celera_Gene';
      $dep{LINKAGE_SOURCE_ID} = $xref->{SOURCE_ID};
      $dep{SOURCE_ID} = $celera_gene_source_id;
      $dep{ACCESSION} = $cg_name;
      push @{$xref->{DEPENDENT_XREFS}}, \%dep; # array of hashrefs
    }

    push @xrefs, $xref;

  }

  $file_io->close();


  XrefParser::BaseParser->upload_xref_object_graphs(\@xrefs);
  print scalar(@xrefs) . " Celera xrefs succesfully parsed\n" if($verbose);

  return 0;
}

1;
