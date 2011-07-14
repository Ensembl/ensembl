package XrefParser::IxodesCAPParser;

use strict;
use File::Basename;

use base qw( XrefParser::BaseParser );

# Ixodes:
#>ISCW800184-PA|GR35|Gustatory receptor|hughrobe|scaffold:IscaW1:DS849364:15863:44694:1
#MTFAYSQFRYSTRLLRWGGVWIVAEATNPGKQSFKTTLKRPYFWYCVLCLSTLVGTEFGN
#IIWALLFSFKHRKVFVSGVYTATQITVLVKTMLSSLMVALAAGRLKKLVARANQFEIIRN
#IKIAPRSKKVTWRDIRIWGRVLFMVLFVSIRNMDNLSILDVENIFGLGALVVVMTASSML

sub run {

  my $self = shift;
  my $source_id = shift;
  my $species_id = shift;
  my $files       = shift;
  my $release_file   = shift;
  my $verbose       = shift;

  my $file = @{$files}[0];

  next if (/^File:/);   # skip header

  my @xrefs;

  local $/ = "\n>";

  my $file_io = $self->get_filehandle($file);

  if ( !defined $file_io ) {
      print STDERR "Could not open $file\n";
      return 1;
  }

  while ( $_ = $file_io->getline() ) {
    my $xref;

    my ($header, $sequence) = $_ =~ /^>?(.+?)\n([^>]*)/s or warn("Can't parse FASTA entry: $_\n");

    # deconstruct header - just use first part
    my ($accession, $symbol, $description, $submitter, $position) = split /\|/, $header;
    if ($symbol eq "") { $symbol = "$accession" ; }

    # make sequence into one long string
    $sequence =~ s/\n//g;

    # build the xref object and store it
    $xref->{ACCESSION}     = $accession;
    $xref->{LABEL}         = $symbol;
    $xref->{DESCRIPTION}   = $description;
    $xref->{SEQUENCE}      = $sequence;
    $xref->{SOURCE_ID}     = $source_id;
    $xref->{SPECIES_ID}    = $species_id;
    $xref->{SEQUENCE_TYPE} = 'peptide';
    $xref->{STATUS}        = 'manual annotation';

    push @xrefs, $xref;

  }

  $file_io->close();


  XrefParser::BaseParser->upload_xref_object_graphs(\@xrefs);

  print scalar(@xrefs) . " Ixodes CAP xrefs succesfully parsed\n" if($verbose);

  return 0;
}

1;
