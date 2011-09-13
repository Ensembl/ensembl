package XrefParser::IxodesCAPParser;

use strict;
use warnings;
use Carp;
use File::Basename;

use base qw( XrefParser::BaseParser );

# Ixodes:
#>ISCW800184-PA|GR35|Gustatory receptor|hughrobe|scaffold:IscaW1:DS849364:15863:44694:1
#MTFAYSQFRYSTRLLRWGGVWIVAEATNPGKQSFKTTLKRPYFWYCVLCLSTLVGTEFGN
#IIWALLFSFKHRKVFVSGVYTATQITVLVKTMLSSLMVALAAGRLKKLVARANQFEIIRN
#IKIAPRSKKVTWRDIRIWGRVLFMVLFVSIRNMDNLSILDVENIFGLGALVVVMTASSML

sub run {

  my ($self, $ref_arg) = @_;
  my $source_id    = $ref_arg->{source_id};
  my $species_id   = $ref_arg->{species_id};
  my $files        = $ref_arg->{files};
  my $release_file = $ref_arg->{rel_file};
  my $verbose      = $ref_arg->{verbose};

  if((!defined $source_id) or (!defined $species_id) or (!defined $files) or (!defined $release_file)){
    croak "Need to pass source_id, species_id, files and rel_file as pairs";
  }
  $verbose |=0;

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
