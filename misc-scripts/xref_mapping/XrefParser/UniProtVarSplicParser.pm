package XrefParser::UniProtVarSplicParser;

# Parse UniProt alternative splice files

use strict;
use File::Basename;

use XrefParser::BaseParser;

use vars qw(@ISA);
@ISA = qw(XrefParser::BaseParser);

# UniProtVarSplic file format: fasta, e.g.

# >P48347-2 (P48347) Splice isoform 2 of P48347
# MENEREKQVYLAKLSEQTERYDEMVEAMKKVAQLDVELTVEERNLVSVGYKNVIGARRAS
# WRILSSIEQKEESKGNDENVKRLKNYRKRVEDELAKVCNDILSVIDKHLIPSSNAVESTV
# FFYKMKGDYYRYLAEFSSGAERKEAADQSLEAYKAAVAAAENGLAPTHPVRLGLALNFSV
# FYYEILNSPESACQLAKQAFDDAIAELDSLNEESYKDSTLIMQLLRDNLTLWTSDLNEEG
# DERTKGADEPQDEV

sub run {

  my ($self, $file, $source_id, $species_id) = @_;

  my @xrefs;

  local $/ = "\n>";

  open(FILE,"<".$file) || die "Could not open $file\n";

  my $species_tax_id = $self->get_taxonomy_from_species_id($species_id);

  while (<FILE>) {

    my $xref;

    my ($header, $sequence) = $_ =~ /^>?(.+?)\n([^>]*)/s or warn("Can't parse FASTA entry: $_\n");

    # deconstruct header
    my ($accession, $original, @description) = split /\s+/, $header;
    my $description = join(" ", @description);

    # make sequence into one long string
    $sequence =~ s/\n//g;

    # build the xref object and store it
    $xref->{ACCESSION}     = $accession;
    $xref->{LABEL}         = $accession;
    $xref->{DESCRIPTION}   = $description;
    $xref->{SEQUENCE}      = $sequence;
    $xref->{SOURCE_ID}     = $source_id;
    $xref->{SPECIES_ID}    = $species_id;
    $xref->{SEQUENCE_TYPE} = 'peptide';
    $xref->{STATUS}        = 'experimental';

    push @xrefs, $xref;

  }

  close (FILE);

  print scalar(@xrefs) . " UniProtVarSplic xrefs succesfully parsed\n";

  XrefParser::BaseParser->upload_xref_object_graphs(\@xrefs);

  print "Done\n";

}


sub new {

  my $self = {};
  bless $self, "XrefParser::UniProtVarSplicParser";
  return $self;

}

1;
