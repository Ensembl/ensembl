package XrefParser::CeleraParser;

use strict;
use File::Basename;

use XrefParser::BaseParser;

use vars qw(@ISA);
@ISA = qw(XrefParser::BaseParser);

# Celera database dump for anopheles - FASTA format
#
# >agCP5429,cg_name=agCG43843,ga_name=GA_x9P1GAV56A9,transcript_name=agCT42178
# MNPNSTGSSSAAGSSISTSSLPGIERLIGRENWETWKFAVQTFLELEDLWCAVKPKKNDD
# GSYESVDTAKDRKARAKIILLLEPVNYVHVKEATTAKEVWSKLEKAFDDSGLTRRVGLLH
#
# This is the parser that provides most functionality, subclasses 
# (CeleraProteinParser, CeleraTranscriptParser) just set sequence type)

sub run {

  my ($self, $file, $source_id, $species_id) = @_;

  my @xrefs;

  local $/ = "\n>";

  open(FILE,"<".$file) || die "Could not open $file\n";

  while (<FILE>) {

    next if (/^File:/);   # skip header

    my $xref;

    my ($header, $sequence) = $_ =~ /^>?(.+?)\n([^>]*)/s or warn("Can't parse FASTA entry: $_\n");

    # deconstruct header - just use first part
    my ($accession, @rest) = split /,/, $header;

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

    push @xrefs, $xref;

  }

  close (FILE);

  print scalar(@xrefs) . " Celera xrefs succesfully parsed\n";

  XrefParser::BaseParser->upload_xref_object_graphs(\@xrefs);

  print "Done\n";

}


sub new {

  my $self = {};
  bless $self, "XrefParser::CeleraParser";
  return $self;

}

1;
