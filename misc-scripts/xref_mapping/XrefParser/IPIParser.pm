package XrefParser::IPIParser;

use strict;
use File::Basename;

use XrefParser::BaseParser;

use vars qw(@ISA);
@ISA = qw(XrefParser::BaseParser);

# IPI file format: fasta, e.g.
# >IPI:IPI00000005.1|SWISS-PROT:P01111|TREMBL:Q15104|REFSEQ_NP:NP_002515|ENSEMBL:ENSP00000261444 Tax_Id=9606 Transforming protein N-Ras
# MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAG
# QEEYSAMRDQYMRTGEGFLCVFAINNSKSFADINLYREQIKRVKDSDDVPMVLVGNKCDL
# PTRTVDTKQAHELAKSYGIPFIETSAKTRQGVEDAFYTLVREIRQYRMKKLNSSDDGTQG
# CMGLPCVVM

sub run {

  my $self = shift;
  my $file = shift;
  my $source_id = shift;
  my $species_id = shift;

  my @xrefs;

  local $/ = "\n>";

  open(IPI,"<".$file) || die "Could not open $file\n";

  my $species_tax_id = $self->get_taxonomy_from_species_id($species_id);

  while (<IPI>) {

    my $xref;

    my ($header, $sequence) = $_ =~ /^>?(.+?)\n([^>]*)/s or warn("Can't parse FASTA entry: $_\n");

    # deconstruct header
    my @header = split /\|/, $header;
    my ($ipi) = $header[0] =~ /^IPI:(IPI(\d)+(\.\d+)?)/ or warn("Can't deduce IPI identifier from " .  $header[0]);
    my ($ipi_ac, $ipi_ver) = $ipi =~ /(IPI\d+)\.(\d+)/;

    my ($tax_id, $description) = $header[-1] =~ /.*Tax_Id=(\d+)\s+(.*)/;

    # only interested in species with the taxonomy ID were looking for
    next if ($tax_id != $species_tax_id);

    # note currently we ignore all the other cross-references in the IPI file

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

    #print "IPI: $ipi tax: $tax_id desc: $description\n";

    push @xrefs, $xref;

  }

  print scalar(@xrefs) . " IPI xrefs succesfully parsed\n";

  XrefParser::BaseParser->upload_xref_object_graphs(\@xrefs);
  print "Done\n";

}


sub new {

  my $self = {};
  bless $self, "XrefParser::IPIParser";
  return $self;

}

1;
