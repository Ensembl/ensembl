# Parse RefSeq files to create xrefs.

package XrefParser::RefSeqParser;

use strict;

use File::Basename;

use XrefParser::BaseParser;

use vars qw(@ISA);
@ISA = qw(XrefParser::BaseParser);

# --------------------------------------------------------------------------------
# Parse command line and run if being run directly

if (!defined(caller())) {

  if (scalar(@ARGV) != 1) {
    print "\nUsage: RefSeqParser.pm file.SPC <source_id> <species_id>\n\n";
    exit(1);
  }

  run($ARGV[0], -1);

}

# --------------------------------------------------------------------------------

sub run {

  my $self = shift if (defined(caller(1)));
  my $file = shift;
  my $source_id = shift;
  my $species_id = shift;

  if (!defined($source_id) or $source_id < 1) {
    $source_id = XrefParser::BaseParser->get_source_id_for_filename(basename($file));
  }
  if(!defined($species_id)){
    $species_id = XrefParser::BaseParser->get_species_id_for_filename($file);
  }

  XrefParser::BaseParser->upload_xref_object_graphs(create_xrefs($source_id, $file, $species_id));

}

# --------------------------------------------------------------------------------
# Parse file into array of xref objects
# There are 2 types of RefSeq files that we are interested in:
# - protein sequence files *.protein.faa
# - mRNA sequence files *.rna.fna
# Slightly different formats

sub create_xrefs {

  my ($source_id, $file, $species_id) = @_;

  my %name2species_id = XrefParser::BaseParser->name2species_id();

  open(REFSEQ, $file) || die "Can't open RefSeq file $file\n";

  my @xrefs;

  local $/ = "\n>";

  while (<REFSEQ>) {

    my $xref;

    my $entry = $_;
    chomp $entry;
    my ($header, $sequence) = split (/\n/, $entry, 2);
    $sequence =~ s/^>//;
    # remove newlines
    my @seq_lines = split (/\n/, $sequence);
    $sequence = join("", @seq_lines);

    (my $gi, my $n, my $ref, my $acc, my $description) = split(/\|/, $header);
    my ($species, $mrna);
    if ($file =~ /\.faa$/) {

      ($mrna, $description, $species) = $description =~ /(\S*)\s+(.*)\s+\[(.*)\]$/;
      $xref->{SEQUENCE_TYPE} = 'peptide';
      $xref->{STATUS} = 'experimental';

    } elsif ($file =~ /\.fna$/) {

      ($species, $description) = $description =~ /\s*(\w+\s+\w+)\s+(.*)$/;
      $xref->{SEQUENCE_TYPE} = 'dna';
      $xref->{STATUS} = 'experimental';

    }

    $species = lc $species;
    $species =~ s/ /_/;

    my $species_id_check = $name2species_id{$species};

    # skip xrefs for species that aren't in the species table
    if (defined($species_id) and $species_id = $species_id_check) {

      my ($acc_no_ver,$ver) = split (/\./,$acc);
      $xref->{ACCESSION} = $acc_no_ver;
      $xref->{VERSION} = $ver;
      $xref->{LABEL} = $acc;
      $xref->{DESCRIPTION} = $description;
      $xref->{SOURCE_ID} = $source_id;
      $xref->{SEQUENCE} = $sequence;
      $xref->{SPECIES_ID} = $species_id;

      # TODO synonyms, dependent xrefs etc

      push @xrefs, $xref;

    }

  }

  close (REFSEQ);

  print "Read " . scalar(@xrefs) ." xrefs from $file\n";

  return \@xrefs;

}

# --------------------------------------------------------------------------------

sub new {

  my $self = {};
  bless $self, "XrefParser::RefSeqParser";
  return $self;

}

# --------------------------------------------------------------------------------

1;
