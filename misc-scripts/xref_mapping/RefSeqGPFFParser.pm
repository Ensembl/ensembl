# Parse RefSeq GPFF files to create xrefs.

package RefSeqGPFFParser;

use strict;

use File::Basename;

use BaseParser;

use vars qw(@ISA);
@ISA = qw(BaseParser);

# --------------------------------------------------------------------------------
# Parse command line and run if being run directly

if (!defined(caller())) {

  if (scalar(@ARGV) != 1) {
    print "\nUsage: RefSeqGPFFParser.pm file.SPC\n\n";
    exit(1);
  }

  run($ARGV[0], -1);

}

# --------------------------------------------------------------------------------

sub run {

  my $self = shift if (defined(caller(1)));
  my $file = shift;
  my $source_id = shift;

  if ($source_id < 1) {
    $source_id = BaseParser->get_source_id_for_filename(basename($file));
    print "Source id for $file: $source_id\n";
  }

  BaseParser->upload_xrefs(create_xrefs($source_id, $file));

}

# --------------------------------------------------------------------------------
# Parse file into array of xref objects
# There are 2 types of RefSeq files that we are interested in:
# - protein sequence files *.protein.faa
# - mRNA sequence files *.rna.fna
# Slightly different formats

sub create_xrefs {

  my ($source_id, $file) = @_;

  my %name2species_id = BaseParser->name2species_id();

  my %dependent_sources = BaseParser->get_dependent_xref_sources();

  open(REFSEQ, $file) || die "Can't open RefSeqGPFF file $file\n";

  my @xrefs;

  local $/ = "\/\/\n";

  while (<REFSEQ>) {
    
    my $xref;
    
    my $entry = $_;
    chomp $entry;
    
    my ($species) = $entry =~ /\s+ORGANISM\s+(.*)\n/;
    $species = lc $species;
    $species =~ s/^\s*//g;
    $species =~ s/\s+/_/g;
    $species =~ s/\n//g;
    my $species_id = $name2species_id{$species};
    
    # skip xrefs for species that aren't in the species table
    if (defined $species_id) {
      
      my ($acc) = $entry =~ /ACCESSION\s+(\S+)/;
      my ($description) = $entry =~ /DEFINITION\s+([^[]*)/s;
      print $entry if (length($description) == 0);
      $description =~ s/\n//g;
      $description =~ s/\s+/ /g;
      $description = substr($description, 0, 255) if (length($description) > 255);
      
      my ($seq) = $_ =~ /ORIGIN\s+(.+)/s; # /s allows . to match newline
      my @seq_lines = split /\n/, $seq;
      my $parsed_seq = "";
      foreach my $x (@seq_lines) {
        my ($seq_only) = $x =~ /\s*\d+\s+(.*)/;
        $parsed_seq .= $seq_only;
      }
      $parsed_seq =~ s/\/\///g;   # remove trailing end-of-record character
      $parsed_seq =~ s/\s//g;     # remove whitespace

      $xref->{ACCESSION} = $acc;
      $xref->{LABEL} = $acc;
      $xref->{DESCRIPTION} = $description;
      $xref->{SOURCE_ID} = $source_id;
      $xref->{SEQUENCE} = $parsed_seq;
      $xref->{SEQUENCE_TYPE} = 'peptide';
      $xref->{SPECIES_ID} = $species_id;

      # TODO experimental/predicted

      # pubmed & medline are simple dependent xrefs; may be several of each
      my @medline = $entry =~ /\s+MEDLINE\s+(\d+)/g;
      my @pubmed = $entry =~ /\s+PUBMED\s+(\d+)/g;
      my @LocusIDline = $entry =~ /db_xref=.LocusID:(\d+)/g;
      my @mimline = $entry =~ /db_xref=.MIM:(\d+)/g;

      foreach my $ll (@LocusIDline) {
	my %dep;
	$dep{SOURCE_ID} = $dependent_sources{LocusLink};
	$dep{ACCESSION} = $ll;
	push @{$xref->{DEPENDENT_XREFS}}, \%dep;
      }
      foreach my $mim (@mimline) {
	my %dep;
	$dep{SOURCE_ID} = $dependent_sources{MIM};
	$dep{ACCESSION} = $mim;
	push @{$xref->{DEPENDENT_XREFS}}, \%dep;
      }
      foreach my $med (@medline) {
	my %dep;
	$dep{SOURCE_ID} = $dependent_sources{MEDLINE};
	$dep{ACCESSION} = $med;
	push @{$xref->{DEPENDENT_XREFS}}, \%dep;
      }
      foreach my $pub (@pubmed) {
	my %dep;
	$dep{SOURCE_ID} = $dependent_sources{PUBMED};
	$dep{ACCESSION} = $pub;
	push @{$xref->{DEPENDENT_XREFS}}, \%dep;
      }

      # Find associated mRNA
      my ($mrna) = $entry =~ /DBSOURCE\s+REFSEQ:\s+accession (.*)\n/;

      my %mrna_dep;
      $mrna_dep{SOURCE_ID} = $source_id; # source is still RefSeq
      $mrna_dep{ACCESSION} = $mrna;
      push @{$xref->{DEPENDENT_XREFS}}, \%mrna_dep;

      push @xrefs, $xref;

    } # if defined species

  } # while <REFSEQ>

  close (REFSEQ);

  print "Read " . scalar(@xrefs) ." xrefs from $file\n";

  return @xrefs;

}

# --------------------------------------------------------------------------------

sub new {

  my $self = {};
  bless $self, "RefSeqGPFFParser";
  return $self;

}

# --------------------------------------------------------------------------------

1;
