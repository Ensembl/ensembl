# Parse UniProt/Swiss-Prot files to create xrefs.

package SwissProtParser;

use strict;
use POSIX qw(strftime);
use File::Basename;

use BaseParser;

use vars qw(@ISA);
@ISA = qw(BaseParser);

# --------------------------------------------------------------------------------
# Parse command line and run if being run directly

if (!defined(caller())) {

  if (scalar(@ARGV) != 1) {
    print "\nUsage: SwissProtParser.pm file.SPC\n\n";
    exit(1);
  }

  run($ARGV[0], -1);

}

# --------------------------------------------------------------------------------

sub run {

  my $self = shift if (defined(caller(1)));
  my $file = shift;
  my $source_id = shift;

  my ($species_id, $species_name) = get_species($file);

  if ($source_id < 1) { # if being run directly
    $source_id = BaseParser->get_source_id_for_filename(basename($file));
    print "Source id for $file: $source_id\n";
  }

  BaseParser->upload_xrefs(create_xrefs($source_id, $species_id, $file));

}

# --------------------------------------------------------------------------------
# Get species (id and name) from file
# For SwissProt files the filename is the taxonomy ID

sub get_species {

  my ($file) = @_;

  my ($taxonomy_id, $extension) = split(/\./, basename($file));

  my $sth = BaseParser->dbi()->prepare("SELECT species_id,name FROM species WHERE taxonomy_id=?");
  $sth->execute($taxonomy_id);
  my ($species_id, $species_name);
  while(my @row = $sth->fetchrow_array()) {
    $species_id = $row[0];
    $species_name = $row[1];
  }
  $sth->finish;

  if (defined $species_name) {

    print "Taxonomy ID " . $taxonomy_id . " corresponds to species ID " . $species_id . " name " . $species_name . "\n";

  } else {

    print "Cannot find species corresponding to taxonomy ID " . $species_id . " - check species table\n";
    exit(1);

  }

  return ($species_id, $species_name);

}

# --------------------------------------------------------------------------------
# Parse file into array of xref objects

sub create_xrefs {

  my ($source_id, $species_id, $file) = @_;

  open(SWISSPROT, $file) || die "Can't open Swissprot file $file\n";

  my @xrefs;

  my $previous_rs = $/;
  $/ = "\/\/\n";

  while (<SWISSPROT>) {

    # if an OX line exists, only store the xref if the taxonomy ID that the OX
    # line refers to is in the species table
    my ($ox) = $_ =~ /OX\s+[a-zA-Z_]+=(\d+);/;
    if (defined $ox) {
      my $taxon = $1;
      my %taxonomy2species_id = BaseParser->taxonomy2species_id();
      if (!exists $taxonomy2species_id{$taxon}) {
	print "Skipping xref for species with taxonomy ID $taxon\n";
	next;
      }
    }

    my $xref;
    my $acc;
    ($acc) = $_ =~ /AC\s+(.+);/; # may catch multiple ; separated accessions
    ($xref->{LABEL})    = $_ =~ /ID\s+(\w+)/;
    ($xref->{SPECIES_ID}) = $species_id;
    ($xref->{SOURCE_ID}) = $source_id;

    $xref->{SEQUENCE_TYPE} = 'peptide';
    $xref->{STATUS} = 'experimental';

    my ($description) = $_ =~ /DE\s(.*)\n/;
    $xref->{DESCRIPTION} = $description;

    # set accession (and synonyms if more than one)
    my @acc = split /;\s*/, $acc;
    $xref->{ACCESSION} = $acc[0];
    for (my $a=1; $a <= $#acc; $a++) {
      push(@{$xref->{"SYNONYMS"} }, $acc[$a]);
    }

    # extract sequence
    my ($seq) = $_ =~ /SQ\s+(.+)/s; # /s allows . to match newline
      my @seq_lines = split /\n/, $seq;
    my $parsed_seq = "";
    foreach my $x (@seq_lines) {
      $parsed_seq .= $x;
    }
    $parsed_seq =~ s/\/\///g;   # remove trailing end-of-record character
    $parsed_seq =~ s/\s//g;     # remove whitespace
    $parsed_seq =~ s/^.*;//g;   # remove everything before last ;

    $xref->{SEQUENCE} = $parsed_seq;
    #print "Adding " . $xref->{ACCESSION} . " " . $xref->{LABEL} ."\n";

    # dependent xrefs - only store those that are from sources listed in the source table
    my ($deps) = $_ =~ /(DR\s+.+)/s; # /s allows . to match newline
      my @dep_lines = split /\n/, $deps;
    foreach my $dep (@dep_lines) {
      if ($dep =~ /^DR\s+(.+)/) {
	my ($source, $acc, @dummy) = split /;\s*/, $1;
	my %dependent_sources = BaseParser->get_dependent_xref_sources(); # name-id hash
	if (exists $dependent_sources{$source}) {
	  # create dependent xref structure & store it
	  my %dep;
	  $dep{SOURCE_NAME} = $source;
	  $dep{SOURCE_ID} = $dependent_sources{$source};
	  $dep{ACCESSION} = $acc;
	  push @{$xref->{DEPENDENT_XREFS}}, \%dep; # array of hashrefs
	}
      }
    }


    push @xrefs, $xref;

  }

  close (SWISSPROT);

  $/ = $previous_rs;

  print "Read " . scalar(@xrefs) ." xrefs from $file\n";

  return @xrefs;

  #TODO - currently include records from other species - filter on OX line??
}

# --------------------------------------------------------------------------------

sub new {

  my $self = {};
  bless $self, "SwissProtParser";
  return $self;

}

# --------------------------------------------------------------------------------

1;
