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

  run();

}

# --------------------------------------------------------------------------------

sub run {

  my $file = $ARGV[0];

  my ($species_id, $species_name) = get_species($file);

  $species_name =~ s/ /_/g;
  my $name = "UniProt_SwissProt_" . $species_name;

  my $source_id = BaseParser->get_source_id($name);

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
# Create source object to be loaded into source table

sub create_source {

  my ($file) = @_;
  my $source;
  my $file_date = POSIX::strftime('%Y%m%d%H%M%S', localtime((stat($file))[9]));
  $source = { NAME => "UniProt/SwissProt",
	      URL  => $file,
	      FILE_MODIFIED_DATE => $file_date
	      # TODO URL? Release?
	    };

  return $source;

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

    my $xref;
    my $acc;
    ($acc) =$_ =~ /AC\s+(.+);/; # may catch multiple ; separated accessions 
    ($xref->{LABEL})    = $_ =~ /ID\s+(\w+)/;
    ($xref->{SPECIES_ID}) = $species_id;
    ($xref->{SOURCE_ID}) = $source_id;

    # set accession (and synonyms if more than one)
    # note synonyms 
    my @acc = split /;/, $acc;
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

    push @xrefs, $xref;

  }

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
