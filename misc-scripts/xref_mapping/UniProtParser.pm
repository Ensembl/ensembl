# Parse UniProt (SwissProt & SPTrEMBL) files to create xrefs.
#
# Files actually contain both types of xref, distinguished by ID line;
#
# ID   143B_HUMAN     STANDARD;      PRT;   245 AA.         SwissProt
# ID   Q9B1S6      PRELIMINARY;      PRT;   260 AA.         SPTrEMBL

package UniProtParser;

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
    print "\nUsage: UniProtParser.pm file.SPC\n\n";
    exit(1);
  }

  run($ARGV[0], -1);

}

# --------------------------------------------------------------------------------

sub run {

  my $self = shift if (defined(caller(1)));
  my $file = shift;
  my $source_id = shift;

  my ($sp_source_id, $sptr_source_id);

  my ($species_id, $species_name) = get_species($file);

  $sp_source_id = BaseParser->get_source_id_for_source_name('UniProtSwissProt');
  $sptr_source_id = BaseParser->get_source_id_for_source_name('UniProtSPTrEMBL');
  print "SwissProt source id for $file: $sp_source_id\n";
  print "SpTREMBL source id for $file: $sptr_source_id\n";

  my @xrefs = create_xrefs($sp_source_id, $sptr_source_id, $species_id, $file);

  # delete previous if running directly rather than via BaseParser
  if (!defined(caller(1))) {
    print "Deleting previous xrefs for these sources\n";
    BaseParser->delete_by_source(\@xrefs);
  }

  # upload
  BaseParser->upload_xrefs(@xrefs);

}

# --------------------------------------------------------------------------------
# Get species (id and name) from file
# For UniProt files the filename is the taxonomy ID

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

  my ($sp_source_id, $sptr_source_id, $species_id, $file) = @_;

  my ($num_sp, $num_sptr) = 0;

  my %dependent_sources = BaseParser->get_dependent_xref_sources(); # name-id hash

  open(UNIPROT, $file) || die "Can't open Swissprot file $file\n";

  my @xrefs;

  local $/ = "\/\/\n";

  while (<UNIPROT>) {

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
    my ($label, $sp_type) = $_ =~ /ID\s+(\w+)\s+(\w+)/;

    # SwissProt/SPTrEMBL are differentiated by having STANDARD/PRELIMINARY here
    if ($sp_type =~ /STANDARD/i) {

      $xref->{SOURCE_ID} = $sp_source_id;
      $num_sp++;

    } elsif ($sp_type =~ /PRELIMINARY/i) {

      $xref->{SOURCE_ID} = $sptr_source_id;
      $num_sptr++;

    } else {

      next; # ignore if it's neither one nor t'other

    }

    $xref->{LABEL} = $label;
    $xref->{SPECIES_ID} = $species_id;
    $xref->{SEQUENCE_TYPE} = 'peptide';
    $xref->{STATUS} = 'experimental';

    my ($description) = $_ =~ /DE\s+(.*)\n/;
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
#	print $dep."\n";
	my ($source, $acc, @extra) = split /;\s*/, $1;
#	print "source is $source \n";
#	print "acc is $acc \n";
	if (exists $dependent_sources{$source}) {
#	  print "EXISTS\n";
	  # create dependent xref structure & store it
	  my %dep;
	  $dep{SOURCE_NAME} = $source;
	  $dep{SOURCE_ID} = $dependent_sources{$source};
	  $dep{ACCESSION} = $acc;
	  push @{$xref->{DEPENDENT_XREFS}}, \%dep; # array of hashrefs
	  if($dep =~ /EMBL/){
#	    print "prtein_id is ".$extra[0]."\n";
	    my ($protein_id) = $extra[0];
	    if($protein_id ne "-"){
	      my %dep2;
	      $dep2{SOURCE_NAME} = $source;
	      $dep2{SOURCE_ID} = $dependent_sources{protein_id};
	      $dep2{ACCESSION} = $protein_id;
	      push @{$xref->{DEPENDENT_XREFS}}, \%dep2; # array of hashrefs
	    }
	  }
	}
      }
    }

    # store PUBMED and MEDLINE dependent xrefs too
    my ($medline) = $_ =~ /RX\s+MEDLINE=(\d+);/;
    if (defined $medline) {

      my %medline_dep;
      $medline_dep{SOURCE_ID} = $dependent_sources{PUBMED};
      $medline_dep{ACCESSION} = $medline;
      push @{$xref->{DEPENDENT_XREFS}}, \%medline_dep;

    }

    my ($pubmed) = $_ =~ /RX\s+PubMed=(\d+);/;
    if (defined $pubmed) {

      my %pubmed_dep;
      $pubmed_dep{SOURCE_ID} = $dependent_sources{PUBMED};
      $pubmed_dep{ACCESSION} = $pubmed;
      push @{$xref->{DEPENDENT_XREFS}}, \%pubmed_dep;

    }

    push @xrefs, $xref;

  }

  close (UNIPROT);

  print "Read $num_sp SwissProt xrefs and $num_sptr SPTrEMBL xrefs from $file\n";

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
