# Parse UniProt (SwissProt & SPTrEMBL) files to create xrefs.
#
# Files actually contain both types of xref, distinguished by ID line;
#
# ID   143B_HUMAN     STANDARD;      PRT;   245 AA.         SwissProt
# ID   Q9B1S6      PRELIMINARY;      PRT;   260 AA.         SPTrEMBL

package XrefParser::UniProtParser;

use strict;
use POSIX qw(strftime);
use File::Basename;

use XrefParser::BaseParser;

use vars qw(@ISA);
@ISA = qw(XrefParser::BaseParser);

# --------------------------------------------------------------------------------
# Parse command line and run if being run directly

if (!defined(caller())) {

  if (scalar(@ARGV) != 3) {
    print "\nUsage: UniProtParser.pm file.SPC <source_id> <species_id>\n\n";
    print scalar(@ARGV);
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
  my $species_name;

  my ($sp_source_id, $sptr_source_id);

  if(!defined($species_id)){
    ($species_id, $species_name) = get_species($file);
  }
  $sp_source_id = XrefParser::BaseParser->get_source_id_for_source_name('Uniprot/SWISSPROT');
  $sptr_source_id = XrefParser::BaseParser->get_source_id_for_source_name('Uniprot/SPTREMBL');
  print "SwissProt source id for $file: $sp_source_id\n";
  print "SpTREMBL source id for $file: $sptr_source_id\n";

  my @xrefs = create_xrefs($sp_source_id, $sptr_source_id, $species_id, $file);

  # delete previous if running directly rather than via BaseParser
  if (!defined(caller(1))) {
    print "Deleting previous xrefs for these sources\n";
    XrefParser::BaseParser->delete_by_source(\@xrefs);
  }

  # upload
  XrefParser::BaseParser->upload_xref_object_graphs(@xrefs);

}

# --------------------------------------------------------------------------------
# Get species (id and name) from file
# For UniProt files the filename is the taxonomy ID

sub get_species {

  my ($file) = @_;

  my ($taxonomy_id, $extension) = split(/\./, basename($file));

  my $sth = XrefParser::BaseParser->dbi()->prepare("SELECT species_id,name FROM species WHERE taxonomy_id=?");
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

  my $num_sp = 0;
  my $num_sptr = 0;

  my %dependent_sources = XrefParser::BaseParser->get_dependent_xref_sources(); # name-id hash

  open(UNIPROT, $file) || die "Can't open Swissprot file $file\n";

  my @xrefs;

  local $/ = "\/\/\n";

  while (<UNIPROT>) {

    # if an OX line exists, only store the xref if the taxonomy ID that the OX
    # line refers to is in the species table
    # due to some records having more than one tax_id, we need to check them 
    # all and only proceed if one of them matches.
    #OX   NCBI_TaxID=158878, 158879;
    #OX   NCBI_TaxID=103690;


    my ($ox) = $_ =~ /OX\s+[a-zA-Z_]+=([0-9 ,]+);/;
#    print "OX  --> $ox\n";
    my @ox = split /\, /, $ox;
    my $found = 0;

    my %taxonomy2species_id = XrefParser::BaseParser->taxonomy2species_id();
    foreach my $taxon_id_from_file (@ox){
#      print "taxon_id= ".$taxon_id_from_file."\n";
      if (exists $taxonomy2species_id{$taxon_id_from_file} 
	  and $taxonomy2species_id{$taxon_id_from_file} eq $species_id) {
#	print "PASS ".$taxon_id_from_file."\n";
	$found = 1;	
      }
#      else{
#	print "FAIL ".$taxon_id_from_file."\n";
#      }
    }
    next if (!$found); # no taxon_id's math, so skip to next record
    my $xref;

    # set accession (and synonyms if more than one)
    # AC line may have primary accession and possibly several ; separated synonyms
    # May also be more than one AC line
    my ($acc) = $_ =~ /(AC\s+.+)/s; # will match first AC line and everything else
    my @all_lines = split /\n/, $acc;

    # extract ^AC lines only & build list of accessions
    my @accessions;
    foreach my $line (@all_lines) {
      my ($accessions_only) = $line =~ /^AC\s+(.+)/;
      push(@accessions, (split /;\s*/, $accessions_only)) if ($accessions_only);
    }

    $xref->{ACCESSION} = $accessions[0];
    for (my $a=1; $a <= $#accessions; $a++) {
      push(@{$xref->{"SYNONYMS"} }, $accessions[$a]);
    }

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

    # May have multi-line descriptions
    my ($description_and_rest) = $_ =~ /(DE\s+.*)/s;
    @all_lines = split /\n/, $description_and_rest;

    # extract ^DE lines only & build cumulative description string
    my $description;
    foreach my $line (@all_lines) {
      my ($description_only) = $line =~ /^DE\s+(.+)/;
      $description .= $description_only if ($description_only);
      $description .= " ";
    }

    $description =~ s/^\s*//g;
    $description =~ s/\s*$//g;

    $xref->{DESCRIPTION} = $description;

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
	my ($source, $acc, @extra) = split /;\s*/, $1;
	if (exists $dependent_sources{$source}) {
	  # create dependent xref structure & store it
	  my %dep;
	  $dep{SOURCE_NAME} = $source;
	  $dep{LINKAGE_SOURCE_ID} = $xref->{SOURCE_ID};
	  $dep{SOURCE_ID} = $dependent_sources{$source};
	  $dep{ACCESSION} = $acc;
	  push @{$xref->{DEPENDENT_XREFS}}, \%dep; # array of hashrefs
	  if($dep =~ /EMBL/){
	    my ($protein_id) = $extra[0];
	    if($protein_id ne "-"){
	      my %dep2;
	      $dep2{SOURCE_NAME} = $source;
	      $dep2{SOURCE_ID} = $dependent_sources{protein_id};
	      $dep2{LINKAGE_SOURCE_ID} = $xref->{SOURCE_ID};
	      $dep2{ACCESSION} = $protein_id;
	      push @{$xref->{DEPENDENT_XREFS}}, \%dep2; # array of hashrefs
	    }
	  }
	}
      }
    }

    # store PUBMED and MEDLINE dependent xrefs too
    #my ($medline) = $_ =~ /RX\s+MEDLINE=(\d+);/;
    #if (defined $medline) {
    #
    #  my %medline_dep;
    #  $medline_dep{SOURCE_ID} = $dependent_sources{PUBMED};
    #  $medline_dep{LINKAGE_SOURCE_ID} = $xref->{SOURCE_ID};
    #  $medline_dep{ACCESSION} = $medline;
    #  push @{$xref->{DEPENDENT_XREFS}}, \%medline_dep;
    #
    #}

    #my ($pubmed) = $_ =~ /RX\s+PubMed=(\d+);/;
    #if (defined $pubmed) {
    #
    #  my %pubmed_dep;
    #  $pubmed_dep{SOURCE_ID} = $dependent_sources{PUBMED};
    #  $pubmed_dep{LINKAGE_SOURCE_ID} = $xref->{SOURCE_ID};
    #  $pubmed_dep{ACCESSION} = $pubmed;
    #  push @{$xref->{DEPENDENT_XREFS}}, \%pubmed_dep;
    #
    #}

    push @xrefs, $xref;

  }

  close (UNIPROT);

  print "Read $num_sp SwissProt xrefs and $num_sptr SPTrEMBL xrefs from $file\n";

  return \@xrefs;

  #TODO - currently include records from other species - filter on OX line??
}

# --------------------------------------------------------------------------------

sub new {

  my $self = {};
  bless $self, "XrefParser::UniProtParser";
  return $self;

}

# --------------------------------------------------------------------------------

1;
