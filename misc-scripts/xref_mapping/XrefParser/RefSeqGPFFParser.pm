# Parse RefSeq GPFF files to create xrefs.

package XrefParser::RefSeqGPFFParser;

use strict;

use File::Basename;

use XrefParser::BaseParser;

use vars qw(@ISA);
@ISA = qw( XrefParser::BaseParser);

# --------------------------------------------------------------------------------
# Parse command line and run if being run directly

if (!defined(caller())) {

  if (scalar(@ARGV) != 1) {
    print "\nUsage: RefSeqGPFFParser.pm file.SPC <source_id>\n\n";
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

  if ($source_id < 1) {
    $source_id =  XrefParser::BaseParser->get_source_id_for_filename(basename($file));
  }
  if(!defined($species_id)){
    $species_id = XrefParser::BaseParser->get_species_id_for_filename($file);
  }

   XrefParser::BaseParser->upload_xrefs(create_xrefs($source_id, $file, $species_id));

}

# --------------------------------------------------------------------------------
# Parse file into array of xref objects
# There are 2 types of RefSeq files that we are interested in:
# - protein sequence files *.protein.faa
# - mRNA sequence files *.rna.fna
# Slightly different formats

sub create_xrefs {

  my ($source_id, $file, $species_id) = @_;

  my %name2species_id =  XrefParser::BaseParser->name2species_id();

  my %dependent_sources =  XrefParser::BaseParser->get_dependent_xref_sources();

  open(REFSEQ, $file) || die "Can't open RefSeqGPFF file $file\n";

  my @xrefs;

  local $/ = "\/\/\n";

  my $type;
  if($file =~ /protein/){
    $type = 'peptide';
  }
  elsif($file =~ /rna/){
    $type = 'dna';
  }
  else{
    die "could not work out type of sequences for fiel $file\n";
  }


  while (<REFSEQ>) {
    
    my $xref;
    
    my $entry = $_;
    chomp $entry;
    
    my ($species) = $entry =~ /\s+ORGANISM\s+(.*)\n/;
    $species = lc $species;
    $species =~ s/^\s*//g;
    $species =~ s/\s+/_/g;
    $species =~ s/\n//g;
    my $species_id_check = $name2species_id{$species};
    
    # skip xrefs for species that aren't in the species table
    if (defined ($species_id) and $species_id = $species_id_check) {
      
      my ($acc) = $entry =~ /ACCESSION\s+(\S+)/;
      my ($ver) = $entry =~ /VERSION\s+(\S+)/;
      my ($description) = $entry =~ /DEFINITION\s+([^[]+)/s;
      print $entry if (length($description) == 0);
      $description =~ s/\nACCESSION.*//s;
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

      my ($acc_no_ver,$ver) = split (/\./,$ver);

      $xref->{ACCESSION} = $acc;
      if($acc eq $acc_no_ver){
         $xref->{VERSION} = $ver;
      }
      else{
         print "$acc NE $acc_no_ver\n";
      }
      $xref->{LABEL} = $acc;
      $xref->{DESCRIPTION} = $description;
      $xref->{SOURCE_ID} = $source_id;
      $xref->{SEQUENCE} = $parsed_seq;
      $xref->{SEQUENCE_TYPE} = $type;
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

      if($mrna){
        my %mrna_dep;
        $mrna_dep{SOURCE_ID} = $source_id; # source is still RefSeq
        my ($mrna_acc,$mrna_ver) = split (/\./,$mrna);

        $mrna_dep{ACCESSION} = $mrna_acc;
        $mrna_dep{VERSION} = $mrna_ver;
        push @{$xref->{DEPENDENT_XREFS}}, \%mrna_dep;
      }
      push @xrefs, $xref;

    }# if defined species

  } # while <REFSEQ>

  close (REFSEQ);

  print "Read " . scalar(@xrefs) ." xrefs from $file\n";

  return @xrefs;

}

# --------------------------------------------------------------------------------

sub new {

  my $self = {};
  bless $self, "XrefParser::RefSeqGPFFParser";
  return $self;

}

# --------------------------------------------------------------------------------

1;
