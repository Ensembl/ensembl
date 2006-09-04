package XrefParser::HUGO_ENSGParser;

use strict;

use DBI;

use XrefParser::BaseParser;

use vars qw(@ISA);
@ISA = qw(XrefParser::BaseParser);

# Parse file of HGNC records and assign direct xrefs
# All assumed to be linked to genes

sub run {

  my ($self, $file, $source_id, $species_id) = @_;

  if(!open(HUGO,"<".$file)){
    print "Could not open $file\n";
    return 1;
  }
  my $line_count = 0;
  my $xref_count = 0;

  my %acc;
  while (<HUGO>) {

    my ($hgnc, $stable_id) = split;

    if(!defined($acc{$hgnc})){
      $acc{$hgnc} = 1;
      my $version ="";
      $line_count++;
      
      my $xref_id = $self->add_xref($hgnc, $version, $hgnc, "", $source_id, $species_id);
      $xref_count++;
      

      $self->add_direct_xref($xref_id, $stable_id, "gene", "");
    }
  }

  print "Parsed $line_count HGNC identifiers from $file, added $xref_count xrefs and $line_count direct_xrefs\n";

  close(HUGO);
  return 0;
}


sub new {

  my $self = {};
  bless $self, "XrefParser::HUGO_ENSGParser";
  return $self;

}

1;
