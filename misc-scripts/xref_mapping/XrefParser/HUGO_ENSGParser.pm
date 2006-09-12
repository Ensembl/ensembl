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



  # becouse the direct mapping have no descriptions etc
  # we have to steal these fromt he previous HUGO parser.
  # This is why the order states this is after the other one.
  # maybe 1091 is not right maybe should use name = HUGO and priority = 2 ??

  my %label;
  my %version;
  my %description;

  my $dbi = dbi();  
  my $sql = "select accession, label, version,  description from xref where source_id = 1091";
  my $sth = $dbi->prepare($sql);
  $sth->execute();
  my ($acc, $lab, $ver, $desc);
  $sth->bind_columns(\$acc, \$lab, \$ver, \$desc);
  while (my @row = $sth->fetchrow_array()) {
    $label{$acc} = $lab;
    $version{$acc} = $ver;
    $description{$acc} = $desc;
  }
  $sth->finish;



  my %acc;
  while (<HUGO>) {

    my ($hgnc, $stable_id) = split;

    if(!defined($acc{$hgnc})){
      $acc{$hgnc} = 1;
      my $version ="";
      $line_count++;
      
      my $xref_id = $self->add_xref($hgnc, $version{$hgnc} , $label{$hgnc}||$hgnc , $description{$hgnc}, $source_id, $species_id);
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
