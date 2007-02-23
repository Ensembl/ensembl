package XrefParser::HUGO_ENSGParser;

use strict;

use DBI;

use base qw( XrefParser::BaseParser );

# Parse file of HGNC records and assign direct xrefs
# All assumed to be linked to genes

sub run {

  my ($self, $file, $source_id, $species_id) = @_;

  my $hugo_io = $self->get_filehandle($file);

  if ( !defined $hugo_io ) {
    print "Could not open $file\n";
    return 1;
  }

  my $line_count = 0;
  my $xref_count = 0;

  # becouse the direct mapping have no descriptions etc
  # we have to steal these fromt he previous HUGO parser.
  # This is why the order states this is after the other one.
  # maybe 1091,1092 is not right maybe should use name = HUGO and priority = 30r4 ??

  my %label;
  my %version;
  my %description;

  my $dbi = $self->dbi();  
  my $sql = "select accession, label, version,  description from xref where source_id in (1091, 1092, 1094)";
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


  my $ignore_count = 0;
  my $ignore_examples ="";
  my %acc;

  while ( $_ = $hugo_io->getline() ) {

    my ($hgnc, $stable_id) = split;

    if(!defined($label{$hgnc})){
      $ignore_count++;
      if($ignore_count < 10){
	$ignore_examples .= " ".$hgnc;
      }
      next;
    }
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
  if($ignore_count){
    print $ignore_count." ignoreed due to numbers no identifiers being no longer valid :- $ignore_examples\n";
  }

  $hugo_io->close();

  return 0;
}

1;
