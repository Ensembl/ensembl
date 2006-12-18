package XrefParser::HUGO_CCDSParser;

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
 
  $sql = 'select x.accession, d.ensembl_stable_id, d.type 
            from xref x, direct_xref d, source s 
             where s.source_id = x.source_id and 
                   x.xref_id = d.general_xref_id and s.name like "CCDS"'; 
 
  $sth = $dbi->prepare($sql);
  $sth->execute();
  my ($access, $stable_id, $type);
  $sth->bind_columns(\$access, \$stable_id, \$type);
  my %ensembl_stable_id;
  my %ensembl_type;
  while (my @row = $sth->fetchrow_array()) {
    $ensembl_stable_id{$access} = $stable_id;
    $ensembl_type{$access} = $type;
  }
  $sth->finish;
  
 
  my $line_count = 0;
  my $xref_count = 0;
  my %seen;
  my $ignore_count = 0;
  my $ignore_examples ="";
  while(<HUGO>){
    chomp;
    my ($ccds,$hgnc) = split;
    
    $line_count++;
    if(!defined($label{$hgnc})){
      $ignore_count++;
      if($ignore_count < 10){
	$ignore_examples .= " ".$hgnc;
      }
      next;
    }
    if(!defined($seen{$hgnc})){
      $seen{$hgnc} = 1;
      my $key = "CCDS".$ccds;
      if(defined($ensembl_stable_id{$key})){
	my $xref_id = $self->add_xref($hgnc, $version{$hgnc} , $label{$hgnc}||$hgnc , 
				      $description{$hgnc}, $source_id, $species_id);
	$self->add_direct_xref($xref_id, $ensembl_stable_id{$key}, $ensembl_type{$key}, "");
	$xref_count++;
      }
    }
  }
  print "Parsed $line_count HGNC identifiers from $file, added $xref_count xrefs and $xref_count direct_xrefs  from $line_count lines.\n";
  if($ignore_count){
    print $ignore_count." ignoreed due to numbers no identifiers being no longer valid :- $ignore_examples\n";
  }

  close(HUGO);
  return 0;

}


sub new {

  my $self = {};
  bless $self, "XrefParser::HUGO_CCDSParser";
  return $self;

}

1;
