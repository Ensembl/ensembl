package XrefParser::IlluminaParser;

use strict;

use XrefParser::BaseParser;

use vars qw(@ISA);
@ISA = qw(XrefParser::BaseParser);

# Parser for Illumina refs
#Index,Target,Transcript,Accession,Symbol,Type,Definition,Ontology,Synony
#1,GI_23097300-A,GI_23097300,NM_021936.1,PLAC3,A,"Homo sapiens placenta-specific 3 (PLAC3), transcript variant 2, mRNA.",,PAPPE;PAPP-A
#2,GI_21070955-A,GI_21070955,NM_015386.1,COG4,A,"Homo sapiens component of oligomeric golgi complex 4 (COG4), mRNA.",,COD1;DKFZp586E1519

sub run {

  my ($self, $file, $source_id, $species_id) = @_;

  my $count = 0;
  my $missed = 0;

  if(!open(FILE,"<".$file)){
    print "Could not open $file\n";
    return 1;
  }

  my (%refseq) = %{XrefParser::BaseParser->get_valid_codes("refseq", $species_id)};

  while (<FILE>) {

    chomp;

    my $xref;

    # strip ^M at end of line
    $_ =~ s/\015//g;

    my ($idx, $target, $transcript, $accession, $symbol, $type, $definition, $ontology, $synonym) = split(/,/);
    next if ($idx == "Index" && $target == "Target"); # skip header

    my ($refseq) = $accession =~ /([^.]+)(\..*)?$/;

    if ($refseq{$refseq}) {
      XrefParser::BaseParser->add_to_xrefs($refseq{$refseq}, $target, '', $target, $definition, '', $source_id, $species_id);
      $count++;
    } else {
      $missed++;
    }


  }

  close(FILE);

  print "$count Illumina xrefs succesfully parsed. Skipped $missed due to non-matching RefSeq IDs.\n";

  return 0;
}


sub new {

  my $self = {};
  bless $self, "XrefParser::IlluminaParser";
  return $self;

}

1;
