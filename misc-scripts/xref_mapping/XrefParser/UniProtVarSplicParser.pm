package XrefParser::UniProtVarSplicParser;

# Parse UniProt alternative splice files

use strict;
use File::Basename;

use XrefParser::BaseParser;

use vars qw(@ISA);
@ISA = qw(XrefParser::BaseParser);

# UniProtVarSplic file format: fasta, e.g.

#>P48347-2|14310_ARATH Isoform 2 of P48347 - Arabidopsis thaliana (Mouse-ear cress)
#MENEREKQVYLAKLSEQTERYDEMVEAMKKVAQLDVELTVEERNLVSVGYKNVIGARRAS
#WRILSSIEQKEESKGNDENVKRLKNYRKRVEDELAKVCNDILSVIDKHLIPSSNAVESTV
#FFYKMKGDYYRYLAEFSSGAERKEAADQSLEAYKAAVAAAENGLAPTHPVRLGLALNFSV
#FYYEILNSPESACQLAKQAFDDAIAELDSLNEESYKDSTLIMQLLRDNLTLWTSDLNEEG
#DERTKGADEPQDEV

sub run {

  my ($self, $file, $source_id, $species_id) = @_;

  my @xrefs;

  local $/ = "\n>";

  if(!open(FILE,"<".$file)){
    print  "ERROR: Could not open $file\n";
    return 1; # 1 error
  }

  my $species_tax_id = $self->get_taxonomy_from_species_id($species_id);
  my (%swiss)  =  %{XrefParser::BaseParser->get_valid_codes("uniprot",$species_id)};
 
  my $missed = 0;
  while (<FILE>) {

    my $xref;

    my ($header, $sequence) = $_ =~ /^>?(.+?)\n([^>]*)/s or warn("Can't parse FASTA entry: $_\n");

    # deconstruct header
    my ($accession, @description) = split /\|/, $header;
    my $description = join(" ", @description);
    
    my ($original, $extension) = split/-/, $accession;

    if(defined($swiss{$original})){
      # make sequence into one long string
      $sequence =~ s/\n//g;
      
      # build the xref object and store it
      $xref->{ACCESSION}     = $accession;
      $xref->{LABEL}         = $accession;
      $xref->{DESCRIPTION}   = $description;
      $xref->{SEQUENCE}      = $sequence;
      $xref->{SOURCE_ID}     = $source_id;
      $xref->{SPECIES_ID}    = $species_id;
      $xref->{SEQUENCE_TYPE} = 'peptide';
      $xref->{STATUS}        = 'experimental';
      
      push @xrefs, $xref;
    }
    else{
     $missed++;
    }
  }

  close (FILE);

  print $missed." ignored as original uniprot not found in database\n";
  print scalar(@xrefs) . " UniProtVarSplic xrefs succesfully parsed\n";

  XrefParser::BaseParser->upload_xref_object_graphs(\@xrefs);

  print "Done\n";
  return 0;
}


sub new {

  my $self = {};
  bless $self, "XrefParser::UniProtVarSplicParser";
  return $self;

}

1;


