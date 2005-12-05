package XrefParser::WormPepParser;

use strict;
use File::Basename;

use XrefParser::BaseParser;

use vars qw(@ISA);
@ISA = qw(XrefParser::BaseParser);
my $xref_sth ;
my $dep_sth;

# wormpep.table file format:
#>B0025.1a	CE24759	vps-34	phosphatidylinositol 3-kinase	Confirmed	SW:Q9TXI7	AAF23184.1
#>B0025.1b	CE24760	vps-34		Confirmed	SW:Q9TXI6	AAF23185.1
#>B0025.1c	CE37691	vps-34		Confirmed	SW:Q5TYK9	AAV34807.1

# Just need direct xref between B0025.1a (=stable ID for C. Elegans) and CE24759

sub run {

  my ($self, $file, $source_id, $species_id) = @_;

  print STDERR "WORMPep source = $source_id\tspecies = $species_id\n";

  my $worm_source_id = XrefParser::BaseParser->get_source_id_for_source_name('wormpep_id');

  my $xref_sth = $self->dbi()->prepare("SELECT xref_id FROM xref WHERE accession=? AND source_id=$worm_source_id AND species_id=$species_id");

  open(PEP,"<".$file) || die "Could not open $file\n";

  my ($x_count, $d_count);

  while (<PEP>) {

    my ($transcript, $wb, $display)  = (split(/\t/,substr($_,1)))[0,1,2];

    #check if the CGC_id is present in wormpep.table
    if ($display !~ m/\w/){ # or if ($display == '')
      $display = $wb;
    }

    # reuse or create xref
    my $xref_id;
    $xref_sth->execute($wb);
    my $xref_id = ($xref_sth->fetchrow_array())[0];
    if (!$xref_id) {
      $xref_id = $self->add_xref($wb, undef, $display, "", $worm_source_id, $species_id);
      $x_count++;
    }

    # and direct xref
    $self->add_direct_xref($xref_id, $transcript, "transcript", "");

    $d_count++;
  }

  close (PEP);

  print "Added $d_count direct xrefs and $x_count xrefs\n";

}

sub new {

  my $self = {};
  bless $self, "XrefParser::WormPepParser";
  return $self;

}

1;

