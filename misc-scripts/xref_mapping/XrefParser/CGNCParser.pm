package XrefParser::CGNCParser;

use strict;
use warnings;
use Carp;
use DBI;

use base qw(XrefParser::BaseParser);

sub run {

  my ($self, $ref_arg) = @_;
  my $source_id    = $ref_arg->{source_id};
  my $species_id   = $ref_arg->{species_id};
  my $files        = $ref_arg->{files};
  my $verbose      = $ref_arg->{verbose};

  if((!defined $source_id) or (!defined $species_id) or (!defined $files) ){
    croak "Need to pass source_id, species_id and files as pairs";
  }
  $verbose |=0;

  my $file = @{$files}[0];
  my $file_io = $self->get_filehandle($file);
  if ( !defined $file_io) {
    print STDERR "ERROR: Could not open $file\n";
    return 1;    # 1 is an error
  }

  my $head = $file_io->getline(); # first record are the headers
  chomp $head;
  # Sanity check on the headers, see if they are what is expected
  my @headers = split(/\t/,$head);
  if ($headers[0] ne 'CGNC id') {
    croak "First column " . $headers[0] . " does not contain CGNC id";
  }
  if ($headers[2] ne 'Ensembl id') {
    croak "Third column " . $headers[2] . " does not contain Ensembl id";
  }
  if ($headers[3] ne 'gene symbol') {
    croak "Fourth column " . $headers[3] . " does not contain gene symbol";
  }
  if ($headers[4] ne 'gene name') {
    croak "Fifth column " . $headers[4] . " does not contain gene name";
  }
  if ($headers[5] ne 'gene synonym') {
    croak "Sixth column " . $headers[5] . " does not contain gene synonym";
  }

  $source_id = $self->get_source_id_for_source_name("CGNC");

  my $dbi = $self->dbi();
  my $sql = "insert ignore into synonym (xref_id, synonym) values (?, ?)";
  my $add_syn_sth = $dbi->prepare($sql);
  my $count = 0;

  while ( $_ = $file_io->getline() ) {
## 2     426934  ENSGALG00000000003      PANX2   pannexin 2              Approved        2011-10-24
    chomp;
    my ($cgnc_id, $entrez_id, $ensid, $symbol, $gene_name, $synonym, $status, $edit_date) = split(/\t/, $_);

    if($ensid =~ /ENSGAL/){
      my $xref_id = $self->add_xref({ acc        => $cgnc_id,
				      version    => 0,
				      label      => $symbol,
				      desc       => $gene_name,
				      source_id  => $source_id,
				      species_id => $species_id,
				      info_type  => "DIRECT"} );

      $self->add_direct_xref( $xref_id, $ensid, "Gene", '');
      $count++;

      if ($synonym) {
        $add_syn_sth->execute($xref_id, $synonym);
      }
    }
  }

  print "$count direct CGNC xrefs added\n";
  return 0;

}

1;
