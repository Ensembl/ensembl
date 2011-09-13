package XrefParser::WormPepParser;

use strict;
use warnings;
use Carp;
use File::Basename;

use base qw( XrefParser::BaseParser );

# wormpep.table file format:
#>B0025.1a	CE24759	vps-34	phosphatidylinositol 3-kinase	Confirmed	SW:Q9TXI7	AAF23184.1
#>B0025.1b	CE24760	vps-34		Confirmed	SW:Q9TXI6	AAF23185.1
#>B0025.1c	CE37691	vps-34		Confirmed	SW:Q5TYK9	AAV34807.1

# Just need direct xref between B0025.1a (=stable ID for C. Elegans) and CE24759

sub run {

  my ($self, $ref_arg) = @_;
  my $source_id    = $ref_arg->{source_id};
  my $species_id   = $ref_arg->{species_id};
  my $files        = $ref_arg->{files};
  my $verbose      = $ref_arg->{verbose};

  if((!defined $source_id) or (!defined $species_id) or (!defined $files)){
    croak "Need to pass source_id, species_id and files as pairs";
  }
  $verbose |=0;

  my $file = @{$files}[0];


  my $worm_source_id = XrefParser::BaseParser->get_source_id_for_source_name('wormpep_id');
  my $worm_locus_id = XrefParser::BaseParser->get_source_id_for_source_name('wormbase_locus');


  my $xref_sth = $self->dbi()->prepare("SELECT xref_id FROM xref WHERE accession=? AND source_id=$worm_source_id AND species_id=$species_id");
  my $xref_sth2 = $self->dbi()->prepare("SELECT xref_id FROM xref WHERE accession=? AND source_id=$worm_locus_id AND species_id=$species_id");

  my $pep_io = $self->get_filehandle($file);

  if ( !defined $pep_io ) {
    print STDERR "ERROR: Could not open $file\n";
    return 1;    # 1 error
  }

  my ($x_count, $d_count);


  while ( $_ = $pep_io->getline() ) {
    my ($transcript, $wb, $display)  = (split(/\t/,substr($_,1)))[0,1,2];

    # reuse or create xref
    $xref_sth->execute($wb);
    my $xref_id = ($xref_sth->fetchrow_array())[0];
    if (!$xref_id) {
      $xref_id = $self->add_xref($wb, undef, $wb, "", $worm_source_id, $species_id);
      $x_count++;
    }
    # and direct xref
    $self->add_direct_xref($xref_id, $transcript, "transcript", "");
    $d_count++;

    if(defined($display) and length($display) > 0 ){
      $xref_sth2->execute($display);
      my $xref_id2 = ($xref_sth2->fetchrow_array())[0];
      if (!$xref_id2) {
	$xref_id2 = $self->add_xref($display, undef, $display, "", $worm_locus_id, $species_id);
	$x_count++;
      }
      # and direct xref
      $self->add_direct_xref($xref_id2, $transcript, "transcript", "");
    }

    $d_count++;
  }

  $pep_io->close();

  print "Added $d_count direct xrefs and $x_count xrefs\n" if($verbose);
  return 0;
}

1;
