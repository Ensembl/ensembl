package XrefParser::XenopusJamboreeParser;

# Parse annotated peptides from Xenopus Jamboree

use strict;
use File::Basename;

use base qw( XrefParser::BaseParser );



sub run {

  my $self = shift;
  my $source_id = shift;
  my $species_id = shift;
  my $files       = shift;
  my $release_file   = shift;
  my $verbose       = shift;

  my $file = @{$files}[0];

  my $file_io = $self->get_filehandle($file);

  if ( !defined $file_io ) {
    print STDERR "ERROR: Could not open $file\n";
    return 1;    # 1 error
  }

  my $count = 0;
  while ( $_ = $file_io->getline() ) {
    chomp;
    my ($acc, $label, $desc, $stable_id) = split /\t/;

    if($label eq "unnamed"){
      $label = $acc;
    }

    XrefParser::BaseParser->add_to_direct_xrefs($stable_id,'gene', $acc, '', $label, $desc, "", $source_id, $species_id);
    $count++;
  }

  $file_io->close();

  print $count . " XenopusJamboreeParser xrefs succesfully parsed\n" if($verbose);

  return 0;
}

1;
