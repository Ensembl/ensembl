package XrefParser::XenopusJamboreeParser;

# Parse annotated peptides from Xenopus Jamboree

use strict;
use warnings;
use Carp;
use File::Basename;

use base qw( XrefParser::BaseParser );



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

    $self->add_to_direct_xrefs({ stable_id  => $stable_id,
				 type       => 'gene',
				 acc        => $acc,
				 label      => $label,
				 desc       => $desc,
				 source_id  => $source_id,
				 species_id => $species_id });
    $count++;
  }

  $file_io->close();

  print $count . " XenopusJamboreeParser xrefs succesfully parsed\n" if($verbose);

  return 0;
}

1;
