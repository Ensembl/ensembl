=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=cut

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
  my $dbi          = $ref_arg->{dbi};
  $dbi = $self->dbi unless defined $dbi;

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
                                 dbi        => $dbi,
				 source_id  => $source_id,
				 species_id => $species_id });
    $count++;
  }

  $file_io->close();

  print $count . " XenopusJamboreeParser xrefs succesfully parsed\n" if($verbose);

  return 0;
}

1;
