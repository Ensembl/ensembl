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

package XrefParser::VGNCParser;

use strict;
use warnings;
use File::Basename;
use Carp;
use base qw( XrefParser::HGNCParser);

sub run {

  my ($self, $ref_arg) = @_;
  my $source_id    = $ref_arg->{source_id};
  my $species_id   = $ref_arg->{species_id};
  my $files        = $ref_arg->{files};
  my $verbose      = $ref_arg->{verbose};
  my $dbi          = $ref_arg->{dbi};
  $dbi = $self->dbi unless defined $dbi;

  if((!defined $source_id) or (!defined $species_id) or (!defined $files) ){
    croak "Need to pass source_id, species_id, files and rel_file as pairs";
  }
  $verbose |=0;

  my $file = @{$files}[0];

  my $count = 0;

  my $file_io = $self->get_filehandle($file);

  if ( !defined $file_io ) {
    print "ERROR: Can't open VGNC file $file\n";
    return 1;
  }

  my $source_name = $self->get_source_name_for_source_id($source_id, $dbi);
  # Create a hash of all valid taxon_ids for this species
  my %species2tax = $self->species_id2taxonomy($dbi);
  my @tax_ids = @{$species2tax{$species_id}};
  my %taxonomy2species_id = map{ $_=>$species_id } @tax_ids;

  # Skip header
  $file_io->getline();

  while ( $_ = $file_io->getline() ) {
    chomp;
    my @array = split /\t/x, $_;

    my $taxon_id         = $array[0];
    my $acc              = $array[1];
    my $symbol           = $array[2];
    my $name             = $array[3];
    my $id               = $array[20];
    my $previous_symbols = $array[9];
    my $synonyms         = $array[11];

    $previous_symbols =~ s/"//g;
    $synonyms =~ s/"//g;

    unless (exists ($taxonomy2species_id{$taxon_id})) { next; }

    if ($id){              # Ensembl direct xref
      $self->add_to_direct_xrefs({ stable_id  => $id,
				   type       => 'gene',
				   acc        => $acc,
				   label      => $symbol,
				   desc       => $name,
                                   dbi        => $dbi,
				   source_id  => $source_id,
				   species_id => $species_id} );

      $self->add_synonyms_for_hgnc( {source_id  => $source_id,
                                     name       => $acc,
                                     species_id => $species_id,
                                     dbi        => $dbi,
                                     dead       => $previous_symbols,
                                     alias      => $synonyms});

      $count++;
    }
  }


  $file_io->close();

  if($verbose){
    print "Loaded a total of $count xrefs\n";
  }
  return 0; # successful
}


1;


