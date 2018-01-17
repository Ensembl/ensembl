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

package XrefParser::ReactomeParser;

use strict;
use warnings;
use Carp;
use DBI;

use base qw( XrefParser::BaseParser );
use XrefParser::Database;

# Parse file of Reactome records and assign direct xrefs


# --------------------------------------------------------------------------------

sub run {

 my ($self, $ref_arg) = @_;
  my $source_id    = $ref_arg->{source_id};
  my $species_id   = $ref_arg->{species_id};
  my $files        = $ref_arg->{files};
  my $release_file = $ref_arg->{rel_file};
  my $verbose      = $ref_arg->{verbose};
  my $dbi          = $ref_arg->{dbi};
  my $species_name = $ref_arg->{species};
  $dbi = $self->dbi unless defined $dbi;

  if((!defined $source_id) or (!defined $species_id) or (!defined $files) ){
    croak "Needs to pass source_id, species_id and  files as pairs";
  }
  $verbose |=0;

  my $file_desc = @{$files}[1];

  if ( defined $release_file ) {
    my $release;
    # Parse and set release information from $release_file.
    my $release_io = $self->get_filehandle($release_file);
    while ( defined( my $line = $release_io->getline() ) ) {
      if ( $line =~ /([0-9]*)/ ) {
        $release = $1;
        print "Reactome release is '$release'\n" if($verbose);
        last;
      }
    }

    if (!$release) {
      croak "Could not find release using $release_file\n";
    }

    $self->set_release( $source_id, $release, $dbi );
  }

  # Create a hash of all valid names for this species
  my %species2alias = $self->species_id2name($dbi);
  if (defined $species_name) { push @{$species2alias{$species_id}}, $species_name; }
  if (!defined $species2alias{$species_id}) { next; }
  my @aliases = @{$species2alias{$species_id}};
  my %alias2species_id = map {$_, 1} @aliases;
 
  my $parsed_count = 0;
  my $err_count = 0;

  my %reactome2ensembl;

  my $reactome_source_id =  $self->get_source_id_for_source_name("reactome", "direct", $dbi);
  my $transcript_reactome_source_id =  $self->get_source_id_for_source_name("reactome_transcript", undef, $dbi);
  my $gene_reactome_source_id =  $self->get_source_id_for_source_name("reactome_gene", undef, $dbi);
  my $reactome_uniprot_source_id = $self->get_source_id_for_source_name("reactome", "uniprot", $dbi);
  if($reactome_source_id < 1 || $transcript_reactome_source_id < 1 || $gene_reactome_source_id < 1){
    die "Could not find source id for reactome sources???\n";
  }
  else{
    print "Source_id = $reactome_source_id\n";
    print "Transcript_source_id = $transcript_reactome_source_id\n";
    print "Gene_source_id = $gene_reactome_source_id\n";
  }

  if($reactome_uniprot_source_id < 1){
    die "Could not find source id for reactome uniprot???\n";
  }
  else{
    print "Source_id = $reactome_uniprot_source_id\n";
  }

  my (%uniprot) = %{$self->get_valid_codes("uniprot/",$species_id, $dbi)};
  my $is_uniprot = 0;

  foreach my $file (@$files) {
    my $reactome_io = $self->get_filehandle($file);
    if ($file =~ /UniProt/) { $is_uniprot = 1; }
  # Example line:
  # ENSG00000138685 REACT_111045    http://www.reactome.org/PathwayBrowser/#REACT_111045    Developmental Biology   TAS     Homo sapiens
    while (my $line = $reactome_io->getline() ) {
      chomp $line;
      my ($ensembl_stable_id, $reactome_id, $url, $description, $evidence, $species) = split /\t+/,$line;
  
      $species =~ s/\s//;
      $species = lc($species);
      if ( $alias2species_id{$species} ){
        $parsed_count++;
  
  # Attempt to guess the object_type based on the stable id
  # Some entries just don't match on stable id, so warn but do not die
  # For example:
  # 00000074047   REACT_268323    http://www.reactome.org/PathwayBrowser/#REACT_268323    Hedgehog 'off' state    TAS     Homo sapiens
        my $type;
        my $info_type = 'DIRECT';
        if ($is_uniprot) {
          if (defined($uniprot{$ensembl_stable_id})) {
    # First check if it is a uniprot id
            foreach my $xref_id (@{$uniprot{$ensembl_stable_id}}){
              $self->add_dependent_xref({ master_xref_id => $xref_id,
                                  acc            => $reactome_id,
                                  label          => $reactome_id,
                                  desc           => $description,
                                  source_id      => $reactome_uniprot_source_id,
                                  dbi            => $dbi,
                                  species_id     => $species_id} );
            }
            $info_type = 'DEPENDENT';
          }
        }
        elsif ($ensembl_stable_id =~ /G[0-9]*$/) { 
          $type = 'gene';
          $reactome_source_id = $gene_reactome_source_id;
        }
        elsif ($ensembl_stable_id =~ /T[0-9]*$/) { 
          $type = 'transcript';
          $reactome_source_id = $transcript_reactome_source_id;
        }
        elsif ($ensembl_stable_id =~ /P[0-9]*$/) { $type = 'translation'; }
        else {
  # Is not in Uniprot and does not match Ensembl stable id format
          print STDERR "Could not find type for $ensembl_stable_id\n";
          $err_count++;
          next;
        }
  
  # Add new entry for reactome xref
  # as well as direct xref to ensembl stable id
        my $xref_id = $self->add_xref({ acc         => $reactome_id,
                          label       => $reactome_id,
                          desc        => $description,
                          info_type   => $info_type,
                          source_id   => $reactome_source_id,
                          dbi         => $dbi,
                          species_id  => $species_id} );
  
        $self->add_direct_xref($xref_id, $ensembl_stable_id, $type, $dbi) if $type;
      }
    }
  }

  print "$parsed_count entries processed\n$err_count not found\n";
  return 0;
}


1;
