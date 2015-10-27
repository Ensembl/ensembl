=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

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

package XrefParser::HGNCParser;

use strict;
use warnings;
use File::Basename;
use Carp;
use base qw( XrefParser::BaseParser );

sub run {

  my ($self, $ref_arg) = @_;
  my $source_id    = $ref_arg->{source_id};
  my $species_id   = $ref_arg->{species_id};
  my $files        = $ref_arg->{files};
  my $verbose      = $ref_arg->{verbose};

  if((!defined $source_id) or (!defined $species_id) or (!defined $files) ){
    croak "Need to pass source_id, species_id, files and rel_file as pairs";
  }
  $verbose |=0;

  my $empty = q{};

  my $file = @{$files}[0];

  my (%swissprot)  =  %{$self->get_valid_codes('Uniprot/SWISSPROT',$species_id)};
  my (%refseq) =  %{$self->get_valid_codes('refseq',$species_id)};
  my @list;
  push @list, 'refseq_peptide';
  push @list, 'refseq_mRNA';
  my (%entrezgene) = %{$self->get_valid_xrefs_for_dependencies('EntrezGene',@list)};


  my %name_count;
  my $mismatch = 0;

  my $hugo_io = $self->get_filehandle($file);

  if ( !defined $hugo_io ) {
    print "ERROR: Can't open HGNC file $file\n";
    return 1;
  }

  my $name_to_source_id = $self->get_hgnc_sources();

  # Skip header
  $hugo_io->getline();

  while ( $_ = $hugo_io->getline() ) {
    chomp;
    my @array = split /\t/x, $_;

    my $seen = 0;

    my $acc              = $array[0];
    my $symbol           = $array[1];
    my $name             = $array[2];
    my $previous_symbols = $array[3];
    my $synonyms         = $array[4];


    my $type = 'lrg';
    my $id = $array[11];
    my $source_id = $name_to_source_id->{$type};
    if($id and $id =~ m/http:\/\/www.lrg-sequence.org\/LRG\/(LRG_\d+)/x){
      my $lrg_stable_id = $1;
      $self->add_to_direct_xrefs({ stable_id   => $lrg_stable_id,
				   type        => 'gene',
                                   acc         => $acc,
				   label       => $symbol,
				   desc        => $name,
				   source_id   => $source_id,
				   species_id  => $species_id} );

      $self->add_synonyms_for_hgnc( {source_id  => $source_id,
				     name       => $acc,
				     species_id => $species_id,
				     dead       => $previous_symbols,
				     alias      => $synonyms} );
      $name_count{$type}++;
    }


    #
    # Direct Ensembl mappings
    #
    $type = 'ensembl_manual';
    $id = $array[9];
    $source_id = $name_to_source_id->{$type};
    if ($id){              # Ensembl direct xref
      $seen = 1;
      $name_count{$type}++;
      $self->add_to_direct_xrefs({ stable_id  => $id,
				   type       => 'gene',
				   acc        => $acc,
				   label      => $symbol,
				   desc       => $name,,
				   source_id  => $source_id,
				   species_id => $species_id} );

      $self->add_synonyms_for_hgnc( {source_id  => $source_id,
				     name       => $acc,
				     species_id => $species_id,
				     dead       => $previous_symbols,
				     alias      => $synonyms});

    }

    #
    # RefSeq
    #
    $type = 'refseq_mapped';
    $id = $array[8];
    $source_id = $name_to_source_id->{$type};
    if ($id) {
      if(defined $refseq{$id} ){
        $seen = 1;
	foreach my $xref_id (@{$refseq{$id}}){
	  $name_count{$type}++;
	  $self->add_dependent_xref({ master_xref_id => $xref_id,
				      acc            => $acc,
				      label          => $symbol,
				      desc           => $name || '',
				      source_id      => $source_id,
				      species_id     => $species_id} );
	}
	$self->add_synonyms_for_hgnc( {source_id  => $source_id,
				       name       => $acc,
				       species_id => $species_id,
				       dead       => $previous_symbols,
				       alias      => $synonyms});
      }
    }

    $type = 'refseq_manual';
    $id = $array[6];
    $source_id = $name_to_source_id->{$type};
    if ($id) {
      if(defined $refseq{$id} ){
        $seen = 1;
        foreach my $xref_id (@{$refseq{$id}}){
          $name_count{$type}++;
          $self->add_dependent_xref({ master_xref_id => $xref_id,
                                      acc            => $acc,
                                      label          => $symbol,
                                      desc           => $name || '',
                                      source_id      => $source_id,
                                      species_id     => $species_id} );
        }
        $self->add_synonyms_for_hgnc( {source_id  => $source_id,
                                       name       => $acc,
                                       species_id => $species_id,
                                       dead       => $previous_symbols,
                                       alias      => $synonyms});
      }
    }

    #
    # Swissprot
    # Skip, Uniprot is protein-centric
    #
    $type = 'swissprot_manual';

    #
    # EntrezGene
    #
    $type = 'entrezgene_manual';
    $id = $array[5];
    $source_id = $name_to_source_id->{$type};
    if(defined $id ){
      if(defined $entrezgene{$id} ){
        $seen = 1;
        $self->add_dependent_xref({ master_xref_id => $entrezgene{$id},
           			    acc            => $acc,
				    label          => $symbol,
				    desc           => $name || '',
				    source_id      => $source_id,
				    species_id     => $species_id} );
        $name_count{$type}++;
        $self->add_synonyms_for_hgnc( {source_id  => $source_id,
				       name       => $acc,
				       species_id => $species_id,
				       dead       => $previous_symbols,
				       alias      => $synonyms});
      }
    }

    $type = 'entrezgene_mapped';
    $id = $array[7];
    if(defined $id ){
      if(defined $entrezgene{$id} ){
        $seen = 1;
        $self->add_dependent_xref({ master_xref_id => $entrezgene{$id},
                                    acc            => $acc,
                                    label          => $symbol,
                                    desc           => $name || '',
                                    source_id      => $source_id,
                                    species_id     => $species_id} );
        $name_count{$type}++;
        $self->add_synonyms_for_hgnc( {source_id  => $source_id,
                                       name       => $acc,
                                       species_id => $species_id,
                                       dead       => $previous_symbols,
                                       alias      => $synonyms});
      }
    }



    if(!$seen){ # Store to keep descriptions etc
      $type = 'desc_only';
      $source_id = $name_to_source_id->{$type};
      $self->add_xref({ acc        => $acc,
			label      => $symbol,
			desc       => $name,
			source_id  => $source_id,
			species_id => $species_id,
			info_type  => "MISC"} );

      $self->add_synonyms_for_hgnc( {source_id  => $source_id,
				     name       => $acc,
				     species_id => $species_id,
				     dead       => $previous_symbols,
				     alias      => $synonyms});
      $mismatch++;
    }
  }


  $hugo_io->close();

  if($verbose){
    print 'Loaded a total of :-';
    foreach my $type (keys %name_count){
      print "\t".$type."\t".$name_count{$type}."\n";
    }
    print "$mismatch xrefs could not be associated via RefSeq, EntrezGene or ensembl\n";
  }
  return 0; # successful
}



sub get_hgnc_sources {
  my $self = shift;
  my %name_to_source_id;

  my @sources = ('entrezgene_manual', 'refseq_manual', 'entrezgene_mapped', 'refseq_mapped', 'ensembl_manual', 'swissprot_manual', 'desc_only');


  foreach my $key (@sources) {
  my $source_id = $self->get_source_id_for_source_name('HGNC', $key);
    if(!(defined $source_id)){
      die 'Could not get source id for HGNC and '. $key ."\n";
    }
    $name_to_source_id{ $key } = $source_id;
  }

  my $source_id = $self->get_source_id_for_source_name('LRG_HGNC_notransfer');
  if(!(defined $source_id) ){
    die 'Could not get source id for LRG_HGNC_notransfer\n';
  }
  $name_to_source_id{'lrg'} = $source_id;

  return \%name_to_source_id;
}

sub add_synonyms_for_hgnc{
  my ($self, $ref_arg) = @_;

  my $source_id  = $ref_arg->{source_id};
  my $name       = $ref_arg->{name};
  my $species_id = $ref_arg->{species_id};
  my $dead_name  = $ref_arg->{dead};
  my $alias      = $ref_arg->{alias};

  if (defined $dead_name ) {     # dead name, add to synonym
    my @array2 = split ',\s*', $dead_name ;
    foreach my $arr (@array2){
      $self->add_to_syn($name, $source_id, $arr, $species_id);
    }
  }

  if (defined $alias ) {     # alias, add to synonym
    my @array2 = split ',\s*', $alias;
    foreach my $arr (@array2){
      $self->add_to_syn($name, $source_id, $arr, $species_id);
    }
  }
  return;
}

1;


