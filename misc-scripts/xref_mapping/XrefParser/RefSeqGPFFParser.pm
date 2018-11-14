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

package XrefParser::RefSeqGPFFParser;

use strict;
use warnings;
use Carp;
use List::MoreUtils qw(uniq);
use Bio::EnsEMBL::IO::Parser::Genbank;

# use Smart::Comments;


use parent qw( XrefParser::BaseParser );


sub run {

  my ($self, $ref_arg) = @_;
  my $source_id    = $ref_arg->{source_id};
  my $species_id   = $ref_arg->{species_id};
  my $species_name = $ref_arg->{species};
  my $files        = $ref_arg->{files};
  my $release_file = $ref_arg->{rel_file};
  my $dbi          = $ref_arg->{dbi} // $self->dbi;
  my $verbose      = $ref_arg->{verbose} // 0;

### $ref_arg

  if((!defined $source_id) or (!defined $species_id) or (!defined $files)){
    croak "Need to pass source_id, species_id, files and rel_file as pairs";
  }

  # set the list of RefSeq source names
  my $refseq_sources = {
    NM => 'RefSeq_mRNA',
    NR => 'RefSeq_ncRNA',
    XM => 'RefSeq_mRNA_predicted',
    XR => 'RefSeq_ncRNA_predicted',
    NP => 'RefSeq_peptide',
    XP => 'RefSeq_peptide_predicted',
  };

  # get RefSeq source ids
  my $source_ids = {};
  while (my ($source_prefix, $source_name) = each %{$refseq_sources}) {
    $source_ids->{$source_name} = $self->get_source_id_for_source_name( $source_name, undef, $dbi )
  }

  # get extra source ids
  $source_ids->{EntrezGene} = $self->get_source_id_for_source_name('EntrezGene', undef, $dbi);
  $source_ids->{WikiGene} = $self->get_source_id_for_source_name('WikiGene', undef, $dbi);

  # Retrieve existing RefSeq mRNA
  my $refseq_ids = { %{$self->get_valid_codes('RefSeq_mRNA', $species_id, $dbi)},
      %{$self->get_valid_codes('RefSeq_mRNA_predicted', $species_id, $dbi)} };
  my $entrez_ids = $self->get_valid_codes('EntrezGene', $species_id, $dbi);
  my $wiki_ids = $self->get_valid_codes('WikiGene', $species_id, $dbi);

  if ($verbose) {
    for my $source_name (sort values %{$refseq_sources}) {
      print "$source_name source ID = $source_ids->{$source_name}\n";
    }
  }

  # populate entrez gene id => label hash
  my $entrez = $self->get_acc_to_label("EntrezGene", $species_id, undef, $dbi);


  # get the species name
  my %id2name = $self->species_id2name($dbi);
  $species_name //= shift @{$id2name{$species_id}};


  # my %species2name = $self->species_id2name($dbi);


  # process the source files
  foreach my $file (@{$files}) {
    ### $file

    # type from the file (peptide or dna)
    my $type = $self->type_from_file($file);

    # get the file handler
    my $refseq_fh = $self->get_filehandle($file);

    if ( !defined $refseq_fh ) {
      warn "WARNING: Can't open RefSeqGPFF file $file\n";
      return;
    }

    # instantiate ensembl-io genbank parser
    my $parser = Bio::EnsEMBL::IO::Parser::Genbank->open($refseq_fh);

    # this will hold the array of xrefs to bulk insert
    my $xrefs;

    # For each record in the file
    while ( $parser->next ) {

      # Get the record species
      my $record_species = lc $parser->get_organism;
      $record_species =~ s/\s+/_/xg;

      # skip if species is not the required
      next unless ( $species_name eq $record_species);

      # get the acc
      my $acc = $parser->get_accession;

      my $taxon = $parser->get_taxon_id;
### $taxon

      # get description and remove the [species] at the end
      my $description = $parser->get_description;
      $description =~ s/\s*\[[^]]*\]\.?\z//x;

      # get refseq pair if available
      my $refseq_pair = $parser->get_dbsource_acc;

      # get the source_id for this acc type, warning and skip if not found
      my $prefix = substr($acc, 0, 2);
      if (!exists $refseq_sources->{$prefix}) {
        warn "WARNING: can't get source ID for $type $acc. Skipping\n";
        next;
      }
      my $acc_source_id = $source_ids->{$refseq_sources->{$prefix}};

      # the "pair" is the id from coded_by, or the dbsource_acc if that does not exist
      my $pair = pop @{$parser->get_coded_by_list};
      $pair //= $parser->get_dbsource_acc;

      # set up the direct xref
      my $xref = {
        ACCESSION     => $acc,
        LABEL         => $parser->get_sequence_name,
        VERSION       => $parser->get_sequence_version,
        DESCRIPTION   => $description,
        INFO_TYPE     => 'SEQUENCE_MATCH',
        PAIR          => $pair,
        SEQUENCE      => $parser->get_sequence,
        SEQUENCE_TYPE => $type,
        SOURCE_ID     => $acc_source_id,
        SPECIES_ID    => $species_id,
        PROTEIN       => pop @{$parser->get_protein_id_list},
      };

      # retrieve and remove duplicates from referenced entrez GeneID
      my $entrez_gene_ids = $parser->get_db_xref_list_for_type('GeneID');
      my @gene_ids = uniq( @{$entrez_gene_ids}); 

      # process existing entrez_gene_ids as dependent xrefs
      foreach my $gene_id (@gene_ids) {

        next unless (defined $entrez->{$gene_id});

        push @{$xref->{DEPENDENT_XREFS}}, {
            SOURCE_ID         => $source_ids->{'EntrezGene'},
            LINKAGE_SOURCE_ID => $acc_source_id,
            ACCESSION         => $gene_id,
            LABEL             => $entrez->{$gene_id}
        };

        push @{$xref->{DEPENDENT_XREFS}}, {
            SOURCE_ID         => $source_ids->{'WikiGene'},
            LINKAGE_SOURCE_ID => $acc_source_id,
            ACCESSION         => $gene_id,
            LABEL             => $entrez->{$gene_id}
        };

        next unless (defined $refseq_pair);

        # discard the version number
        $refseq_pair =~ s/\.\d+//x;

        # Add xrefs for RefSeq mRNA as well where available
        foreach my $refseq_id (@{ $refseq_ids->{$refseq_pair} }) {
          foreach my $entrez_id (@{ $entrez_ids->{$gene_id} }) {
            $self->add_dependent_xref_maponly(
                $entrez_id,
                $source_ids->{'EntrezGene'},
                $refseq_id,
                $source_ids->{$refseq_sources->{substr($refseq_id, 0, 2)}},
                $dbi
            );
          }
          foreach my $wiki_id (@{ $wiki_ids->{$gene_id} }) {
            $self->add_dependent_xref_maponly(
                $wiki_id,
                $source_ids->{'WikiGene'},
                $refseq_id,
                $source_ids->{$refseq_sources->{substr($refseq_id, 0, 2)}},
                $dbi
            );
          }
        }

      }

### $xref
      push @{$xrefs}, $xref;

    }

    $refseq_fh->close();

    # no xrefs in this record...
    if ( !defined( $xrefs ) ) {
      next;
    }

    # upload the xrefs
    $self->upload_xref_object_graphs( $xrefs, $dbi );

  }


  # process the release file
  if ( defined $release_file ) {
    # get filehandle
    my $release_io = $self->get_filehandle($release_file);

    # get file header
    my $release = do { local $/ = "\n*"; <$release_io> };
    $release =~ s/\s+/ /xg;

    if ( $release =~ m/(NCBI.*Release\s\d+)\s(.*)\sDistribution/x ) {
      my ($rel_number, $rel_date) = ($1, $2);
      my $release_string = "$rel_number, $rel_date";
 
      # set release info
      $self->set_release( $source_id, $release_string, $dbi );
      for my $source_name (sort values %{$refseq_sources}) {
        $self->set_release( $source_ids->{$source_name}, $release_string, $dbi );
      }
      if ($verbose) {
        print "RefSeq release: '$release_string'\n";
      }
    } else {
      warn "WARNING: Could not set release info from release file '$release_file'\n";
    }

  } else {
    warn "WARNING: No release_file available\n";
  }

  return 0;
}



# get type from filename path. this includes the source name and that's enough to extract it
sub type_from_file {
  my ($self, $file) = @_;

  my ($type) = $file =~ /RefSeq_(\w+)\//x;

  warn "WARNING: Could not work out sequence type for $file\n" unless $type;

  return $type;
}

1;
