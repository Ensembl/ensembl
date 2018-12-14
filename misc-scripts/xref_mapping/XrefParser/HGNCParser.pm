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

package XrefParser::HGNCParser;

use strict;
use warnings;
use Carp;
use Text::CSV;
use utf8;

use parent qw( XrefParser::BaseParser );


# HGNC sources to be processed
my @SOURCES = (
  'ccds',
  'entrezgene_manual',
  'refseq_manual',
  'ensembl_manual',
  'desc_only'
);



sub run_script {
  my ($self, $ref_arg) = @_;

  my $source_id    = $ref_arg->{source_id};
  my $species_id   = $ref_arg->{species_id};
  my $file         = $ref_arg->{file};
  my $db           = $ref_arg->{dba};
  my $verbose      = $ref_arg->{verbose} // 0;
  my $dbi          = $ref_arg->{dbi} // $self->dbi;

  if ((!defined $source_id) or (!defined $species_id) or (!defined $file) ){
    confess "Need to pass source_id, species_id, file as pairs";
  }

  # parse the file string and set default user
  my $file_params = $self->parse_file_string($file);
  $file_params->{user} //= 'ensro';

  # Prepare lookup lists
  my (%swissprot) = %{$self->get_valid_codes('Uniprot/SWISSPROT',$species_id, $dbi)};
  my (%refseq) = %{$self->get_valid_codes('refseq',$species_id, $dbi)};
  my @list = ('refseq_peptide', 'refseq_mRNA');
  my (%entrezgene) = %{$self->get_valid_xrefs_for_dependencies('EntrezGene', $dbi, @list)};

  # Prepare sources
  my $self_source_name = $self->get_source_name_for_source_id($source_id, $dbi);

  # get RefSeq source ids
  foreach my $source_name (@SOURCES) {
    $self->{source_ids}->{$source_name} = $self->get_source_id_for_source_name( $self_source_name, $source_name , $dbi );
  }
  $self->{source_ids}->{'lrg'} = $self->get_source_id_for_source_name( 'LRG_HGNC_notransfer', undef, $dbi );

  # statistics counts
  my %name_count;
  my $mismatch = 0;

  # Get CCDS data from core db
  my $core_db;
  if (defined $db) {
    $core_db = $db->dbc();
  } elsif (defined $file_params->{host}) {
    $core_db = XrefParser::Database->new({
        host   => $file_params->{host},
        port   => $file_params->{port},
        user   => $file_params->{user},
        dbname => $file_params->{dbname},
        pass   => $file_params->{pass}
    })->dbi;
  } else {
    confess "No ensembl core database provided\n";
  }

  if (!defined $core_db) {
    confess "No ensembl core database!\n";
  }

  my $sql =(<<'CCDS');
  SELECT ta.value, t.stable_id
  FROM transcript t
  INNER JOIN transcript_attrib ta ON t.transcript_id = ta.transcript_id
  INNER JOIN attrib_type a ON ta.attrib_type_id = a.attrib_type_id
  WHERE a.code = 'ccds_transcript';
CCDS

  my %ccds_to_ens;
  my $sth = $core_db->prepare($sql);
  $sth->execute() or croak( $core_db->errstr() );
  while ( my ($ccds_id, $ens_id) = $sth->fetchrow_array() ) {
    # Remove version
    $ccds_id =~ s/\.\d+//x;
    $ccds_to_ens{$ccds_id} = $ens_id;
  }
  $sth->finish;

  # in memory HGNC file
  my $mem_file;

  # use wget link to get file
  if (defined $file_params->{wget}) {
    my $ua = LWP::UserAgent->new();
    $ua->timeout(10);
    $ua->env_proxy();
    my $request = HTTP::Request->new(
      GET => $file_params->{wget}
    );
    my $response = $ua->request($request);

    if ( !$response->is_success() ) {
      confess $response->status_line;
    }

    $mem_file = $response->decoded_content;

  # else get file from disk
  } else {
    my $disk_fh = $self->get_filehandle($file);
    if ( !defined $disk_fh ) {
      confess "Can't open HGNC file '$file'\n";
    }
    $mem_file = do { local $/; <$disk_fh> };
  }

  my $input_file = Text::CSV->new({
    sep_char       => "\t",
    empty_is_undef => 1,
    binary         => 1,
    auto_diag      => 1
  }) or croak "Cannot use file $file: ".Text::CSV->error_diag ();


## COLUMNS OF FILE ARE AS FOLLOW:

# HGNC ID
# Approved symbol
# Approved name
# Previous symbols
# Synonyms
# NCBI Gene ID
# Ensembl gene ID
# RefSeq IDs
# CCDS IDs
# Locus specific databases


  # make sure it's utf8
  utf8::encode($mem_file);
  # get rid of non-conventional " used in the Locus specific databases field
  $mem_file =~ s/"//xg;

  open my $fh, '<', \$mem_file or confess "Can't open HGNC in-memory file: $!\n";

  $input_file->column_names( @{ $input_file->getline( $fh ) } );


  # loop through each row
  while ( my $data = $input_file->getline_hr( $fh ) ) {

    my $acc              = $data->{'HGNC ID'};
    my $symbol           = $data->{'Approved symbol'};
    my $name             = $data->{'Approved name'};
    my $previous_symbols = $data->{'Previous symbols'};
    my $synonyms         = $data->{'Synonyms'};

    my $seen = 0;

    # Direct CCDS to ENST mappings
    my $ccds = $data->{'CCDS IDs'};
    my @ccds_list;

    if ( defined $ccds ) {
      @ccds_list = split( /,\s/x, $ccds );
    }

    CCDS:
    foreach my $ccds (@ccds_list) {
      my $enst_id = $ccds_to_ens{$ccds};

      if (!defined $enst_id) {
        next CCDS;
      }

      $self->add_to_direct_xrefs({
          stable_id  => $enst_id,
          type       => 'gene',
          acc        => $acc,
          label      => $symbol,
          desc       => $name,
          source_id  => $self->{source_ids}->{'ccds'},
          dbi        => $dbi,
          species_id => $species_id
      });

      $self->add_synonyms_for_hgnc({
          source_id  => $self->{source_ids}->{'ccds'},
          name       => $acc,
          species_id => $species_id,
          dbi        => $dbi,
          dead       => $previous_symbols,
          alias      => $synonyms
      });
      $name_count{'ccds'}++;
    }

    # Direct LRG to ENST mappings
    my $lrg_id = $data->{'Locus specific databases'};
    if ( defined $lrg_id && $lrg_id =~ m/(LRG_\d+)|/x ){
      $lrg_id = $1;
    }

    if ( defined $lrg_id ){
      $self->add_to_direct_xrefs({
          stable_id   => $lrg_id,
          type        => 'gene',
          acc         => $acc,
          label       => $symbol,
          desc        => $name,
          source_id   => $self->{source_ids}->{'lrg'},
          dbi         => $dbi,
          species_id  => $species_id
      });

      $self->add_synonyms_for_hgnc({
          source_id  => $self->{source_ids}->{'lrg'},
          name       => $acc,
          species_id => $species_id,
          dbi        => $dbi,
          dead       => $previous_symbols,
          alias      => $synonyms
      });
      $name_count{'lrg'}++;
    }

    # Direct Ensembl mappings
    my $ensg_id = $data->{'Ensembl gene ID'};
    if ( defined $ensg_id ){
      $seen = 1;

      $self->add_to_direct_xrefs({
          stable_id  => $ensg_id,
          type       => 'gene',
          acc        => $acc,
          label      => $symbol,
          desc       => $name,
          dbi        => $dbi,
          source_id  => $self->{source_ids}->{'ensembl_manual'},
          species_id => $species_id
      });

      $self->add_synonyms_for_hgnc({
          source_id  => $self->{source_ids}->{'ensembl_manual'},
          name       => $acc,
          species_id => $species_id,
          dead       => $previous_symbols,
          dbi        => $dbi,
          alias      => $synonyms
      });
      $name_count{'ensembl_manual'}++;
    }

    # RefSeq
    my $refseq_id = $data->{'RefSeq IDs'};
    if ($refseq_id) {
      if ( defined $refseq{$refseq_id} ){
        $seen = 1;
        foreach my $xref_id ( @{$refseq{$refseq_id}} ){
          $self->add_dependent_xref({
              master_xref_id => $xref_id,
              acc            => $acc,
              label          => $symbol,
              desc           => $name,
              source_id      => $self->{source_ids}->{'refseq_manual'},
              dbi            => $dbi,
              species_id     => $species_id
          });
          $name_count{'refseq_manual'}++;
        }

        $self->add_synonyms_for_hgnc({
            source_id  => $self->{source_ids}->{'refseq_manual'},
            name       => $acc,
            species_id => $species_id,
            dbi        => $dbi,
            dead       => $previous_symbols,
            alias      => $synonyms
        });
      }
    }

    # EntrezGene
    my $entrez_id = $data->{'NCBI Gene ID'};
    if ( defined $entrez_id ){
      if ( defined $entrezgene{$entrez_id} ){
        $seen = 1;
        $self->add_dependent_xref({
            master_xref_id => $entrezgene{$entrez_id},
            acc            => $acc,
            label          => $symbol,
            desc           => $name,
            source_id      => $self->{source_ids}->{'entrezgene_manual'},
            dbi            => $dbi,
            species_id     => $species_id
        });

        $self->add_synonyms_for_hgnc({
            source_id  => $self->{source_ids}->{'entrezgene_manual'},
            name       => $acc,
            species_id => $species_id,
            dead       => $previous_symbols,
            dbi        => $dbi,
            alias      => $synonyms
        });
        $name_count{'entrezgene_manual'}++;
      }
    }

    # Store to keep descriptions if stored yet
    if ( !$seen ){
      $self->add_xref({
          acc        => $acc,
          label      => $symbol,
          desc       => $name,
          source_id  => $self->{source_ids}->{'desc_only'},
          species_id => $species_id,
          dbi        => $dbi,
          info_type  => "MISC"
      });

      $self->add_synonyms_for_hgnc({
          source_id  => $self->{source_ids}->{'desc_only'},
          name       => $acc,
          species_id => $species_id,
          dbi        => $dbi,
          dead       => $previous_symbols,
          alias      => $synonyms
      });
      $mismatch++;
    }

  }

  close $fh;

  if ( $verbose ){
    print "HGNC xrefs loaded:\n";
    foreach my $type (sort keys %name_count){
      print "\t$type\t$name_count{$type}\n";
    }
    print "$mismatch HGNC ids could not be associated in xrefs\n";
  }
  return 0; # successful
}



sub add_synonyms_for_hgnc {
  my ($self, $ref_arg) = @_;

  my $source_id    = $ref_arg->{source_id};
  my $name         = $ref_arg->{name};
  my $species_id   = $ref_arg->{species_id};
  my $dead_string  = $ref_arg->{dead};
  my $alias_string = $ref_arg->{alias};
  my $dbi          = $ref_arg->{dbi};

  # dead name, add to synonym
  if (defined $dead_string) {
    $dead_string =~ s/"//xg;
    my @dead_array = split( ',\s', $dead_string );
    foreach my $dead (@dead_array){
      $self->add_to_syn($name, $source_id, $dead, $species_id, $dbi);
    }
  }

  # alias name, add to synonym
  if (defined $alias_string) {
    $alias_string =~ s/"//xg;
    my @alias_array = split( ',\s', $alias_string );
    foreach my $alias (@alias_array){
      $self->add_to_syn($name, $source_id, $alias, $species_id, $dbi);
    }
  }

  return;
}



# parses the input string $file into an hash
# string $file is in the format as the example:
# script:project=>ensembl,host=>ens-staging1,dbname=>homo_sapiens_core_70_37,ofhost=>ens-staging1,...
# string until : is ignored, hash is built with keys=>values provided
sub parse_file_string {
  my ($self, $file_string) = @_;

  $file_string =~ s/\A\w+://x;

  my @param_pairs = split( /,/x, $file_string );

  my $params;

  # set provided values
  foreach my $pair ( @param_pairs ) {
    my ($key, $value) = split( /=>/x, $pair );
    $params->{$key} = $value;
  }

  return $params;
}



1;
