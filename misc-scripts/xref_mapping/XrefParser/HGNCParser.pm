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
use File::Basename;
use Carp;
use base qw( XrefParser::BaseParser );

sub run_script {

  my ($self, $ref_arg) = @_;
  my $source_id    = $ref_arg->{source_id};
  my $species_id   = $ref_arg->{species_id};
  my $file         = $ref_arg->{file};
  my $verbose      = $ref_arg->{verbose};
  my $db           = $ref_arg->{dba};
  my $dbi          = $ref_arg->{dbi};
  $dbi = $self->dbi unless defined $dbi;

  if((!defined $source_id) or (!defined $species_id) or (!defined $file) ){
    croak "Need to pass source_id, species_id, file as pairs";
  }
  $verbose |=0;

  my $user = "ensro";
  my $host;
  my $port;
  my $dbname;
  my $pass;
  my $wget;

  if($file =~ /host[=][>](\S+?)[,]/){
    $host = $1;
  }
  if($file =~ /port[=][>](\S+?)[,]/){
    $port =  $1;
  }
  if($file =~ /dbname[=][>](\S+?)[,]/){
    $dbname = $1;
  }
  if($file =~ /pass[=][>](\S+?)[,]/){
    $pass = $1;
  }
  if($file =~ /wget[=][>](\S+?)[,]/){
    $wget = $1;
  }

  my @lines;
  if (defined $wget) {
    my $ua = LWP::UserAgent->new();
    $ua->timeout(10);
    $ua->env_proxy();
    my $request = HTTP::Request->new(GET => $wget);
    my $response = $ua->request($request);

    if ( !$response->is_success() ) {
      warn($response->status_line);
      return 1;
    }
    @lines = split(/\n/, $response->decoded_content);
  } else {
    my $file_io = $self->get_filehandle($file);
    if ( !defined $file_io ) {
      print "ERROR: Can't open HGNC file $file\n";
      return 1;
    }
    while (my $line = $file_io->getline()) {
      push(@lines, $line);
    }
  }


  my (%swissprot)  =  %{$self->get_valid_codes('Uniprot/SWISSPROT',$species_id, $dbi)};
  my (%refseq) =  %{$self->get_valid_codes('refseq',$species_id, $dbi)};
  my @list;
  push @list, 'refseq_peptide';
  push @list, 'refseq_mRNA';
  my (%entrezgene) = %{$self->get_valid_xrefs_for_dependencies('EntrezGene', $dbi, @list)};
  my $source_name = $self->get_source_name_for_source_id($source_id, $dbi);
  my $name_to_source_id = $self->get_sources($source_name, $dbi);


  my %name_count;
  my $mismatch = 0;

  # Get CCDS data
  my ($ccds_db, $dbi2);
  if (defined $host) {
    $ccds_db =  XrefParser::Database->new({ host   => $host,
                                             port   => $port,
                                             user   => $user,
                                             dbname => $dbname,
                                             pass   => $pass});
    $dbi2 = $ccds_db->dbi();
  } else {
    $dbi2 = $db->dbc();
  }

  if(!defined($dbi2)){
    return 1;
  }

  my $sql = "select t.stable_id, x.dbprimary_acc from transcript t, object_xref ox, xref x, external_db e where t.transcript_id = ox.ensembl_id and ox.ensembl_object_type = 'Transcript' and ox.xref_id = x.xref_id and x.external_db_id = e.external_db_id and e.db_name like 'Ens_%_transcript'";
  my (%ccds_to_ens, $ccds);
  my $sth = $dbi2->prepare($sql);
  $sth->execute() or croak( $dbi2->errstr() );
  while ( my @row = $sth->fetchrow_array() ) {
    $ccds = $row[0];
    $ccds =~ s/\.[0-9]//; # Remove version
    $ccds_to_ens{$ccds} = $row[1];
  }
  $sth->finish;

  # Skip header
  pop @lines;

  foreach my $line (@lines){
    chomp $line;
    my @array = split /\t/x, $line;

    my $seen = 0;

    my $acc              = $array[0];
    my $symbol           = $array[1];
    my $name             = $array[2];
    my $status           = $array[5];
    my $previous_symbols = $array[8];
    my $synonyms         = $array[10];

    if ($status ne 'Approved') { next; }

    my $type = 'lrg';
    my $id = $array[29];
    my $source_id = $name_to_source_id->{$type};
    if($id and $id =~ m/http:\/\/www.lrg-sequence.org\/LRG\/(LRG_\d+)/x){
      my $lrg_stable_id = $1;
      $self->add_to_direct_xrefs({ stable_id   => $lrg_stable_id,
				   type        => 'gene',
                                   acc         => $acc,
				   label       => $symbol,
				   desc        => $name,
				   source_id   => $source_id,
                                   dbi         => $dbi,
				   species_id  => $species_id} );

      $self->add_synonyms_for_hgnc( {source_id  => $source_id,
				     name       => $acc,
				     species_id => $species_id,
                                     dbi        => $dbi,
				     dead       => $previous_symbols,
				     alias      => $synonyms} );
      $name_count{$type}++;
    }


    # 
    # Direct CCDS to ENST mappings
    #
    $type = 'ccds';
    $source_id = $name_to_source_id->{$type};

    my $ccds = $array[24];
    $ccds =~ s/"//g;
    my @ccds_list = split(/\|/,$ccds);

    foreach my $ccds (@ccds_list) {
      $id = $ccds_to_ens{$ccds};
      if (!defined $id) { next; }
      $self->add_to_direct_xrefs({ stable_id   => $id,
                                   type        => 'gene',
                                   acc         => $acc,
                                   label       => $symbol,
                                   desc        => $name,
                                   source_id   => $source_id,
                                   dbi         => $dbi,
                                   species_id  => $species_id} );

      $self->add_synonyms_for_hgnc( {source_id  => $source_id,
                                     name       => $acc,
                                     species_id => $species_id,
                                     dbi        => $dbi,
                                     dead       => $previous_symbols,
                                     alias      => $synonyms} );
      $name_count{$type}++;
    }


    #
    # Direct Ensembl mappings
    #
    $type = 'ensembl_manual';
    $id = $array[19];
    $source_id = $name_to_source_id->{$type};
    if ($id){              # Ensembl direct xref
      $seen = 1;
      $name_count{$type}++;
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
				     dead       => $previous_symbols,
                                     dbi        => $dbi,
				     alias      => $synonyms});

    }

    $type = 'refseq_manual';
    $id = $array[23];
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
                                      dbi            => $dbi,
                                      species_id     => $species_id} );
        }
        $self->add_synonyms_for_hgnc( {source_id  => $source_id,
                                       name       => $acc,
                                       species_id => $species_id,
                                       dbi        => $dbi,
                                       dead       => $previous_symbols,
                                       alias      => $synonyms});
      }
    }

    #
    # EntrezGene
    #
    $type = 'entrezgene_manual';
    $id = $array[18];
    $source_id = $name_to_source_id->{$type};
    if(defined $id ){
      if(defined $entrezgene{$id} ){
        $seen = 1;
        $self->add_dependent_xref({ master_xref_id => $entrezgene{$id},
           			    acc            => $acc,
				    label          => $symbol,
				    desc           => $name || '',
				    source_id      => $source_id,
                                    dbi            => $dbi,
				    species_id     => $species_id} );
        $name_count{$type}++;
        $self->add_synonyms_for_hgnc( {source_id  => $source_id,
				       name       => $acc,
				       species_id => $species_id,
				       dead       => $previous_symbols,
                                       dbi        => $dbi,
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
                        dbi        => $dbi,
			info_type  => "MISC"} );

      $self->add_synonyms_for_hgnc( {source_id  => $source_id,
				     name       => $acc,
				     species_id => $species_id,
                                     dbi        => $dbi,
				     dead       => $previous_symbols,
				     alias      => $synonyms});
      $mismatch++;
    }
  }


  if($verbose){
    print 'Loaded a total of :-';
    foreach my $type (keys %name_count){
      print "\t".$type."\t".$name_count{$type}."\n";
    }
    print "$mismatch xrefs could not be associated via RefSeq, EntrezGene or ensembl\n";
  }
  return 0; # successful
}



sub get_sources {
  my $self = shift;
  my $source_name = shift;
  my $dbi = shift;
  my %name_to_source_id;

  my @sources = ('entrezgene_manual', 'refseq_manual', 'entrezgene_mapped', 'refseq_mapped', 'ensembl_manual', 'swissprot_manual', 'desc_only', 'ccds');


  foreach my $key (@sources) {
  my $source_id = $self->get_source_id_for_source_name($source_name, $key, $dbi);
    if(!(defined $source_id)){
      die 'Could not get source id for HGNC and '. $key ."\n";
    }
    $name_to_source_id{ $key } = $source_id;
  }

  my $source_id = $self->get_source_id_for_source_name('LRG_HGNC_notransfer', undef, $dbi);
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
  my $dbi        = $ref_arg->{dbi};

  if (defined $dead_name ) {     # dead name, add to synonym
    my @array2 = split ',\s*', $dead_name ;
    foreach my $arr (@array2){
      $self->add_to_syn($name, $source_id, $arr, $species_id, $dbi);
    }
  }

  if (defined $alias ) {     # alias, add to synonym
    my @array2 = split ',\s*', $alias;
    foreach my $arr (@array2){
      $self->add_to_syn($name, $source_id, $arr, $species_id, $dbi);
    }
  }
  return;
}

1;


