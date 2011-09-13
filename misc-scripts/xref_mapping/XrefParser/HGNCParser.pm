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

  if(!defined $species_id){
    $species_id = $self->get_species_id_for_filename($file);
  }

  my (%swissprot)  =  %{$self->get_valid_codes('Uniprot/SWISSPROT',$species_id)};
  my (%refseq) =  %{$self->get_valid_codes('refseq',$species_id)};
  my @list;
  push @list, 'refseq_peptide';
  push @list, 'refseq_dna';
  my (%entrezgene) = %{$self->get_valid_xrefs_for_dependencies('EntrezGene',@list)};


  my %name_count;
  my $mismatch = 0;

  my $hugo_io = $self->get_filehandle($file);

  if ( !defined $hugo_io ) {
    print "ERROR: Can't open HGNC file $file\n";
    return 1;
  }

  my ($name_to_source_id, $name_to_array_index)
    = $self->process_header($hugo_io->getline());

  while ( $_ = $hugo_io->getline() ) {
    chomp;
    my @array = split /\t/x, $_;

    my $seen = 0;

    my $acc    = $array[$name_to_array_index->{'desc_only'}];
    my $symbol = $array[$name_to_array_index->{'Approved Symbol'}];
    my $name   = $array[$name_to_array_index->{'Approved Name'}];
    my $previous_symbols = $array[$name_to_array_index->{'Previous Symbols'}];
    my $synonyms = $array[$name_to_array_index->{'Synonyms'}];


    my $type = 'Locus Specific Databases';
    my $id = $array[$name_to_array_index->{$type}];
    my $source_id = $name_to_source_id->{$type};
    if($id and $id =~ m/http:\/\/www.lrg-sequence.org\/LRG\/(LRG_\d+)\'/x){
      my $lrg_stable_id = $1;
      $self->add_to_direct_xrefs($lrg_stable_id, 'gene', $acc, $empty, $symbol,
						  $name, $empty, $source_id, $species_id);

      $self->add_synonyms_for_hgnc( {source_id  => $source_id,
				     name       => $acc,
				     species_id => $species_id,
				     dead       => $previous_symbols,
				     alias      => $synonyms});
      $name_count{$type}++;
    }


    #
    # Direct Ensembl mappings
    #
    $type = 'ensembl_manual';
    $id = $array[$name_to_array_index->{$type}];
    $source_id = $name_to_source_id->{$type};
    if ($id){              # Ensembl direct xref
      $seen = 1;
      $name_count{$type}++;
      $self->add_to_direct_xrefs($id,'gene', $acc, $empty,
						  $symbol, $name, $empty,
						  $source_id, $species_id);
      $self->add_synonyms_for_hgnc( {source_id  => $source_id,
				     name       => $acc,
				     species_id => $species_id,
				     dead       => $previous_symbols,
				     alias      => $synonyms});

    }

    #
    # RefSeq
    #
    foreach my $type (qw(refseq_mapped refseq_manual)){
      $id = $array[$name_to_array_index->{$type}];
      $source_id = $name_to_source_id->{$type};
      if ($id) {
	if(defined $refseq{$id} ){
	  $seen = 1;
	  foreach my $xref_id (@{$refseq{$id}}){
	    $name_count{$type}++;
	    $self->add_to_xrefs($xref_id, $acc, $empty, $symbol, $name, $empty,
						 $source_id, $species_id);
	  }
	  $self->add_synonyms_for_hgnc( {source_id  => $source_id,
					 name       => $acc,
					 species_id => $species_id,
					 dead       => $previous_symbols,
					 alias      => $synonyms});
	}
      }
    }

    #
    # Swissprot
    #
    $type = 'swissprot_manual';
    $id = $array[$name_to_array_index->{$type}];
    $source_id = $name_to_source_id->{$type};
    if ($id) {             # Swissprot
      if(defined $swissprot{$id} ){
	$seen = 1;
	foreach my $xref_id (@{$swissprot{$id}}){
	  $name_count{$type}++;
	  $self->add_to_xrefs($xref_id, $acc, $empty, $symbol, $name, $empty,
					       $source_id, $species_id);
	}
	  $self->add_synonyms_for_hgnc( {source_id  => $source_id,
					 name       => $acc,
					 species_id => $species_id,
					 dead       => $previous_symbols,
					 alias      => $synonyms});
      }
    }


    #
    # EntrezGene
    #
    foreach my $type (qw(entrezgene_manual entrezgene_mapped)){
      $id = $array[$name_to_array_index->{$type}];
      $source_id = $name_to_source_id->{$type};
      if(defined $id ){
	if(defined $entrezgene{$id} ){
	  $seen = 1;
	  $self->add_to_xrefs($entrezgene{$id}, $acc, $empty, $symbol, $name, $empty,
					       $source_id, $species_id);
	  $name_count{$type}++;
	  $self->add_synonyms_for_hgnc( {source_id  => $source_id,
					 name       => $acc,
					 species_id => $species_id,
					 dead       => $previous_symbols,
					 alias      => $synonyms});
	}
      }
    }


    if(!$seen){ # Store to keep descriptions etc
      $type = 'desc_only';
      $source_id = $name_to_source_id->{$type};
      $self->add_xref($acc, $empty, $symbol, $name, $source_id, $species_id, "MISC");
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



sub process_header{
  my $self = shift;
  my $header = shift;
  my %name_to_source_id;
  my %name_to_array_index;

  my %header_to_name = ( 'HGNC ID' => 'desc_only',
			 'Approved Symbol' => 0,
			 'Approved Name'   => 0,
			 'Previous Symbols' => 0,
			 'Synonyms' => 0,
			 'Entrez Gene ID' => 'entrezgene_manual',
			 'RefSeq IDs' => 'refseq_manual',
			 'Entrez Gene ID ' => 'entrezgene_mapped', #note space needed
			 'RefSeq ' => 'refseq_mapped',
			 'Ensembl Gene ID' => 'ensembl_manual',
			 'UniProt ID ' => 'swissprot_manual',
			 'Locus Specific Databases' => 0);


  foreach my $key (keys %header_to_name){
    if($header_to_name{$key}){
      my $source_id =
	$self->get_source_id_for_source_name('HGNC',$header_to_name{$key});
      if(!(defined $source_id)){
	die 'Could not get source id for '.$header_to_name{$key}."\n";
      }
      $name_to_source_id{ $header_to_name{$key} } = $source_id;
    }
  }

  my $source_id = $self->get_source_id_for_source_name('LRG_HGNC_notransfer');
  if(!(defined $source_id) ){
    die 'Could not get source id for LRG_HGNC_notransfer\n';
  }
  $name_to_source_id{'Locus Specific Databases'} = $source_id;

  chomp $header;
  my @items = split /\t/x, $header;
  my $i=0;
  foreach my $h (@items){
    $h =~ s/\(.*\)//x;
    if((defined $header_to_name{$h}) and $header_to_name{$h}){
      $name_to_array_index{$header_to_name{$h}} = $i++;
    }
    elsif(defined $header_to_name{$h}){
      $name_to_array_index{$h} = $i++;
    }
    else{
      die "Problem with $h not listed\n";
    }
  }

  return \%name_to_source_id, \%name_to_array_index;
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


