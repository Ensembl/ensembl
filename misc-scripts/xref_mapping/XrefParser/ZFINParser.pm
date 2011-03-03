package XrefParser::ZFINParser;

use strict;
use POSIX qw(strftime);
use File::Basename;
use File::Spec::Functions;

use base qw( XrefParser::BaseParser );

# --------------------------------------------------------------------------------
# Parse command line and run if being run directly

if (!defined(caller())) {

  if (scalar(@ARGV) != 1) {
    print STDERR "\nUsage: ZFINParser.pm file <source_id> <species_id>\n\n";
    exit(1);
  }

  run($ARGV[0]);

}


sub run {

  my $self = shift if (defined(caller(1)));

  my $source_id = shift;
  my $species_id = shift;
  my $files       = shift;
  my $release_file   = shift;
  my $verbose       = shift;

  my $file = @{$files}[0];

  if(!defined($source_id)){
    $source_id = XrefParser::BaseParser->get_source_id_for_filename($file);
  }
  if(!defined($species_id)){
    $species_id = XrefParser::BaseParser->get_species_id_for_filename($file);
  }

  my $dir = dirname($file);
  

  my (%swiss) = %{XrefParser::BaseParser->get_valid_codes("uniprot",$species_id)};
  my (%refseq) = %{XrefParser::BaseParser->get_valid_codes("refseq",$species_id)};

  my $swissprot_io =
    $self->get_filehandle( catfile( $dir, 'uniprot.txt' ) );

  if ( !defined $swissprot_io ) {
    print STDERR "ERROR: Could not open " . catfile( $dir, 'uniprot.txt' ). "\n" ;
    return 1;    # 1 error
  }

#e.g.
#ZDB-GENE-000112-30      couptf2 O42532
#ZDB-GENE-000112-32      couptf3 O42533
#ZDB-GENE-000112-34      couptf4 O42534


  my %description;

  my $dbi = $self->dbi();


  my $sql = "insert into synonym (xref_id, synonym) values (?, ?)";
  my $add_syn_sth = $dbi->prepare($sql);    
  
#  my $syn_hash = $self->get_zfin_synonyms();
  


  #get the source ids for HGNC refseq, entrezgene and unitprot
  $sql = 'select source_id, priority_description from source where name like "ZFIN_ID"';
  my $sth = $dbi->prepare($sql);
  
  $sth->execute();


  my ($hgnc_source_id, $desc);
  $sth->bind_columns(\$hgnc_source_id, \$desc);
  my @arr;
  while($sth->fetch()){
    push @arr, $hgnc_source_id;
  }
  $sth->finish;
  
  $sql = "select accession, label, version,  description from xref where source_id in (".join(", ",@arr).")";
#  print "$sql\n";;
  $sth = $dbi->prepare($sql);
  $sth->execute();
  my ($acc, $lab, $ver);
  my $hgnc_loaded_count = 0;
  $sth->bind_columns(\$acc, \$lab, \$ver, \$desc);
  while (my @row = $sth->fetchrow_array()) {
    $description{$acc} = $desc if(defined($desc));
    $hgnc_loaded_count++;
  }
  $sth->finish;

  my $spcount =0;
  my $rscount =0;
  my $mismatch=0;

  while ( $_ = $swissprot_io->getline() ) {
    chomp;
    my ($zfin, $label, $acc) = split (/\s+/,$_);
    if(defined($swiss{$acc})){
      XrefParser::BaseParser->add_to_xrefs($swiss{$acc},$zfin,'',$label,$description{$zfin},'',$source_id,$species_id);
      $spcount++;
    }
    else{
      $mismatch++;
    }
  }

  $swissprot_io->close();

  my $refseq_io = $self->get_filehandle( catfile( $dir, 'refseq.txt' ) );

  if ( !defined $refseq_io ) {
    print STDERR "ERROR: Could not open " . catfile( $dir, 'refseq.txt' ),"\n" ;
    return 1;
  }

#ZDB-GENE-000125-12      igfbp2  NM_131458
#ZDB-GENE-000125-12      igfbp2  NP_571533
#ZDB-GENE-000125-4       dlc     NP_571019

  while ( $_ = $refseq_io->getline() ) {
    chomp;
    my ($zfin, $label, $acc) = split (/\s+/,$_);
    if(defined($refseq{$acc})){
      XrefParser::BaseParser->add_to_xrefs($refseq{$acc},$zfin,'',$label,$description{$zfin},'',$source_id,$species_id);
      $rscount++;
    }
    else{
      $mismatch++;
    }
  }

  $refseq_io->close();

  my (%zfin) = %{XrefParser::BaseParser->get_valid_codes("zfin",$species_id)};

  my $zfin_io = $self->get_filehandle( catfile( $dir, 'aliases.txt' ) );

  if ( !defined $zfin_io ) {
    print STDERR  "ERROR: Could not open " . catfile( $dir, 'aliases.txt' ), "\n" ;
    return 1;
  }

#ZDB-GENE-000125-4       deltaC  dlc     bea
#ZDB-GENE-000125-4       deltaC  dlc     beamter

  my $syncount = 0;

#  my $dbi = $self-> dbi();
#  my $sth =  
  $sth = $dbi->prepare('SELECT source_id from source where name like "ZFIN_ID"');

  $sth->execute;
  my $s1;
  $sth->bind_columns(\$s1);
  my $sources;
  while($sth->fetch()){
    push @$sources, $s1;
  }
  $sth->finish;

  while ( $_ = $zfin_io->getline() ) {
    chomp;
    my ($acc, undef, undef, $syn) = split (/\t/,$_);
    if(defined($zfin{$acc})){
      XrefParser::BaseParser->add_to_syn_for_mult_sources($acc, $sources, $syn, $species_id);
      $syncount++;
    }
  }

  $zfin_io->close();

  if($verbose){
    print "\t$spcount xrefs from UniProt and\n";
    print "\t$rscount xrefs from RefSeq succesfully loaded\n";
    print "\t$syncount synonyms loaded\n";
    print "\t$mismatch xrefs ignored\n";
  }
  return 0;
}

1;
