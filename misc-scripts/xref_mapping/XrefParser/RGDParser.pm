package XrefParser::RGDParser;

use strict;
use warnings;
use Carp;
use POSIX qw(strftime);
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

  my $dbi = $self->dbi();

  my (%refseq) = %{XrefParser::BaseParser->get_valid_codes("refseq",$species_id)};

  my $rgd_io = $self->get_filehandle($file);

  if ( !defined $rgd_io ) {
    print "ERROR: Could not open $file\n";
    return 1;
  }
  my $line;
  my $found =0;
  while((!$found) and ($line = $rgd_io->getline())){ # ignore comments
    if(!($line =~ /^#/)){
      $found = 1;
    }
  };
  chomp $line;
  my @linearr = split(/\t/,$line);

  #
  #warn if sanity check fails
  #

  if($linearr[0] =~ /GENE_RDB_ID/){
   die ($linearr[0]."!= GENE_RDB_ID is not the first element in the header\n$line\n");
  }
  if($linearr[1] ne "SYMBOL"){
    die ("SYMBOL is not the second element in the header\n$line\n");
  }
  if($linearr[2] ne "NAME"){
    die ("NAME is not the third element in the header\n$line\n");
  }
  if($linearr[23] ne "GENBANK_NUCLEOTIDE"){
    die ("GENBANK_NUCLEOTIDE is not the twentysixth element in the header but ".$linearr[23]." is.\\n");
  }
  if($linearr[29] ne "OLD_SYMBOL"){
    die ("NAME is not the third element in the header\n$line\n");
  }  

  my $sql = "insert into synonym (xref_id, synonym) values (?, ?)";
  my $add_syn_sth = $dbi->prepare($sql);    
  
  my $count= 0;
  my $mismatch = 0;
  my $syn_count = 0;
  while ( $line = $rgd_io->getline() ) {
    chomp $line;
    my ($rgd, $symbol, $name, $refseq,$old_name) = (split (/\t/,$line))[0,1,2,23,29];
    my @nucs = split(/\;/,$refseq);
    my $done = 0;
    my $failed_list ="";
    foreach my $nuc (reverse @nucs){
      if(!$done){
	my $xref=undef; 
	if(defined($refseq{$nuc})){
	  foreach my $xref (@{$refseq{$nuc}}){
	    $done = 1;
	    my $xref_id = XrefParser::BaseParser->add_to_xrefs($xref,$rgd,"",$symbol,$name,"",$source_id,$species_id);
	    $count++;
	    my @syns  = split(/\;/,$old_name);
	    foreach my $syn(@syns){
	      $add_syn_sth->execute($xref_id, $syn);
	      $syn_count++;
	    }
	  }
	}
	else{
	  $failed_list .= " $nuc";
	}
      }
    }

    if(!$done){
#      print STDERR "$rgd FAILED for $failed_list\n";
      $self->add_xref($rgd,"",$symbol,$name,$source_id,$species_id,"MISC");
      $mismatch++;
    }

  }

  $rgd_io->close();

  print "\t$count xrefs succesfully loaded and dependent on refseq\n" if($verbose);
  print "\t$mismatch xrefs added but with NO dependencies\n" if($verbose);
  print "added $syn_count synonyms\n" if($verbose);
  return 0;
}

1;
