package XrefParser::MGDParser;

use strict;
use warnings;
use POSIX qw(strftime);
use File::Basename;

use base qw( XrefParser::BaseParser );

my $xref_sth ;
my $dep_sth;

sub run {
  my $self = shift;
  my $source_id = shift;
  my $species_id = shift;
  my $files       = shift;
  my $release_file   = shift;
  $verbose       = shift;

  my $file = @{$files}[0];

  die "No longer used. MGI is taken form the uniprot file\n";
}
1;

__END__

  if(!defined($source_id)){
    $source_id = XrefParser::BaseParser->get_source_id_for_filename($file);
  }
  if(!defined($species_id)){
    $species_id = XrefParser::BaseParser->get_species_id_for_filename($file);
  }

  my (%swiss) = %{XrefParser::BaseParser->get_valid_codes("uniprot",$species_id)};
  my (%refseq) = %{XrefParser::BaseParser->get_valid_codes("refseq",$species_id)};


  my $count = 0;
  my $mismatch = 0;
  my %mgi_good;

  my $file_io = $self->get_filehandle($file);

  if ( !defined $file_io ) {
    print "ERROR: Could not open file $file";
    return 1;    # 1 is an error
  }

  while ( my $line = $file_io->getline() ) {
    chomp $line;
    my ($key,$label,$desc,$sps) = (split("\t",$line))[0,1,3,6];
    my @sp = split(/\s/,$sps); 
    foreach my $value (@sp){
      if(defined($value) and $value and defined($swiss{$value})){
	XrefParser::BaseParser->add_to_xrefs($swiss{$value},$key,'',$label,$desc,"",$source_id,$species_id);
	$mgi_good{$key} = 1;
	$count++;
      }
      elsif(defined($value) and $value and defined($refseq{$value})){
	$mismatch++;
      }
    }

  }
  $file_io->close();

  my $dir = dirname($file);
  my $syn_file = $dir."/MRK_Synonym.sql.rpt";

  $file_io = $self->get_filehandle($syn_file);

  if ( !defined $file_io ) {
    print "ERROR: Could not open file $syn_file";
    return 1;
  }

  my $synonyms=0;

  while ( $_ = $file_io->getline() ) {
    if(/MGI:/){
      chomp ;
      my ($key,$syn) = (split)[0,4];
      if(defined($mgi_good{$key})){
	$self->add_to_syn($key, $source_id, $syn, $species_id);
	$synonyms++;
      }
    }
  }

  $file_io->close();

  print "\t$count xrefs succesfully loaded\n";
  print "\t$synonyms synonyms successfully loaded\n";
  print "\t$mismatch xrefs failed to load\n";
     
  return 0;
}

1;
    
