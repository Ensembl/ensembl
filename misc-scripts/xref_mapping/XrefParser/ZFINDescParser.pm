package XrefParser::ZFINDescParser;

use strict;
use POSIX qw(strftime);
use File::Basename;
use File::Spec::Functions;

use base qw( XrefParser::BaseParser );


sub run {

  my $self = shift;
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


  open( my $FH, "<", $file) || die "could not open file $file";
  

#e.g.
#ZDB-GENE-050102-6       WITHDRAWN:zgc:92147     WITHDRAWN:zgc:92147     0
#ZDB-GENE-060824-3       apobec1 complementation factor  a1cf    0
#ZDB-GENE-090212-1       alpha-2-macroglobulin-like      a2ml    15      ZDB-PUB-030703-1


  my $count =0;
  my $withdrawn = 0;
  while ( <$FH> ) {
    chomp;
    my ($zfin, $desc, $label) = split (/\t/,$_);

    if($label =~ /^WITHDRAWN/){
      $withdrawn++;
    }
    else{
      $self->add_xref($zfin,"",$label,$desc,$source_id,$species_id,"MISC");
      $count++;
    }
  }
  close($FH);

  if($verbose){
    print "\t$count xrefs added, $withdrawn withdrawn entries ignored\n";
  }
  return 0;
}

1;
