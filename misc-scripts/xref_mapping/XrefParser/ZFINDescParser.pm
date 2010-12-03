package XrefParser::ZFINDescParser;

use strict;
use POSIX qw(strftime);
use File::Basename;
use File::Spec::Functions;

use base qw( XrefParser::BaseParser );

# --------------------------------------------------------------------------------
# Parse command line and run if being run directly

if (!defined(caller())) {

  if (scalar(@ARGV) != 1) {
    print STDERR "\nUsage: ZFINDescParser.pm file <source_id> <species_id>\n\n";
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


  open(FH,"<$file") || die "could not open file $file";
  

#e.g.
#ZDB-GENE-000112-30      couptf2 long_desc

  my $count =0;
  while ( <FH> ) {
    chomp;
    my ($zfin, $label, $desc) = split (/\t/,$_);

    $self->add_xref($zfin,"",$label,$desc,$source_id,$species_id,"MISC");
    $count++;
  }

  if($verbose){
    print "\t$count xrefs added\n";
  }
  return 0;
}

1;
