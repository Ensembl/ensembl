package XrefParser::FantomParser;

use strict;
use POSIX qw(strftime);
use File::Basename;
use File::Spec::Functions;

use base qw( XrefParser::BaseParser );

# --------------------------------------------------------------------------------
# Parse command line and run if being run directly

if (!defined(caller())) {

  if (scalar(@ARGV) != 1) {
    print STDERR "\nUsage: FantomParser.pm file <source_id> <species_id>\n\n";
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
  

  my (%embl) = %{XrefParser::BaseParser->get_valid_codes("embl",$species_id)};
  my (%ddbj) = %{XrefParser::BaseParser->get_valid_codes("ddbj",$species_id)};
  my (%genbank) = %{XrefParser::BaseParser->get_valid_codes("genbank",$species_id)};

  my $fantom_io =
    $self->get_filehandle( $file  );

  if ( !defined $fantom_io ) {
    print STDERR "ERROR: Could not open " . $file . "\n" ;
    return 1;    # 1 error
  }

#e.g.


  my $ecount =0;
  my $dcount =0;
  my $gcount =0;
  
  my $mismatch=0;

  $fantom_io->getline(); # remove header

  while ( $_ = $fantom_io->getline() ) {
    chomp;
    my ($master, $label, $acc) = split (/\s+/,$_);
    if(defined($embl{$master})){
      XrefParser::BaseParser->add_to_xrefs($embl{$master},$label,'',$label,'','',$source_id,$species_id);
      $ecount++;
    }
    elsif(defined($ddbj{$master})){
      XrefParser::BaseParser->add_to_xrefs($ddbj{$master},$label,'',$label,'','',$source_id,$species_id);
      $dcount++;
    }
    elsif(defined($genbank{$master})){
      XrefParser::BaseParser->add_to_xrefs($genbank{$master},$label,'',$label,'','',$source_id,$species_id);
      $gcount++;
    }
    else{
      if($mismatch < 10){
	print STDERR "Could not find master $master\n";
      }
      $mismatch++;
    }
  }

  $fantom_io->close();

  if($verbose){
    print "\t$ecount xrefs from EMBL and\n";
    print "\t$dcount xrefs from DDBJ succesfully loaded\n";
    print "\t$gcount xrefs from GenBank\n";
    print "\t$mismatch xrefs ignored as no master found\n";
  }
  return 0;
}

1;
