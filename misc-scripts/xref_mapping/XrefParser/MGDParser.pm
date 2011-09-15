package XrefParser::MGDParser;

use strict;
use warnings;
use POSIX qw(strftime);
use File::Basename;

use base qw( XrefParser::BaseParser );

sub run {
  my $self = shift;
  my $source_id = shift;
  my $species_id = shift;
  my $files       = shift;
  my $release_file   = shift;
  my $verbose       = shift;

  my $file = @{$files}[0];

  die "No longer used. MGI is taken form the uniprot file\n";
}

1;
