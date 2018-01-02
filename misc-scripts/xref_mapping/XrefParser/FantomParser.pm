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

package XrefParser::FantomParser;

use strict;
use warnings;
use Carp;
use POSIX qw(strftime);
use File::Basename;
use File::Spec::Functions;

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

  my $dir = dirname($file);


  my (%embl) = %{$self->get_valid_codes("embl",$species_id)};

  my $fantom_io =
    $self->get_filehandle( $file  );

  if ( !defined $fantom_io ) {
    print STDERR "ERROR: Could not open " . $file . "\n" ;
    return 1;    # 1 error
  }

  my $ecount =0;

  my $mismatch=0;

  $fantom_io->getline(); # remove header

  while ( $_ = $fantom_io->getline() ) {
    chomp;
    my ($master, $label, $acc) = split (/\s+/,$_);
    if(defined($embl{$master})){
      foreach my $xref_id (@{$embl{$master}}){
	$self->add_dependent_xref({ master_xref_id => $xref_id,
				    acc            => $label,
				    label          => $label,
				    source_id      => $source_id,
				    species_id     => $species_id} );
	$ecount++;
      }
    }
    else{
      if($mismatch < 10){
	print "Could not find master $master\n";
      }
      $mismatch++;
    }
  }

  $fantom_io->close();

  if($verbose){
    print "\t$ecount xrefs from EMBL and\n";
    print "\t$mismatch xrefs ignored as no master found\n";
  }
  return 0;
}

1;
