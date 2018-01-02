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

package XrefParser::PhytozomeGmaxParser;

use strict;
use warnings;
use Carp;
use POSIX qw(strftime);
use File::Basename;

use base qw( XrefParser::BaseParser );

# --------------------------------------------------------------------------------
# Parse command line and run if being run directly

sub run {

 my ($self, $ref_arg) = @_;
  my $source_id    = $ref_arg->{source_id};
  my $species_id   = $ref_arg->{species_id};
  my $files        = $ref_arg->{files};
  my $verbose      = $ref_arg->{verbose};

  if((!defined $source_id) or (!defined $species_id) or (!defined $files) ){
    croak "Need to pass source_id, species_id and  files as pairs";
  }
  $verbose |=0;

  my $file = @{$files}[0];

  my $gene_source_id = $self->get_source_id_for_source_name("PHYTOZOME_GMAX_GENE");

  my $gmax_io = $self->get_filehandle($file);

  if ( !defined $gmax_io ) {
    print STDERR "ERROR: Could not open $file\n";
    return 1;    # 1 is an error
  }

  my $xref_count =0;
  my $syn_count =0;

  while ( $_ = $gmax_io->getline() ) {

    chomp;

    if ($_ =~ /^([^\t]+)\t[^\t]*\t[^\t]*\t[^\t]*\t[^\t]*\t[^\t]*\t[^\t]*\t[^\t]*\t([^\t]*)/) {
	my ($gmax_id, $desc) = ($1,$2);
	
	if ($verbose) {
	    #print STDERR "gmax_id, $gmax_id\n";
	    #print STDERR "desc: $desc\n";
	}
	
	my $locus_name = $gmax_id;
	my $gene_xref_id = undef;
	
	if ((defined $desc) && ($desc ne "")) {
	    $gene_xref_id = $self->add_xref({ acc        => $gmax_id,
						 label      => $locus_name,
						 desc       => $desc,
						 source_id  => $gene_source_id,
						 species_id => $species_id,
						 info_type  => "DIRECT"} );
	}
	else {
	    
	    # no description given
	    
	    $gene_xref_id = $self->add_xref({ acc        => $gmax_id,
						 label      => $locus_name,
						 source_id  => $gene_source_id,
						 species_id => $species_id,
						 info_type  => "DIRECT"} );
	}
	$self->add_direct_xref($gene_xref_id, $gmax_id, "Gene", "DIRECT");
	
	$xref_count++;
	
    }
    else {
	if ($verbose) {
	    print STDERR "failed to parse line, $_\n\n";
	}
    }
  }
 
 $gmax_io->close();
 
 print $xref_count." Phytozome_GMAX Xrefs added\n" if($verbose);
 return 0; #successful
}

1;
