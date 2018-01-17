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

package XrefParser::SGDParser;

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
  my $dbi          = $ref_arg->{dbi};
  $dbi = $self->dbi unless defined $dbi;

  if((!defined $source_id) or (!defined $species_id) or (!defined $files) ){
    croak "Need to pass source_id, species_id and  files as pairs";
  }
  $verbose |=0;

  my $file = @{$files}[0];

  my $gene_source_id = $self->get_source_id_for_source_name("SGD_GENE", undef, $dbi);
  #my $transcript_source_id = $self->get_source_id_for_source_name("SGD_TRANSCRIPT");
  my $translation_source_id = $self->get_source_id_for_source_name("SGD_TRANSLATION", undef, $dbi);

  my $sgd_io = $self->get_filehandle($file);

  if ( !defined $sgd_io ) {
    print STDERR "ERROR: Could not open $file\n";
    return 1;    # 1 is an error
  }

  my $xref_count =0;
  my $syn_count =0;

  while ( $_ = $sgd_io->getline() ) {

    chomp;

    if ($_ =~ /^([^\t]+)\t([^\t]+)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t[^\t]+\t[^\t]*\t[^\t]+\t[^\t]+\t[^\t]+\t[^\t]+\t[^\t]*\t[^\t]+\t[^\t]+\t([^\t]*)$/) {
	my ($sgd_id, $biotype, $status, $orf_name, $locus_name, $alias_name, $desc) = ($1,$2,$3,$4,$5,$6,$7);
	
	# parse the lines corresponding to the gene entries
	# and filter out lines corresponding to the CDS for example
	
	if ($biotype =~ /ORF|.+RNA|transposable_element_gene|pseudogene/) {

	    if ($verbose) {
		#print STDERR "parsing line for biotype, $biotype\n";
		#print STDERR "sgd_id, biotype, status, orf_name, locus_name, alias_name, $sgd_id, $biotype, $status, $orf_name, $locus_name, $alias_name\n";
		#print STDERR "desc: $desc\n";
	    }

	    if (!defined $locus_name || ($locus_name eq "")) {
		
		if (!defined $orf_name || ($orf_name eq "")) {
		    print STDERR "can't assign the orf_name as the locus_name!\n";
		}
		else {
		    if ($verbose) {
			#print STDERR "assigning the orf_name as the locus_name(ie the gene_name)\n";
		    }
		    $locus_name = $orf_name;
		}
	    }
	    
	    my (@syn) = split(/\|/,$alias_name);

	    my $gene_xref_id = $self->add_xref({ acc        => $sgd_id,
						 label      => $locus_name,
						 desc       => $desc,
						 source_id  => $gene_source_id,
						 species_id => $species_id,
                                                 dbi        => $dbi,
						 info_type  => "DIRECT"} );
	    $self->add_direct_xref($gene_xref_id, $orf_name, "Gene", "DIRECT", $dbi);

	    my $translation_xref_id = $self->add_xref({ acc        => $sgd_id,
			 			        label      => $locus_name,
						        desc       => $desc,
						        source_id  => $translation_source_id,
						        species_id => $species_id,
                                                        dbi        => $dbi,
						        info_type  => "DIRECT"} );
	    $self->add_direct_xref($translation_xref_id, $orf_name, "Translation", "DIRECT", $dbi);

	    $xref_count++;

	    foreach my $synonym (@syn){
		if ($verbose) {
		    # print STDERR "adding synonym, $synonym\n";
		}
		$self->add_to_syn($sgd_id, $gene_source_id, $synonym, $species_id, $dbi);
		$syn_count++;
	    }
	    
	}
	else {
	    if ($verbose) {
		#print STDERR "filtering biotype, $biotype\n";
	    }
	}
    }
    else {
	if ($verbose) {
	    print STDERR "failed to parse line, $_\n\n";
	}
    }
  }

  $sgd_io->close();

  print $xref_count." SGD Xrefs added with $syn_count synonyms\n" if($verbose);
  return 0; #successful
}

1;
