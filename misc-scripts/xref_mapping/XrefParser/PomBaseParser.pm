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

package XrefParser::PomBaseParser;

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
    croak "Need to pass source_id, species_id and files as pairs";
  }
  $verbose |=0;

  my $file = @{$files}[0];

  my $gene_source_id = $self->get_source_id_for_source_name("PomBase_GENE", undef, $dbi);
  my $transcript_source_id = $self->get_source_id_for_source_name("PomBase_TRANSCRIPT", undef, $dbi);

  my $pombase_io = $self->get_filehandle($file);

  if ( !defined $pombase_io ) {
    print STDERR "ERROR: Could not open $file\n";
    return 1;    # 1 is an error
  }

  my $xref_count =0;
  my $syn_count =0;

  while ( $_ = $pombase_io->getline() ) {

    chomp;

    if ($_ =~ /^([^\t]+)\t([^\t]+)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)$/) {

	    my @line = split(m/\t/ms, $_);
	    my ($pombase_id, $name, $info_type, $biotype, $external_db_source, $desc, $ensembl_object_type, $synonyms) = undef;

	    $pombase_id          = $line[0];
	    $name                = $line[1];
	    $info_type           = $line[2];
	    $biotype             = $line[3];
            $external_db_source  = $line[4];
	    $desc                = $line[5];
	    $ensembl_object_type = $line[6];

	    if (scalar @line == 8) {
	        $synonyms = $line[7];
	    }
	    # parse the lines corresponding to the gene entries
	    # and filter out lines corresponding to the CDS for example

            #print "$ensembl_object_type\n"; 
	    if ($ensembl_object_type eq 'Gene') {
	        my $ensembl_xref_id = $self->add_xref({ acc        => $pombase_id,
							label      => $name,
							desc       => $desc,
							source_id  => $gene_source_id,
							species_id => $species_id,
                                                        dbi        => $dbi,
							info_type  => $info_type} );

	        $self->add_direct_xref($ensembl_xref_id, $pombase_id, $ensembl_object_type, $info_type, $dbi);
	    } elsif ($ensembl_object_type eq 'Transcript') {
	        my $ensembl_xref_id = $self->add_xref({ acc        => $pombase_id,
							label      => $name,
							desc       => $desc,
                                                        dbi        => $dbi,
							source_id  => $transcript_source_id,
							species_id => $species_id,
							info_type  => $info_type} );

	        $self->add_direct_xref($ensembl_xref_id, $pombase_id, $ensembl_object_type, $info_type, $dbi);
	    }

	    $xref_count++;
	    if ($synonyms) {
	   	 my (@syn) = split(/,/,$synonyms);
	    	foreach my $synonym (@syn){
			    if ($verbose) {
			        print STDERR "adding synonym, $synonym\n";
			    }
			    $self->add_to_syn($pombase_id, $gene_source_id, $synonym, $species_id, $dbi);
			    $syn_count++;
		    }
	    }
    } else {
	    if ($verbose) {
	        print STDERR "failed to parse line, $_\n\n";
	    }
    }
  }

  $pombase_io->close();

  print $xref_count." PomBase Xrefs added with $syn_count synonyms\n" if($verbose);
  return 0; #successful
}

1;
