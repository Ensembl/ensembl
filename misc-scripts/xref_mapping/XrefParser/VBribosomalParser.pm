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

package XrefParser::VBribosomalParser;

use strict;
use warnings;
use Carp;
use POSIX qw(strftime);
use File::Basename;
use base qw( XrefParser::BaseParser );

# Parse the external description file
#
# RPG_DB ID     stable_id       feature Name    Description
# MOS10001	AGAP011777	gene	RpsA	40S ribosomal protein SA
# MOS10002	AGAP003768	gene	Rps2	40S ribosomal protein S2
# MOS10003	AGAP001910	gene	Rps3	40S ribosomal protein S3
# ...

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

  print "source_id = $source_id, species= $species_id, file = $file\n" if($verbose);

  my $added = 0;
  my $count = 0;

  my $file_io = $self->get_filehandle($file);

  if ( !defined $file_io ) {
    print STDERR "ERROR: Could not open file $file\n";
    return 1;
  }

  while ( my $line = $file_io->getline() ) {
    if ($line !~ /^#/) {
      chomp $line;
      my ($DB_id, $gene_id, $rien, $name, $description) = split("\t",$line);  #and use the gene_id as accession
      my $descr_full = "$DB_id - $description ($name)\n" ;

      my $xref_id = $self->get_xref($gene_id,$source_id, $species_id);
      if(!defined($xref_id)){
	$xref_id = $self->add_xref({ acc        => $gene_id,
				     label      => $gene_id,
				     desc       => $descr_full,
				     source_id  => $source_id,
				     species_id => $species_id,
				     info_type  => "DIRECT"} );
	$count++;
      }
      if(defined($gene_id) and $gene_id ne "-"){
	$self->add_direct_xref($xref_id, $gene_id, "Gene", "") ;
	$added++;
      }
    }
  }
  $file_io->close();

  print "Added $count xrefs and $added Direct xrefs to genes for VBribosomal\n" if($verbose);
  return 0;
}
1;
