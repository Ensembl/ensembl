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

package XrefParser::XenopusJamboreeParser;

# Parse annotated peptides from Xenopus Jamboree

use strict;
use warnings;
use Carp;
use File::Basename;
use Text::CSV;

use parent qw( XrefParser::BaseParser );



sub run {
 my ($self, $ref_arg) = @_;
  my $source_id    = $ref_arg->{source_id};
  my $species_id   = $ref_arg->{species_id};
  my $files        = $ref_arg->{files};
  my $verbose      = $ref_arg->{verbose} || 0;
  my $dbi          = $ref_arg->{dbi} || $self->dbi;

  if((!defined $source_id) or (!defined $species_id) or (!defined $files) ){
    croak "Need to pass source_id, species_id and files as pairs";
  }

  my $file = shift @{$files};

  my $file_io = $self->get_filehandle($file);

  if ( !defined $file_io ) {
    print STDERR "ERROR: Could not open $file\n";
    return 1;    # 1 error
  }

  my $input_file = Text::CSV->new({
                                    sep_char           => "\t",
				    empty_is_undef     => 1,
				   }) or croak "Cannot use file $file: ".Text::CSV->error_diag ();

  $input_file->column_names([qw(acc label desc stable_id)] );
  my $count = 0;
  while ( my $data = $input_file->getline_hr( $file_io ) ) {

    my $desc;
    if(defined $data->{'desc'}){
      $desc = $self->parse_description($data->{'desc'});
    }

    if($data->{'label'} eq "unnamed"){
      $data->{'label'} = $data->{'acc'};
    }

    $self->add_to_direct_xrefs({ stable_id  => $data->{'stable_id'},
				 type       => 'gene',
				 acc        => $data->{'acc'},
				 label      => $data->{'label'},
				 desc       => $desc,
                                 dbi        => $dbi,
				 source_id  => $source_id,
				 species_id => $species_id });
    $count++;
  }

  $file_io->close();

  print $count . " XenopusJamboreeParser xrefs succesfully parsed\n" if($verbose);

  return 0;
}


=begin comment
Regex handles lines in the following desc formats

XB-GENE-940410	unnamed	Putative ortholog of g2/mitotic-specific cyclin B3, 3 of 14	ENSXETG00000007206
XB-GENE-956173	hba4	alpha-T4 globin, Putative ortholog of hemoglobin alpha chain. [Source:Uniprot/SWISSPROT;Acc:P01922], 2 of 3	ENSXETG00000001141

=end comment
=cut
sub parse_description{
  my ($self, $desc) = @_;

  # Remove some provenance information encoded in the description
  $desc =~ s/\[.*\]//;
  # Remove labels of type 5 of 14 from the description
  $desc =~ s/ , [0-9]+ of [0-9]+//;
  return $desc;
}


1;
