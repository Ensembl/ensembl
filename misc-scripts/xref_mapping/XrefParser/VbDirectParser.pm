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

package XrefParser::VbDirectParser;

use strict;
use warnings;
use Carp;
use File::Basename;
use Bio::SeqIO;
use XrefParser::BaseParser;

use base qw(XrefParser::BaseParser);

# OParser for FASTA-format probe mappings from Agilent
# >A_23_P253586
# CTGTCAGGATTCTAGAACTTCTAAAATTAAAGTTTGGGGAAATCAGTAGCTCTGATGAGA
# >A_23_P217507
# AGAAAGACGTTTTCCAACATGTAGAACTGCTTTTTAACTGGAGGAAAAATACTTCAGGAG

sub run {
  my ($self, $ref_arg) = @_;
  my $source_id    = $ref_arg->{source_id};
  my $species_id   = $ref_arg->{species_id};
  my $files        = $ref_arg->{files};
  my $verbose      = $ref_arg->{verbose};

  if((!defined $source_id) or (!defined $species_id) or (!defined $files)){
    croak "Need to pass source_id, species_id, files and rel_file as pairs";
  }
  $verbose |=0;

  my $file = @{$files}[0];

  open my $FH, "<", $file || return;
  my $i=0;
  print STDERR "source: ".$source_id."\tspecies: ".$species_id."\n";
  my $type = "transcript";

  while (my $ln = <$FH>) {
    chomp($ln);
    my ($probe,$id, $version, $description, $ensembl_id) = split("\t",$ln);      
    $i++;

    my $xref_id = $self->get_xref($probe, $source_id, $species_id);
    if (!defined($xref_id) || $xref_id eq "") {
      $xref_id = $self->add_xref({ acc        => $probe,
				   version    => 1,
				   label      => $probe,
				   desc       => $description,
				   source_id  => $source_id,
				   species_id => $species_id,
				   info_type  =>"DIRECT"} );
    }
    $self->add_direct_xref($xref_id, $ensembl_id, $type, $probe);
  }

  print $i." VB direct xrefs succesfully parsed\n" if($verbose);

  close($FH);

  return 0;

}

1;
