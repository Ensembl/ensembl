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

package XrefParser::ncRNAParser;
 
use strict;
use warnings;
use Carp;
use POSIX qw(strftime);
use File::Basename;
 
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

  print "source_id = $source_id, species= $species_id, file = $file\n";

  my %name_2_source_id=();
  my $added=0;

  my $file_io = $self->get_filehandle($file);

  if ( !defined $file_io ) {
    print STDERR "ERROR: Could not open file $file\n";
    return 1;
  }

  while ( my $line = $file_io->getline() ) {
    chomp $line;
    my ($gene_id,$transcript_id,$source_name,$acc,$display_label,$full_description, $status)
      = split("\t",$line);

    #trim the description.
    my ($description,$junk) = split("[[]Source:",$full_description);
    if($source_name eq "miRNA_Registry"){
      if($status eq "KNOWN"){
	$source_name = "miRBase";
      }
      else{
	$source_name = "miRBase_predicted";
      }
    }
    if(!defined($name_2_source_id{$source_name})){
      my $tmp = $self->get_source_id_for_source_name($source_name);
      if(!$tmp){
	die("Could not get source_id for $source_name\n");
      }
      $name_2_source_id{$source_name} = $tmp;
    }
    my $xref_id = $self->get_xref($acc,$name_2_source_id{$source_name}, $species_id);
    if(!defined($xref_id)){
      $xref_id = $self->add_xref({ acc        => $acc,
				   label      => $display_label,
				   desc       => $description,
				   source_id  => $name_2_source_id{$source_name},
				   species_id => $species_id,
				   info_type  => "DIRECT"} );
      $added++;
    }
    $self->add_direct_xref($xref_id, $transcript_id, "Transcript", "") if (defined($transcript_id));    

    #just add to the transcript ONLY as the check at the end will move all
    #the those mapped to the transcript to the genes anyway due to the
    #biomart check
    #    $self->add_direct_xref($xref_id, $gene_id, "Gene", "")             if (defined($gene_id)); 
  }

  $file_io->close();

  print "Added $added Xrefs for ncRNAs\n" if($verbose);
  return 0;
}

1;
