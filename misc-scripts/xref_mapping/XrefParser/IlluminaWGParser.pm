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

package XrefParser::IlluminaWGParser;

use strict;
use Carp;
use base qw( XrefParser::BaseParser );


sub run {

  my ($self, $ref_arg) = @_;
  my $source_id    = $ref_arg->{source_id};
  my $species_id   = $ref_arg->{species_id};
  my $files        = $ref_arg->{files};
  my $release_file = $ref_arg->{rel_file};
  my $verbose      = $ref_arg->{verbose};

  if((!defined $source_id) or (!defined $species_id) or (!defined $files) or (!defined $release_file)){
    croak "Need to pass source_id, species_id, files and rel_file as pairs";
  }
  $verbose |=0;

  my $file = @{$files}[0];


  my @xrefs;

  my $file_io = $self->get_filehandle($file);

  if ( !defined $file_io ) {
    print STDERR "Could not open $file\n";
    return 1;
  }

  my $read = 0;
  my $name_index;
  my $seq_index;
  my $defin_index;
  while ( $_ = $file_io->getline() ) {
    chomp;

    my $xref;

    # strip ^M at end of line
    $_ =~ s/\015//g;

    if(/^\[/){
#      print $_."\n";
      if(/^\[Probes/){
	my $header = $file_io->getline();
#	print $header."\n";
	$read =1;
	my @bits = split("\t", $header);
	my $index =0;
	foreach my $head (@bits){
	  if($head eq "Search_Key"){
	    $name_index = $index;
	  }
	  elsif($head eq "Probe_Sequence"){
	    $seq_index = $index;
	  }
	  elsif($head eq "Definition"){
	    $defin_index = $index;
	  }
	  $index++;
	}
	if(!defined($name_index) or !defined($seq_index) or !defined($defin_index)){
	  die "Could not find index for search_key->$name_index, seq->$seq_index, definition->$defin_index";
	}
      
	next;
      }
      else{
	$read = 0;
      }	
    }
    if($read){
#      print $_."\n";
      my @bits = split("\t", $_);
      my $sequence = $bits[$seq_index];

      my $description = $bits[$defin_index];
      my $illumina_id = $bits[$name_index];

      # build the xref object and store it
      $xref->{ACCESSION}     = $illumina_id;
      $xref->{LABEL}         = $illumina_id;
      $xref->{SEQUENCE}      = $sequence;
      $xref->{SOURCE_ID}     = $source_id;
      $xref->{SPECIES_ID}    = $species_id;
      $xref->{DESCRIPTION}   = $description;
      $xref->{SEQUENCE_TYPE} = 'dna';
      $xref->{STATUS}        = 'experimental';

      push @xrefs, $xref;
      
    }
  }

  $file_io->close();

  $self->upload_xref_object_graphs(\@xrefs);


  print scalar(@xrefs) . " Illumina V2 xrefs succesfully parsed\n" if($verbose);

  return 0;
}

1;
