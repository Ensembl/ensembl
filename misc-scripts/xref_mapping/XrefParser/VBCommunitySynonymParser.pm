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

package XrefParser::VBCommunitySynonymParser;
 
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
  my $release_file = $ref_arg->{rel_file};
  my $verbose      = $ref_arg->{verbose};

  if ((!defined $source_id) or (!defined $species_id) or (!defined $files)) {
    croak "Need to pass source_id, species_id, files and rel_file as pairs";
  }
  $verbose |=0;

  my $file = @{$files}[0];
  #The Synonyms are for gene symboles. Therefore using the source_id of symboels
  my $symboel_source_id;
  open IN,"<","source_id_file" or die "could not open source_id_file";
  while(my $line = <IN>){
  	  chomp $line;
  	  my($organism_id,$source_id) = split /\t/ , $line;
  	  if($organism_id eq $species_id){
  	  	$symboel_source_id = $source_id;
  	  }
  }
  close IN;
  #my $symboel_source_id = $self->get_source_id_for_source_name("VB_Community_Annotation");
  print "source_id = $symboel_source_id, species= $species_id, file = $file\n" if $verbose;

  my $added = 0;
  my $count = 0;

  my $file_io = $self->get_filehandle($file);

  if ( !defined $file_io ) {
    print STDERR "ERROR: Could not open file $file\n";

    return 1;
  }

  while ( my $line = $file_io->getline() ) {
    next unless $line =~ /^\w+/;

    chomp $line;
    my ($stable_id,$synonym) = split "\t", $line;
    $count++;

    
    #Add synonym for an xref given by accession and source_id
    $self->add_to_syn($stable_id, $symboel_source_id,$synonym,$species_id);

  }

  $file_io->close();

  print "Added $count synonyms to synonym  for VBCommunitySynonyms\n" if $verbose;

  return 0;
}

1;
