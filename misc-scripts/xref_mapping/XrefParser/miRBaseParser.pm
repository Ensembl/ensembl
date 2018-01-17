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

package XrefParser::miRBaseParser;

use strict;
use warnings;
use Carp;
use DBI;

use base qw(XrefParser::BaseParser);


sub run {

  my ($self, $ref_arg) = @_;
  my $source_id    = $ref_arg->{source_id};
  my $species_id   = $ref_arg->{species_id};
  my $species_name = $ref_arg->{species};
  my $files        = $ref_arg->{files};
  my $verbose      = $ref_arg->{verbose};
  my $dbi          = $ref_arg->{dbi};
  $dbi = $self->dbi unless defined $dbi;

  if((!defined $source_id) or (!defined $species_id) or (!defined $files) ){
    croak "Need to pass source_id, species_id and files as pairs";
  }
  $verbose |=0;

  my $file = @{$files}[0];

  if(!defined($species_id)){
    $species_id = $self->get_species_id_for_filename($file);
  }

  my $xrefs = $self->create_xrefs($source_id, $file, $species_id, $dbi, $species_name);
  if(!defined($xrefs)){
    return 1; #error
  }
  # upload
  if(!defined($self->upload_xref_object_graphs($xrefs, $dbi))){
    return 1; 
  }
  return 0; # successfull

}

# --------------------------------------------------------------------------------
# Parse file into array of xref objects

sub create_xrefs {

  my ($self, $source_id, $file, $species_id, $dbi, $species_name) = @_;

  my %species2name = $self->species_id2name($dbi);
  if (defined $species_name) { push @{$species2name{$species_id}}, $species_name; }
  if (!defined $species2name{$species_id}) { next; }
  my @names   = @{$species2name{$species_id}};

  my %name2species_id     = map{ $_=>$species_id } @names;

  my $file_io = $self->get_filehandle($file);
  if ( !defined $file_io ) {
    print STDERR "ERROR: Could not open $file\n";
    return 1;    # 1 is an error
  }

  my @xrefs;

  local $/ = "\n\/\/";

  while ($_ = $file_io->getline()) {

    my $xref;

    my $entry = $_;
    chomp $entry;

    next if (!$entry);

    my ($header, $sequence) = split (/\nSQ/, $entry, 2);
    # remove newlines
    my @seq_lines = split (/\n/, $sequence) if ($sequence);
    # drop the information line
    shift @seq_lines;
    # put onto one line
    $sequence = join("", @seq_lines);
    # make uppercase
    $sequence = uc($sequence);
    # replace Ts for Us
    $sequence =~ s/U/T/g;
    # remove numbers and whitespace
    $sequence =~ s/[\d+,\s+]//g;

#    print "$header\n";
    my ($name) = $header =~ /\nID\s+(\S+)\s+/;
    my ($acc) = $header =~ /\nAC\s+(\S+);\s+/;
    my ($description) = $header =~ /\nDE\s+(.+)\s+stem-loop/;
    my @description_parts = split (/\s+/, $description) if ($description);
    # remove the miRNA identifier
    pop @description_parts;
    my $species =  join(" ", @description_parts);
    $xref->{SEQUENCE_TYPE} = 'dna';
    $xref->{STATUS} = 'experimental';
    $xref->{SOURCE_ID} = $source_id;
    $species = lc $species;
    $species =~ s/ /_/;
    
    my $species_id_check = $name2species_id{$species};


    next if (!defined($species_id_check));
    
    # skip xrefs for species that aren't in the species table
    if (defined($species_id) and $species_id == $species_id_check) {
      
      $xref->{ACCESSION} = $acc;
      $xref->{LABEL} = $name;
      $xref->{DESCRIPTION} = $name;
      $xref->{SEQUENCE} = $sequence;
      $xref->{SPECIES_ID} = $species_id;

      # TODO synonyms, dependent xrefs etc
      push @xrefs, $xref;
    }
  }

  $file_io->close();

  print "Read " . scalar(@xrefs) ." xrefs from $file\n";
 
  return \@xrefs;

}


1;
