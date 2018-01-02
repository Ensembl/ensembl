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

package XrefParser::AedesGenBankParser;

use strict;
use File::Basename;
use Carp;
use base qw( XrefParser::BaseParser );

#Aedes GenBank protein - because not yet in UniProt
#>EAT48991.1
#MGKSKAHRIKGLTGPKMSLGDQITEGRVSKKPKAPKIRLRAEEEEFVDSRTTKKILQQAR
#KQQAELNLLDDSFGPSLAESAAAASVGKRRHRLGDAASSDESDEEYREEADVDGQDFFDD
#IKINEEDERALEMFQNKDGVKTRTLADLIMDKITEKQTEIQTQFSDTGSLKMEEIDPRVR

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

  my $cpt = 0 ;

  next if (/^File:/);   # skip header

  my @xrefs;

  local $/ = "\n>";

  my $file_io = $self->get_filehandle($file);
  if ( !defined $file_io ) {
      print "Could not open $file\n";
      return 1;
  }

  while ( $_ = $file_io->getline() ) {

    my $xref;

    my ($header, $sequence) = $_ =~ /^>?(.+?)\n([^>]*)/s or warn("Can't parse FASTA entry: $_\n");
    #print "My header is -$header-\n" ;
    #print "My sequence is -$sequence-\n" ;

    if ($header eq "") {
      $header = "Aedes_GenBank".$cpt ;
      print STDERR "One sequence with a random name ... \n" if($verbose);
      $cpt++ ;
    }

    # deconstruct header - just use first part
    #my ($accession, $description) = split /\|/, $header;  #if description
    my $accession = $header;                               #if no description



    # make sequence into one long string
    $sequence =~ s/\n//g;

    # build the xref object and store it
    #print "ACCESSION & LABEL are $accession\n" ;
    #print "SEQUENCE is $sequence\n" ;
    #print "SOURCE_ID is $source_id\n" ;
    #print "SPECIES_ID is $species_id\n" ;
    #print "SEQUENCE_TYPE is peptide!\n";
    #print "STATUS is experimental!\n" ;

    $xref->{ACCESSION}     = $accession;
    $xref->{LABEL}         = $accession;
    #$xref->{DESCRIPTION}   = $description;
    $xref->{SEQUENCE}      = $sequence;
    $xref->{SOURCE_ID}     = $source_id;
    $xref->{SPECIES_ID}    = $species_id;
    $xref->{SEQUENCE_TYPE} = 'peptide';
    $xref->{STATUS}        = 'experimental';

    push @xrefs, $xref;

  }

  $file_io->close();

  print scalar(@xrefs) . " AedesGenBank xrefs succesfully parsed\n" if($verbose);

  $self->upload_xref_object_graphs(\@xrefs);

  print "Done\n";
  return 0;
}

1;
