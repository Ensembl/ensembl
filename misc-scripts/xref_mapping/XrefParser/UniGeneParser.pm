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

# Parse UniGene Hs.seq.uniq files to create xrefs.

package XrefParser::UniGeneParser;

use strict;
use warnings;
use Carp;
use File::Basename;

use base qw( XrefParser::BaseParser );

sub run {

  my ($self, $ref_arg) = @_;
  my $source_id    = $ref_arg->{source_id};
  my $species_id   = $ref_arg->{species_id};
  my $files        = $ref_arg->{files};
  my $release_file = $ref_arg->{rel_file};
  my $verbose      = $ref_arg->{verbose};

  if((!defined $source_id) or (!defined $species_id) or (!defined $files)){
    croak "Need to pass source_id, species_id and  files as pairs";
  }
  $verbose |=0;

  my $uniq_file = @{$files}[0];
  my $data_file = @{$files}[1];

  my $unigene_source_id = $self->get_source_id_for_source_name('UniGene');

  print "UniGene source ID = $unigene_source_id.\n" if($verbose);

  my $xrefs =
    $self->create_xrefs( $unigene_source_id, $unigene_source_id,
      $uniq_file, $data_file, $species_id );

  if(!defined($xrefs)){
    return 1; #error
  }
  if(!defined($self->upload_xref_object_graphs($xrefs))){
    print "Read " . scalar(@$xrefs) ." xrefs from $uniq_file\n" if($verbose);
    return 1; # error
  }

  if ( defined $release_file ) {
    # Get species name from species ID.
    my $species_name;

    my $sth =
      $self->dbi()
	->prepare("SELECT name FROM species WHERE species_id = ?");

    $sth->execute($species_id);
    $sth->bind_columns( \$species_name );
    $sth->fetchrow_array();

    $species_name =~ tr/_/ /;

    # Parse and set release info.
    my $release;
    my $release_io = $self->get_filehandle($release_file);

    while ( defined( my $line = $release_io->getline() ) ) {
      if ( $line =~ /^(.*$species_name)/i ) {
	$release = $1;
      }
    }
    $release_io->close();

    if ( defined $release ) {
      $release =~ s/\s{2,}/ /g;
      $release =~ s/^(.*) UniGene/$1, UniGene/;

      print "UniGene release: '$release'\n" if($verbose);
      $self->set_release( $unigene_source_id, $release );
    }
  }

  return 0; # successfull

}


my %geneid_2_desc;

sub get_desc{
  my $self = shift;
  my $data_file = shift;

  my $dir = dirname($data_file);

  local $/ = "//";

  my $desc_io = $self->get_filehandle( $data_file );

  if ( !defined $desc_io ) {
    print STDERR "ERROR: Can't open $data_file\n";
    return;
  }

  while ( $_ = $desc_io->getline() ) {
    #ID          Hs.159356
    #TITLE       Hypothetical LOC388277
    
    (my $id) = $_ =~ /ID\s+(\S+)/;
    (my $descrip) = $_ =~ /TITLE\s+(.+)\n/;

    if ( defined $id && defined $descrip ) {
        $geneid_2_desc{$id} = $descrip;
    }
    
  }

  $desc_io->close();

  return 1;
}


sub create_xrefs {
  my ($self, $peptide_source_id, $unigene_source_id, $uniq_file, $data_file, $species_id ) = @_;

  # Create a hash of all valid names for this species. Not used...
  # my %species2name = $self->species_id2name();
  # my @names   = @{$species2name{$species_id}};
  # my %name2species_id     = map{ $_=>$species_id } @names;

  if ( !defined( $self->get_desc($data_file) ) ) {
    return;
  }

  my $unigene_io = $self->get_filehandle($uniq_file);

  if ( !defined $unigene_io ) {
    print STDERR "Can't open RefSeq file $uniq_file\n";
    return;
  }

#>gnl|UG|Hs#S19185843 Homo sapiens N-acetyltransferase 2 (arylamine N-acetyltransferase)
  # , mRNA (cDNA clone MGC:71963 IMAGE:4722596), complete cds /cds=(105,977) /gb=BC067218 /gi=45501306 /ug=Hs.2 /len=1344
#GGGGACTTCCCTTGCAGACTTTGGAAGGGAGAGCACTTTATTACAGACCTTGGAAGCAAG


  my @xrefs;

  local $/ = "\n>";

  while ( $_ = $unigene_io->getline() ) {

    my $xref;

    my $entry = $_;
    chomp $entry;
    my ($header, $sequence) = split (/\n/, $entry, 2);
    $sequence =~ s/^>//;
    # remove newlines
    my @seq_lines = split (/\n/, $sequence);
    $sequence = join("", @seq_lines);

#    (my $gnl, my $n, my $rest) = split(/\|/, $header,3);

    (my $acc_no_ver) = $header =~ /\/ug=(\S*)/;

    if(!defined($geneid_2_desc{$acc_no_ver})){
      print "****$_\n";
      $geneid_2_desc{$acc_no_ver} = "";
      warn "No desc for $acc_no_ver\n";
    }


    $xref->{SEQUENCE_TYPE} = 'dna';
    $xref->{STATUS} = 'experimental';
    $xref->{SOURCE_ID} = $unigene_source_id;
   

    ##No species check as files contain data  fro only one species.
    
    $xref->{ACCESSION} = $acc_no_ver;
    $xref->{LABEL} = $acc_no_ver;
    $xref->{DESCRIPTION} = $geneid_2_desc{$acc_no_ver};
    $xref->{SEQUENCE} = $sequence;
    $xref->{SPECIES_ID} = $species_id;
    $xref->{INFO_TYPE} = "SEQUENCE_MATCH";
    
    push @xrefs, $xref;
    
  }

  $unigene_io->close();

  %geneid_2_desc=();

  return \@xrefs;

}

1;
