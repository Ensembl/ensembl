=head1 LICENSE

See the NOTICE file distributed with this work for additional information
regarding copyright ownership.

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

package XrefParser::UniProtVarSplicParser;

# Parse UniProt alternative splice files

use strict;
use warnings;
use Carp;
use File::Basename;

use base qw( XrefParser::BaseParser );

# UniProtVarSplic file format: fasta, e.g.

#>P48347-2|14310_ARATH Isoform 2 of P48347 - Arabidopsis thaliana (Mouse-ear cress)
#MENEREKQVYLAKLSEQTERYDEMVEAMKKVAQLDVELTVEERNLVSVGYKNVIGARRAS
#WRILSSIEQKEESKGNDENVKRLKNYRKRVEDELAKVCNDILSVIDKHLIPSSNAVESTV
#FFYKMKGDYYRYLAEFSSGAERKEAADQSLEAYKAAVAAAENGLAPTHPVRLGLALNFSV
#FYYEILNSPESACQLAKQAFDDAIAELDSLNEESYKDSTLIMQLLRDNLTLWTSDLNEEG
#DERTKGADEPQDEV

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

  local $/ = "\n>";

  my $file_io = $self->get_filehandle($file);

  if ( !defined $file_io ) {
    print STDERR "ERROR: Could not open $file\n";
    return 1;    # 1 error
  }

  my %swiss = %{ $self->get_valid_codes( "uniprot", $species_id ) };

  print scalar(%swiss)." uniprot entries will be used as tests\n" if($verbose);
  my $missed = 0;
  while ( $_ = $file_io->getline() ) {
    my $xref;

    my ($header, $sequence) = $_ =~ /^>?(.+?)\n([^>]*)/s or warn("Can't parse FASTA entry: $_\n");

    # deconstruct header
    my ($accession, @description) = split /\|/, $header;
    my $description = join(" ", @description);
    
    my ($original, $extension) = split/-/, $accession;

    if(defined($swiss{$original})){
      # make sequence into one long string
      $sequence =~ s/\n//g;
      
      # build the xref object and store it
      $xref->{ACCESSION}     = $accession;
      $xref->{LABEL}         = $accession;
      $xref->{DESCRIPTION}   = $description;
      $xref->{SEQUENCE}      = $sequence;
      $xref->{SOURCE_ID}     = $source_id;
      $xref->{SPECIES_ID}    = $species_id;
      $xref->{SEQUENCE_TYPE} = 'peptide';
      $xref->{STATUS}        = 'experimental';
      
      push @xrefs, $xref;
    }
    else{
     $missed++;
    }
  }

  $file_io->close();

  print $missed." ignored as original uniprot not found in database\n" if($verbose);
  print scalar(@xrefs) . " UniProtVarSplic xrefs succesfully parsed\n" if($verbose);

  $self->upload_xref_object_graphs(\@xrefs);

    if ( defined $release_file ) {
        # Parse and apply the Swiss-Prot release info
        # from $release_file.
        my $release_io = $self->get_filehandle($release_file);
        while ( defined( my $line = $release_io->getline() ) ) {
            if ( $line =~ m#(UniProtKB/Swiss-Prot Release .*)# ) {
                print "Swiss-Prot release is '$1'\n" if($verbose);
                $self->set_release( $source_id, $1 );
            }
        }
        $release_io->close();
    }


  return 0;
}

1;
